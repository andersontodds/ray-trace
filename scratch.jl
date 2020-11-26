# scratch.jl
# scratch work ray-trace
#

##
#= Functions:
	phase_refractive_index!():
		inputs:	μ, dμdψ, r, θ, χ, freq
		outputs: mutates μ, dμdψ

	ddr!():
		inputs:	μ, dμdψ, r, θ, χ, freq, dμdr
		outputs: mutates dμdr
		calls: phase_refractive_index!()

	ddθ!():
		inputs:	μ, dμdψ, r, θ, χ, freq, dμdθ
		outputs: mutates dμdθ
		calls: phase_refractive_index!()

	ddχ!():
		inputs:	μ, dμdψ, r, θ, χ, freq, dμdχ
		outputs: mutates dμdχ
		calls: phase_refractive_index!()

	ddf!():
		inputs:	μ, dμdψ, r, θ, χ, freq, dμdf
		outputs: mutates dμdf
		calls: phase_refractive_index!()

	haselgrove!():
		inputs: μ, dμdψ, r, θ, χ, freq, t, dμ, ddt
				dμdr,dμdθ,dμdχ,dμdf = dμ
				ddt = [drdt, dθdt, dχdt]
		outputs: mutates ddt


=#

#=TODO:
	1. Get functions and argument-passing working with single input (i.e. non-vectorized, non-looping)
	2. Vectorize!
	3. Run and time solver
	4. Switch to in-place functions where possible
	5. Run and time solver, optimize
=#
#const rₑ = 6378100
c = 2.99792458e8
re = 6.3712e6      	# radius of Earth in meters
B0 = 3.0696381e-5   # magnitude of B field at Earth's equator surface in Tesla
e = 1.602e-19 		# electric charge in Coulombs
me = 9.1093e-31		# electron rest mass
mp = 1.6726219e-27	# proton rest mass
eps = 8.854e-12 	# permittivity of free space in mks

using DifferentialEquations

# calculate phase refractive index μ
# DONE! todo: convert this function to mutating form
# 	QUESTION: is this possible, given that phase_refractive_index!() must be called multiple times in a single ddr!() (e.g.) call, before using the result?
# 	ANSWER: yes! see call/assignment in ddr!()

u = function phase_refractive_index(r,θ,χ,freq)
# u = [μ, dμdψ]
    # convert from radial angle to wave normal angle
    dip = atan(2.0*cot(θ))     			# dip angle: angle between horizontal and B field;
	# 11/9: fixed error tan(pi/2.0-θ) ➡ tan(θ)
	# 11/13: no! dip = 0 when θ = π/2 ➡ tan(θ - π/2) is correct!
	# 11/18: wrong again! Dip is defined as positive toward the North pole, i.e. should be positive in the Northern hemisphere!
	# 11/19: really need to check this; looks like factor of 2 was probably not correct
#	ψ = 2*π - (3π/2 - dip - χ) 		# this should be ≡ to the below, but the plotted ray paths are not all that similar!
	ψ = π/2.0 + dip + χ
	#ψ = (3.0/2.0)*pi - dip - χ
#	ψ = χ - (3.0/2.0)*pi + dip		# wave normal angle: angle between wave direction and B; 11/9: fixed error (χ - ϕ) ➡ (ϕ - χ)
	# NOTE 11/24/2020: Rice 1997 uses ψ = χ - 3π/2 + dip, as well as a strange method of calculating χ (see FORTRAN code on p. 112)

    # convert from frequency to angular frequency
    ω = 2.0*pi*freq

    # convert radius to fraction of Earth radius
    rE = r/re

    # find magnetic field at (r,θ) from dipole field model
    # Bmag = B0*((re^3)/(r^3))*sqrt(4.0-3.0*cos((pi/2.)-θ)*cos((pi/2.0)-θ))
	# 	   = B0*((re^3)/(r^3))*sqrt(4 - 3sin²θ)
	# 	   ≡ B0*((re^3)/(r^3))*sqrt(1 + 3cos²θ)
	# dipole field magnitude from Parks p. 61:
	# B(r,λ) = μ₀M/(4πr³)*(1+3sin²λ)^(1/2)
	# 		 = μ₀M/(4πr³)*(1+3sin²(π/2 - θ))^(1/2)
	# 		 = μ₀M/(4πr³)*(1+3cos²θ)^(1/2)
	# here B0*re³ = μ₀M/(4π)
	Bmag = B0*(re^3/(r^3))*sqrt(1+3*cos(θ)*cos(θ))

    # calculate electron and proton density profiles
    n_e = 1.e6*(1.8e5*exp(-4.183119*(rE-1.0471)))
	n_p = 1.e6*(1.8e5*exp(-4.183119*(rE-1.0471)))

	# electron plasma angular frequency squared
    ω_e2 = (n_e*(e^2))/(eps*me)
	# proton plasma angular frequency squared
    ω_p2 = (n_p*(e^2))/(eps*mp)

	# electron angular gyrofrequency
    Ω_e = (e*Bmag)/(me)
	# proton angular gyrofrequency
    Ω_p = (e*Bmag)/(mp)

	# define R, L, P, D, S,
	# R ≡ 1 - Σ(ωₖ²/ω²)(ω/(ω+ϵₖΩₖ))
	R = 1.0 - (ω_e2/ω^2.0)*(ω/(ω - Ω_e)) - (ω_p2/ω^2.0)*(ω/(ω + Ω_p))
	#println("R = ",R)

	# L ≡ 1 - Σ(ωₖ²/ω²)(ω/(ω-ϵₖΩₖ))
	L = 1.0 - (ω_e2/ω^2.0)*(ω/(ω + Ω_e)) - (ω_p2/ω^2.0)*(ω/(ω - Ω_p))
	#println("L = ",L)

	# P ≡ 1 - Σ(ωₖ²/ω²)
	P = 1.0 - (ω_e2/ω^2.0) - (ω_p2/ω^2.0)
	#println("P = ",P)

	# D ≡ ½(R-L); S ≡ ½(R+L)
	D = (R - L)/2.0
	S = (R + L)/2.0

	# define parts of dispersion relation
	# Aμ⁴ - Bμ² + C = 0
	# A = Ssin²ψ + Pcos²ψ
	A = S*sin(ψ)^2.0 + P*cos(ψ)^2.0

	# B = RLsin²ψ + PS(1 + cos²ψ)
	B = R*L*sin(ψ)^2.0 + P*S*(1.0+cos(ψ)^2.0)

	# C  = PRL
	C = P*R*L

	# solve the dispersion relation for μ:
	# μ² = (B +- F)/2A
	# where F² = (RL-PS)²sin⁴ψ + 4P²D²cos²ψ
	#QUESTION: is it faster to (a) solve for F², then sqrt(); or (b) solve for F directly? (a) requires fewer sqrt() calls
	F2 = (R*L - P*S)^2.0*sin(ψ)^4.0 + 4.0*(P*D*cos(ψ))^2.0
	F = sqrt(F2)

	# Rice 1997: typical solution to dispersion relation with quadratic formula
	μ2_minus = (B - F)/(2.0*A)
	μ2_plus = (B + F)/(2.0*A)

	# Bortnik 2004:
	# if B > 0
	# 	μ2_minus = (B - F)/(2.0*A)
	# else
	# 	μ2_plus = (2.0*C)/(B + F)
	# end


	μ_minus = try
		sqrt(abs(μ2_minus))	# abs() is not physical! for test only

	catch err
		if isa(err,DomainError)
			println("DomainError in μ_minus!")
			println("B=", B)
			println("F=", F)
			println("A=", A)
		else
			# println(err)
		end
	end

	μ_plus = try
		sqrt(abs(μ2_plus)) # abs() is not physical! for test only

	catch err
		if isa(err,DomainError)
			# println("DomainError in μ_plus!")
			# println("B=", B)
			# println("F=", F)
			# println("A=", A)
		else
			# println(err)
		end
	end



	# Electron whistler case: ψ = 0, μ² = R; this is the μ_plus case
	#μ = μ_minus		# EMIC case
	μ = μ_plus		# electron whistler case

	# Find dA/dψ, dB/dψ, dC/dψ, dμ/dψ
	dAdψ = 2.0*(S-P)*sin(ψ)*cos(ψ)
	dBdψ = 2.0*(R*L-P*S)*sin(ψ)*cos(ψ)
    dCdψ = 0.0
    #dμdψ = ((μ^4.0)*dAdψ-(μ^2.0)*dBdψ+dCdψ)/(4.0*A*(μ^3.0)-2.0*B*μ)

	dFdψ = 1/(2.0*F)*((R*L-P*S)^2 * 4*sin(ψ)^3*cos(ψ) - 8*(P*D)^2*sin(ψ)*cos(ψ))
	#dFdψ = sqrt(abs((R*L-P*S)^2 * 4*sin(ψ)^3*cos(ψ) - 8*(P*D)^2*sin(ψ)*cos(ψ)))
	dμdψ = 1/(2.0*μ)*((dBdψ + dFdψ)/(2*A) - 2*dAdψ*(B + F)/(2*A^2))
	# NOTE 11/24/2020: from Rice 1997; choose '+' for (B +- F) and (dBdψ +- dFdψ)

	#DEBUG check values
	#println("μ = ",μ)
	#println("dμdψ = ", dμdψ)

	# non-mutating output
	u = [μ,dμdψ]

end

# Derivatives of phase refractive index w.r.t. r, θ, χ,
# calculate derivative of μ w.r.t. r
dμdr = function ddr(r,θ,χ,freq)
	#μ,dμdψ = u			# unpack
	dr = 1.0e-11		# r step size

	μ_l = phase_refractive_index(r-dr/2.0,θ,χ,freq)
	#μ_l = μ
	μ_r = phase_refractive_index(r+dr/2.0,θ,χ,freq)
	#μ_r = μ

	dμdr = (μ_r[1] - μ_l[1])/dr
end

dμdθ = function ddθ(r,θ,χ,freq)
	#μ,dμdψ = u			# unpack
	dθ = 1.0e-11		# θ step size

	μ_l = phase_refractive_index(r,θ-dθ/2.0,χ,freq)
	#μ_l = μ
	#println("μ_l = ", μ)
	μ_r = phase_refractive_index(r,θ+dθ/2.0,χ,freq)
	#μ_r = μ
	#println("μ_r = ", μ)

	dμdθ = (μ_r[1] - μ_l[1])/dθ
end

dμdχ = function ddχ(r,θ,χ,freq)
	#μ,dμdψ = u			# unpack
	dχ = 1.0e-11		# χ step size

	μ_l = phase_refractive_index(r,θ,χ-dχ/2.0,freq)
	#μ_l = μ
	μ_r = phase_refractive_index(r,θ,χ+dχ/2.0,freq) # NOTE Nov 3: corrected error 'χ-dχ/2.0'➡'χ+dχ/2.0'
	#μ_r = μ

	dμdχ = (μ_r[1] - μ_l[1])/dχ
end

dμdf = function ddf(r,θ,χ,freq)
	#μ,dμdψ = u			# unpack
	df = 1.0e-5		# f step size

	μ_l = phase_refractive_index(r,θ,χ,freq-df/2.0)
	#μ_l = μ
	μ_r = phase_refractive_index(r,θ,χ,freq+df/2.0)
	#μ_r = μ

	dμdf = (μ_r[1] - μ_l[1])/df
end

# calculate derivatives w.r.t. time
function haselgrove!(du,u,p,t)
# 	du = drdt, dθdt, dχdt, dfdt
# 	u = r, θ, χ, f
#  	p = freq, μ, dμdψ, dμdr, dμdθ, dμdχ, dμdf
# 		^ this compiles and runs, but parameters
# 	t = time

	r, θ, χ, freq = u 						# unpack u
	#freq = p[1] 	# unpack p

	dμdr = ddr(r,θ,χ,freq)			# mutating form: ddr!([μ,dμdψ],r,θ,χ,freq,dμdr)
	dμdθ = ddθ(r,θ,χ,freq)			# calls phase_refractive_index(r,θ,χ,freq)
	dμdχ = ddχ(r,θ,χ,freq)
	dμdf = ddf(r,θ,χ,freq)

	v = phase_refractive_index(r,θ,χ,freq)
	μ = v[1]
	dμdψ = v[2]

	# 11/20: try dμdχ --> -1*dμdψ; these should be equal but are calcualted differently
    du[1] = 1/(μ^2)*(μ*cos(χ)-dμdψ*sin(χ))
    du[2] = 1/(r*μ^2)*(μ*sin(χ)+dμdψ*cos(χ))
    du[3] = 1/(r*μ^2)*(dμdθ*cos(χ)-(r*dμdr + μ)*sin(χ))
	du[4] = 1/c*(1+(freq/μ)*dμdf)

	#du = [drdt, dθdt, dχdt, dfdt]

	#println("drdt = ", drdt)
	#println("dθdt = ", dθdt)
	#println("dχdt = ", dχdt)


	#DEBUG check values
	#println("drdt = ",drdt)
	#println("dθdt = ",dθdt)
	#println("dχdt = ",dχdt)


    # dμdψ = 1/(2*μ)*((dBdψ + dFdψ)/(2*A) - 2*dAdψ*(B+F)/(2*A^2))
    # dμdω = 1/(2*μ)*((dBdω + dFdω)/(2*A) - 2*dAdω*(B+F)/(2*A^2))
	#
    # dμdr = (μ*(r+δr/2) - μ*(r-δr/2))/δr
    # dμdθ = (μ*(θ+δθ/2) - μ*(θ-δθ/2))/δθ



end

#DEBUG check values
#haselgrove!(μ,dμdψ,r,θ,χ,freq,dμ,ddt,t)
#haselgrove!(1.0,1.0,10000.0,pi/4.0,0.0,5000.0,[1.0,1.0,1.0,1.0],[1.0,1.0,1.0],0.0)

# du0 = [1.0, 1.0, 1.0, 1.0]
# u0 = [8000.0, pi/4.0, 0.0, 5000.0]
# p = []
# t0 = 0.0
# haselgrove!(du0, u0, p, t0)
# println(du0)

## ODE solver!
#QUESTION: how does ODEProblem work?

using BenchmarkTools
using LSODA
using Sundials

u0 = [re+1.0e+6, 1.0*pi/4, 0.0, 1000.0]					# r0, θ0, χ0
p = []	# f0, dμdψ, dμdr, dμdθ, dμdχ, dμdf
tspan = (0.0,5.0e+9)

hasel_prob = ODEProblem(haselgrove!,u0,tspan,p)
hasel_soln = solve(hasel_prob, CVODE_BDF(), reltol=1e-7) #reltol=1e-4
# LSODA_BDF()?

using Plots
plotly()
plot(hasel_soln)
plot(hasel_soln, vars=(1,2,3))

t = hasel_soln.t
r = hasel_soln[1,:]
θ = hasel_soln[2,:]
χ = hasel_soln[3,:]
f = hasel_soln[4,:]

x = r.*sin.(θ)
y = r.*cos.(θ)

plot(x,y, aspect_ratio=1, legend=:none)


#earthx = re/1.0e3.*cos([0:0.01:2*π;])
plot!(re.*cos.([0:0.01:2*π;]),re.*sin.([0:0.01:2*π;]), aspect_ratio = 1)
# plot!(re.*[0:0.01:2;],re.*[0:0.01:2;], aspect_ratio = 1)
# plot!(re.*[0:0.01:4;],re.*[0:0.005:2;], aspect_ratio = 1)


gridrange = [-pi/2.0:0.01:pi/2.0;]
#thtgrid = np.arange(0.,np.pi,0.1)
gridmag = pi/2.0.-gridrange
rmag = re.*(cos.(gridmag)).^2

#rw = 2.*ray_cmks.Re
xmag1 = 1.0.*rmag.*cos.(gridmag)
ymag1 = 1.0.*rmag.*sin.(gridmag)

xmag2 = 2.0.*rmag.*cos.(gridmag)
ymag2 = 2.0.*rmag.*sin.(gridmag)

xmag3 = 3.0.*rmag.*cos.(gridmag)
ymag3 = 3.0.*rmag.*sin.(gridmag)

xmag4 = 4.0.*rmag.*cos.(gridmag)
ymag4 = 4.0.*rmag.*sin.(gridmag)


plot!(xmag4,ymag4,color = :red)
plot!(xmag3,ymag3,color = :red)
plot!(xmag2,ymag2,color = :red)
plot!(xmag1,ymag1,color = :red)

plot!(xmag4,-ymag4,color = :red)
plot!(xmag3,-ymag3,color = :red)
plot!(xmag2,-ymag2,color = :red)
plot!(xmag1,-ymag1,color = :red)

θ_test = pi/2
plot!([-2e7:1e5:2e7;].*sin(θ_test),[-2e7:1e5:2e7;].*cos(θ_test))
x_dip = [0:1e5:1e6;].*sin(dip_test +pi/2 - θ_test)
y_dip = [0:1e5:1e6;].*cos(dip_test +pi/2 - θ_test)
x_dip4 = x_dip .+ 4.0*re*sin(θ_test)^2*sin(θ_test)
y_dip4 = y_dip .+ 4.0*re*sin(θ_test)^2*cos(θ_test)
plot!(x_dip4,y_dip4,aspect_ratio=1)

## calculate and plot refractive index surface
# initial conditions
θ_test = pi/3
r_test = 2.0*re*(sin(θ_test))^2
ψ_test = [0:0.01:2π;]
dip_test = atan(2.0*cot(θ_test))
χ_test = ψ_test .- 3π/2 .+ dip_test
f_test = 5000.0

# calculate μ, dμdψ
u_test = phase_refractive_index.(r_test, θ_test, χ_test, f_test)
u_test_r = reduce(hcat,u_test)
u_test_r = u_test_r'
μ_test =  u_test_r[:,1]
dμdψ_test = u_test_r[:,2]

# xμ_test = μ_test.*sin.(ψ_test .+ π/2)
# yμ_test = μ_test.*cos.(ψ_test .+ π/2)
# # plot surface without rotation into B|| frame
xμ_test = μ_test.*sin.(ψ_test)
yμ_test = μ_test.*cos.(ψ_test)

xdμ_test = dμdψ_test.*sin.(ψ_test)
ydμ_test = dμdψ_test.*cos.(ψ_test)

xμ_test_fieldline = 4.0*(re*(sin(θ_test))^2)*sin(θ_test) .+ xμ_test.*1e5
yμ_test_fieldline = 4.0*(re*(sin(θ_test))^2)*cos(θ_test) .+ yμ_test.*1e5

plot(xμ_test, yμ_test, aspect_ratio = 1)
#plot(xdμ_test, ydμ_test, aspect_ratio = 1)

plot!(xμ_test_fieldline, yμ_test_fieldline, aspect_ratio = 1)


plot!([-30:1:30;].*cos(dip_test),[-30:1:30;].*sin(dip_test))
plot!([-30:1:30;].*sin(θ_test+π/2),[-30:1:30;].*cos(θ_test+π/2))


##
function lorenz!(du, u, p, t)
    du[1] = 10.0*(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end
# call this function in a problem
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
sol = solve(prob)

plot(sol,vars=(1,2,3))

## B-field check

re = 6.3712e6      	# radius of Earth in meters
B0 = 3.0696381e-5	# magnetic field strentgh (where?) in Tesla
r = 5000.0
θ = pi/4.0
Bmag = B0*((re^3)/(r^3))*sqrt(4.0-3.0*cos((pi/2.)-θ)*cos((pi/2.0)-θ))

##

out = function fun(in1,in2)
	out = in1*in2
end

a = fun(1.0,2.0)

function mutefun!(out::Array{Int64},in1,in2)
	out .= in1*in2		# don't forget '.' notation for vectorized i/o!
end

function right1_fill_with_twos!(v::Vector{Int64})
	v[:] =[2 for ii in 1:length(v)]
end

a = [1:5;]

function varargs(args...)
	return args
end
varargs(1,2,3)

a = [1:5;]
sum(a)

## type/value check
r = re + 1.0e6
θ = pi/4.0
χ = 0.0
freq = 5000.0
dip = atan(2.0*cot(θ))
dip*180/pi	# dip angle: angle between horizontal and B field
ϕ = (3.0/2.0)*pi - dip              # intermediate angle -- NOT azimuth
ψ = ϕ - χ							# wave normal angle: angle between wave direction and B
ψ*180/π

# convert from frequency to angular frequency
ω = 2.0*pi*freq

# convert radius to fraction of Earth radius
rE = r/re
re

# find magnetic field at (r,θ) from dipole field model
B0
Bmag = B0*((re^3)/(r^3))*sqrt(4.0-3.0*cos((pi/2.0)-θ)*cos((pi/2.0)-θ))
B0*((re^3)/(r^3))
sqrt(4.0-3.0*cos((pi/2.0)-θ)*cos((pi/2.0)-θ))
# calculate electron and proton density profiles
n_e = 1.e6*(1.8e5*exp(-4.183119*(rE-1.0471)))
n_p = 1.e6*(1.8e5*exp(-4.183119*(rE-1.0471)))

# electron plasma angular frequency squared
ω_e2 = (n_e*(e^2))/(eps*me)
# proton plasma angular frequency squared
ω_p2 = (n_p*(e^2))/(eps*mp)

# electron angular gyrofrequency
Ω_e = (e*Bmag)/(me)
# proton angular gyrofrequency
Ω_p = (e*Bmag)/(mp)

# define R, L, P, D, S,
# R ≡ 1 - Σ(ωₖ²/ω²)(ω/(ω+ϵₖΩₖ))
R = 1.0 - (ω_e2/ω^2.0)*(ω/(ω - Ω_e)) - (ω_p2/ω^2.0)*(ω/(ω + Ω_p))
(ω_e2/ω^2.0)*(ω/(ω - Ω_e))
(ω_p2/ω^2.0)*(ω/(ω + Ω_p))
#println("R = ",R)

# L ≡ 1 - Σ(ωₖ²/ω²)(ω/(ω-ϵₖΩₖ))
L = 1.0 - (ω_e2/ω^2.0)*(ω/(ω + Ω_e)) - (ω_p2/ω^2.0)*(ω/(ω - Ω_p))
(ω_e2/ω^2.0)*(ω/(ω + Ω_e))
(ω_p2/ω^2.0)*(ω/(ω - Ω_p))
#println("L = ",L)

# P ≡ 1 - Σ(ωₖ²/ω²)
P = 1.0 - (ω_e2/ω^2.0) - (ω_p2/ω^2.0)
#println("P = ",P)

# D ≡ ½(R-L); S ≡ ½(R+L)
D = (R - L)/2.0
S = (R + L)/2.0

# define parts of dispersion relation
# Aμ⁴ - Bμ² + C = 0
# A = Ssin²ψ + Pcos²ψ
A = S*sin(ψ)^2.0 + P*cos(ψ)^2.0

# B = RLsin²ψ + PS(1 + cos²ψ)
B = R*L*sin(ψ)^2.0 + P*S*(1.0+cos(ψ)^2.0)

# C  = PRL
C = P*R*L

# solve the dispersion relation for μ:
# μ² = (B +- F)/2A
# where F² = (RL-PS)²sin⁴ψ + 4P²D²cos²ψ
#QUESTION: is it faster to (a) solve for F², then sqrt(); or (b) solve for F directly? (a) requires fewer sqrt() calls
F2 = (R*L - P*S)^2.0*sin(ψ)^4.0 + 4.0*(P*D*cos(ψ))^2.0
F = sqrt(F2)

μ2_minus = (B - F)/(2.0*A)
μ2_plus = (B + F)/(2.0*A)


##
hbi = 1.0e6
b = -π/3
re
α = acos((re*sin(b)^2.0 + cos(b)*((re + hbi)^2.0 - re^2.0 * sin(b)^2.0)^(0.5))/(re + hbi))
α*180/π
chi = b - α
chi*180/π


## plot magnetic field strength
# both plots: need to figure out how plot a log scale surface
histogram2d(randn(10000),randn(10000),nbins = 20)

plot(contourf(randn(10,20)))

# r,θ: need to figure out how to convert grid to x,y
r_gridmag = re:1e5:5e7
θ_gridmag = 0:0.01:2*π
f_bmag(r_gridmag, θ_gridmag) = begin
	B0*(re^3/(r_gridmag^3))*sqrt(1+3*cos(θ_gridmag)*cos(θ_gridmag))
end
R_gridmag = repeat(reshape(r_gridmag,1,:), length(θ_gridmag), 1)
T_gridmag = repeat(θ_gridmag, 1, length(r_gridmag))
Z_gridmag = map(f_bmag, R_gridmag, T_gridmag)
p1 = contour(r_gridmag, θ_gridmag, f_bmag, fill = true)
p2 = contour(r_gridmag, θ_gridmag, Z_gridmag)
plot(p1,p2)

# x,y: need to figure out how to exclude area inside Earth from grid
x_gridmag = -5e7:1e5:5e7
y_gridmag = x_gridmag
f_bmag_xy(x_gridmag, y_gridmag) = begin
	B0*(re^3/(sqrt(x_gridmag^2 + y_gridmag^2)^3))*sqrt(1+3*cos(atan(x_gridmag/y_gridmag))^2)
end
X_gridmag = repeat(reshape(x_gridmag,1,:), length(y_gridmag), 1)
Y_gridmag = repeat(y_gridmag, 1, length(x_gridmag))
Z_gridmag_xy = map(f_bmag_xy, X_gridmag, Y_gridmag)
p1_xy = contour(x_gridmag, y_gridmag, f_bmag_xy, fill = (true,cgrad(scale = :log)))
p2_xy = contour(x_gridmag, y_gridmag, Z_gridmag_xy)
plot(p1,p2)
