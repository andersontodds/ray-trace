##RayTrace_lat.jl
# Todd Anderson
# 5 October 2020
#
# Ray-tracing of electromagnetic waves from upper ionosphere to plasmasphere using Haselgrove equations, following ray-tracing equations in Rice 1997 combined with plasmasphere model adapted from Carpenter and Anderson 1992.  Includes main ODE solver function calls.
#

using DifferentialEquations
using BenchmarkTools
using LSODA
using Sundials

c = 2.99792458e8
re = 6.3712e6      	# radius of Earth in meters
B0 = 3.0696381e-5   # magnitude of B field at Earth's equator surface in Tesla
e = 1.602e-19 		# electric charge in Coulombs
me = 9.1093e-31		# electron rest mass
mp = 1.6726219e-27	# proton rest mass
eps = 8.854e-12 	# permittivity of free space in mks

# pre-solve plasmasphere
r_range = [re:1000:(10*re);]
λ_range = 0
L = r_range./(re*cos(λ_range)^2)

Kp_max = 3                             # for example
Lppi = 5.6 - 0.46*Kp_max                # plasmapause inner limit in L
#Lppo = 5.75                             # placeholder
d = 0                                   # day number
R̄ = 90                                  # 13-month average sunspot number
mlt = 2                                   # magnetic local time

ne_Lppi = 10^((-0.3145*Lppi + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2-Lppi)/1.5))

# pre-solve individual components in order to find Lppo
ne_plasma_1 = @. 10 .^((-0.3145.*L + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2 .-L)./1.5))
ne_plasma_2 = ne_Lppi*10 .^(-1.0.*(L .- Lppi)./0.1)
ne_plasma_3 = (5800 + 300*mlt).*L.^(-4.5) + (1 .- exp.((2 .-L)./10))

min_index = findmin(abs.(ne_plasma_2 - ne_plasma_3))[2]
Lppo = L[min_index]

# calculate phase refractive index μ
u = function phase_refractive_index(r,λ,χ,freq)
	# u = [μ, dμdψ]
    # convert from radial angle to wave normal angle
    dip = atan(2.0*tan(λ)) # dip angle: angle between horizontal and B field;

	#	ψ = 2*π - (3π/2 - dip - χ) 		# this should be ≡ to the below, but the plotted ray paths are not all that similar!
	ψ = π/2.0 + dip + χ
	#	ψ = χ - (3.0/2.0)*pi + dip		# wave normal angle: angle between wave direction and B; 11/9: fixed error (χ - ϕ) ➡ (ϕ - χ)
	# NOTE 11/24/2020: Rice 1997 uses ψ = χ - 3π/2 + dip, as well as a strange method of calculating χ (see FORTRAN code on p. 112)

    # convert from frequency to angular frequency
    ω = 2.0*pi*freq

    # find magnetic field at (r,θ) from dipole field model
    # Bmag = B0*((re^3)/(r^3))*sqrt(4.0-3.0*cos((pi/2.)-θ)*cos((pi/2.0)-θ))
	# 	   = B0*((re^3)/(r^3))*sqrt(4 - 3sin²θ)
	# 	   ≡ B0*((re^3)/(r^3))*sqrt(1 + 3cos²θ)
	# dipole field magnitude from Parks p. 61:
	# B(r,λ) = μ₀M/(4πr³)*(1+3sin²λ)^(1/2)
	# 		 = μ₀M/(4πr³)*(1+3sin²(π/2 - θ))^(1/2)
	# 		 = μ₀M/(4πr³)*(1+3cos²θ)^(1/2)
	# here B0*re³ = μ₀M/(4π)
	Bmag = B0*(re^3/(r^3))*sqrt(1+3*sin(λ)*sin(λ))

    # calculate electron and proton density profiles
	Lshell = r./(re*cos(λ)^2)
	ne_iono = (1.8e5.*exp.(-4.183119.*((r./re).-1.0471))) # cm^-3

	if Lshell <= Lppi
	     log_ne = (-0.3145*Lshell + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2-Lshell)/1.5)
	     ne_plasma = 10^(log_ne)
	elseif Lppi < Lshell <= Lppo
	     ne_plasma = ne_Lppi*10^((Lppi-Lshell)/0.1)
	elseif Lppo < Lshell
	     ne_plasma = (5800 + 300*mlt)*Lshell^(-4.5) + (1 - exp((2-Lshell)/10))
	else
	     ne_plasma = 0.0
	end

	n_e = (ne_iono + ne_plasma)*1e6
	n_p = n_e

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
dμdr = function ddr(r,λ,χ,freq)
	#μ,dμdψ = u			# unpack
	dr = 1.0e-11		# r step size

	μ_l = phase_refractive_index(r-dr/2.0,λ,χ,freq)
	#μ_l = μ
	μ_r = phase_refractive_index(r+dr/2.0,λ,χ,freq)
	#μ_r = μ

	dμdr = (μ_r[1] - μ_l[1])/dr
end

dμdλ = function ddλ(r,λ,χ,freq)
	#μ,dμdψ = u			# unpack
	dλ = 1.0e-11		# θ step size

	μ_l = phase_refractive_index(r,λ-dλ/2.0,χ,freq)
	#μ_l = μ
	#println("μ_l = ", μ)
	μ_r = phase_refractive_index(r,λ+dλ/2.0,χ,freq)
	#μ_r = μ
	#println("μ_r = ", μ)

	dμdλ = (μ_r[1] - μ_l[1])/dλ
end

dμdχ = function ddχ(r,λ,χ,freq)
	#μ,dμdψ = u			# unpack
	dχ = 1.0e-11		# χ step size

	μ_l = phase_refractive_index(r,λ,χ-dχ/2.0,freq)
	#μ_l = μ
	μ_r = phase_refractive_index(r,λ,χ+dχ/2.0,freq) # NOTE Nov 3: corrected error 'χ-dχ/2.0'➡'χ+dχ/2.0'
	#μ_r = μ

	dμdχ = (μ_r[1] - μ_l[1])/dχ
end

dμdf = function ddf(r,λ,χ,freq)
	#μ,dμdψ = u			# unpack
	df = 1.0e-5		# f step size

	μ_l = phase_refractive_index(r,λ,χ,freq-df/2.0)
	#μ_l = μ
	μ_r = phase_refractive_index(r,λ,χ,freq+df/2.0)
	#μ_r = μ

	dμdf = (μ_r[1] - μ_l[1])/df
end

# calculate derivatives w.r.t. time

# calculate derivatives w.r.t. time
function haselgrove!(du,u,p,t)
# 	du = drdt, dλdt, dχdt, dfdt
# 	u = r, λ, χ, f
#  	p = []
# 	t = time

	r, λ, χ, freq = u 						# unpack u
	#freq = p[1] 	# unpack p

	dμdr = ddr(r,λ,χ,freq)			# mutating form: ddr!([μ,dμdψ],r,λ,χ,freq,dμdr)
	dμdλ = ddλ(r,λ,χ,freq)			# calls phase_refractive_index(r,λ,χ,freq)
	dμdχ = ddχ(r,λ,χ,freq)
	dμdf = ddf(r,λ,χ,freq)

	v = phase_refractive_index(r,λ,χ,freq)
	μ = v[1]
	dμdψ = v[2]

	# 11/20: try dμdχ --> -1*dμdψ; these should be equal but are calcualted differently
    du[1] = 1/(μ^2)*(μ*cos(χ)+dμdψ*sin(χ))
    du[2] = 1/(r*μ^2)*(μ*sin(χ)-dμdψ*cos(χ))
    du[3] = 1/(r*μ^2)*(dμdλ*cos(χ)-(r*dμdr + μ)*sin(χ))
	du[4] = 1/c*(1+(freq/μ)*dμdf)

	#du = [drdt, dλdt, dχdt, dfdt]

	#println("drdt = ", drdt)
	#println("dλdt = ", dθdt)
	#println("dχdt = ", dχdt)


	#DEBUG check values
	#println("drdt = ",drdt)
	#println("dλdt = ",dθdt)
	#println("dχdt = ",dχdt)


    # dμdψ = 1/(2*μ)*((dBdψ + dFdψ)/(2*A) - 2*dAdψ*(B+F)/(2*A^2))
    # dμdω = 1/(2*μ)*((dBdω + dFdω)/(2*A) - 2*dAdω*(B+F)/(2*A^2))
	#
    # dμdr = (μ*(r+δr/2) - μ*(r-δr/2))/δr
    # dμdθ = (μ*(θ+δθ/2) - μ*(θ-δθ/2))/δθ



end

#DEBUG check values
#haselgrove!(1.0,1.0,10000.0,pi/4.0,0.0,5000.0,0.0,[1.0,1.0,1.0,1.0],[1.0,1.0,1.0])

## callbacks!
# re_cb: terminate when ray intersects Earth surface
function re_term_condition(u,t,integrator)
	u[1] - re
end

# function μ2_term_condition(u,t,integrator)
# 	body
# end

function terminate_affect!(integrator)
	terminate!(integrator)
end

re_cb = ContinuousCallback(re_term_condition,terminate_affect!)

# saveμ_cb: save μ, dμ values at each timestep
function save_func(u,t,integrator)
	r, λ, χ, f = u
	v = phase_refractive_index(r, λ, χ, f)
	dip = atan(2.0*tan(λ))
	ψ = π/2.0 + dip + χ
	vcat(v,dip,ψ)
end

saved_μ = SavedValues(Float64, Array{Float64,1})
saveμ_cb = SavingCallback(save_func, saved_μ)

cb = CallbackSet(re_cb, saveμ_cb)

## Plot results
u0 = [re+1.0e+6, 1.0*pi/4, 0.0, 5000.0]					# r0, λ0, χ0, f0
p = []	# f0, dμdψ, dμdr, dμdθ, dμdχ, dμdf
tspan = (0.0,5.0e+9)

hasel_prob = ODEProblem(haselgrove!,u0,tspan,p)
hasel_soln = solve(hasel_prob, CVODE_BDF(), reltol=1e-7, callback=cb)

using Plots
plotly()
plot(hasel_soln)
plot(hasel_soln, vars=(1,2,3))

t = hasel_soln.t
r = hasel_soln[1,:]
λ = hasel_soln[2,:]
χ = hasel_soln[3,:]
f = hasel_soln[4,:]

x = r.*cos.(λ)
y = r.*sin.(λ)

plot(re.*cos.([0:0.01:2*π;]),re.*sin.([0:0.01:2*π;]), aspect_ratio = 1, legend=:none)
plot!(x,y, aspect_ratio=1)

## saved values
t_s = saved_μ.t
val_s = saved_μ.saveval
val_s_r = reduce(hcat, val_s)
val_s_r = val_s_r'
μ_s = val_s_r[:,1]
dμdψ_s = val_s_r[:,2]
dip_s = val_s_r[:,3]
ψ_s = val_s_r[:,4]

plot(t_s, μ_s, title="refractive index v. time")
plot(t_s, dμdψ_s, title="d/dpsi refractive index v. time")
plot(t_s, dip_s.*(180/pi), title="dip angle v. time")
plot(t_s, ψ_s.*(180/pi), title="wave normal angle v. time")

## refractive index surface
# initial conditions
x_test = 1*9.031E+6
y_test = 1*1.537E+6
λ_test = atan(abs(y_test/x_test))
r_test = sqrt(x_test^2 + y_test^2) #2.0*re*(sin(θ_test))^2 # on plotted B field line
#θ_test = θ_pw
#r_test = r_pw
ψ_test = [0:0.001:2π;]
dip_test = atan(2.0*tan(λ_test))
#χ_test = ψ_test .+ 3π/2 .- dip_test
χ_test = -1.0.*ψ_test .+ 3.0*π/2 .- dip_test
f_test = 5000.0

# calculate μ, dμdψ
u_test = phase_refractive_index.(r_test, λ_test, χ_test, f_test)
u_test_r = reduce(hcat,u_test)
u_test_r = u_test_r'
μ_test =  u_test_r[:,1]
dμdψ_test = u_test_r[:,2]

# xμ_test = μ_test.*sin.(ψ_test .+ π/2)
# yμ_test = μ_test.*cos.(ψ_test .+ π/2)
# # plot surface rotated into B|| frame
xμ_test_B = μ_test.*sin.(ψ_test)
yμ_test_B = μ_test.*cos.(ψ_test)

plot(xμ_test_B, yμ_test_B, aspect_ratio = 1)

# plot surface in x, y frame
xμ_test_xy = μ_test.*sin.(χ_test .- (π/2 - λ_test))
yμ_test_xy = μ_test.*cos.(χ_test .- (π/2 - λ_test))

plot(xμ_test_xy, yμ_test_xy, aspect_ratio = 1)

xdμ_test = dμdψ_test.*sin.(ψ_test)
ydμ_test = dμdψ_test.*cos.(ψ_test)
