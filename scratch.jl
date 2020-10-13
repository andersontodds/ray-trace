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
re = 6.3781e6      	# radius of Earth in meters
B0 = 3.0696381e-5   # magnitude of B field (where?) in Tesla
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
    dip = atan(2.0*tan(pi/2.0-θ))     	# dip angle: angle between horizontal and B field
    ϕ = (3.0/2.0)*pi - dip              # intermediate angle -- NOT azimuth
    ψ = χ - ϕ							# wave normal angle: angle between wave direction and B

    # convert from frequency to angular frequency
    ω = 2.0*pi*freq

    # convert radius to fraction of Earth radius
    rE = r/re

    # find magnetic field at (r,θ) from dipole field model
    Bmag = B0*((re^3)/(r^3))*sqrt(4.0-3.0*cos((pi/2.)-θ)*cos((pi/2.0)-θ))

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

	μ2_minus = (B - F)/(2.0*A)
	μ_minus = sqrt(μ2_minus)
	μ2_plus = (B + F)/(2.0*A)
	μ_plus = sqrt(μ2_plus)

	# Electron whistler case: ψ = 0, μ² = R; this is the μ_minus case
	μ = μ_minus

	# Find dA/dψ, dB/dψ, dC/dψ, dμ/dψ
	dAdψ = 2.0*(S-P)*sin(ψ)*cos(ψ)
	dBdψ = 2.0*(R*L-P*S)*sin(ψ)*cos(ψ)
    dCdψ = 0.0
    dμdψ = ((μ^4.0)*dAdψ-(μ^2.0)*dBdψ+dCdψ)/(4.0*A*(μ^3.0)-2.0*B*μ)

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
	μ_r = phase_refractive_index(r,θ,χ-dχ/2.0,freq)
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
# 	du = drdt, dθdt, dχdt
# 	u = r, θ, χ
#  	p = freq, μ, dμdψ, dμdr, dμdθ, dμdχ, dμdf
# 		^ this compiles and runs, but parameters
# 	t = time

	r, θ, χ = u 						# unpack u
	freq = p[1] 	# unpack p

	dμdr = ddr(r,θ,χ,freq)			# mutating form: ddr!([μ,dμdψ],r,θ,χ,freq,dμdr)
	dμdθ = ddθ(r,θ,χ,freq)			# calls phase_refractive_index(r,θ,χ,freq)
	dμdχ = ddχ(r,θ,χ,freq)
	dμdf = ddf(r,θ,χ,freq)

	v = phase_refractive_index(r,θ,χ,freq)
	μ = v[1]
	dμdψ = v[2]

    du[1] = 1/(μ^2)*(μ*cos(χ)+dμdχ*sin(χ))
    du[2] = 1/(r*μ^2)*(μ*sin(χ)-dμdχ*cos(χ))
    du[3] = 1/(r*μ^2)*(dμdθ*cos(χ)-(r*dμdr + μ)*cos(χ))

	#du = [drdt, dθdt, dχdt]

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

du0 = [1.0, 1.0, 1.0]
u0 = [8000.0, pi/4.0, 0.0]
p = [5000.0]
t0 = 0.0
haselgrove!(du0, u0, p, t0)
println(du0)

## ODE solver!
#QUESTION: how does ODEProblem work?

using BenchmarkTools

u0 = [8000.0, pi/4.0, 0.0]					# r0, θ0, χ0
p = [5000.0]	# f0, dμdψ, dμdr, dμdθ, dμdχ, dμdf
tspan = (0.0,1.0)

hasel_prob = ODEProblem(haselgrove!,u0,tspan,p)
@benchmark hasel_soln = solve(hasel_prob)

using Plots
plot(hasel_soln)

##
function lorenz!(du, u, p, t)
    du[1] = 10.0*(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end
# call this function in a problem
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob2 = ODEProblem(lorenz!,u0,tspan)
sol2 = solve(prob2)

plot(sol2,vars=(1,2,3))

## try a nested ODE
function someODE!(du,u,p,t)
	du

end

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
# r = 8000.0
# θ = pi/4.0
# χ = 0.0
# freq = 5000.0
# dip = atan(2.0*tan(pi/2.0-θ))
# ϕ = (3.0/2.0)*pi - dip
# ψ = χ - ϕ
# ω = 2.0*pi*freq
# rE = r/re
# Bmag = B0*((re^3)/(r^3))*sqrt(4.0-3.0*cos((pi/2.)-θ)*cos((pi/2.0)-θ))
# n_e = 1.e6*(1.8e5*exp(-4.183119*(rE-1.0471)))
# n_p = 1.e6*(1.8e5*exp(-4.183119*(rE-1.0471)))
# ω_e2 = (n_e*(e^2))/(eps*me)
# ω_p2 = (n_p*(e^2))/(eps*mp)
# Ω_e = (e*Bmag)/(me)
# Ω_p = (e*Bmag)/(mp)
