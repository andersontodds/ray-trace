## BortnikRayTrace.jl
# Todd Anderson
# 5 October 2020
#
# Ray-tracing of electromagnetic waves from upper ionosphere to plasmasphere using Haselgrove equations, following method in Jacob Bortnik's thesis.  Includes main ODE solver function calls.
#
# Ray-tracing functions:
# a.  dr/dt   = (1/μ²)*(ρᵣ - μ*(δμ/δρᵣ))
# b.  dθ/dt   = (1/rμ²)*(ρ_θ - μ*(δμ/δρ_θ))
# c.  dϕ/dt   = (1/rμ²sinθ)*(ρ_ϕ - μ*(δμ/δρ_ϕ))
# d.  dρᵣ/dt  = (1/μ)*(δμ/δr) + ρ_θ*dθ/dt + ρ_ϕ*dϕ/dt*sinθ
# e.  dρ_θ/dt = (1/r)*((1/μ)*(δμ/δθ) - ρ_θ*dr/dt + r*ρ_ϕ*dϕ/dt*cosθ)
# f.  dρ_ϕ/dt = (1/r*sinθ)*((1/μ)*(δμ/δϕ) - ρ_ϕ*dr/dt*sinθ - r*ρ_ϕ*dθ/dt*cosθ)
#
# Assume all rays fixed to meridional plane -> dϕ/dt, ρ_ϕ, dρ_ϕ/dt = 0
#   -> ignore eqns (c) and (f), new forms for (d) and (e):
#   d2. dρᵣ/dt  = (1/μ)*(δμ/δr) + ρ_θ*dθ/dt
#   e2. dρ_θ/dt = (1/r)*((1/μ)*(δμ/δθ) - ρ_θ*dr/dt)
# Functions:
# 1. refractive_index: calculate
#
# TODO April 16 2021
# 1. Calculate ψ, μ from components ρ_r, ρ_λ, ρ_ϕ
# 	-> most likely do this in refractive_index, or in separate magnetic field/plasma environment functions
# 	-> calculation of ψ requires magnetic field direction, from dip calculation or lookup
# 	--> need three-component B field, not just dip angle w.r.t ground
# 2. Check for latitude/colatitude consistency: should be using colatitude (θ) for ray-tracing

using DifferentialEquations
using LinearAlgebra
using BenchmarkTools
using LSODA
using Sundials
using CairoMakie
#using Plots #--> replace Plots with Makie

# physical constants
c = 2.99792458e8	# speed of light in m/s
re = 6.3712e6      	# radius of Earth in meters
B0 = 3.12e-5   		# magnitude of B field at Earth's equator surface in Tesla
e = 1.602e-19 		# electric charge in Coulombs
me = 9.1093e-31		# electron rest mass
mp = 1.6726219e-27	# proton rest mass
eps = 8.854e-12 	# permittivity of free space in mks

# plasmasphere parameters
Kp_max = 3                             	# maximum Kp index
Lppi = 5.6 - 0.46*Kp_max                # plasmapause inner limit in L
d = 0                              		# day number
R̄ = 90                                  # 13-month average sunspot number
mlt = 2                             	# magnetic local time
Lppo = initialize_plasmasphere(Kp_max, Lppi, d, R̄, mlt)

b = function magnetic_field(r,θ,ϕ)
	# get magnetic field components at any location in magnetic local time
	# centered dipole field model (i.e. no correction for solar forcing or higher-order geomagnetic model terms): valid for L ≈ 2 - 5
	λ = π/2 - θ								# convert colatitude to latitude
	Br = -2*B0*(re/r)^3*sin(λ)              # radial compoenent of B-field
    Bλ = B0*(re/r)^3*cos(λ)                 # latitude component of B-field
	Bθ = -1.0*Bλ							# colatitude component of B-field
	Bϕ = 0.0
    Bmag = B0*(re/r)^3*sqrt(1 + 3*sin(λ)^2)  # magnitude of B-field == sqrt(Br^2 + Bλ^2)

	b = [Br,Bθ,Bϕ]

end

Lppo = function initialize_plasmasphere(Kp_max, Lppi, d, R̄, mlt)
	# pre-solve plasmasphere from Carpenter and Anderson 1992 model
	r_range = [re:1000:(10*re);]
	λ_range = 0
	L = r_range./(re*cos(λ_range)^2) # range of L-values where Lppo will be searched for

	Kp_max = 3              		# for example
	Lppi = 5.6 - 0.46*Kp_max        # plasmapause inner limit in L
	#Lppo = 5.75                    # placeholder
	d = 0                           # day number
	R̄ = 90							 # 13-month average sunspot number
	mlt = 2                         # magnetic local time

	ne_Lppi = 10^((-0.3145*Lppi + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2-Lppi)/1.5))

	# pre-solve individual components in order to find Lppo
	ne_plasma_1 = @. 10 .^((-0.3145.*L + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2 .-L)./1.5))
	ne_plasma_2 = ne_Lppi*10 .^(-1.0.*(L .- Lppi)./0.1)
	ne_plasma_3 = (5800 + 300*mlt).*L.^(-4.5) + (1 .- exp.((2 .-L)./10))

	min_index = findmin(abs.(ne_plasma_2 - ne_plasma_3))[2]
	Lppo = L[min_index]

end

u = function refractive_index(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
    # 1. set up medium
    # 1a. magnetic field

	Bvec = magnetic_field(r,θ,ϕ)
	Br, Bθ, Bϕ = Bvec
	Bmag = sqrt(Br^2 + Bθ^2 + Bϕ^2)

    # 1b. plasma
    # assumptions: cold plasma, 2.25 < L < 8, midnight sector, based on Carpenter and Anderson 1992
	#Lppo = initialize_plasmasphere(Kp_max, Lppi, d, R̄, mlt) # include this line if any Lppo variables change during integration (they shouldn't)

    L = r/(re*cos(π/2 - θ)^2)
	ne_iono = (1.8e5.*exp.(-4.183119.*((r./re).-1.0471))) # cm^-3

    if L <= Lppi
        log_ne = (-0.3145*L + 3.9043) + (0.15*(cos((2*π*(d+9)))/365 - 0.5*cos((4*π*(d+9)))/365) + 0.00127*R̄ - 0.0635)*exp((2-L)/1.5)
		ne_plasma = 10^(log_ne)
    elseif Lppi < L <= Lppo
        ne_plasma = ne_Lppi*10^((Lppi-L)/0.1)
    elseif Lppo < L
        ne_plasma = (5800 + 300*mlt)*L^(-4.5) + (1 - exp((2-L)/10))
    else
        ne_plasma = 0.0
    end

	n_e = (ne_plasma + ne_iono)*1e6 	# add plasmasphere + ionosphere, convert from cm^-3 to m^-3
    n_p = n_e

    # electron plasma angular frequency squared
    ω_e2 = (n_e*(e^2))/(eps*me)
	# proton plasma angular frequency squared
    ω_p2 = (n_p*(e^2))/(eps*mp)

	# electron angular gyrofrequency
    Ω_e = (e*Bmag)/(me)
	# proton angular gyrofrequency
    Ω_p = (e*Bmag)/(mp)

	# convert wave frequency to gyrofrequency
	ω = 2*π*freq

	# 2. define dispersion relation
	# 2a. calculate ψ from ρ_r, ρ_θ, ρ_ϕ

	μvec = [ρ_r, ρ_θ, ρ_ϕ]
	μmag = sqrt(ρ_r^2 + ρ_θ^2 + ρ_ϕ^2)
	cos_ψ = (Bvec ⋅ μvec)/(Bmag*μmag)
	ψ = acos(cos_ψ)

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
	u = [μ,dμdψ,ψ,Bvec]

end

# Derivatives of phase refractive index w.r.t. r, λ, ϕ, ρ_r, ρ_λ, ρ_ϕ, freq
# calculate derivative of μ w.r.t. r
dμdr = function ddr(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ,dμdψ = u			# unpack
	dr = 1.0e-11		# r step size

	μ_l = refractive_index(r-dr/2.0,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ_l = μ
	μ_r = refractive_index(r+dr/2.0,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ_r = μ

	dμdr = (μ_r[1] - μ_l[1])/dr
end

dμdθ = function ddθ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ,dμdψ = u			# unpack
	dθ = 1.0e-11		# θ step size

	μ_l = refractive_index(r,θ-dθ/2.0,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ_l = μ
	#println("μ_l = ", μ)
	μ_r = refractive_index(r,θ+dθ/2.0,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ_r = μ
	#println("μ_r = ", μ)

	dμdθ = (μ_r[1] - μ_l[1])/dθ
end

dμdϕ = function ddϕ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ,dμdψ = u			# unpack
	dϕ = 1.0e-11		# χ step size

	μ_l = refractive_index(r,θ,ϕ-dϕ/2.0,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ_l = μ
	μ_r = refractive_index(r,θ,ϕ+dϕ/2.0,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ_r = μ

	dμdϕ = (μ_r[1] - μ_l[1])/dϕ
end

dμdρr = function ddρr(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dρr = 1.0e-11		# ρ_r step size

	# Kimura 1966:
	# dμ/dρ_k = (dμ/dψ)(dψ/dρ_k) = (dμ/dψ)((ρ_k*cosψ - μcos(α_Bk)/(μ²sinψ))
	# where k = r, λ, ϕ; and α_Bk = angle between the magnetic field vector and k

	μ, dμdψ, ψ, Bvec = refractive_index(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq) #DEBUG: BoundsError at this location 04/21/2021

	ρrvec = [ρ_r, 0.0, 0.0]
	Bmag = sqrt(Bvec[1]^2 + Bvec[2]^2 + Bvec[3]^2)
	ρrmag = abs(ρ_r)

	cos_αr = (Bvec ⋅ ρrvec)/(Bmag*ρrmag)

	dμdρr = dμdψ*((ρ_r*cos(ψ) - μ*cos_αr)/(μ^2 * sin(ψ)))
end

dμdρθ = function ddρλ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dρθ = 1.0e-11		# ρ_λ step size

	# Kimura 1966:
	# dμ/dρ_k = (dμ/dψ)(dψ/dρ_k) = (dμ/dψ)((ρ_k*cosψ - μcos(α_Bk)/(μ²sinψ))
	# where k = r, λ, ϕ; and α_Bk = angle between the magnetic field vector and k
	μ, dμdψ, ψ, Bvec = refractive_index(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq) # unpack

	ρθvec = [0.0, ρ_θ, 0.0]
	Bmag = sqrt(Bvec[1]^2 + Bvec[2]^2 + Bvec[3]^2)
	ρθmag = abs(ρ_θ)

	cos_αθ = (Bvec ⋅ ρθvec)/(Bmag*ρθmag)

	dμdρθ = dμdψ*((ρ_θ*cos(ψ) - μ*cos_αθ)/(μ^2 * sin(ψ)))
end

dμdρϕ = function ddρϕ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq) #TODO λ → θ
	dρϕ = 1.0e-11		# ρ_ϕ step size

	# Kimura 1966:
	# dμ/dρ_k = (dμ/dψ)(dψ/dρ_k) = (dμ/dψ)((ρ_k*cosψ - μcos(α_Bk)/(μ²sinψ))
	# where k = r, λ, ϕ; and α_Bk = angle between the magnetic field vector and k
	μ, dμdψ, ψ, Bvec = refractive_index(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq) # unpack

	ρϕvec = [0.0, 0.0, ρ_ϕ]
	Bmag = sqrt(Bvec[1]^2 + Bvec[2]^2 + Bvec[3]^2)
	ρϕmag = abs(ρ_ϕ)

	cos_αϕ = (Bvec ⋅ ρϕvec)/(Bmag*ρϕmag)

	dμdρϕ = dμdψ*((ρ_ϕ*cos(ψ) - μ*cos_αϕ)/(μ^2 * sin(ψ)))
end

dμdf = function ddf(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	#μ,dμdψ = u			# unpack
	df = 1.0e-5			# f step size

	μ_l = refractive_index(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq-df/2.0)
	#μ_l = μ
	μ_r = refractive_index(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq+df/2.0)
	#μ_r = μ

	dμdf = (μ_r[1] - μ_l[1])/df
end


# calculate derivatives w.r.t. time
function haselgrove!(du,u,p,t)
# 	du = drdt, dθdt, dϕdt, dρrdt, dρθdt, dρϕdt
# 	u = r, θ, ϕ, ρ_r, ρ_θ, ρ_ϕ, T
#  	p = freq
# 	t = time

	r, θ, ϕ, ρ_r, ρ_θ, ρ_ϕ, T = u 	# unpack u
	freq = p[1] 						# unpack p

	dμdr = ddr(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dμdθ = ddθ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dμdϕ = ddϕ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dμdρr = ddρr(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dμdρθ = ddρλ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dμdρϕ = ddρϕ(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	dμdf = ddf(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)

	v = refractive_index(r,θ,ϕ,ρ_r,ρ_θ,ρ_ϕ,freq)
	μ = v[1]
	dμdψ = v[2]

	# TODO: ensure ρ_k have appropriate initial conditions!

    du[1] = 1/(μ^2)*(ρ_r-μ*dμdρr) 							# r: radial distance (m)
    du[2] = 1/(r*μ^2)*(ρ_θ-μ*dμdρθ)							# θ: colatitude (radians)
    du[3] = 1/(r*μ^2*sin(θ))*(ρ_ϕ-μ*dμdρϕ)					# ϕ: longitude (radians)
	du[4] = (1/μ)*dμdr + ρ_θ*du[2] + ρ_ϕ*du[3]*sin(θ)		# ρ_r: r-component of μ
	du[5] = (1/r)*((1/μ)*dμdθ - ρ_θ*du[1] + r*ρ_ϕ*du[3]*cos(θ)) # ρ_θ: θ-component of μ
	du[6] = (1/(r*sin(θ)))*((1/μ)*dμdϕ - ρ_ϕ*du[1]*sin(θ) - r*ρ_ϕ*du[2]*cos(θ)) # ρ_ϕ: ϕ-component of μ
	du[7] = (1/c)*(1 + (freq/μ)*dμdf) # T: group delay time (s)


end


## callbacks!
# re_cb: terminate when ray intersects Earth surface
function re_term_condition(u,t,integrator)
	u[1] - re
end

function terminate_affect!(integrator)
	terminate!(integrator)
end

re_cb = ContinuousCallback(re_term_condition,terminate_affect!)

# TODO: adapt saving callback to new ray-tracing equations
# #saveμ_cb: save μ, dμ values at each timestep
# function save_func(u,t,integrator)
# 	r, θ, ϕ, freq = u
# 	v = refractive_index(r, θ, χ, f)
# 	dip = atan(2.0*tan(λ))
# 	ψ = π/2.0 + dip + χ
# 	vcat(v,dip,ψ)#,r,λ)
# end

# saved_μ = SavedValues(Float64, Array{Float64,1})
# saveμ_cb = SavingCallback(save_func, saved_μ)

cb = CallbackSet(re_cb) #, saveμ_cb)

## solve ODEs and plot
u0 = [re+1.0e+6, 1.0*pi/4, 0.0, 1.0, 1.0, 0.0, 0.0]	# r0, θ0, ϕ0, ρ_r0, ρ_θ0, ρ_ϕ0, T0
p = [1000.0]										# f0
tspan = (0.0,5.0e+9)

hasel_prob = ODEProblem(haselgrove!,u0,tspan,p)
hasel_soln = solve(hasel_prob,  callback=cb)#, reltol=1e-7, dtmax = 1e6, dtmin = 1e-8)


## scratch
# plasma model:
# superposition of ionosphere and plasmasphere. ionosphere is modeled as isotropic exponential decay; plasmasphere model from Carpenter and Anderson 1992, extened to L = 1
r = [re:1000:(10*re);]
λ = 0
L = r./(re*cos(λ)^2)

Kp_max = 3                             # for example
Lppi = 5.6 - 0.46*Kp_max                # plasmapause inner limit in L
#Lppo = 5.75                             # placeholder
d = 0                                   # day number
R̄ = 90                                  # 13-month average sunspot number
t = 2                                   # magnetic local time

ne_iono = (1.8e5.*exp.(-4.183119.*((r./re).-1.0471))) # cm^-3

ne_Lppi = 10^((-0.3145*Lppi + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2-Lppi)/1.5))

# individual components
ne_plasma_1 = @. 10 .^((-0.3145.*L + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2 .-L)./1.5))
ne_plasma_2 = ne_Lppi*10 .^(-1.0.*(L .- Lppi)./0.1)
ne_plasma_3 = (5800 + 300*t).*L.^(-4.5) + (1 .- exp.((2 .-L)./10))

plot(L, ne_plasma_1, yaxis=:log, label="plasmasphere")
plot!(L, ne_plasma_2, yaxis=:log, label="plasmapause inner limit")
plot!(L, ne_plasma_3, yaxis=:log, label="plasmasphere outer limit")
plot!(L,ne_iono, yaxis=:log, label="ionosphere")

f = Figure()
ax = Axis(f[1,1], xlabel = "x", ylabel = "f(x)", yscale = log10)
ylims!(10^-1, 10^4)
lines!(L, ne_plasma_1, label="plasmasphere")
lines!(L, ne_plasma_2, label="plasmapause inner limit")
lines!(L, ne_plasma_3, label="plasmapause outer limit")
f


min_index = findmin(abs.(ne_plasma_2 - ne_plasma_3))[2]
Lppo = L[min_index]

ne_plasma = ones(length(L))

for i = 1:length(L)
    if L[i] <= Lppi
        log_ne = (-0.3145*L[i] + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2-L[i])/1.5)
        ne_plasma[i] = 10^(log_ne)
    elseif Lppi < L[i] <= Lppo
        ne_plasma[i] = ne_Lppi*10^((Lppi-L[i])/0.1)
    elseif Lppo < L[i]
        ne_plasma[i] = (5800 + 300*t)*L[i]^(-4.5) + (1 - exp((2-L[i])/10))
    else
        ne_plasma[i] = 0.0
    end
end

ne_plas = function plasmasphere(L,Lppi,Lppo)
	if L <= Lppi
		log_ne = (-0.3145*L + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2-L)/1.5)
        ne_plasma = 10^(log_ne)
	elseif Lppi < L <= Lppo
        ne_plasma = ne_Lppi*10^((Lppi-L)/0.1)
    elseif Lppo < L
        ne_plasma = (5800 + 300*t)*L^(-4.5) + (1 - exp((2-L)/10))
    else
        ne_plasma = 0.0
    end
end

ne_ion(r) = (1.8e5.*exp.(-4.183119.*((r).-1.0471))) # cm^-3; where r is in units of Re

ne_tot(x,y) = ne_ion(r_xy) + ne_plas(L_xy)

#using Plots
#plot(L, ne_iono, yaxis=:log, ylims=[1,10^6], label="ionosphere")
#plot!(L, ne_plasma, yaxis=:log, label="plasmasphere")
#plot!(L, ne_total, yaxis=:log, label="aggregate n_e")

f2 = Figure()
ax = Axis(f2[1,1], xlabel = "L(Rₑ)", ylabel = "n (m⁻³)", yscale = log10)
ylims!(10^-1, 10^5)
lines!(L, ne_iono)
lines!(L, ne_plasma)
lines!(L, ne_total)
f2

# plot magnetic field lines, strength
#using CairoMakie

let
    x = -1:0.02:1
    y = -1.5:0.02:2
    egg(x,y) = x^2 + y^2/(1.4 + y/5)^2
    segg = [egg(x,y) for x in x, y in y]
    fig = Figure(resolution = (470, 550))
    ax = Axis(fig, xlabel = "x", ylabel = "y", backgroundcolor = :black,
    xgridstyle=:dash, ygridstyle=:dash, xgridcolor = :grey, ygridcolor = :grey)
    cl =contour!(x, y, segg, linewidth = 0.85,colormap = :viridis,
                levels = 0:0.02:0.5)
    cbar = Colorbar(fig, cl, label ="egg-l", labelpadding = 0, width = 15,
                ticksize=15, tickalign = 1, height = Relative(1))
    fig[1, 1] = ax
    fig[1, 2] = cbar
    colgap!(fig.layout, 7)
    fig
end

let 
	x = -2:0.01:2
	y = -2:0.01:2
	gridrange = -pi:0.01:0
	r_xy(x,y) = (x^2 + y^2)^(1/2);
	θ_xy(x,y) = atan(y/x);
	
	L_xy(x,y) = r_xy(x,y)/cos(θ_xy(x,y))^2
	sL = [L_xy(x,y) for x in x, y in y]

	Bmagnitude(x,y) = B0*(1/r_xy(x,y))^3*(1 + 3*cos(θ_xy(x,y))*cos(θ_xy(x,y)))
	log_Bmag(x,y) = log10(Bmagnitude(x,y))
	sB = [Bmagnitude(x,y) for x in x, y in y]
	slB = [log_Bmag(x,y) for x in x, y in y]

	fig = Figure()
	ax = Axis(fig, xlabel = "Lₓ", ylabel = "Ly", backgroundcolor = :black, 
		xgridstyle = :dash, ygridstyle = :dash, xgridcolor = :grey, ygridcolor = :grey,
		aspect = DataAspect())
	#c1 = contour!(x, y, sB, linewidth = 0.85, colormap = :viridis, levels = 10^-6:5*10^-7:2*10^-4)
	#cbar1 = Colorbar(fig, c1, label = "|B| (T)", labelpadding = 0, width = 15, 
	#	ticksize = 15, tickalign = 1, height = Relative(1))
	
	c1 = heatmap!(x, y, sB, colormap = :viridis, colorrange = (10^-6,10^-4))
	cbar1 = Colorbar(fig, c1, label = "|B| (T)", labelpadding = 0, width = 15, 
		ticksize = 15, tickalign = 1, height = Relative(1), scale = log10)
		

	c2 = contour!(x, y, sL, linewidth = 0.85, color = :red, levels = 1:0.5:6)
	#cbar2 = Colorbar(fig, c1, label = "B magnitude", labelpadding = 0, width = 15, 
	#	ticksize = 15, tickalign = 1, height = Relative(1))
	
	
	poly!(Circle(Point2f(0,0), 1f0), color = :black)
	fig[1,1] = ax
	fig[1,2] = cbar1
	colgap!(fig.layout, 7)
	fig
end

# next: plot B field lines, heatmap plasma density
let 
	x = -4:0.01:4
	y = -4:0.01:4
	r_xy(x,y) = (x^2 + y^2)^(1/2);
	θ_xy(x,y) = atan(y/x);
	
	# B-field lines: aka contours of constant L
	L_xy(x,y) = r_xy(x,y)/cos(θ_xy(x,y))^2
	sL = [L_xy(x,y) for x in x, y in y]

	# plasma density: ne_total is a combination of ne_iono(r) and ne_plasma(L).  
	# Additionally, ne_plasma is piecewise based on plasmasphere location (Lppi, Lppo).
	ne_tot(x,y) = ne_ion(r_xy(x,y)) + ne_plas(L_xy(x,y),Lppi, Lppo)
	#ne_tot(x,y) = ne_plas(L_xy(x,y),Lppi, Lppo)
	sn = [ne_tot(x,y) for x in x, y in y]

	fig = Figure()
	ax = Axis(fig, xlabel = "Lₓ", ylabel = "Ly", backgroundcolor = :black, 
		xgridstyle = :dash, ygridstyle = :dash, xgridcolor = :grey, ygridcolor = :grey,
		aspect = DataAspect())
	
	cmap = cgrad(:magma, scale = :log10)
	c1 = heatmap!(x, y, sn, colormap = cmap, colorrange = (10^2,10^3))
	cbar1 = Colorbar(fig, c1, label = "n (cm⁻³)", labelpadding = 0, width = 15, 
		ticksize = 15, tickalign = 1, height = Relative(1), scale = log10)
		

	c2 = contour!(x, y, sL, linewidth = 0.85, color = :red, levels = 1:0.5:6)
	#cbar2 = Colorbar(fig, c1, label = "B magnitude", labelpadding = 0, width = 15, 
	#	ticksize = 15, tickalign = 1, height = Relative(1))
	
	
	poly!(Circle(Point2f(0,0), 1f0), color = :black)
	fig[1,1] = ax
	fig[1,2] = cbar1
	colgap!(fig.layout, 7)
	fig
end


##
r = re*2
θ = π/4
ϕ = 0.0
Bvec = magnetic_field(r,θ,ϕ)
Br, Bθ, Bϕ = Bvec
Bmag = sqrt(Br^2 + Bθ^2 + Bϕ^2)

# 1b. plasma
# assumptions: cold plasma, 2.25 < L < 8, midnight sector, based on Carpenter and Anderson 1992
#Lppo = initialize_plasmasphere(Kp_max, Lppi, d, R̄, mlt) # include this line if any Lppo variables change during integration (they shouldn't)

L = r/(re*cos(π/2 - θ)^2)
ne_iono = (1.8e5.*exp.(-4.183119.*((r./re).-1.0471))) # cm^-3

if L <= Lppi
	log_ne = (-0.3145*L + 3.9043) + (0.15*(cos((2*π*(d+9)))/365 - 0.5*cos((4*π*(d+9)))/365) + 0.00127*R̄ - 0.0635)*exp((2-L)/1.5)
	ne_plasma = 10^(log_ne)
elseif Lppi < L <= Lppo
	ne_plasma = ne_Lppi*10^((Lppi-L)/0.1)
elseif Lppo < L
	ne_plasma = (5800 + 300*mlt)*L^(-4.5) + (1 - exp((2-L)/10))
else
	ne_plasma = 0.0
end

n_e = (ne_plasma + ne_iono)*1e6 	# add plasmasphere + ionosphere, convert from cm^-3 to m^-3
n_p = n_e

# electron plasma angular frequency squared
ω_e2 = (n_e*(e^2))/(eps*me)
# proton plasma angular frequency squared
ω_p2 = (n_p*(e^2))/(eps*mp)

# electron angular gyrofrequency
Ω_e = (e*Bmag)/(me)
# proton angular gyrofrequency
Ω_p = (e*Bmag)/(mp)

# convert wave frequency to gyrofrequency
ω = 2*π*freq

# 2. define dispersion relation
# 2a. calculate ψ from ρ_r, ρ_θ, ρ_ϕ
ρ_r = 1.0
ρ_θ = 1.0
ρ_ϕ = 0.0
μvec = [ρ_r, ρ_θ, ρ_ϕ]
μmag = sqrt(ρ_r^2 + ρ_θ^2 + ρ_ϕ^2)
cos_ψ = (B ⋅ μvec)/(Bmag*μmag)
ψ = acos(cos_ψ)

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
u = [μ,dμdψ,ψ,Bvec]

# try unpacking u

μ1, dμdψ1, ψ1, B1 = u

B1[1]
