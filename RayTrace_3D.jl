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
# 	->
# 2. Check for latitude/colatitude consistency: should be using colatitude (θ) for ray-tracing

using DifferentialEquations
using BenchmarkTools
using LSODA
using Sundials
using Plots

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
mlt = 2                             		# magnetic local time

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

end

u = function refractive_index(r,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
    # 1. set up medium
    # 1a. magnetic field

	b = magnetic_field(r,θ,ϕ)
	Br, Bθ, Bϕ = b
	Bmag = sqrt(Br^2 + Bθ^2 + Bϕ^2)

    # 1b. plasma
    # assumptions: cold plasma, 2.25 < L < 8, midnight sector, based on Carpenter and Anderson 1992

	Lppo = initialize_plasmasphere(Kp_max, Lppi, d, R̄, mlt)

    L = r/(re*cos(λ)^2) #TODO λ → θ

    if L <= Lppi
        log_ne = (-0.3145*L + 3.9043) + [0.15*(cos((2*π*(d+9)))/365 - 0.5*cos((4*π*(d+9)))/365) + 0.00127*R̄ - 0.0635]*exp((2-L)/1.5)
    elseif Lppi < L <= Lppo
        ne = ne_Lppi*10^((Lppi-L)/0.1)
    elseif Lppo < L
        ne = (5800 + 300*mlt)*L^(-4.5) + (1 - exp((2-L)/10))
    else
        ne = 0.0
    end

    np = ne

    # electron plasma angular frequency squared
    ω_e2 = (n_e*(e^2))/(eps*me)
	# proton plasma angular frequency squared
    ω_p2 = (n_p*(e^2))/(eps*mp)

	# electron angular gyrofrequency
    Ω_e = (e*Bmag)/(me)
	# proton angular gyrofrequency
    Ω_p = (e*Bmag)/(mp)

	# 2. define dispersion relation
	# 2a. calculate ψ from ρ_r, ρ_λ, ρ_ϕ
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
	u = [μ,dμdψ]

end

# Derivatives of phase refractive index w.r.t. r, θ, χ,
# calculate derivative of μ w.r.t. r
dμdr = function ddr(r,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ,dμdψ = u			# unpack
	dr = 1.0e-11		# r step size

	μ_l = refractive_index(r-dr/2.0,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ_l = μ
	μ_r = refractive_index(r+dr/2.0,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ_r = μ

	dμdr = (μ_r[1] - μ_l[1])/dr
end

dμdλ = function ddλ(r,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ,dμdψ = u			# unpack
	dλ = 1.0e-11		# θ step size #TODO λ → θ

	μ_l = refractive_index(r,λ-dλ/2.0,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ_l = μ
	#println("μ_l = ", μ)
	μ_r = refractive_index(r,λ+dλ/2.0,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ_r = μ
	#println("μ_r = ", μ)

	dμdλ = (μ_r[1] - μ_l[1])/dλ #TODO λ → θ
end

dμdϕ = function ddϕ(r,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ,dμdψ = u			# unpack
	dϕ = 1.0e-11		# χ step size

	μ_l = refractive_index(r,λ,ϕ-dϕ/2.0,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ_l = μ
	μ_r = refractive_index(r,λ,ϕ+dϕ/2.0,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ_r = μ

	dμdϕ = (μ_r[1] - μ_l[1])/dϕ
end

dμdρr = function ddρr(r,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	#μ,dμdψ = u			# unpack
	dρr = 1.0e-11		# χ step size

	# Kimura 1966:
	# dμ/dρ_k = (dμ/dψ)(dψ/dρ_k) = (dμ/dψ)((ρ_k*cosψ - μcos(α_Bk)/(μ²sinψ))
	# where k = r, λ, ϕ; and α_Bk = angle between the magnetic field vector and k

	# μ_l = refractive_index(r,λ,ϕ,ρ_r-dρr/2.0,ρ_λ,ρ_ϕ,freq)
	# #μ_l = μ
	# μ_r = refractive_index(r,λ,ϕ,ρ_r+dρr/2.0,ρ_λ,ρ_ϕ,freq)
	# #μ_r = μ

	dμdρr = (μ_r[1] - μ_l[1])/dρr
end


# calculate derivatives w.r.t. time
function haselgrove!(du,u,p,t)
# 	du = drdt, dλdt, dϕdt, dρrdt, dρλdt, dρϕdt
# 	u = r, λ, ϕ, ρ_r, ρ_λ, ρ_ϕ
#  	p = freq
# 	t = time

	r, λ, ϕ, ρ_r, ρ_λ, ρ_ϕ = u 		# unpack u #TODO λ → θ
	freq = p 						# unpack p

	dμdr = ddr(r,λ,ϕ,freq) #TODO λ → θ


	v = refractive_index(r,λ,ϕ,ρ_r,ρ_λ,ρ_ϕ,freq) #TODO λ → θ
	μ = v[1]
	dμdψ = v[2]

	#TODO finish changing these to Kimura 1966 eqns (Bortnik thesis 2.3, group time 2.4)
    du[1] = 1/(μ^2)*(ρ_r - μ*dμdρr)
    du[2] = 1/(r*μ^2)*(μ*sin(χ)-dμdψ*cos(χ))
    du[3] = 1/(r*μ^2)*(dμdλ*cos(χ)-(r*dμdr + μ)*sin(χ))
	du[4] = 1/c*(1+(freq/μ)*dμdf)


end

## scratch
# plasma model:
# superposition of ionosphere and plasmasphere. ionosphere is modeled as isotropic exponential decay; plasmasphere model from Carpenter and Anderson 1992, extened to L = 1
r = [re:1000:(10*re);]
λ = 0
L = r./(re*cos(λ)^2)

Kp_max = 3                             # for example
Lppi = 5.6 - 0.46*Kp_max                # plasmapause inner limit in L
Lppo = 5.75                             # placeholder
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

ne_total = ne_iono + ne_plasma

using Plots
plot(L, ne_iono, yaxis=:log, ylims=[1,10^6], label="ionosphere")
plot!(L, ne_plasma, yaxis=:log, label="plasmasphere")
plot!(L, ne_total, yaxis=:log, label="aggregate n_e")
