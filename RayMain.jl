## RayMain.jl
# Todd Anderson
# 5 October 2020
#
# Ray-tracing of electromagnetic waves from upper ionosphere to plasmasphere using Haselgrove equations, following method in Rice 1997.  Includes main ODE solver function calls.
#
# Dependencies:
#   Haselgrove.jl:      Ray-tracing equations derived by Haselgrove (1957)
#   RTConstants.jl:     Global constants
#   RTParameters.jl:    Parameter choices for ODE solver
#   RTPlot.jl:          Plotting function calls and options

#= Outline of Python implementation:
    1. Import packages
    2. Read input parameters and initial conditions
        includes parameters in ray_param.dat
        a. load parameters
        b. load constants
        c. bundle initial conditions
        d. select and initialize ionosphere
    3. Calculate RHS of ray-tracing equations
        a. f(t, rtcG, freq)
            rtcG: current values of r, θ, χ, G
        b. calculate μ with call of phase_refractive_index(r,θ,χ,freq,ion[])
        c. calculate dμ/dx, where x = r, θ, χ, freq
            i) note that dμ/dχ = -1*(dμ/dψ)
        d. calculate derivatives w.r.t:
            derivs = [dr/dt, dθ/dt, dχ/dt, dfreq/dt]
    4. ODE solver calls
        a. select solver (e.g. LSODA)
        b. set plot parameters -- see ray_plot.py
        c. loop over initial frequencies to calculate ray paths:
            i.  determine number of evenly-spaced frequencies
            ii. initialize freq, group delay time, nphase arrays
            iii.for loop:
                1. set frequency to new value in range
                2. call ode solver: soln = solve f (see 3a)
                    a. set integration method
                    b. set initial values of radius, latitude, gdt, nphase
                3. while: solver is successful & time < tstop & 1Re < altitude < max alt
                    a.  soln.integrate(soln.t+dt)
                    b.  radius.append(soln.y[0])
                    c.  latitude.apppend(soln.y[1])
                    d.  gdt.append(soln.y[3])
                    e.  nphase_single =
                        phase_refractive_index(radius,latitude,chi,freq,[ion])  [from soln]
                    f.  nphase.append(nphase_single[0]) [= μ; nphase_single[1] is dμ/dψ]
                4. calculate x = radius[:]*cos(latitude[:]); y = radius[:]*sin(latitude[:])
                5. append freq, gdt, nphase to arrays
                    QUESTION: this program uses the initialization of empty (unsized) arrays/lists and subsequent appending to fill them in, resulting in unknown list size and differing size between list elements (i.e. different frequency rays).  How slow is this compared to a pre-allocated implementation?
                    TODO: optimize this ODE solve process.
                6.  plot x,y for this frequency

    5. Show plots

=#
## Read parameters, constants, and initial conditions




## Ray-tracing equations calculation

#=TODO 10/6/2020:
	1. Get functions and argument-passing working with single input (i.e. non-vectorized, non-looping)
	2. Vectorize!
	3. Run and time solver
	4. Switch to in-place functions where possible
	5. Run and time solver, seek out optimization
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
u = function phase_refractive_index!(dμdψ,μ,r,θ,χ,freq)

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
	R = 1.0 - (ω_e2/ω^2)*(ω/(ω - Ω_e)) - (ω_p2/ω^2)*(ω/(ω + Ω_p))

	# L ≡ 1 - Σ(ωₖ²/ω²)(ω/(ω-ϵₖΩₖ))
	L = 1.0 - (ω_e2/ω^2)*(ω/(ω + Ω_e)) - (ω_p2/ω^2)*(ω/(ω - Ω_p))

	# P ≡ 1 - Σ(ωₖ²/ω²)
	P = 1.0 - (ω_e2/ω^2) - (ω_p2/ω^2)

	# D ≡ ½(R-L); S ≡ ½(R+L)
	D = (R - L)/2
	S = (R + L)/2

	# define parts of dispersion relation
	# Aμ⁴ - Bμ² + C = 0
	# A = Ssin²ψ + Pcos²ψ
	A = S*sin(ψ)^2 + P*cos(ψ)^2

	# B = RLsin²ψ + PS(1 + cos²ψ)
	B = R*L*sin(ψ)^2 + P*S*(1+cos(ψ)^2)

	# C  = PRL
	C = P*R*L

	# solve the dispersion relation for μ:
	# μ² = (B +- F)/2A
	# where F² = (RL-PS)²sin⁴ψ + 4P²D²cos²ψ
	#QUESTION: is it faster to (a) solve for F², then sqrt(); or (b) solve for F directly? (a) requires fewer sqrt() calls
	F2 = (R*L - P*S)^2*sin(ψ)^4 + 4*(P*D*cos(ψ))^2
	F = sqrt(F2)

	μ2_minus = (B - F)/(2*A)
	μ_minus = sqrt(μ2_minus)
	μ2_plus = (B + F)/(2*A)
	μ_plus = sqrt(μ2_plus)

	# Electron whistler case: ψ = 0, μ² = R; this is the μ_minus case
	μ = μ_minus

	# Find dA/dψ, dB/dψ, dC/dψ, dμ/dψ
	dAdψ = 2.0*(S-P)*sin(ψ)*cos(ψ)
	dBdψ = 2.0*(R*L-P*S)*sin(ψ)*cos(ψ)
    dCdψ = 0.0
    dμdψ = ((μ^4)*dAdψ-(μ^2)*dBdψ+dCdψ)/(4.0*A*(μ^3)-2.0*B*μ)

	#DEBUG check values
	print("μ = ",μ,"\n")
	print("dμdψ = ", dμdψ)

	u = [μ,dμdψ]

end

# Derivatives of phase refractive index w.r.t. r, θ, χ,
# calculate derivative of μ w.r.t. r
dμdr = function ddr(μ,r,θ,χ,freq)
	dr = 1.0e-11		# r step size

	μ_l = phase_refractive_index!(dμdψ,μ,r-dr/2.0,θ,χ,freq)
	μ_r = phase_refractive_index!(dμdψ,μ,r+dr/2.0,θ,χ,freq)

	dμdr = (μ_r[1] - μ_l[1])/dr
end

dμdθ = function ddθ(μ,r,θ,χ,freq)
	dθ = 1.0e-11		# θ step size

	μ_l = phase_refractive_index!(dμdψ,μ,r,θ-dθ/2.0,χ,freq)
	μ_r = phase_refractive_index!(dμdψ,μ,r,θ+dθ/2.0,χ,freq)

	dμdθ = (μ_r[1] - μ_l[1])/dθ
end

dμdχ = function ddχ(μ,r,θ,χ,freq)
	dχ = 1.0e-11		# χ step size

	μ_l = phase_refractive_index!(dμdψ,μ,r,θ,χ-dχ/2.0,freq)
	μ_r = phase_refractive_index!(dμdψ,μ,r,θ,χ-dχ/2.0,freq)

	dμdχ = (μ_r[1] - μ_l[1])/dχ
end

dμdf = function ddf(μ,r,θ,χ,freq)
	df = 1.0e-5		# f step size

	μ_l = phase_refractive_index!(dμdψ,μ,r,θ,χ,freq-df/2.0)
	μ_r = phase_refractive_index!(dμdψ,μ,r,θ,χ,freq+df/2.0)

	dμdf = (μ_r[1] - μ_l[1])/df
end

# calculate derivatives w.r.t. time
function haselgrove!(dμ,μ,r,θ,χ,freq,t)

	dμdr,dμdθ,dμdχ,dμdf = dμ		# unpack dμ

	dμdr = ddr(μ,r,θ,χ,freq)
	dμdθ = ddθ(μ,r,θ,χ,freq)
	dμdχ = ddχ(μ,r,θ,χ,freq)
	dμdf = ddf(μ,r,θ,χ,freq)


    drdt = 1/(μ^2)*(μ*cos(χ)+dμdχ*sin(χ))
    dθdt = 1/(r*μ^2)*(μ*sin(χ)-dμdχ*cos(χ))
    dχdt = 1/(r*μ^2)*(dμdθ*cos(χ)-(r*dμdr + μ)*cos(χ))

	#DEBUG check values
	print("drdt = ",drdt,"\n")
	print("drdθ = ",drdθ,"\n")
	print("drdχ = ",drdχ)


    # dμdψ = 1/(2*μ)*((dBdψ + dFdψ)/(2*A) - 2*dAdψ*(B+F)/(2*A^2))
    # dμdω = 1/(2*μ)*((dBdω + dFdω)/(2*A) - 2*dAdω*(B+F)/(2*A^2))
	#
    # dμdr = (μ*(r+δr/2) - μ*(r-δr/2))/δr
    # dμdθ = (μ*(θ+δθ/2) - μ*(θ-δθ/2))/δθ



end

#DEBUG check values
haselgrove!([1,1,1,1],1,10000,pi/4,0,5000,0)


## ODE solver calls



## Plot results
