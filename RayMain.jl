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
	0. Rewrite particle model for realistic plasmasphere! See Sousa thesis p. 39 for plots of electron density models.
		a. see GCPM at plasmasphere.nasa.gov/models/
		b. the simplified GCPM looks like it could be approximated analytically as an exponential(?) decay multiplied by a dipole field equation
		i.e. n_e = a*exp(-b*r)*c*sqrt(1 + 3cos²θ)
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
    dip = atan(2.0*cot(θ)) # dip angle: angle between horizontal and B field;

#	ψ = 2*π - (3π/2 - dip - χ) 		# this should be ≡ to the below, but the plotted ray paths are not all that similar!
	ψ = π/2.0 + dip + χ
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
#haselgrove!(1.0,1.0,10000.0,pi/4.0,0.0,5000.0,0.0,[1.0,1.0,1.0,1.0],[1.0,1.0,1.0])


## ODE solver calls



## Plot results
using BenchmarkTools
using LSODA
using Sundials

u0 = [re+1.0e+6, 1.0*pi/4, 0.0, 5000.0]					# r0, θ0, χ0
p = []	# f0, dμdψ, dμdr, dμdθ, dμdχ, dμdf
tspan = (0.0,5.0e+9)

hasel_prob = ODEProblem(haselgrove!,u0,tspan,p)
hasel_soln = solve(hasel_prob, CVODE_BDF(), reltol=1e-7)

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

plot(re.*cos.([0:0.01:2*π;]),re.*sin.([0:0.01:2*π;]), aspect_ratio = 1, legend=:none)
plot!(x,y, aspect_ratio=1)
