## plasmasphere.jl
# Todd Anderson
# 11 January 2022
#
# Specify plasma density profile in the Earth's plasmasphere.
# Model options:
#   1.  Carpenter and Anderson 1992
#       from https://doi.org/10.1029/91JA01548
#       Empirical model based on parameters:
#           Kp_max: maximum Kp index in preceding 24 hours; also MLT dependant
#           Lppi:   plasmapause inner limit in L
#           Lppo:   plasmapause outer limit in L (determined from initializaton function)
#           d:      day number
#           R̄:      13-month average sunspot number
#           mlt:    magnetic local time
#       see also Ngo model (e.g. Bortnik thesis)
#   2.  Global Core Plasmasphere Model (GCPM)
#       from https://doi.org/10.1029/1999JA000241
#       TODO: add GCPM, Sousa simplified GCPM (Sousa thesis? other source?)
#
# TODO:
#   0. incorporate diffusive_equilibrium() !!!
#   1. clean up functions for inclusion in Geospace.jl
#       - argument defaults: possibly include physical constants in own file rather than at start of this file
#   2. separate ionosphere into ionosphere.jl
#       - work on ionosphere profile fit, see note in ionosphere_eq()
#   3. separate test lines into file that calls plasmasphere

# physical constants
c = 2.99792458e8	# speed of light in m/s
re = 6.3712e6      	# radius of Earth in meters
B0 = 3.12e-5   		# magnitude of B field at Earth's equator surface in Tesla
e = 1.602e-19 		# electric charge in Coulombs
me = 9.1093e-31		# electron rest mass in kg
mp = 1.6726219e-27	# proton rest mass in kg
eps = 8.854e-12 	# permittivity of free space in mks

# plasmasphere parameters
Kp_max = 3                             	# maximum Kp index in preceding 24 hours
Lppi = 5.6 - 0.46*Kp_max                # plasmapause inner limit in L
d = 0                              		# day number
R̄ = 90                                  # 13-month average sunspot number
mlt = 2                             	# magnetic local time

Lppo = function initialize_plasmasphere(Lppi, d, R̄, mlt)
	# pre-solve plasmasphere from Carpenter and Anderson 1992 model
	r_range = [re:1000:(10*re);]
	λ_range = 0
	L = r_range./(re*cos(λ_range)^2) # range of L-values where Lppo will be searched for

	#Kp_max = 3              		 # for example
	#Lppi = 5.6 - 0.46*Kp_max        # plasmapause inner limit in L
	##Lppo = 5.75                    # placeholder
	#d = 0                           # day number
	#R̄ = 90						  # 13-month average sunspot number
	#mlt = 2                         # magnetic local time

	ne_Lppi = 10^((-0.3145*Lppi + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2-Lppi)/1.5))

	# pre-solve individual components in order to find Lppo
	ne_plasma_1 = @. 10 .^((-0.3145.*L + 3.9043) + (0.15*(cos((2*π*(d+9))/365) - 0.5*cos((4*π*(d+9))/365)) + 0.00127*R̄ - 0.0635)*exp((2 .-L)./1.5))
	ne_plasma_2 = ne_Lppi*10 .^(-1.0.*(L .- Lppi)./0.1)
	ne_plasma_3 = (5800 + 300*mlt).*L.^(-4.5) + (1 .- exp.((2 .-L)./10))

	min_index = findmin(abs.(ne_plasma_2 - ne_plasma_3))[2]
	Lppo = L[min_index]

end

ne_plasma_eq = function plasmasphere_eq(L, Lppi, Lppo, d, R̄, mlt)
    # Calculates plasma density at a location in the equatorial plane using Carpenter and Anderson 1992 model.
    # Inputs:
    #   L: L-shell parameter, i.e. distance between Earth center and where magnetic field line at location crosses equator
    #   Lppi: L-value of plasmapause inner limit
    #   d: day number
    #   R̄: 13-month average sunspot number
    #   mlt: magnetic local time at location
    
    #Lppo = initialize_plasmasphere(Kp_max, Lppi, d, R̄, mlt)

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

ne_plasma_de = function diffusive_equilibrium(r, ne_eq)
    
    rb = 7370000 # geocentric distance to base of diffusive equilibrium model, in m
    G = rb*(1-rb/r)
    S = 1.506*T*(rb/7370)^2*(1/4^(i-1)) # parameter i specifies ion species
    αᵢ = 1 # relative concentration of each 

    ne_plasma_de = ne_eq*(αᵢ*exp(-G/S)) # for single ion species
end

ne_iono = function ionosphere_eq(r)
    # currently this is an eyeballed fit to dayside ionosphere profile up to 
    # 2000km in Sousa dissertation p32.  Next steps:
    #   1. More careful fits to dayside and nightside
    #   2. Interpolate between these based on input MLT
    #   3. Alternatively, find simplified IRI to use (which is where this comes from) 
    #ne_iono = (1.8e5*exp(-4.183119*(((r/re)-1.0471))))
    ne_iono = (1e5*exp(-10*(((r/re)-1.0471))))
    
end


# test

L_rλ(r,λ) = r./(re*cos(λ)^2)

Lppo = initialize_plasmasphere(Lppi, d, R̄, mlt)
r_test = [re:1000:10*re;]
λ_test = zeros(size(r_test))
L_test = L_rλ.(r_test,λ_test)
ne_plasma_eq = plasmasphere_eq.(L_test,Lppi,Lppo,d,R̄,mlt)
ne_iono = ionosphere_eq.(r_test)
ne_total = ne_plasma_eq .+ ne_iono

ionosphere_eq(1000000 + re)

let
    fig = Figure()
	ax = Axis(fig, xlabel = "Lₓ", ylabel = "nₑ", backgroundcolor = :black, 
		xgridstyle = :dash, ygridstyle = :dash, xgridcolor = :grey, ygridcolor = :grey,
		yscale = log10, aspect = AxisAspect(1))
    xlims!(1,7)
    ylims!(10^-1, 10^6)
	#c1 = contour!(x, y, sB, linewidth = 0.85, colormap = :viridis, levels = 10^-6:5*10^-7:2*10^-4)
	#cbar1 = Colorbar(fig, c1, label = "|B| (T)", labelpadding = 0, width = 15, 
	#	ticksize = 15, tickalign = 1, height = Relative(1))
	
	#c1 = heatmap!(x, y, sB, colormap = :viridis, colorrange = (10^-6,10^-4))
	#cbar1 = Colorbar(fig, c1, label = "|B| (T)", labelpadding = 0, width = 15, 
	#	ticksize = 15, tickalign = 1, height = Relative(1), scale = log10)
    
    lines!(L_test, ne_iono)
    lines!(L_test, ne_plasma_eq)
    lines!(L_test, ne_total)

	fig[1,1] = ax
	fig
end

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

	log_ne_tot(x,y) = log10(ne_tot(x,y))
	sln = [log_ne_tot(x,y) for x in x, y in y]

	fig = Figure()
	ax = Axis(fig, xlabel = "Lₓ", ylabel = "Ly", backgroundcolor = :black, 
		xgridstyle = :dash, ygridstyle = :dash, xgridcolor = :grey, ygridcolor = :grey,
		aspect = DataAspect())
	
	#cmap = cgrad(:magma, scale = :log10)
	#c1 = heatmap!(x, y, sn, colormap = cmap, colorrange = (10^2,10^3))
	c1 = heatmap!(x, y, sln, colormap = :magma, colorrange = (0,5))
	cbar1 = Colorbar(fig, c1, label = "log₁₀n (cm⁻³)", labelpadding = 0, width = 15, 
		ticksize = 15, tickalign = 1, height = Relative(1))
		

	c2 = contour!(x, y, sL, linewidth = 0.85, colormap = :bilbao, levels = 1:0.5:6)
	#cbar2 = Colorbar(fig, c1, label = "B magnitude", labelpadding = 0, width = 15, 
	#	ticksize = 15, tickalign = 1, height = Relative(1))
	
	
	poly!(Circle(Point2f(0,0), 1f0), color = :black)
	fig[1,1] = ax
	fig[1,2] = cbar1
	colgap!(fig.layout, 7)
	fig
end