## plasmasphere.jl
# Todd Anderson
# 11 January 2022
#
# Specify plasma density profile.
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

# physical constants
c = 2.99792458e8	# speed of light in m/s
re = 6.3712e6      	# radius of Earth in meters
B0 = 3.12e-5   		# magnitude of B field at Earth's equator surface in Tesla
e = 1.602e-19 		# electric charge in Coulombs
me = 9.1093e-31		# electron rest mass in kg
mp = 1.6726219e-27	# proton rest mass in kg
eps = 8.854e-12 	# permittivity of free space in mks

# plasmasphere parameters
#Kp_max = 3                             	# maximum Kp index in preceding 24 hours
#Lppi = 5.6 - 0.46*Kp_max                # plasmapause inner limit in L
#d = 0                              		# day number
#R̄ = 90                                  # 13-month average sunspot number
#mlt = 2                             	# magnetic local time

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

ne_plasmasphere = function plasmasphere(L, Lppi, Lppo, d, R̄, mlt)
    # Calculates plasma density at a location in geospace using Carpenter and Anderson 1992 model.
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
