# ray-trace
Julia implementation of electromagnetic plasma wave ray-tracing from terrestrial sources into the magnetosphere.

This repository is in its early stages!  Header documentation in individual files may not match the actual workings of the files (see RayMain.jl, for example).  Documentation in individual files will be revised before the v1.0 release.

The current version (v0.13) uses only RayMain.jl (plus dependencies outside this repo).  The eventual plan the structure of this package is:

 * ODE solver calls, ray-tracing function calls, and plotting, will be handled in RayMain.jl.
 * the "environment" variables will either be outsourced from an existing Julia geospace package, or I will start developing such a package (Geospace.jl).  This includes physical constants, as well as various particle and magnetic field models.
 * ray-tracing equations will be contained in appropriate function modules.  For now, the Haselgrove equations will be contained in Haselgrove.jl; the phase refractive index function, and dx/dt equations, will be contained in RayFunctions.jl
 * Inputs, i.e. a list of (r,θ,χ,f,t) that are used by the ODE solver, will be in ray_start.dat
 * If needed, solver parameters will be contained in RayParam.jl

#RayMain.jl
This module traces electron whistler waves from the upper ionosphere to the plasmasphere.  The ray-tracing equations used are the Haselgrove equations.

#scratch.jl
This is a scratch pad for the rest of the project.
