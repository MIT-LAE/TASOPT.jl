"""
`aerodynamics` is a module that contains all aerodynamic calculations
required in the aircraft sizing. Only functions used in other routines
are exported other functions are kept within the `aerodynamics` namespace.
"""
module aerodynamics

using StaticArrays
using LinearAlgebra
using ..atmosphere
import ..TASOPT: __TASOPTindices__, __TASOPTroot__

export airfoil, cdsum!, surfcm, wingsc, wingpo, wingcl, fusebl!

# Define the __init__ function
#This function gets executed automatically when the module is loaded
function __init__()
    BLAS.set_num_threads(1) #This sets the number of threads in BLAS to be equal to 1. 
    #It prevents multithreading but ensures consistent speed across CPU families. Without it,
    #the LU calculation in blax() can take up to 1000x longer.
    #TODO this may cause issues if parallelization is attempted in the future. Other approaches are to 
    #match the number of BLAS threads to the CPUs available on the machine or server
end

#include index to access arrays
include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils/spline.jl"))

idim::Int = 360
jdim::Int = 360
 t     = zeros(Float64, jdim)
 y     = zeros(Float64, jdim)
 yp    = zeros(Float64, jdim)
 z     = zeros(Float64, jdim)
 zp    = zeros(Float64, jdim)
 gw    = zeros(Float64, jdim)

 yc    = zeros(Float64, idim)
 ycp   = zeros(Float64, idim)
 zc    = zeros(Float64, idim)
 zcp   = zeros(Float64, idim)
 gc    = zeros(Float64, idim)
 vc    = zeros(Float64, idim)
 wc    = zeros(Float64, idim)
 vnc   = zeros(Float64, idim)
 
# Aerofoil calculations
include("airfoil.jl")
include("airtable.jl")
include("airfun.jl")

airfoil_data = joinpath(__TASOPTroot__,"airfoil_data/C.air")
airsection = airtable(airfoil_data);

include("surfcd.jl")
include("surfcm.jl")
include("wingpo.jl")
include("wingsc.jl")

# Fuselage IBLT calculations
include("fusebl.jl")
include("axisol.jl")
include("blax.jl")
include("blsys.jl")

# Trefftz plane CDi calcs
include("trefftz.jl")

# Total CD calculations 
include("cdsum.jl")


end