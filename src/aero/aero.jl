"""
`aerodynamics` is a module that contains all aerodynamic calculations
required in the aircraft sizing. Only functions used in other routines
are exported other functions are kept within the `aerodynamics` namespace.
"""
module aerodynamics

using StaticArrays
using ..atmosphere

export cdsum!, surfcm, wingsc, wingpo, wingcl, fusebl!

#include index to access arrays
include("../misc/index.inc")

include("../utils/spline.jl")

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
include("airtable2.jl")
include("airfun2.jl")
include("surfcd.jl")
include("surfcm.jl")
include("wingpo.jl")
include("wingsc.jl")

# Fuselage IBLT calculations
include("fusebl.jl")
include("axisol.jl")
include("blax.jl")
include("blax2.jl")
include("blsys.jl")

# Trefftz plane CDi calcs
include("trefftz.jl")

# Total CD calculations 
include("cdsum.jl")


end