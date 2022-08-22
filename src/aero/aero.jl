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

include("../../utils/spline.jl")

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
include("blsys.jl")

# Trefftz plane CDi calcs
include("trefftz.jl")

# Total CD calculations 
include("cdsum.jl")


end