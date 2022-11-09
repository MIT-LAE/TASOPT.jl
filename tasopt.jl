"""
TASOPT
"""

include("src/misc/constants.jl")
include("atmos/atmos.jl")
include("src/aero/aero.jl")
include("src/structures/structures.jl")
include("src/propsys/propsys.jl")
using .atmosphere
using .aerodynamics
using .structures
using .propsys
