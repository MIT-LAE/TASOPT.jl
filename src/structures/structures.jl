"""
`structures` is a module that contains all the structural modules
required to size an aircraft
"""

module structures
using NLopt

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, place_cabin_seats, find_cabin_width, find_floor_angles, optimize_double_decker_cabin, find_double_decker_cabin_length

export surfw, surfdx, fusew, tailpo, update_fuse!, update_fuse_for_pax!

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"misc/constants.jl"))
#include fuselage sizing
include("fuseW.jl")

#include sizing of surfaces
include("surfdx.jl")
include("surfw.jl")
include("tailpo.jl")

#Hydrogen tank related code
include("update_fuse.jl")

end
