"""
`structures` is a module that contains all the structural modules
required to size an aircraft
"""

module structures

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, place_cabin_seats

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
