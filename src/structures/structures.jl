"""
`structures` is a module that contains all the structural modules
required to size an aircraft
"""

module structures
using NLopt

using ..atmosphere
using ..materials

using NLsolve
using Roots
using NLopt
import ..TASOPT: __TASOPTindices__, __TASOPTroot__

export surfw, surfdx, fusew!, tailpo,
 update_fuse!, update_fuse_for_pax!, place_cabin_seats, find_cabin_width, find_floor_angles, arrange_seats


include("../misc/index.inc")
include("../misc/constants.jl")
include("loads.jl")
export î, ĵ, k̂, WORLD, Weight

#include fuselage sizing
include("../misc/layout.jl")
export SingleBubble, MultiBubble, scaled_cross_section
include("../misc/structuralMember.jl")
include("../misc/fuselage.jl")
export Fuselage
include("fuseW.jl")
include("../misc/fuselage_geometry.jl")

#include sizing of surfaces
include("surfdx.jl")
include("surfw.jl")
include("tailpo.jl")

include("size_cabin.jl") #Seat layouts and cabin length

#Hydrogen tank related code
include("update_fuse.jl")


end
