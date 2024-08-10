"""
`structures` is a module that contains all the structural modules
required to size an aircraft
"""

module structures

using ..atmosphere
using ..materials

using NLsolve
using Roots
using NLopt

export surfw, surfdx, fusew!, tailpo, tanksize!,
 update_fuse!, update_fuse_for_pax!


include("../misc/index.inc")
include("../misc/constants.jl")
include("loads.jl")
export î, ĵ, k̂, WORLD, Weight

#include fuselage sizing
include("../misc/layout.jl")
include("cross_sections.jl")
export SingleBubble, MultiBubble
include("../misc/structuralMember.jl")
include("../misc/fuselage.jl")
export Fuselage
include("fuseW.jl")
include("../misc/fuselage_geometry.jl")

#include sizing of surfaces
include("surfdx.jl")
include("surfw.jl")
include("tailpo.jl")

#Hydrogen tank related code
include("update_fuse.jl")

end
