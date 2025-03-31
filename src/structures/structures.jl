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
import ..TASOPT: __TASOPTindices__, __TASOPTroot__, unpack_ac, unpack_ac_components

export get_wing_weights!, calculate_centroid_offset!, calculate_centroid_offset, fusew!,
 update_fuse!, update_fuse_for_pax!, place_cabin_seats, find_cabin_width, find_floor_angles, arrange_seats,
size_landing_gear!

include("../data_structs/index.inc")
include("../data_structs/constants.jl")
include("loads.jl")
export î, ĵ, k̂, WORLD, Weight

#include fuselage sizing
include("../data_structs/layout.jl")
export SingleBubble, MultiBubble, scaled_cross_section
include("../data_structs/structuralMember.jl")
export StructuralMember
include("../data_structs/fuselage.jl")
export Fuselage
export AbstractCrossSection
include("fuseW.jl")
include("../data_structs/fuselage_geometry.jl")

include("../data_structs/wingSections.jl")
include("../data_structs/wing.jl")
export WingSection,TailSection,Wing,Tail,wing_additional_weight 

#include sizing of surfaces
include("calculate_centroid_offset.jl")
include("get_wing_weights.jl")

include("size_cabin.jl") #Seat layouts and cabin length

#Hydrogen tank related code
include("update_fuse.jl")

#Landing gear 
include("../misc/landing_gear.jl")# type definitions for the landing gear
include("size_landing_gear.jl") 

end
