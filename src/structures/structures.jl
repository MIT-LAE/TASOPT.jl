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
import ..TASOPT: __TASOPTindices__, __TASOPTroot__, unpack_ac, unpack_ac_components, TailSizing

export wing_weights!, calculate_centroid_offset!, calculate_centroid_offset, fusew!,
 update_fuse!, update_fuse_for_pax!, place_cabin_seats, find_cabin_width, find_floor_angles, arrange_seats,
size_landing_gear!

include(joinpath(__TASOPTroot__,"data_structs/index.inc"))
include(joinpath(__TASOPTroot__,"utils/constants.jl"))
include(joinpath(__TASOPTroot__,"structures/loads.jl"))
export î, ĵ, k̂, WORLD, Weight

#include fuselage sizing
include(joinpath(__TASOPTroot__,"data_structs/layout.jl"))
export SingleBubble, MultiBubble, scaled_cross_section
include(joinpath(__TASOPTroot__,"data_structs/structuralMember.jl"))
export StructuralMember
include(joinpath(__TASOPTroot__,"data_structs/fuselage.jl"))
export Fuselage
export AbstractCrossSection
include(joinpath(__TASOPTroot__,"structures/fuseW.jl"))
include(joinpath(__TASOPTroot__,"data_structs/fuselage_geometry.jl"))

include(joinpath(__TASOPTroot__,"data_structs/wingSections.jl"))
include(joinpath(__TASOPTroot__,"data_structs/wing.jl"))
export WingSection,TailSection,Wing,Tail,wing_additional_weight 

#include sizing of surfaces
include(joinpath(__TASOPTroot__,"structures/calculate_centroid_offset.jl"))
include(joinpath(__TASOPTroot__,"structures/wing_weights.jl"))
export WingSectionDimensions

include(joinpath(__TASOPTroot__,"structures/size_cabin.jl")) #Seat layouts and cabin length

#Hydrogen tank related code
include(joinpath(__TASOPTroot__,"structures/update_fuse.jl"))

#Landing gear 
include("../data_structs/landing_gear.jl")# type definitions for the landing gear
include("size_landing_gear.jl") 

end
