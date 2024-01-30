"""
`structures` is a module that contains all the structural modules
required to size an aircraft
"""

module structures

export surfw, surfdx, fusew, tailpo, tanksize

include("../misc/index.inc")
include("../misc/constants.jl")
#include fuselage sizing
include("fuseW.jl")

#include sizing of surfaces
include("surfdx.jl")
include("surfw.jl")
include("tailpo.jl")

#Hydrogen tank related code
include("tankWmech.jl")
include("tankWthermal.jl")
include("tanksize.jl")

end
