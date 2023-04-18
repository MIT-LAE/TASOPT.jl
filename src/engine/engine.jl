"""
`engine` is a module that contains all low-fidelity (NPSS is included in a different directory) calculations
required in the aircraft sizing. 
"""

module engine

using StaticArrays

export tfcalc!, mcool, tfweight

include("gasfun.jl")
include("gascalc.jl")
# include("tfan.jl")
include("tfmap.jl")
include("tfcool.jl")
include("tfsize.jl")
include("gaussn.jl")
include("compare.jl")
include("tfoper.jl")
include("tfcalc.jl")
include("tfweight.jl")

end