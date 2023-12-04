"""
`engine` is a module that contains all low-fidelity (NPSS is included in a different directory) calculations
required in the aircraft sizing. 
"""
module engine

using StaticArrays

export tfcalc!, mcool, Tmcalc, gas_tset, gaschem
export tfweight, ddct, ddat, gct, gat, tfsize!, Ncmap, ecmap, Ncmap1, ecmap1, etmap, Pimap, tfoper!
export gassum, gassumd, gas_prat, gas_delh, gas_delhd, gas_burn, gas_burnd, gas_mach, gas_machd, gas_mass, gasfuel, gasPr
export hxdesign!, hxweight

include("../misc/index.inc")
include("../misc/constants.jl")
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
include("hxfun.jl")

end