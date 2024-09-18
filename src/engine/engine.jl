"""
`engine` is a module that contains all low-fidelity (NPSS is included in a different directory) calculations
required in the aircraft sizing. 
"""
module engine

using NLopt
using Roots
using NLsolve
using LinearAlgebra

export tfcalc!, mcool, Tmcalc, gas_tset, gaschem, enginecalc!
export tfweight, ddct, ddat, gct, gat, tfsize!, Ncmap, ecmap, Ncmap1, ecmap1, etmap, Pimap, tfoper!

export gassum, gassumd, gas_prat, gas_delh, gas_delhd, gas_burn, gas_burnd, gas_mach, gas_machd, gas_mass, gasfuel, fuelLHV, gasPr
export hxdesign!, hxweight, resetHXs

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, StructuralAlloy

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"misc/constants.jl"))
include("gasfun.jl")
include("gascalc.jl")
# include("tfan.jl")
include("enginecalc.jl")
include("turbofan/tfmap.jl")
include("turbofan/tfcool.jl")
include("turbofan/tfsize.jl")
include("gaussn.jl")
include("compare.jl")
include("turbofan/tfoper.jl")
include("turbofan/tfcalc.jl")
include("turbofan/tfweight.jl")
include("hxfun.jl")
include("PEMfuelcell.jl")
include("ducted_fan/ductedfancalc.jl")
include("ducted_fan/ductedfansize.jl")
include("ducted_fan/ductedfanoper.jl")

end