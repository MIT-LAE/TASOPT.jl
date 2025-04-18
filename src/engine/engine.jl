"""
`engine` is a module that contains all low-fidelity calculations
required in the aircraft sizing. 
"""
module engine

using NLopt
using Roots
using LinearAlgebra

export Engine

export tfwrap!, tfcalc!, mcool, Tmcalc, gas_tset, gaschem
export tfweightwrap!, tfweight, ddct, ddat, gct, gat, tfsize!, Ncmap, ecmap, Ncmap1, ecmap1, etmap, Pimap, tfoper!

export gassum, gassumd, gas_prat, gas_delh, gas_delhd, gas_burn, gas_burnd, gas_mach, gas_machd, gas_mass, gasfuel, fuelLHV, gasPr
export hxdesign!, radiator_design!, hxweight, resetHXs, HXOffDesign!

export check_engine_convergence_failure

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, StructuralAlloy, unpack_ac, compare_strings

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"utils/constants.jl"))
include("gasfun.jl")
include("gascalc.jl")
# include("tfan.jl")
include("turbomachinery/tfmap.jl")
include("turbomachinery/maps.jl")
include("turbofan/tfcool.jl")
include("turbofan/tfsize.jl")
include("gaussn.jl")
include("compare.jl")
include("turbofan/tfoper.jl")
include("turbofan/tfcalc.jl")
include("turbofan/tfweight.jl")
include("turbofan/tfwrap.jl")
include("turbofan/tfweightwrap.jl")
include("hxfun.jl")
include("PEMfuelcell.jl")
include(joinpath(__TASOPTroot__,"data_structs/engine.jl"))

end