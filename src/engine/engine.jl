"""
`engine` is a module that contains all low-fidelity (NPSS is included in a different directory) calculations
required in the aircraft sizing. 
"""
module engine

using NLopt
using Roots
using LinearAlgebra

export Engine, store_engine_model!, extract_engine_model

export tfwrap!, tfcalc!, mcool, Tmcalc, gas_tset, gaschem
export tfweightwrap!, tfweight, ddct, ddat, gct, gat, tfsize!, Ncmap, ecmap, Ncmap1, ecmap1, etmap, Pimap, tfoper!

export gassum, gassumd, gas_prat, gas_delh, gas_delhd, gas_burn, gas_burnd, gas_mach, gas_machd, gas_mass, gasfuel, fuelLHV, gasPr
export hxdesign!, radiator_design!, hxweight, resetHXs

import ..TASOPT: __TASOPTindices__, __TASOPTroot__, StructuralAlloy

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"misc/constants.jl"))
include("../misc/engine.jl")
include("gasfun.jl")
include("gascalc.jl")
# include("tfan.jl")
include("tfmap.jl")
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

end