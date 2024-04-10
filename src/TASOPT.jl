"""
TASOPT
"""
module TASOPT

# Add basic pacakges required by TASOPT
using Base: SignedMultiplicativeInverse, @kwdef
using NLopt: G_MLSL_LDS, GN_MLSL_LDS, GN_CRS2_LM, GN_DIRECT_L

using BenchmarkTools
using Printf

using StaticArrays
# using PyCall
using PyPlot
# pygui(true)
using Dates
using ForwardDiff
using CSV, Tables
using DocStringExtensions

const __TASOPTroot__ = @__DIR__

# Constants and array indices
include("./misc/constants.jl")
export ft_to_m, in_to_m, nmi_to_m, deg_to_rad, 
       lbf_to_N, kts_to_mps, hp_to_W, lb_N
export gee, gamSL, cpSL, μAir, pref, Tref

include("./misc/units.jl")
export convertMass, convertForce, convertDist, 
       convertSpeed, convertPower, convertAngle

include("./misc/materials.jl")
export StructuralAlloy, Conductor, Insulator

include("./misc/index.inc")
include("./misc/aircraft.jl")
export aircraft, fuselage_tank

#Load modules
include(joinpath(__TASOPTroot__,"atmos/atmos.jl"))
include(joinpath(__TASOPTroot__,"sizing/wsize.jl"))
include(joinpath(__TASOPTroot__,"mission/mission.jl"))
include(joinpath(__TASOPTroot__,"mission/takeoff.jl"))
include(joinpath(__TASOPTroot__,"aero/aero.jl"))
include(joinpath(__TASOPTroot__,"propsys/propsys.jl"))
include(joinpath(__TASOPTroot__,"balance/balance.jl"))
include(joinpath(__TASOPTroot__,"engine/engine.jl"))
include(joinpath(__TASOPTroot__,"structures/structures.jl"))
include(joinpath(__TASOPTroot__,"cryo_tank/CryoTank.jl"))


# Off-design performance via BADA file like output
#  and LTO output for EDB points for use in AEIC
include(joinpath(__TASOPTroot__,"mission/odperformance.jl"))
include(joinpath(__TASOPTroot__,"mission/woper.jl"))
include(joinpath(__TASOPTroot__,"mission/LTO.jl"))
include(joinpath(__TASOPTroot__,"mission/AircraftDeck.jl"))

include(joinpath(__TASOPTroot__,"engine/PT.inc"))

# Input and output functions
include(joinpath(__TASOPTroot__,"IO/read_input.jl"))
include(joinpath(__TASOPTroot__,"IO/outputs.jl"))
include(joinpath(__TASOPTroot__,"IO/save_model.jl"))

include(joinpath(__TASOPTroot__,"IO/quicksave_load.jl"))
include(joinpath(__TASOPTroot__,"IO/par_array_opers.jl"))
include(joinpath(__TASOPTroot__,"IO/read_externals.jl"))
include(joinpath(__TASOPTroot__,"IO/output_csv.jl"))

include(joinpath(__TASOPTroot__,"cost/cost_est.jl"))
include(joinpath(__TASOPTroot__,"cost/cost_val.jl"))
include(joinpath(__TASOPTroot__,"utils/printBADA.jl"))

export size_aircraft!


using .atmosphere
using .aerodynamics
using .structures
using .propsys
using .engine
using .CryoTank

#------------------------------------------------------
#End imports/loading files
#------------------------------------------------------

# Derived constants
TSL, pSL, ρSL, aSL, μSL = atmos(0.0)
RSL = pSL / (ρSL * TSL)
ρAir = ρSL

# ----------------------
# Sizing function
# ----------------------

"""
    size_aircraft(ac::aircraft; iter=35, initwgt=false, Ldebug=false,
        printiter=true, saveOD=false)

sizes the given `aircraft` instance
"""
function size_aircraft!(ac::aircraft; iter=35, initwgt=false, Ldebug=false,
        printiter=true, saveOD=false)

    Ldebug && println("Max weight iterations = $iter")
    wsize(ac, itermax = iter, initwgt = initwgt,
        Ldebug = Ldebug, printiter = printiter,
        saveODperf = saveOD)

    #if sized properly, mark as such
    #TODO: apply logic and exit codes to make check more robust
    ac.sized .= true
end
end
