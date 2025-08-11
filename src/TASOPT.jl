"""
TASOPT - Transport Aircraft System Optimization
"""
module TASOPT

const TASOPT = @__MODULE__
# Add basic packages required by TASOPT
using Base: SignedMultiplicativeInverse, @kwdef
using NLopt: G_MLSL_LDS, GN_MLSL_LDS, GN_CRS2_LM, GN_DIRECT_L

using Printf

using StaticArrays

using Dates
using ForwardDiff
using CSV, Tables
using DocStringExtensions
using Plots, StatsPlots, Plots.PlotMeasures

#convenient directories
const __TASOPTroot__ = @__DIR__
const __TASOPTindices__ = joinpath(__TASOPTroot__,"data_structs/index.inc") #include(__TASOPTindices__) in REPL
export __TASOPTroot__, __TASOPTindices__

# Constants and array indices
include(__TASOPTindices__)

include(joinpath(__TASOPTroot__,"utils/constants.jl"))
export ft_to_m, in_to_m, nmi_to_m, deg_to_rad, 
       lbf_to_N, kts_to_mps, hp_to_W, lb_N
export gee, gamSL, cpSL, μAir, pref, Tref

include(joinpath(__TASOPTroot__,"utils/units.jl"))
export convertMass, convertForce, convertDist, 
       convertSpeed, convertPower, convertAngle

include(joinpath(__TASOPTroot__,"data_structs/materials.jl"))
using .materials
export StructuralAlloy, Conductor, Insulator, ThermalInsulator

#miscellaneous auxiliary fxns
include("./utils/helper_functions.jl")


#Load modules
include(joinpath(__TASOPTroot__,"utils/aircraft_utils.jl"))
include(joinpath(__TASOPTroot__,"atmos/atmos.jl"))
include(joinpath(__TASOPTroot__,"sizing/size_aircraft.jl"))
include(joinpath(__TASOPTroot__,"mission/mission_iteration.jl"))
include(joinpath(__TASOPTroot__,"mission/fly_mission.jl"))
include(joinpath(__TASOPTroot__,"mission/takeoff.jl"))
include(joinpath(__TASOPTroot__,"aero/aero.jl"))
export plot_airf
include(joinpath(__TASOPTroot__,"structures/structures.jl"))
include(joinpath(__TASOPTroot__,"balance/balance.jl"))
include(joinpath(__TASOPTroot__,"engine/engine.jl"))
include(joinpath(__TASOPTroot__,"propsys/propsys.jl"))

include(joinpath(__TASOPTroot__,"data_structs/fuselage_tank.jl"))

export fuselage_tank

#Use above modules
using .atmosphere
using .aerodynamics
using .structures
using .engine

# Load primary aircraft structure 
include(joinpath(__TASOPTroot__,"data_structs/landing_gear.jl"))
include(joinpath(__TASOPTroot__,"data_structs/options.jl"))
include(joinpath(__TASOPTroot__,"data_structs/aircraft.jl"))
export aircraft, fuselage_tank

#Include cryogenic tanks after loading Fuselage and fuselage_tank
include(joinpath(__TASOPTroot__,"cryo_tank/CryoTank.jl"))
using .CryoTank

# Off-design performance via BADA file like output
#  and LTO output for EDB points for use in AEIC
include(joinpath(__TASOPTroot__,"mission/odperformance.jl"))
export fly_mission!
include(joinpath(__TASOPTroot__,"mission/AircraftDeck.jl"))

# Input and output functions
include(joinpath(__TASOPTroot__,"IO/read_input.jl"))
include(joinpath(__TASOPTroot__,"IO/output_texts.jl"))
include(joinpath(__TASOPTroot__,"IO/output_plots.jl"))
export stickfig, plot_details, PayloadRange
include(joinpath(__TASOPTroot__,"IO/save_model.jl"))

include(joinpath(__TASOPTroot__,"IO/quicksave_load.jl"))
include(joinpath(__TASOPTroot__,"IO/par_array_opers.jl"))
include(joinpath(__TASOPTroot__,"IO/output_csv.jl"))

include(joinpath(__TASOPTroot__,"cost/cost_est.jl"))
include(joinpath(__TASOPTroot__,"cost/cost_val.jl"))
include(joinpath(__TASOPTroot__,"utils/printBADA.jl"))
include(joinpath(__TASOPTroot__,"utils/sensitivity.jl"))

export size_aircraft!


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

sizes the given `aircraft` instance. A light wrapper around the `_size_aircraft!` function, which does the actual work.
"""
function size_aircraft!(ac::aircraft; iter=35, initwgt=false, Ldebug=false,
        printiter=true, saveOD=false)

    Ldebug && println("Max weight iterations = $iter")
    _size_aircraft!(ac, itermax = iter, initwgt = initwgt,
        Ldebug = Ldebug, printiter = printiter,
        saveODperf = saveOD)

    #if sized properly, mark as such
    #TODO: apply logic and exit codes to make check more robust
    ac.is_sized .= true
    ;
end # size_aircraft!

end # module TASOPT
