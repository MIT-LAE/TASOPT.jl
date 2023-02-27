"""
TASOPT
"""

# Add basic pacakges required by TASOPT
using Base: SignedMultiplicativeInverse
using NLopt: G_MLSL_LDS, GN_MLSL_LDS, GN_CRS2_LM, GN_DIRECT_L

using BenchmarkTools
using Printf

using StaticArrays
using Profile, UnicodePlots
using PyPlot
using Dates

# Constants and array indices
include("src/misc/constants.jl")
include("src/misc/index.inc")

#Load modules
include("src/atmos/atmos.jl")
include("src/sizing/wsize.jl")
include("src/mission/mission.jl")
include("src/aero/aero.jl")
include("src/structures/structures.jl")
include("src/propsys/propsys.jl")
include("src/balance/balance.jl")
include("src/engine/engine.jl")

# Off-design performance via BADA file like output
#  and LTO output for EDB points for use in AEIC
include("src/mission/odperformance.jl")
include("src/mission/woper.jl")
include("src/mission/LTO.jl")
include("src/mission/AircraftDeck.jl")

include("src/fuel/hydrogen.jl")
include("src/engine/PT.inc")

# Input and output functions
include("src/IO/outputs.jl")
include("src/IO/savemodel.jl")

include("src/cost/cost_est.jl")
include("src/utils/printBADA.jl")

using .atmosphere
using .aerodynamics
using .structures
using .propsys
using .engine

#------------------------------------------------------
#End imports/loading files
#------------------------------------------------------

# Derived constants
TSL, pSL, ρSL, aSL, μSL = atmos(0.0)
RSL =  pSL/(ρSL * TSL)
ρAir = ρSL

# ----------------------
# Sizing function
# ----------------------
time_writing = 0.0
time_run_NPSS = 0.0

initwgt = 0
saveOD = false
track_fig = nothing
opt_iter_counter = 0

function size_aircraft(iter, initwgt, Ldebug, printiter, saveOD)
    global time_writing = 0.0
    global time_run_NPSS = 0.0
    parpt[ipt_time_NPSS] = 0.0
    parpt[ipt_calls_NPSS] = 0
    global opt_iter_counter += 1

    Ldebug && println("Max weight iterations = $iter")
    wsize(pari, parg, parm, view(para,:,:,1), view(pare, :,:,1),
                    iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, saveOD)
    
    # global track_fig = stickfig(parg, pari,  parm; ax = track_fig)
    if (opt_iter_counter%5 == 0) || (opt_iter_counter == 1)
        global track_fig = plot_details(parg, pari, para, parm; ax = track_fig)
    end
end
