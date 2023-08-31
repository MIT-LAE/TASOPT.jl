"""
TASOPT
"""
module TASOPT

export atmos, size_aircraft

# Add basic pacakges required by TASOPT
using Base: SignedMultiplicativeInverse
using NLopt: G_MLSL_LDS, GN_MLSL_LDS, GN_CRS2_LM, GN_DIRECT_L

using BenchmarkTools
using Printf

using StaticArrays
using Profile, UnicodePlots
using PyPlot
using Dates
using ForwardDiff

const __TASOPTroot__ = @__DIR__

# Constants and array indices
include("./misc/constants.jl")
export ft_to_m, in_to_m, nmi_to_m, deg_to_rad, 
       lbf_to_N, kts_to_mps, hp_to_W, lb_N
export gee, gamSL, cpSL, μAir, pref, Tref

include("./misc/units.jl")
export convertMass, convertForce, convertDist, 
       convertSpeed, convertPower, convertAngle

include("./misc/index.inc")
include("./misc/aircraft.jl")
export aircraft

#Load modules
include(joinpath(__TASOPTroot__,"atmos/atmos.jl"))
include(joinpath(__TASOPTroot__,"sizing/wsize.jl"))
include(joinpath(__TASOPTroot__,"mission/mission.jl"))
include(joinpath(__TASOPTroot__,"mission/takeoff.jl"))
include(joinpath(__TASOPTroot__,"aero/aero.jl"))
include(joinpath(__TASOPTroot__,"structures/structures.jl"))
include(joinpath(__TASOPTroot__,"propsys/propsys.jl"))
include(joinpath(__TASOPTroot__,"balance/balance.jl"))
include(joinpath(__TASOPTroot__,"engine/engine.jl"))

# Off-design performance via BADA file like output
#  and LTO output for EDB points for use in AEIC
include(joinpath(__TASOPTroot__,"mission/odperformance.jl"))
include(joinpath(__TASOPTroot__,"mission/woper.jl"))
include(joinpath(__TASOPTroot__,"mission/LTO.jl"))
include(joinpath(__TASOPTroot__,"mission/AircraftDeck.jl"))

include(joinpath(__TASOPTroot__,"fuel/hydrogen.jl"))
include(joinpath(__TASOPTroot__,"engine/PT.inc"))

# Input and output functions
include(joinpath(__TASOPTroot__,"IO/outputs.jl"))
include(joinpath(__TASOPTroot__,"IO/savemodel.jl"))

include(joinpath(__TASOPTroot__,"cost/cost_est.jl"))
include(joinpath(__TASOPTroot__,"cost/cost_val.jl"))
include(joinpath(__TASOPTroot__,"utils/printBADA.jl"))

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
RSL = pSL / (ρSL * TSL)
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

function size_aircraft(ac::aircraft, iter, initwgt, Ldebug, printiter, saveOD)

    Ldebug && println("Max weight iterations = $iter")
    wsize(ac.pari, ac.parg, ac.parm, 
        view(ac.para, :, :, 1), view(ac.pare, :, :, 1),
        iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, saveOD)

end

end