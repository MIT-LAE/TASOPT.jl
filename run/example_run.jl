using Base: SignedMultiplicativeInverse
using NLopt: G_MLSL_LDS, GN_MLSL_LDS, GN_CRS2_LM, GN_DIRECT_L

using BenchmarkTools
using Printf

using StaticArrays
using Profile, UnicodePlots
using PyPlot
using Dates

include("../src/misc/index.inc")
include("../tasopt.jl")

include("../balance/balance.jl")

include("../engine/PMSM.jl")  # Motor/generator functions
include("../engine/PMSM.inc") # Motor/generator properties array
include("../engine/NPSS_functions.jl") # NPSS functions
include("../engine/PT.inc")
include("../engine/propsys.jl")

include("../fuel/hydrogen.jl")

include("../sizing/wsize.jl")
include("../src/mission/mission.jl")

include("input.jl")
include("../IO/outputs.jl")
include("../IO/savemodel.jl")
include("../cost/cost_est.jl")
include("../utils/printBADA.jl")
include("../src/mission/odperformance.jl")
include("../contrail/AircraftDeck.jl")
include("../src/mission/LTO.jl")

const gee = 9.81
TSL, pSL, ρSL, aSL, μSL = atmos(0.0)
RSL =  pSL/(ρSL * TSL)
const gamSL = 1.4
const cpSL  = 1004.0

# temporary
μAir = 1.78e-5 
ρAir = ρSL

time_writing = 0.0
time_run_NPSS = 0.0
f = open("temp.results", "w")
# f = stdout

initwgt = 0
saveOD = false
track_fig = nothing
opt_iter_counter = 0
function run_wsize(iter, initwgt, Ldebug, printiter, saveOD)
    global time_writing = 0.0
    global time_run_NPSS = 0.0
    parpt[ipt_time_NPSS] = 0.0
    parpt[ipt_calls_NPSS] = 0
    global opt_iter_counter += 1

    Ldebug && println("Max weight iterations = $iter")
    wsize(pari, parg, parm, view(para,:,:,1), view(pare, :,:,1),
                    iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, saveOD)

    ## Write outputs to io stream f
    printoutput = false
    if printoutput
        println(f, "AR  = $(parg[igAR])\n")
        println(f, "CL  = $(para[iaCL, ipcruise1])\n")
        println(f, "Alt  = $(para[iaalt, ipcruise1])\n")
        @printf(f,"--------------------\n")
        println(f, "Fe    = np.array(",pare[ieFe,:], ")")

        println(f, "ηmot     = np.array(", pare[ieemot    , :], ")")
        println(f, "ηinv     = np.array(", pare[ieeinv    , :], ")")
        println(f, "ηcable   = np.array(", pare[ieecable  , :], ")")
        println(f, "ηgen     = np.array(", pare[ieegen    , :], ")")
        println(f, "ηthermal = np.array(", pare[ieethermal, :], ")")

        println(f, "Hmot     = np.array(", pare[ieHrejmot , :], ")")
        println(f, "Hinv     = np.array(", pare[ieHrejinv , :], ")")
        println(f, "Hcable   = np.array(", pare[ieHrejcab , :], ")")
        println(f, "Hgen     = np.array(", pare[ieHrejgen , :], ")")
        println(f, "Htot     = np.array(", pare[ieHrejtot , :], ")")
        println(f, "EINOx1 = np.array(", pare[ieEINOx1, :], ")")
        println(f, "EINOx2 = np.array(", pare[ieEINOx2, :], ")")
        println(f, "FAR = np.array(", pare[ieFAR, :], ")")

        println(f, "h     = np.array(",para[iaalt,:],")")
        println(f, "R     = np.array(",para[iaRange,:],")")
        println(f, "deNOx = np.array(",pare[iedeNOx, :],")")
        println(f, "fracW = np.array(",para[iafracW, :],")")
        println(f, "mdotf = np.array(",pare[iemdotf, :],")")
        println(f, "mdotH2O = np.array(",pare[iemdotf, :].* 9.0,")")
        println(f, "Ptank = np.array(",pare[iePLH2, :],")")
        println(f, "CL = np.array(",para[iaCL, :],")")
        println(f, "CD = np.array(",para[iaCD, :],")")
        println(f, "CLh = np.array(",para[iaCLh, :],")\n")


        @printf(f, "\nPFEI = %.5f J/Nm\n", parm[imPFEI])

        weight_buildup(parg, io = f)
        aero(parg, para, io = f)
        geometry(parg, io = f)
    end

    # global track_fig = stickfig(parg, pari,  parm; ax = track_fig)
    if (opt_iter_counter%5 == 0) || (opt_iter_counter == 1)
        global track_fig = plot_details(parg, pari, para, parm; ax = track_fig)
    end
end

close(f)
time_wsize = @elapsed run_wsize(35, 0, false, true, saveOD)
# @profview run_wsize(25, 0, false, true)

geometry(parg)
weight_buildup(parg)

