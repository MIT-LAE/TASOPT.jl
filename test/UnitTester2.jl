using Base: SignedMultiplicativeInverse
using NLopt: G_MLSL_LDS, GN_MLSL_LDS, GN_CRS2_LM, GN_DIRECT_L

using BenchmarkTools
using Printf


using StaticArrays
using Profile, UnicodePlots
using PyPlot
# pygui(true)
# plt.style.use(["./prash.mplstyle", "tableau-colorblind10"])
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
# include("../input_fwdTank.jl")
include("../IO/outputs.jl")
include("../IO/savemodel.jl")
# include("../fuselage/update_fuse.jl")
include("../cost/cost_est.jl")
include("../utils/printBADA.jl")
include("../src/mission/odperformance.jl")
include("../src/mission/woper.jl")
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

    # println(parg[igWwing])
    Ldebug && println("Max weight iterations = $iter")
    wsize(pari, parg, parm, view(para,:,:,1), view(pare, :,:,1),
                    iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, saveOD)
    # @benchmark PowerTrain(0.0, 0.8, 25.0e3*2,0.0, 0.0, parg, parpt, parmot, pargen)

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


        # println("Time netNPSS       = $(parpt[ipt_time_NPSS])")
        # println("Time writing       = $time_writing")
        # println("Time runnning NPSS = $time_run_NPSS")
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

time_wsize = @elapsed run_wsize(35, 0, false, true, saveOD)
# @profview run_wsize(25, 0, false, true)
# println("Wsize time = $time_wsize s")
# println("NPSS Time = $(parpt[ipt_time_NPSS]) s")
# println("NPSS Calls = $(parpt[ipt_calls_NPSS])")


# include("../LoadModel.jl")

# Final aircraft:
# include(".././ZIHA/ZIA_BLI_10_8_0.760_15.3.mdl")
# include(".././ZISA/ZIA_SAF_BLI_10_8_0.647_17.4.mdl")

# saveOD = false
# time_wsize = @elapsed run_wsize(30, 1, false, false, saveOD)
# println(para[iagamV, ipclimbn])

# load_aircraft("ZIA_BLI_10_8_0.811_15.8.mdl", "./Models", true)
CostEst(pari, parg, pare, parm, parpt, 1000)


# # ------ Optimzation ---------------
# using NLopt
# xarray = []
# farray = []
# PFEIarray = []
# CDarray = []
# WMTOarray = []

# function obj(x, grad)
#     parg[igAR] = x[1]
#     para[iaalt, ipcruise1, :] .=  x[2] *ft_to_m
#     para[iaCL  , ipclimb1+1:ipdescentn-1, :] .= x[3]
#     parpt[ipt_pifan] = x[4]
#     parg[igsweep   ] = x[5]
    
#     # para[iaMach, ipclimbn:ipdescent1  , :] .= x[6]
#     pare[ieTt4, 1:iptotal, :] .= x[6] # [R]
#     parpt[ipt_Tt41] = pare[ieTt4,ipcruise1, 1]

#     pare[ieTt4, ipstatic:iptakeoff, :] .= x[7] #[R]

#     parg[iglambdas]  = x[8]
#     parg[iglambdat]  = x[9]

#     parg[ighboxo   ] = x[10]
#     parg[ighboxs   ] = x[11]

#     parpt[ipt_piHPC] = x[12]

#     para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= x[13]   #  rcls  
#     para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= x[14]   #  rclt  

#     parg[igARh] = x[15]
#     parg[igsweeph] = x[16]
#     parpt[ipt_Fnsplit] = x[17]
#     # parg[igsweeph] = parg[igsweep]

#     # parg[igRfuse] = x[16]
#     # parg[iglftankin] = x[17]

#     wsize_time = @elapsed run_wsize(50, 1, false, false, saveOD)

#     f = parm[imPFEI]
#     push!(PFEIarray, parm[imPFEI])
#     push!(xarray, x)
#     push!(CDarray, para[iaCD, ipcruise1, 1])
#     push!(WMTOarray, parg[igWMTO])
    
#     # Max span constriant
#     bmax = parg[igbmax]
#     b    = parg[igb]
#     constraint  = b/bmax - 1.0
#     penfac  = 25.0* parg[igWpay]
#     # f = f + penfac*max(0.0, constraint)^2

#     # Min climb gradient
#     gtocmin = parg[iggtocmin]
#     gtoc    = para[iagamV, ipclimbn]
#     constraint = 1.0 - gtoc/gtocmin
#     penfac = 1.0*parg[igWpay]
#     f = f + penfac*max(0.0, constraint)^2

#     # Max Tt3 at TOC 
#     Tt3max = 900.0
#     Tt3    = maximum(pare[ieTt3, :, 1])
#     constraint = Tt3/Tt3max - 1
#     penfac = 5.0*parg[igWpay]
#     f = f + penfac*max(0.0, constraint)^2
    
#     # Max Tmetal at Rotation
#     Tvanemax = 2400 #Rankine
#     Tvane    = maximum(pare[ieTmet1, :, 1]) #Rankine
#     constraint = Tvane/Tvanemax - 1
#     penfac = 5.0*parg[igWpay]
#     # f = f + penfac*max(0.0, constraint)^2


#     # Ensure fuel volume makes sense
#     Wfmax = parg[iglftankin]  #parg[igWfmax]
#     Wf    = parg[iglftank]    #parg[igWfuel]
#     if (pari[iifwing] ==1)
#         Wfmax = parg[igWfmax]
#         Wf    = parg[igWfuel]
#     end
#     constraint = Wf/Wfmax - 1.0
#     penfac = 10*parg[igWpay]
#     f = f + penfac*max(0.0, constraint)^2
#     end

#     # Max Fan diameters
#     dfanmax = 2.0
#     dfan = parg[igdfan]
#     constraint = dfan/dfanmax - 1.0
#     penfac = parg[igWpay]
#     f = f + penfac*max(0.0, constraint)^2

#     # Ensure fans will fit within fuselage
#     daftfanmax = min(3.0, 0.9*parg[igRfuse])
#     daftfan = parg[igdaftfan]
#     constraint = daftfan/daftfanmax - 1.0
#     penfac = parg[igWpay]
#     f = f + penfac*max(0.0, constraint)^2

#     # Ensure Fans will fit on wing
#     bo = parg[igbo]
#     b  = parg[igb]
#     lmax = b/2*4/5 - (bo/2+2*parg[igdfan])
#     lfans = parg[igneng]/2*parg[igdfan]*1.25
#     constraint = lfans/lmax - 1.0
#     penfac = parg[igWpay]
#     f = f + penfac*max(0.0, constraint)^2
    
#     if (opt_iter_counter%10 == 0) || (opt_iter_counter == 2)
#     printstyled(@sprintf("%6s %4s %4s %5s %5s %5s  %4s  %4s  %5s %5s %5s  %5s  %5s %5s %5s %5s %5s  %5s | %-8s %10s %5s %6s %5s  %6s  %6s %6s %6s %5s %6s %5s  %5s \n",
#     "t", "AR", "h", "CLcr", "FPR", "λ", "T4c", "T4r", "λs", "λt", "hbo", "hbs", "πHPC", "rcls", "rclt", "λh", "ARh", "Fnrat",
#            "f", "PFEI[J/Nm]", "Wf[tons]", "Wwing", "CDwing", "Wf/Wfmax", "γ", "Tt3", "Tvane", "b", "L/D", "dfan", "dfanaft" ); color = :light_green)
#     end
    
#     printstyled(@sprintf("%5.1fs %4.1f %4.1f %5.4f %5.2f %5.2fᵒ %4.0fR %4.0fR %5.4f %5.4f %5.4f  %5.4f  %5.2f %5.3f %5.3f %5.2f %5.2fᵒ %5.2f | %8.6e %10.6f %5.2f %6.3f %5.4f  %6.4f  %6.4f %6.1f %6.1f %5.1fm %6.1f %5.3fm  %5.3fm\n",
#                     wsize_time, x[1], x[2]/1e3, x[3],x[4],x[5],  x[6],  x[7], x[8], x[9],  x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17], 
#     f, parm[imPFEI], parm[imWfuel]/9.81/1000, parg[igWwing]/9.81/1000, para[iaCDwing, ipcruise1], parg[igWfuel]/parg[igWfmax], para[iagamV, ipclimbn, 1], Tt3, Wf, Wfmax, dfan, daftfan); color = :light_green)

#     # println("X̄ = $x  ⇨  PFEI = $(parm[imPFEI]) f = $f, MTOW = $(parg[igWMTO]), Wtank = $(parg[igWftank]), Wfan = $(parg[igWfan])")
#     push!(farray, f)

#     return f
# end
# # Des. vars:  AR    Alt      Cl    FPR   Λ     Tt4 CR  Tt4 TO  λs   λt   hboxo  hboxs  πHPC  rcls  rclt  ARh     Λh  Fnsplit Rf   lftank   Λh
# lower      = [7.0 , 20000.0, 0.40, 1.20, 10.0, 3200.0, 3200.0, 0.1, 0.1,  0.10,  0.10, 10.0, 0.1,  0.1,  4.0,   5.0, 0.45 ]#, 2.83,  6.0]#,  5.0] 
# upper      = [12.0, 60000.0, 0.65, 1.60, 40.0, 3400.0, 3400.0, 1.0, 1.0,  0.15,  0.15, 20.0, 1.2,  1.0,  8.0,  30.0, 0.7 ]#, 2.87, 20.0]#, 30.0] 

# # got from prior optimization run 
# initial = [7.773081944713115, 40181.603962553054, 0.5527508383478996, 1.230675138085258, 28.62473706060067, 3203.5726595187757, 3598.340134020528, 0.1964082835266504, 0.14368441782959446, 0.11825599806307492, 17.045255858009348, 1.1116365925303677, 0.922992859170234, 5.0445229907657145, 24.70786719622403, 0.5499037770914472]
# initial =  [11.0, 34937.0, 0.500, 1.2126591103869193, 26.686300660179352, 3200.0030186046956, 3400.0, 0.1752945276646451, 0.14, 0.14, 16.71, 1.1030458855952907, 0.9, 6.0, 0.3, 25.0]
# # [9.997431096436161, 38856.75360556306, 0.6045675203477272, 1.2235783865690488, 27.657069015636147, 3218.0635400626825, 3594.8081821795568, 0.13775246610016734, 0.12342780723694166, 0.11848599767102688, 7.071610375877828, 1.199171364736705, 0.9976217470483848, 5.769453459306792, 0.55, 30.0]
# initial = [9.800808554263122, 39209.734944039534, 0.5435602263841047, 1.2187326499096856, 26.415133350952154, 3203.546649444036, 3599.4481633997107, 0.9, 0.10009111870159837, 0.1398701097353828, 0.12337129359538086, 18.881719264268169, 1.0369829799562946, 0.9734146535956018, 7.668225570367148, 29.515164166878304, 0.5490999469170115]
# initial = [8.228454349920403, 39100.86744271755, 0.5852077089083737, 1.2524597900932806, 26.45720878711301, 3201.530593286635, 3599.9974291129274, 0.9655671908068356, 0.10164563435148659, 0.14922775264607838, 0.10853595139952857, 19.364634746000252, 1.0478356751809224, 0.9735662117689496, 5.787416863799637, 29.56757563957012, 0.5308090351667302]

# initial = [parg[igAR], 33000.0, 0.57,
# 1.20, parg[igsweep], 3200, 3300,
# parg[iglambdas], parg[iglambdat], parg[ighboxo   ], parg[ighboxs   ],
# 15.0, para[iarcls, ipcruise1], para[iarclt, ipcruise1], 
# parg[igARh], parg[igsweeph], 0.5]

# initial_dx = [ 0.5  ,  1000.0, 0.005 , 0.05,  0.1,   15.0,   15.0, 0.01, 0.01,  0.001,  0.001,  2.0, 0.01,  0.01, 2.0, 0.1, 0.1]#, 0.2, 2.0]#, 5.0]
# # initial_dx = [ 0.1  ,  10.0, 0.001 , 0.01,  0.01,   1.0,   1.0, 0.001, 0.001,  0.005,  0.005,  0.5, 0.001,  0.001, 0.10, 0.01, 0.05]#, 0.2, 2.0]#, 5.0]

# # if pari[iifuel] == 2
# #     initial = [8.706823321294936, 39324.133600477624, 0.6334358715728483, 1.1932525623004469, 26.344302712899477, 3200.3360554740043, 3459.6712394659494, 0.9, 0.1094619379183989, 0.14110572372270128, 0.11737767098627176, 13.208536486517602, 1.0628009618249974, 0.9933437602824484, 6.542082548728352, 29.246724429560665, 0.5241858300077095]
# #     initial = [8.498746697625975, 36107.236835519856, 0.6412478193423614, 1.2196594875246485, 25.55061885049667, 3200.061100009788, 3300.622159663861, 0.9446505022053928, 0.10000058746957677, 0.12809079514244112, 0.10249631007922386, 15.0, 1.0068733417143467, 0.9546250624086424, 6.4065976827230315, 29.182176695467813, 0.5499999999999989]
# # end
# # include("../ZIA.mdl") # note this will override constraints too!!
# # include("../ZIA_343_918.mdl")

# # x_tol_abs = [0.01, 50.0, 0.0001, 0.001, 0.05, 1.0, 1.0, 0.0001, 0.001, 0.001, 0.01, 0.001, 0.05]
# # f_tol_rel = 1e-6
# f_tol_rel = 1e-5

# opt = NLopt.Opt(:LN_NELDERMEAD, length(initial))
# # opt = NLopt.Opt(:LN_BOBYQA, length(initial))
# # opt = NLopt.Opt(:LN_COBYLA, length(initial))
# # opt = NLopt.Opt(:LN_SBPLX, length(initial))
# # opt = NLopt.Opt(:GN_DIRECT_L, length(initial))

# # opt = NLopt.Opt(:GN_MLSL, length(initial))
# # local_opt = NLopt.Opt(:LN_NELDERMEAD, length(initial))
# # local_opt.ftol_rel = 1e-4
# # local_opt.lower_bounds = lower
# # local_opt.upper_bounds = upper
# # local_opt.min_objective = obj
# # local_opt.initial_step = initial_dx

# # opt.local_optimizer = local_opt

# opt.lower_bounds = lower
# opt.upper_bounds = upper
# opt.min_objective = obj
# opt.initial_step = initial_dx

# # opt.xtol_abs = x_tol_abs
# opt.ftol_rel = f_tol_rel
# # opt.maxeval = 10
# # # nprop = [4, 6, 8, 10, 12, 14, 16]
# # include(".././Models/ZIA_BLI_10_ 8_0.781_14.5.mdl")
# # include(".././Models/ZIA_BLI_10_6_0.823_15.2.mdl")
# # include(".././Models/ZIA_BLI_10_12_0.847_15.2.mdl")    
# nprop = [6, 8, 10, 12]
# nprop = [8]
# # para[iaMach, ipclimbn:ipdescent1  , :] .= 0.8
# Wpay = parg[igWpay]
# frac = [1.0, 0.9, 0.8, 0.6]
# # frac = [0.6]
# Alts = []
# Weight = []
# PaxWeight = []
# S = []
# CL = []
# CD = []
# b = []
# # for n in nprop
# #     # parm[imWpay] = Wpay*n
# #     parpt[ipt_nfan] = n
# #     parg[igneng] =  parpt[ipt_nfan]

# # #     printstyled(@sprintf("\t%10s  %5s  %10s  %6s  %5s  %5s  %10s  %10s  %5s  %5s  %5s  %10s  %5s  %5s %5s  %5s| %10s  %10s  %10s  %10s  %10s  %10s  %10s \n",
# # #     "time2size", "AR", "hcr[ft]", "CLcr", "FPR", "λ[deg]", "Tt4cr[R]", "Tt4ro[R]", "λt", "hbo", "hbs", "πHPC", "rcls", "rclt", "λh", "ARh",
# # #             "PFEI[J/Nm]", "Wf[tons]", "Wf/Wfmax", "γ", "Tt3", "ltank", "ltankin" ); color = :light_green)
# # println(pari[iiVTsize])
#     # opt_time = @elapsed (optf, optx, ret) = NLopt.optimize(opt, initial)
#     # numevals = opt.numevals # the number of function evaluations

#     # # global initial = optx 

#     # println("got $optf at $optx after $numevals iterations which took $(opt_time/60) min (returned $ret)")

#     # savedir = "./Figures/"
#     # if pari[iifuel] == 1
#     #     figname = @sprintf("ZIA_BLI_%d_%d_%.3f_%.1f", seats_per_row, parg[igneng], parm[imPFEI],  para[iaCL, ipcruise1]/para[iaCD, ipcruise1])
#     # elseif pari[iifuel] == 2
#     #     figname = @sprintf("ZIA_SAF_BLI_%d_%d_%.3f_%.1f", seats_per_row, parg[igneng], parm[imPFEI],  para[iaCL, ipcruise1]/para[iaCD, ipcruise1])
#     # end
#     # global track_fig = plot_details(parg, pari, para, parm; ax = track_fig)
#     # plt.savefig(savedir*figname*".png")
#     # savemodel("./Models/"*figname*".mdl", pari, parg, parm, para, pare, parpt, parmot, pargen)
# # push!(Alts, para[iaalt, ipcruise1])
# # push!(Weight, parg[igWMTO])
# # push!(PaxWeight, parg[igWpay])
# # push!(S, parg[igS])
# # push!(CL, para[iaCL, ipcruise1])
# # push!(CD, para[iaCD, ipcruise1])
# # push!(b, parg[igb])

# # #     fig, ax = plt.subplots()
# # #     ax.plot(farray)
# # #     ax.set_xlabel("Iterations")
# # #     ax.set_ylabel("Obj")

# # end

# # fig, ax = plt.subplots(3,1, dpi = 100, sharex = true)
# # ax[1].plot(Weight/9.81/1000, Alts/ft_to_m/1000, ".-", ms = 10)
# # ax[1].set_ylabel("Cruise altitude [kft]")
# # ax[2].plot(Weight/9.81/1000, CL, ".-", ms = 10, label = "CL")
# # ax[2].plot(Weight/9.81/1000, CD*10, ".-", ms = 10, label = "CD*10")
# # ax[2].legend()
# # ax[3].plot(Weight/9.81/1000, b, ".-", ms = 10)

# # ax[end].set_xlabel("MTOW [tonnes]")

# # # # Prepare comparison 
# # models = readdir("./OldModels")
# # PFEI  = zeros(length(models))
# # nProp = zero(PFEI)
# # Ssmax = zero(PFEI)
# # Somax = zero(PFEI)
# # Msmax = zero(PFEI)
# # Momax = zero(PFEI)
# # LD    = zero(PFEI)
# # CD    = zero(PFEI)
# # CDfuse    = zero(PFEI)
# # CDwing    = zero(PFEI)
# # CDi    = zero(PFEI)
# # CDnace    = zero(PFEI)
# # CL    = zero(PFEI)
# # WMTO  = zero(PFEI)
# # Wfuel = zero(PFEI)
# # Wwing = zero(PFEI)
# # Wtesys = zero(PFEI)
# # Wdfans = zero(PFEI)
# # Sref = zero(PFEI)
# # hcr = zero(PFEI)
# # Cost = zero(PFEI)

# # for (i,model) in enumerate(models)
# #     include(".././OldModels/"*model)
# #     Cost[i] = CostEst(parg, pare, parm, parpt, 500)/500

# #     PFEI[i]  = parm[imPFEI]
# #     nProp[i] = parg[igneng]
# #     Ssmax[i] = parg[igSsmax]
# #     Somax[i] = parg[igSomax]
# #     Msmax[i] = parg[igMsmax]
# #     Momax[i] = parg[igMomax]
# #     LD[i]    = para[iaCL, ipcruise1]/para[iaCD, ipcruise1]   
# #     CD[i]    = para[iaCD, ipcruise1]  
# #     CDfuse[i]    = para[iaCDfuse, ipcruise1]  
# #     CDwing[i]    = para[iaCDwing, ipcruise1]  
# #     CDi[i]    = para[iaCDi, ipcruise1]  
# #     CDnace[i]    = para[iaCDnace, ipcruise1]  
# #     CL[i]    = para[iaCL, ipcruise1]   
# #     WMTO[i]  = parg[igWMTO] 
# #     Wfuel[i] = parg[igWfuel]
# #     Wwing[i] = parg[igWwing]
# #     Sref[i] = parg[igS]
# #     hcr[i] = para[iaalt, ipcruise1]
# #     Wtesys[i] = parg[igWtesys]
# #     Wdfans[i] = parg[igneng]*(parg[igWfan] + parg[igWmot])
# # end

# # fig, ax = plt.subplots(4,1, figsize = (6, 8), dpi = 100)
# # ax[1].plot(nProp, PFEI, "o-k")
# # ax[1].set_ylabel("PFEI")
# # ax[2].plot(nProp, WMTO/9.81/1000, "o-k")
# # ax[2].set_ylabel("Max. TO weight [tonnes]")
# # ax[3].plot(nProp, Wfuel/9.81/1000, "o-k")
# # ax[3].set_ylabel("Fuel weight [tonnes]")
# # ax[4].plot(nProp, Cost/1e6,    "o-k")
# # ax[4].set_ylabel("Cost (million \$)")
# # ax[4].set_xlabel("Number of ducted fans")
# # plt.tight_layout()

# # fig, ax = plt.subplots(2,1, figsize = (8,5.5), dpi = 100)
# # ax[1].plot(nProp, Ssmax/1000, "o-k")
# # ax[1].plot(nProp, Somax/1000, "o-b")
# # ax[2].plot(nProp, Msmax/1000, "o-k")
# # ax[2].plot(nProp, Momax/1000, "o-b")

# # fig, ax = plt.subplots(5,1, figsize = (8,5.5), dpi = 100)
# # ax[1].plot(nProp, Sref.*CD, "o-k")
# # ax[2].plot(nProp, Sref.*CDfuse, "o-k")
# # ax[3].plot(nProp, Sref.*CDwing, "o-k")
# # ax[4].plot(nProp, Sref.*CDnace,    "o-k")
# # ax[5].plot(nProp, Sref.*CDi,    "o-k")

# # fig, ax = plt.subplots(5,1, figsize = (8,5.5), dpi = 100)
# # ax[1].plot(nProp, Sref.*CL, "o-k")
# # ax[1].set_ylabel("\$C_L\\times S_{ref}\$")
# # ax[2].plot(nProp, Sref.*CD, "o-k")
# # ax[2].set_ylabel("\$C_D\\times S_{ref}\$")
# # ax[3].plot(nProp, LD, "o-k")
# # ax[3].set_ylabel("\$ \\frac{L}{D} \$")
# # ax[4].plot(nProp, Sref, "o-k")
# # ax[4].set_ylabel("\$ S_{ref}\$")
# # ax[5].plot(nProp, hcr, "o-k")

# # fig, ax = plt.subplots(3,1, figsize = (8,5.5), dpi = 100)
# # ax[1].plot(nProp, Wwing./9.81./1000, "o-k")
# # ax[2].plot(nProp, Wdfans./9.81./1000, "o-k")
# # ax[3].plot(nProp, Wtesys./9.81./1000, "o-k")
