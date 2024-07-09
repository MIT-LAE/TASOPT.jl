using TASOPT
using Printf
include(joinpath(TASOPT.__TASOPTroot__, "./misc/index.inc"))
# using Base: SignedMultiplicativeInverse
# using NLopt: G_MLSL_LDS, GN_MLSL_LDS, GN_CRS2_LM, GN_DIRECT_L
# using BenchmarkTools
# using Printf

# ## Testing build up of wsize
# using StaticArrays
# using Profile, UnicodePlots
# using PyPlot

# if Sys.iswindows()
#     pygui(true)
# end

# using Dates

# include("index.inc")
# include("atmos.jl")
# include("fuseW.jl")

# include("fusebl.jl")
# include("axisol.jl")
# include("blax.jl")
# include("blsys.jl")
# include("wingsc.jl")
# include("surfcm.jl")
# include("surfdx.jl")
# include("surfw.jl")
# include("wingpo.jl")
# include("tailpo.jl")

# include("airtable2.jl")
# include("airfun2.jl")

# include("balance.jl")
# include("surfcd.jl")
# include("trefftz.jl")
# include("cdsum.jl")

# include("PMSM.jl")  # Motor/generator functions
# include("PMSM.inc") # Motor/generator properties array

# include("TF_NPSS_functions.jl") 

# # include("TE_NPSS_functions.jl") 
# include("PT.inc")
# include("propsys.jl")

# include("hydrogen.jl")

# include("tankWmech.jl")
# include("tankWthermal.jl")
# include("tanksize.jl")

# include("wsize.jl")
# include("woper.jl")
# include("mission.jl")
# include("LTO_extra_pts.jl")

# # include("737.jl")

# # WE STILL NEED THE CONSTANTS DEFINED IN 737.jl
# nmisx = 5
# pari = zeros(Int64, iitotal)
# parg = zeros(Float64, igtotal)
# parm = zeros(Float64, (imtotal, nmisx))
# para = zeros(Float64, (iatotal, iptotal, nmisx))
# pare = zeros(Float64, (ietotal, iptotal, nmisx))

# ft_to_m = 0.3048
# in_to_m = 0.0254
# nmi_to_m = 1852.0
# deg_to_rad = π/180.0
# lbf_to_N = 4.448222
# kts_to_mps = 0.51444
# hp_to_W    = 745.7

# pax = 180
# seat_pitch = 30.0 * in_to_m  
# seat_width = 19.0 * in_to_m
# aisle_halfwidth = 10.0 * in_to_m # per CFR § 25.815 

# # include("input_fwdTank.jl")
# include("outputs.jl")
# # include("outputs1.jl")
# include("savemodel.jl")
# include("update_fuse.jl")
# include("cost_est.jl")
# include("printBADA.jl")
# include("odperformance.jl")
# include("AircraftDeck.jl")
# include("LTO.jl")

# const gee = 9.81
# TSL, pSL, ρSL, aSL, μSL = atmos(0.0)
# RSL =  pSL/(ρSL * TSL)
# const gamSL = 1.4
# const cpSL  = 1004.0
# # temporary
# μAir = 1.78e-5 
# ρAir = ρSL

# time_writing = 0.0
# time_run_NPSS = 0.0

# # f = stdout

# initwgt = 0
# saveOD = false # Optimize first, and call mdl run odperf then
# optimize = false
# track_fig = nothing
# modify_thrustreq = false

function run_wsize(ac, iter, initwgt, Ldebug, printiter, saveOD, optimize, nmisx)
    # global time_writing = 0.0
    # global time_run_NPSS = 0.0
    # parpt[ipt_time_NPSS] = 0.0
    # parpt[ipt_calls_NPSS] = 0

    kmdes = 1
    # println(parg[igWwing])
    # Ldebug && println("Max weight iterations = $iter")
    # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS = wsize(pari, parg, view(parm,:,kmdes), view(para,:,:,kmdes), view(pare, :,:,kmdes),
    #                 iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, saveOD, modify_thrustreq, false)
    # @benchmark PowerTrain(0.0, 0.8, 25.0e3*2,0.0, 0.0, parg, parpt, parmot, pargen)

    # println("Time netNPSS       = $(parpt[ipt_time_NPSS])")
    # println("Time writing       = $time_writing")
    # println("Time runnning NPSS = $time_run_NPSS")
    size_aircraft!(ac)

    local RMS_error
    
    for km = 2 : nmisx 

        for ip = 1: iptotal
            for ia = 1:iatotal
                ac.para[ia,ip,km] = ac.para[ia,ip,kmdes]
            end
            for ie = 1 : ietotal
                ac.pare[ie,ip,km] = ac.pare[ie,ip,kmdes]
            end
        end
        
        initwgt = 1
        initeng = 1

        # if km == nmisx
        #     endNPSS_flag = true
        # else
        #     endNPSS_flag = false
        # end

        # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS, RMS_error = woper(pari, parg,view(parm, :,km),view(para,:,:,km),view(pare, :,:,km), view(para,:,:,kmdes),view(pare, :,:,kmdes),
        #     iter, initeng,  NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS, endNPSS_flag)

        TASOPT.woper(ac, kmdes, itermax = iter, initeng = initeng, saveOffDesign = saveOD)

    end

    f = open("temp.results", "w")

    for km = 1:nmisx
        @printf(f, "\nPFEI = %.5f J/Nm for mission range %5.3f NM\n", ac.parm[imPFEI,km], 0.000539957*ac.parm[imRange,km])
    end

    @printf(f, "\nPFEI_weighted = %.5f J/Nm", sum(ac.parm[imwOpt,:].*ac.parm[imPFEI,:]))

    # weight_buildup(parg, io = f)
    # aero(parg, para[:,:,kmdes], io = f)
    # geometry(parg, io = f)

    TASOPT.weight_buildup(ac, io = f) 
    TASOPT.aero(ac, io = f)
    TASOPT.geometry(ac, io = f)

    close(f)

    if (ac.parg[8] >= ac.parg[7] && optimize == false)
        ############################################################
        # sleep(10) # Wait 10 second for NPSS to shutdown and restart

        kmdes = 1
        # Run once more for saveOD
        println("Verified that Wfmax is enough to hold Wfuel, RUNNIGN ONCE MORE FOR saveOD")
        # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS = wsize(pari, parg, view(parm, :,kmdes), view(para,:,:,kmdes), view(pare, :,:,kmdes),
        # iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, true, modify_thrustreq, false)
        TASOPT.woper(ac, kmdes, itermax = iter, initeng = initeng, saveOffDesign = saveOD)
        for km = 2 : nmisx 

            for ip = 1: iptotal
                for ia = 1:iatotal
                    ac.para[ia,ip,km] = ac.para[ia,ip,kmdes]
                end
                for ie = 1 : ietotal
                    ac.pare[ie,ip,km] = ac.pare[ie,ip,kmdes]
                end
            end
            
            initwgt = 1
            initeng = 1
    
            # if km == nmisx
            #     endNPSS_flag = true
            # else
            #     endNPSS_flag = false
            # end
    
            # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS = woper(pari, parg,view(parm, :,km),view(para,:,:,km),view(pare, :,:,km), view(para,:,:,kmdes),view(pare, :,:,kmdes),
            #     iter, initeng,  NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS, endNPSS_flag)
            TASOPT.woper(ac, kmdes, itermax = iter, initeng = initeng, saveOffDesign = saveOD)
    
        end

    end

    # reduce range if it Wfmax < Wfuel
    while ac.parg[8] < ac.parg[7]
        ############################################################
        # sleep(10) # Wait 10 second for NPSS to shutdown and restart
        
        println("Wfmax:" *string(ac.parg[8]) *", Wfuel:" *string(ac.parg[7]))
        #f = open("temp_FPRchanged" * "_range_" * string(parm[2]*0.99/1000*0.539957) * ".results", "w")
        
        println("Wfmax < Wfuel\nReducing mission range to", string(ac.parm[2,kmdes]*0.998/1000*0.539957), " nm, RERUNNING!")
        ac.parm[2,kmdes] = ac.parm[2,kmdes] * 0.998 # reduce range
        
        # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS = wsize(pari, parg, view(parm, :,kmdes), view(para,:,:,kmdes), view(pare, :,:,kmdes),
        #                 iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, saveOD, modify_thrustreq, false)
        TASOPT.wsize(ac, itermax = iter, initwgt = initwgt,
        Ldebug = Ldebug, printiter = printiter,
        saveODperf = saveOD)
        for km = 2 : nmisx 

            for ip = 1: iptotal
                for ia = 1:iatotal
                    ac.para[ia,ip,km] = ac.para[ia,ip,kmdes]
                end
                for ie = 1 : ietotal
                    ac.pare[ie,ip,km] = ac.pare[ie,ip,kmdes]
                end
            end
            
            initwgt = 1
            initeng = 1
    
            # if km == nmisx
            #     endNPSS_flag = true
            # else
            #     endNPSS_flag = false
            # end
    
            # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS = woper(pari, parg,view(parm, :,km),view(para,:,:,km),view(pare, :,:,km), view(para,:,:,kmdes),view(pare, :,:,kmdes),
            #     iter, initeng,  NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS, endNPSS_flag)
            TASOPT.woper(ac, kmdes, itermax = iter, initeng = initeng, saveOffDesign = saveOD)
    
        end

        if ac.parg[7] < ac.parg[8]
            # Run once more for saveOD
            # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS = wsize(pari, parg, view(parm, :,kmdes), view(para,:,:,kmdes), view(pare, :,:,kmdes),
            # iter, 0.5, 0.9, 0.5, initwgt, 1, 1, Ldebug, printiter, true, modify_thrustreq, false)
            TASOPT.wsize(ac, itermax = iter, initwgt = initwgt,
            Ldebug = Ldebug, printiter = printiter,
            saveODperf = saveOD)
        
            for km = 2 : nmisx 

                for ip = 1: iptotal
                    for ia = 1:iatotal
                        ac.para[ia,ip,km] = ac.para[ia,ip,kmdes]
                    end
                    for ie = 1 : ietotal
                        ac.pare[ie,ip,km] = ac.pare[ie,ip,kmdes]
                    end
                end
                
                initwgt = 1
                initeng = 1
        
                # if km == nmisx
                #     endNPSS_flag = true
                # else
                #     endNPSS_flag = false
                # end
        
                # NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS = woper(pari, parg,view(parm, :,km),view(para,:,:,km),view(pare, :,:,km), view(para,:,:,kmdes),view(pare, :,:,kmdes),
                #     iter, initeng,  NPSS_TS, NPSS_Fan, NPSS_AftFan, NPSS_PT, NPSS, endNPSS_flag)
                TASOPT.woper(ac, kmdes, itermax = iter, initeng = initeng, saveOffDesign = saveOD)
        
            end

            f = open("temp.results", "w")

            for km = 1:nmisx
                @printf(f, "\nPFEI = %.5f J/Nm for mission range %5.3f NM\n", ac.parm[imPFEI,km], 0.000539957*ac.parm[imRange,km])
            end
        
            @printf(f, "\nPFEI_weighted = %.5f J/Nm", sum(ac.parm[imwOpt,:].*ac.parm[imPFEI,:]))
        
            # weight_buildup(parg, io = f)
            # aero(parg, para[:,:,kmdes], io = f)
            # geometry(parg, io = f)

            TASOPT.weight_buildup(ac, io = f) 
            TASOPT.aero(ac, io = f, mi = kmdes)
            TASOPT.geometry(ac, io = f)
        
            close(f)

        end

    end

    # global track_fig = stickfig(parg, pari,  parm; ax = track_fig)
    # global track_fig = plot_details(parg, pari, para, parm; ax = track_fig)

    # figname = "B737-800_temp"
    # savemodel("./Models/"*figname*".mdl", pari, parg, parm, para, pare, parpt, parmot, pargen)

end

# time_wsize = @elapsed run_wsize(100, 0, false, true, saveOD, false)

# @profview run_wsize(25, 0, false, true)
# println("Wsize time = $time_wsize s")
# println("NPSS Time = $(parpt[ipt_time_NPSS]) s")
# println("NPSS Calls = $(parpt[ipt_calls_NPSS])")

# include("LoadModel.jl")

# # # Final aircraft:
# include("./Models/B737-800_optimized.mdl")

# time_wsize = @elapsed run_wsize(100, 1, false, false, saveOD)

# savedir = "./Figures/"
# figname = "B737-800"

# global track_fig = plot_details(parg, pari, para[:,:,1], pare[:,:,1], parm[:,1]; ax = track_fig)

# plt.savefig(savedir * figname * ".jpg", dpi = 300)
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/EEI_input.toml"))
size_aircraft!(ac)

mofftpax = 0.00633
mofftmMTO = 0.0 * 0.001
Pofftpax  = 200.0
PofftmMTO = 1.8

ac.parg[igmofWpay] = mofftpax / (230.0 * lbf_to_N)  #This later gets multiplied by Wpay, resulting in mofftpax*pax
ac.parg[igmofWMTO] = mofftmMTO  / 9.81

ac.parg[igPofWpay] = Pofftpax / (230.0 * lbf_to_N)
ac.parg[igPofWMTO] = PofftmMTO / 9.81

saveOD = true
nmisx = 5

run_wsize(ac, 200, 1, false, true, saveOD, false, nmisx)

# load_aircraft("B737-800_optimized.mdl", "./Models", saveOD, modify_thrustreq, run = true)
#load_aircraft("B737-800_Mach.mdl", "./Models", saveOD, modify_thrustreq, run = true)

# if optimize == false
#     savedir = "./Figures/"
#     figname = "B737-800"

#     global track_fig = plot_details(parg, pari, para[:,:,1], pare[:,:,1], parm[:,1]; ax = track_fig)

#     plt.savefig(savedir * figname * ".jpg", dpi = 300)
#     high_res_airplane_plot(parg, pari, parm; ax = nothing, label_fs = 11, save_name = figname * "_highres")

# end