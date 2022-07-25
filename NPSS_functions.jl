"""
A trivial wrapper that runs a NPSS batch file `bat_file` at the specifed `dir`

>Note: 
>Remember to use `/` as separators
"""
function NPSS_run(dir, bat_file)
    global time_run_NPSS += @elapsed run(`cmd /c cd $dir '&&' $bat_file `)
    return nothing
end

"""
This function starts up and returns an NPSS process that can then be written to
"""
function startNPSS(dir, bat_file)
    if Sys.iswindows()
        NPSS = open(`cmd /c cd $dir '&&' $bat_file `, "w+")
    elseif Sys.islinux()
        cd(dir) # quick fix cant figure out how to cd in open()
        NPSS = open(`bash $bat_file `, "w+")
        cd("../")
    end
    return NPSS
end

"""
Ends NPSS process that can then be written to
"""
function endNPSS(NPSS)
    write(NPSS, "999 \\n")
    if(NPSS.exitcode == 0)
        close(NPSS)
    end
end


function NPSS_TEsys(NPSS, alt_in, MN_in, Fn, Tt41,
     πaftfan, πfan, Kinlaft, Φinlaft, Kinl, Φinl, first, parg, parpt)

    nfans   = parpt[ipt_nfan]
    ngen    = parpt[ipt_ngen]
    nTshaft = parpt[ipt_nTshaft]
    FnSplit = parpt[ipt_Fnsplit]

    πLPC = parpt[ipt_piLPC]
    πHPC = parpt[ipt_piHPC]
    cpsi = parpt[ipt_cpsi]
    w    = parpt[ipt_wcat]
    lcat = parpt[ipt_lcat]
    deNOx= parpt[ipt_deNOx]
    LHV  = parg[igLHVfuel]

    rSnace = parg[igrSnace]
    fpylon = parg[igfpylon]
    feadd  = parg[igfeadd]

    input_string = "111 "*
    "first = $first;" *
    "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
    "Eng.Amb.MN_in  = $MN_in ;" *

    "Fn_target         = $(Fn/nTshaft);"*
    "Tt41              = $Tt41 ;" *
    "FnSplitTarget     = $FnSplit;" *
    "Eng.n_poddedprops = $(nfans/nTshaft);" *
    "deNOx_target      = $deNOx ;"*

    "Eng.FusEng.LHV = $(LHV*429.923); "*
    "Eng.CmpF.PRdes = $(πaftfan) ;" *
    "Eng.CmpL.PRdes = $(πLPC) ;" *
    "Eng.CmpH.PRdes = $(πHPC) ;" *

    "Eng.ShL.Nmech = 13000.0;"*
    "Eng.ShF.Nmech = 3400.0;"*
    "Eng.PodProp.MotorShaft.Nmech = 12000.0;"*

    "Eng.PodProp.Fan.PRdes  = $πfan;" *
    
    "Eng.InEng.Kinl         = $Kinlaft;" *
    "Eng.PodProp.InEng.Kinl = $Kinl;" *
    "Eng.Phiinl             = $Φinlaft;" *
    "Eng.PodProp.Phiinl     = $Φinl;" *

    "Eng.PCEC.l = $lcat ;"*
    "Eng.PCEC.w = $w ;"*
    "Eng.PCEC.cpsi = $cpsi ;"*
    
    "Eng.rSnace = $(rSnace*3/4) ;" *
    "Eng.PodProp.rSnace = $(rSnace) ;" *
    "Eng.PodProp.fpylon = $(fpylon) ;" *
    "Eng.feadd = $(feadd) ;" *
    "\n"

    write(NPSS, input_string)
    # write(stdout, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        mdotf      = parse(Float64, out[2] )
        deNOx      = parse(Float64, out[3] )
        mcat       = parse(Float64, out[4] )
        EINOx1     = parse(Float64, out[5] )
        EINOx2     = EINOx1*(1-deNOx)
        mdot3      = parse(Float64, out[6] )
        Tt3        = parse(Float64, out[7] )
        OPR        = parse(Float64, out[8])
        Wc3        = parse(Float64, out[9])
        LHV        = parse(Float64, out[10])       
        W_in       = parse(Float64, out[11])

        SPmot      = parse(Float64, out[12])
        SPinv      = parse(Float64, out[13])
        SPgen      = parse(Float64, out[14])
        SPtshaft   = parse(Float64, out[15])   
    
        Pshaft_fan = parse(Float64, out[16])
        Pelec_mot  = parse(Float64, out[17])
        Pelec_gen  = parse(Float64, out[18])
        Pshaft_gen = parse(Float64, out[19])

        mgen       = parse(Float64, out[20])
        mmot       = parse(Float64, out[21])

        Wtshaft    = parse(Float64, out[22])
        Waftfan    = parse(Float64, out[23])
        Wfan       = parse(Float64, out[24])
        
        Daftfan    = parse(Float64, out[25])
        Dfan       = parse(Float64, out[26])

        Wgb        = parse(Float64, out[27])

        FAR        = parse(Float64, out[28])

        Saftnace1  = parse(Float64, out[29])
        Snace1     = parse(Float64, out[30])

        Tt41       = parse(Float64, out[31])
        EGT        = parse(Float64, out[32])

        ηmot      = parse(Float64, out[33])
        ηinv      = parse(Float64, out[34])
        ηcable    = 0.0
        ηgen      = parse(Float64, out[35])
        ηtshaft   = parse(Float64, out[36]) 

        Wgbpods   = parse(Float64, out[37])

        Vcable    = parse(Float64, out[38])

        Pshaft_aftfan = parse(Float64, out[39])

    end
    mdotf_tot = mdotf*nTshaft
    Wcat = mcat*gee*1.5
    Wgen = mgen*gee
    Wmot = mmot*gee

    parpt[ipt_Ptshaft] = Pshaft_aftfan + ngen*Pshaft_gen/nTshaft #Add this for cost model


    Wpowertrain = 0.0
    xWpowertrain = 0.0
    parg[igdaftfan] = Daftfan
    parg[igWaftfan] = Waftfan
  
     naftfan = nTshaft
     Wpowertrain += (Waftfan + Wgb)*naftfan #Aft fan and gearbox per fan
    xWpowertrain += (Waftfan + Wgb)*naftfan*parg[igxtshaft]
     
     parg[igdfan] = Dfan
     parg[igWfan] = Wfan
  
     Wpowertrain += Wfan*nfans
    xWpowertrain += Wfan*nfans*parg[igxfan]
  
     parg[igWmot] = Wmot
     parpt[ipt_Wmot]    = Wmot
     parg[igWfanGB]     = Wgbpods
     parg[igWaftfanGB]  = Wgb
  
     parg[igWgen] = Wgen
     parpt[ipt_Wgen] = Wgen
     Wpowertrain += Wgen*ngen
    xWpowertrain += Wgen*ngen*parg[igxgen]  

    Wrect = Pelec_gen/SPinv * gee
    parg[igWrect] = Wrect
    Wpowertrain += Wrect*ngen # Add to total powertrain weight
   xWpowertrain += Wrect*ngen*parg[igxgen]    
  
     Wpowertrain += (Wmot+Wgbpods)*nfans  # Add to total powertrain weight
     xWpowertrain += (Wmot+Wgbpods)*nfans*parg[igxmot]  
  
  
     Winv = Pelec_mot/SPinv * gee
     parg[igWinv] = Winv
     parpt[ipt_Winv] = Winv
     Wpowertrain += Winv*nfans # Add to total powertrain weight
    xWpowertrain += Winv*nfans*parg[igxinv]    
   
    lcable = (parg[igxgen] - parg[igxwbox]) + parg[igyeout]
    parpt[ipt_lcable] = lcable
    ηcable, Wcable = cable(Pelec_mot, 2*Vcable, lcable, parpt) 
     parg[igWcables] = Wcable*nfans # multiply by number of fans
     parpt[ipt_Wcables] = Wcable*nfans
     Wpowertrain += Wcable*nfans # Add to total powertrain weight
#    xWpowertrain += Wcable*parg[igxcable] #need to add x of cable 


   parg[igWcat] = Wcat
   parpt[ipt_Wcatalyst] = Wcat
   Wpowertrain += Wcat*nTshaft
  xWpowertrain += Wcat*nTshaft*parg[igxtshaft]

  feadd = parg[igfeadd] # Added weight fraction for engine accessories, fuel systems etc.
  Wtshaft = Wtshaft*(1+feadd)
        
  parg[igWtshaft] = Wtshaft
  parpt[ipt_Wtshaft] = Wtshaft
  Wpowertrain += Wtshaft*nTshaft
 xWpowertrain += Wtshaft*nTshaft*parg[igxtshaft]

    Hrejmot = nfans*(Pelec_mot - Pshaft_fan)
    Hrejinv = nfans*(Pelec_mot*(1-ηinv))
    Hrejcab = ngen*(Pelec_gen*(1-ηcable))
    Hrejgen = ngen*(Pshaft_gen - Pelec_gen)
    Hrejrect= ngen*(Pelec_gen*(1-ηinv))
    Hrejtot = Hrejmot + Hrejinv + Hrejcab + Hrejgen + Hrejrect

    if LHV<100
        Qtms = Hrejtot
    else
        Qtms = Hrejmot + Hrejinv + Hrejcab
    end
    P_mTMS = 8.0*hp_to_W/(1/2.205)  # convert 8 hp/lbm to W/kg. Value based on MIT NASA LEARN report Table A8
    Wtms = Qtms/P_mTMS * 9.81
    parg[igWtms] = Wtms
    Wpowertrain += Wtms
   xWpowertrain += Wtms*(parg[igxtshaft] + parg[igxmot])/2
    
    parg[ igWtesys] = Wpowertrain
    parg[igxWtesys] = xWpowertrain



    ηpt = [0.0, ηmot, ηinv, ηcable, ηgen, ηtshaft]
    Hpt = [Hrejmot, Hrejinv, Hrejcab, Hrejgen+Hrejrect, Hrejtot]
    Ppt = [Pshaft_fan, Pelec_mot, Pelec_gen, Pshaft_gen]
    SPpt = [SPmot,SPinv, SPgen, SPtshaft]

    heatexcess = 0.0

    return NPSS_success,ηpt, SPpt, Ppt, Hpt, heatexcess, 
            mdotf_tot, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, EGT,
            Snace1, Saftnace1

end

""" 
Off Des NPSS_TEsys
"""
function NPSS_TEsysOD(NPSS, alt_in, MN_in, Fn, Tt41,
     Kinlaft, Φinlaft, Kinl, Φinl, first, parg::Array{Float64, 1}, 
     parpt::Array{Union{Float64, Int64},1}, pare, ip::Int64)
    
    #  println("Starting Off des at ip = $ip")
   nTshaft = parpt[ipt_nTshaft]
   FnSplit = parpt[ipt_Fnsplit]
   nfans   = parpt[ipt_nfan]
   ngen    = parpt[ipt_ngen]

   if Fn == 0.0
    mode = "222" #Tt41 specified
   else
    mode = "333" #Fn specified
   end

   input_string = "$(mode) "*
   "first = $first;" *
   "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
   "Eng.Amb.MN_in  = $MN_in ;" *

   "Fn_target         = $(Fn/nTshaft);"*
   "Tt41              = $Tt41 ;" *
   "FnSplitTarget     = $FnSplit;" *

   "Eng.InEng.Kinl         = $Kinlaft;" *
   "Eng.PodProp.InEng.Kinl = $Kinl;" *
   "Eng.Phiinl             = $Φinlaft;" *
   "Eng.PodProp.Phiinl     = $Φinl;" *

#    "solver.independentValues = {$(pare[iegsFnsplit , ip]) ,"*   
#    "$(pare[iegsMotShP  , ip]) ,"*     
#    "$(pare[iegsmdotf   , ip]) ,"*       
#    "$(pare[iegsWin     , ip]) ,"*   
#    "$(pare[iegsRlineF  , ip]) ,"*               
#    "$(pare[iegsBPR     , ip]) ,"*    
#    "$(pare[iegsRlinelc , ip]) ,"*               
#    "$(pare[iegsRlinehc , ip]) ,"*               
#    "$(pare[iegsPRhtrb  , ip]) ,"*             
#    "$(pare[iegsPRltrb  , ip]) ,"*             
#    "$(pare[iegsNmechH  , ip]) ,"*     
#    "$(pare[iegsGBtrq   , ip]) ,"*  
#    "$(pare[iegsNmechL  , ip]) ,"*     
#    "$(pare[iegsNmechF  , ip]) ,"*     
#    "$(pare[iegsPodWin  , ip]) ,"*           
#    "$(pare[iegsPodRlineF  , ip]) ,"*                    
#    "$(pare[iegsPodGBtrq   , ip]) ,"*        
#    "$(pare[iegsPodMotNmech, ip]) ,"*                  
#    "$(pare[iegsPodFanNmech, ip]) };"*       

   "\n"

   write(NPSS, input_string)
#    write(stdout, input_string)

   out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
   NPSS_success = parse(Float64, out[1] )
   
   if length(out) > 1
       mdotf      = parse(Float64, out[2] )
       deNOx      = parse(Float64, out[3] )
       Ftotal     = parse(Float64, out[4] )
       EINOx1     = parse(Float64, out[5] )
       EINOx2     = EINOx1*(1-deNOx)
       mdot3      = parse(Float64, out[6] )
       Tt3        = parse(Float64, out[7] )
       OPR        = parse(Float64, out[8])
       Wc3        = parse(Float64, out[9])
       W_in       = parse(Float64, out[10])
   
       Pshaft_fan = parse(Float64, out[11])
       Pelec_mot  = parse(Float64, out[12])
       Pelec_gen  = parse(Float64, out[13])
       Pshaft_gen = parse(Float64, out[14])

       FAR        = parse(Float64, out[15])

       Tt41       = parse(Float64, out[16])
       EGT        = parse(Float64, out[17])

       ηmot      = parse(Float64, out[18])
       ηinv      = parse(Float64, out[19])
       ηcable    = 0.0
       ηgen      = parse(Float64, out[20])
       ηtshaft   = parse(Float64, out[21]) 

       pare[iegsFnsplit , ip] = parse(Float64, out[22])
       pare[iegsMotShP  , ip] = parse(Float64, out[23])
       pare[iegsmdotf   , ip] = parse(Float64, out[24])
       pare[iegsWin     , ip] = parse(Float64, out[25])
       pare[iegsRlineF  , ip] = parse(Float64, out[26])
       pare[iegsBPR     , ip] = parse(Float64, out[27])
       pare[iegsRlinelc , ip] = parse(Float64, out[28])
       pare[iegsRlinehc , ip] = parse(Float64, out[29])
       pare[iegsPRhtrb  , ip] = parse(Float64, out[30])
       pare[iegsPRltrb  , ip] = parse(Float64, out[31])
       pare[iegsNmechH  , ip] = parse(Float64, out[32])
       pare[iegsGBtrq   , ip] = parse(Float64, out[33])
       pare[iegsNmechL  , ip] = parse(Float64, out[34])
       pare[iegsNmechF  , ip] = parse(Float64, out[35])
       pare[iegsPodWin  , ip] = parse(Float64, out[36])
       pare[iegsPodRlineF  , ip] = parse(Float64, out[37])
       pare[iegsPodGBtrq   , ip] = parse(Float64, out[38])
       pare[iegsPodMotNmech, ip] = parse(Float64, out[39])
       pare[iegsPodFanNmech, ip] = parse(Float64, out[40])
       pare[ieTmet1, ip] = parse(Float64, out[41])

   end
   mdotf_tot = mdotf*nTshaft
   Ftotal    = Ftotal*nTshaft
   ηcable = cable(Pelec_gen, parpt) 

   Hrejmot = nfans*(Pelec_mot - Pshaft_fan)
   Hrejinv = nfans*(Pelec_mot*(1-ηinv))
   Hrejcab = nfans*(Pelec_mot*(1-ηcable))
   Hrejgen = ngen*(Pshaft_gen - Pelec_gen)
   Hrejrect= ngen*(Pelec_gen*(1-ηinv))
   Hrejtot = Hrejmot + Hrejinv + Hrejcab + Hrejgen + Hrejrect

   ηpt = [ηmot, ηinv, ηcable, ηgen, ηtshaft]
   Hpt = [Hrejmot, Hrejinv, Hrejcab, Hrejgen+Hrejrect, Hrejtot]
   Ppt = [ Pshaft_fan, Pelec_mot, Pelec_gen, Pshaft_gen]
   heatexcess = 0.0

   return NPSS_success, Ftotal, ηpt, Ppt, Hpt, heatexcess, 
           mdotf_tot, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT

end

"""
# Design mode

Writes an input file for NPSS Turboshaft model 

    - Abmient altitude and mach number
    - Specifies πᵢ ∀ i ∈ {HPC, LPC}
    - Shaft power demand, Tt41 
    - DeNOx target
    - SCR paramters -> CPSI, w, l

"""
function NPSS_TShaft_input(NPSS, alt_in, MN_in, 
    SHP_dmd, Tt41, 
    LPC_PR, HPC_PR,
    cpsi, w, lcat, deNOx, first, LHV)

    input_string = "111 "*
    "first = $first;" *
    "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
    "Eng.Amb.MN_in  = $MN_in ;" *
   
    "Eng.ShP.HPX = $(SHP_dmd/hp_to_W) ;" * #convert to hp before sending to NPSS
    "Tt41 = $Tt41 ;" *
    "deNOx_target = $deNOx ;"*

    "Eng.FusEng.LHV = $(LHV*429.923); "*
    "Eng.CmpL.PRdes = $(LPC_PR) ;" *
    "Eng.CmpH.PRdes = $(HPC_PR) ;" *

    "Eng.PCEC.l = $lcat ;"*
    "Eng.PCEC.w = $w ;"*
    "Eng.PCEC.cpsi = $cpsi ;"*
    "\n"

    write(NPSS, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        ηtherm = parse(Float64, out[2] )
        mdotf  = parse(Float64, out[3] )
        BSFC   = parse(Float64, out[4] )
        deNOx  = parse(Float64, out[5] )
        mcat   = parse(Float64, out[6] )
        EINOx1 = parse(Float64, out[7] )
        EINOx2 = EINOx1*(1-deNOx)
        mdot   = parse(Float64, out[8] )
        Tt3    = parse(Float64, out[9] )
        OPR    = parse(Float64, out[10])
        Wc3    = parse(Float64, out[11])
        LHV    = parse(Float64, out[12])
        W_in    = parse(Float64, out[13])
    end

    # println(W_in, "lbm/s; If m/mdot^1.2 = 38.9   -> Wtshaft = ", (W_in/2.205)^1.2 * 38.9   * gee, " SP [kW/kg] = ", SHP_dmd/1000/((W_in/2.205)^1.2 * 38.9) ) # 38.9 kg/(kg/s) from Hall et al. https://arc.aiaa.org/doi/pdf/10.2514/6.2018-3973
    # println(W_in, "lbm/s; If m/mdot^1.2 = 45.605 -> Wtshaft = ", (W_in/2.205)^1.2 * 45.605 * gee, " SP [kW/kg] = ", SHP_dmd/1000/((W_in/2.205)^1.2 * 45.605) ) # 38.9 kg/(kg/s) from Dowdle et al
    Wtshaft = (W_in/2.205)^1.2 * 45.605 * gee
    
    return NPSS_success, ηtherm, mdotf, BSFC, deNOx, mcat, EINOx1, EINOx2, mdot, Tt3, OPR, Wc3, Wtshaft

end

"""
OffdesMode

Writes an input file for NPSS Turboshaft model in off-des conditions

    - Abmient altitude [in m]
    - Flight mach number
    - Tt41 
    - the flag first sets whether NPSS should restart or use previous value

"""
function NPSS_TShaft_run(NPSS, alt_in, MN_in, 
                            Tt41, N2_dmd, first)

    write(NPSS, "222 Tt41=$Tt41; N2_dmd = $N2_dmd; Eng.Amb.alt_in=$(alt_in/ft_to_m); Eng.Amb.MN_in = $MN_in;first=$first;\n")
    # write(stdout, "222 Tt41=$Tt41; N2_dmd = $N2_dmd; Eng.Amb.alt_in=$(alt_in/ft_to_m); Eng.Amb.MN_in = $MN_in;first=$first;\n")
    
    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        ShP    = parse(Float64, out[2] ) #NPSS returns in W
        ηtherm = parse(Float64, out[3] )
        mdotf  = parse(Float64, out[4] )
        BSFC   = parse(Float64, out[5] )
        deNOx  = parse(Float64, out[6] )
        EGT    = parse(Float64, out[7] )
        Tt3    = parse(Float64, out[8] )
        W3     = parse(Float64, out[9] )
        EINOx1 = parse(Float64, out[10] )
        EINOx2 = EINOx1*(1-deNOx)
        FAR    = parse(Float64, out[11])
    end
    
    return NPSS_success, ShP, ηtherm, mdotf, BSFC, deNOx, EGT, Tt3, W3, EINOx1, EINOx2, FAR
end

function NPSS_TShaft_run2(NPSS, alt_in, MN_in, 
    SHP_dmd, N2_dmd, first)

write(NPSS, "333 Eng.ShP.HPX=$(SHP_dmd/hp_to_W); N2_dmd = $N2_dmd; Eng.Amb.alt_in=$(alt_in/ft_to_m); Eng.Amb.MN_in = $MN_in; first=$first;\n")

out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
NPSS_success = parse(Float64, out[1] )

if length(out) > 1
ShP    = parse(Float64, out[2] )
ηtherm = parse(Float64, out[3] )
mdotf  = parse(Float64, out[4] )
BSFC   = parse(Float64, out[5] )
deNOx  = parse(Float64, out[6] )
end

return NPSS_success, ShP, ηtherm, mdotf, BSFC, deNOx
end
"""
Fan - Design mode:

Writes the desired conditions to the external NPSS process
"""
function runNPSS_Fan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Fn::Float64,
                    Kinl::Float64, Φinl::Float64, 
                    πfan::Float64, first)
    input_string = "111 "*
                "first = $first;" *
                "DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
                "DuctedFan.Amb.MN_in  = $MN_in ;" *
                "Fn_target = $Fn ;" *
                "DuctedFan.InEng.Kinl = $Kinl;" *
                "DuctedFan.Phiinl     = $Φinl;" *
                "DuctedFan.Fan.PRdes = $πfan;" *
                "\n"

    write(NPSS, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        Dfan = parse(Float64, out[2] )
        Power  = parse(Float64, out[3] )
        Torque = parse(Float64, out[4] )
        Nmech  = parse(Float64, out[5] )
        Mtip   = parse(Float64, out[6] )
        ηprop  = parse(Float64, out[7] )
        ηDF    = parse(Float64, out[8] )
        ANoz   = parse(Float64, out[9] )
    end

    return NPSS_success, Dfan, Power, Torque, Nmech, Mtip, ηprop, ηDF, ANoz
end

"""
Fan - Off-design mode:

Writes the desired conditions to the external NPSS process
"""
function runNPSS_Fan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Pin::Float64,
                    Kinl::Float64, Φinl::Float64, first)
    input_string = "222 "*
                "first = $first;" *
                "DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
                "DuctedFan.Amb.MN_in  = $MN_in ;" *
                "ShP_input = $(-Pin/745.7) ;" *
                "DuctedFan.InEng.Kinl = $Kinl;" *
                "DuctedFan.Phiinl     = $Φinl;" *
                "\n"

    write(NPSS, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        Fn_N = parse(Float64, out[2] )
        Power  = parse(Float64, out[3] )
        Torque = parse(Float64, out[4] )
        Nmech  = parse(Float64, out[5] )
        Mtip   = parse(Float64, out[6] )
        ηprop  = parse(Float64, out[7] )
        ηDF    = parse(Float64, out[8] )
    end

    return NPSS_success, Fn_N, Power, Torque, Nmech, Mtip, ηprop, ηDF
end
function runNPSS_Fan2(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Fn::Float64,
    Kinl::Float64, Φinl::Float64, first)

    input_string = "333 "*
"first = $first;" *
"DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
"DuctedFan.Amb.MN_in  = $MN_in ;" *
"Fn_target = $Fn ;" *
"DuctedFan.InEng.Kinl = $Kinl;" *
"DuctedFan.Phiinl     = $Φinl;" *
"\n"

write(NPSS, input_string)

out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
NPSS_success = parse(Float64, out[1] )

if length(out) > 1
Fn_N = parse(Float64, out[2] )
Power  = parse(Float64, out[3] )
Torque = parse(Float64, out[4] )
Nmech  = parse(Float64, out[5] )
Mtip   = parse(Float64, out[6] )
ηprop  = parse(Float64, out[7] )
ηDF    = parse(Float64, out[8] )
end

return NPSS_success, Fn_N, Power, Torque, Nmech, Mtip, ηprop, ηDF
end
