
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


function NPSS_TFsys(NPSS, alt_in, MN_in, Fn, Tt41,
    πfan, first, parg, parpt)

    nfans   = parpt[ipt_nfan]
    neng = parg[igneng]

    πLPC = parpt[ipt_piLPC]
    πHPC = parpt[ipt_piHPC]
    LHV  = parg[igLHVfuel]

    rSnace = parg[igrSnace]
    fpylon = parg[igfpylon]
    feadd  = parg[igfeadd]

    input_string = "111 "*
    "first = $first;" *
    "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
    "Eng.Amb.MN_in  = $MN_in ;" *

    "Fn_target         = $(Fn/neng);"*
    "Tt41              = $Tt41 ;" *

    "Eng.FusEng.LHV = $(LHV*429.923); "*
    "Eng.CmpF.PRdes = $(πfan) ;" *
    "Eng.CmpL.PRdes = $(πLPC) ;" *
    "Eng.CmpH.PRdes = $(πHPC) ;" *

    # "Eng.ShL.Nmech = 3500.0;"*
    # "Eng.ShF.Nmech = 3400.0;"*

    
    
    "Eng.rSnace = $(rSnace) ;" *
    "Eng.feadd = $(feadd) ;" *
    "Eng.fpylon = $(fpylon) ;" *
    "\n"

    write(NPSS, input_string)
    # write(stdout, input_string)

    out = split(String(readavailable(NPSS.out)), "_") # `readavailable(stream)` is blocking only if no data is available
    NPSS_success = parse(Float64, out[1] )
    
    if length(out) > 1
        mdotf      = parse(Float64, out[2] )
        EINOx1     = parse(Float64, out[3] )
        mdot3      = parse(Float64, out[4] )
        Tt3        = parse(Float64, out[5] )
        OPR        = parse(Float64, out[6])
        Wc3        = parse(Float64, out[7])
        LHV        = parse(Float64, out[8])       
        W_in       = parse(Float64, out[9])


        Weng    = parse(Float64, out[10])
        
        Dfan    = parse(Float64, out[11])


        FAR        = parse(Float64, out[12])

        Snace1     = parse(Float64, out[13])

        Tt41       = parse(Float64, out[14])
        EGT        = parse(Float64, out[15])

        ηtshaft   = parse(Float64, out[16]) 

        Pshaft_fan = parse(Float64, out[17])

    end

    neng = parg[igneng    ]

    mdotf_tot = mdotf*neng

    parpt[ipt_Ptshaft] = Pshaft_fan #Add this for cost model


    Wpowertrain = 0.0
    xWpowertrain = 0.0
    parg[igdfan] = Dfan
    parg[igWfan] = 0.0

    
  
      
    
#    ASK: Do I need this for TF??
#     lcable = (parg[igxgen] - parg[igxwbox]) + parg[igyeout]
#     parpt[ipt_lcable] = lcable
#     ηcable, Wcable = cable(Pelec_mot, 2*Vcable, lcable, parpt) 
#      parg[igWcables] = Wcable*nfans # multiply by number of fans
#      parpt[ipt_Wcables] = Wcable*nfans
#      Wpowertrain += Wcable*nfans # Add to total powertrain weight
# #    xWpowertrain += Wcable*parg[igxcable] #need to add x of cable 


        


    parg[igWeng  ] = Weng

    # parg[igWeng  ] = 12280.0 * 4.44822

    # ASK: Is this necessary??
    # Hrejmot = nfans*(Pelec_mot - Pshaft_fan)
    # Hrejinv = nfans*(Pelec_mot*(1-ηinv))
    # Hrejcab = ngen*(Pelec_gen*(1-ηcable))
    # Hrejgen = ngen*(Pshaft_gen - Pelec_gen)
    # Hrejrect= ngen*(Pelec_gen*(1-ηinv))
    # Hrejtot = Hrejmot + Hrejinv + Hrejcab + Hrejgen + Hrejrect

    # Qtms = Hrejtot

    # ASK: what about this stuff?? 
#     P_mTMS = 8.0*hp_to_W/(1/2.205)  # convert 8 hp/lbm to W/kg. Value based on MIT NASA LEARN report Table A8
#     Wtms = Qtms/P_mTMS * 9.81
#     parg[igWtms] = Wtms
#     Wpowertrain += Wtms
#    xWpowertrain += Wtms*(parg[igxtshaft] + parg[igxmot])/2

    Wpowertrain = 0.0
    xWpowertrain = 0.0
    
    parg[ igWtesys] = Wpowertrain
    parg[igxWtesys] = xWpowertrain
    



    # ηpt = [0.0, ηmot, ηinv, ηcable, ηgen, ηtshaft]
    # Hpt = [Hrejmot, Hrejinv, Hrejcab, Hrejgen+Hrejrect, Hrejtot]
    # Ppt = [Pshaft_fan, Pelec_mot, Pelec_gen, Pshaft_gen]
    # SPpt = [SPmot,SPinv, SPgen, SPtshaft]

    heatexcess = 0.0

    # return NPSS_success,ηpt, SPpt, Ppt, Hpt, heatexcess, 
    #         mdotf_tot, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, EGT,
    #         Snace1, Saftnace1

     return NPSS_success, heatexcess, 
    mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, EGT,
    Snace1

end

""" 
Off Des NPSS_TEsys
"""
function NPSS_TFsysOD(NPSS, alt_in, MN_in, Fn, Tt41,
     first, parg::Array{Float64, 1}, parpt::Array{Union{Float64, Int64},1}, 
     pare, ip::Int64)
    
    #  println("Starting Off des at ip = $ip")
   nTshaft = parpt[ipt_nTshaft]
   nfans   = parpt[ipt_nfan]

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
       Ftotal     = parse(Float64, out[3] )
       EINOx1     = parse(Float64, out[4] )
       mdot3      = parse(Float64, out[5] )
       Tt3        = parse(Float64, out[6] )
       OPR        = parse(Float64, out[7])
       Wc3        = parse(Float64, out[8])
       W_in       = parse(Float64, out[9])
   

       FAR        = parse(Float64, out[10])
       Tt41       = parse(Float64, out[11])
       EGT        = parse(Float64, out[12])


       ηtshaft   = parse(Float64, out[13]) 

       pare[iegsmdotf   , ip] = parse(Float64, out[14])
       pare[iegsWin     , ip] = parse(Float64, out[15])
       pare[iegsRlineF  , ip] = parse(Float64, out[16])
       pare[iegsBPR     , ip] = parse(Float64, out[17])
       pare[iegsRlinelc , ip] = parse(Float64, out[18])
       pare[iegsRlinehc , ip] = parse(Float64, out[19])
       pare[iegsPRhtrb  , ip] = parse(Float64, out[20])
       pare[iegsPRltrb  , ip] = parse(Float64, out[21])
       pare[iegsNmechH  , ip] = parse(Float64, out[22])
       pare[iegsNmechL  , ip] = parse(Float64, out[23])
    #    pare[iegsNmechF  , ip] = parse(Float64, out[24])

   end
   mdotf_tot = mdotf*nTshaft
   Ftotal    = Ftotal*nTshaft

#    Hrejmot = nfans*(Pelec_mot - Pshaft_fan)
#    Hrejinv = nfans*(Pelec_mot*(1-ηinv))
#    Hrejcab = nfans*(Pelec_mot*(1-ηcable))
#    Hrejgen = ngen*(Pshaft_gen - Pelec_gen)
#    Hrejrect= ngen*(Pelec_gen*(1-ηinv))
#    Hrejtot = Hrejmot + Hrejinv + Hrejcab + Hrejgen + Hrejrect

#    ηpt = [ηmot, ηinv, ηcable, ηgen, ηtshaft]
#    Hpt = [Hrejmot, Hrejinv, Hrejcab, Hrejgen+Hrejrect, Hrejtot]
#    Ppt = [ Pshaft_fan, Pelec_mot, Pelec_gen, Pshaft_gen]
   heatexcess = 0.0

#    return NPSS_success, Ftotal, ηpt, Ppt, Hpt, heatexcess, 
#            mdotf_tot, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT

    return NPSS_success, Ftotal, heatexcess, 
    mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT

end
