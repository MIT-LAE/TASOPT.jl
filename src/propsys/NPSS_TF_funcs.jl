"""
NPSS_TFsys

Turbofan sizing code for design point sizing
"""
function NPSS_TFsys(NPSS, alt_in, MN_in, Fn, Tt41,
    πfan, first, parg, parpt, mofft, Pofft, Elec_PowerTot)

    neng = parg[igneng]

    πLPC = parpt[ipt_piLPC]
    πHPC = parpt[ipt_piHPC]

    LHV  = parg[igLHVfuel]

    rSnace = parg[igrSnace]
    fpylon = parg[igfpylon]
    feadd  = parg[igfeadd]

    fanPCT = parg[igfanPCT]

    ConvFromJkgToBtuLbm = 0.9478171*0.4535924/1000

    input_string = "111 "*
    "first = $first;" *
    "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
    "Eng.Amb.MN_in  = $MN_in ;" *

    "Fn_target              = $(Fn/neng);"*
    "Eng.CmpH.Fl_B027.Wbld  = $(mofft*2.20462);"*
    "Eng.Gen.ShP            = $(Pofft*0.00134102);"*

    "Tt41    = $Tt41 ;" *
    "fanPCT  = $fanPCT;"*

    "Eng.FusEng.LHV = $(LHV*429.923); "*

    "Eng.CmpF.PRdes = $(πfan) ;" *
    "Eng.CmpL.PRdes = $(πLPC) ;" *
    "Eng.CmpH.PRdes = $(πHPC) ;" *

    "Eng.rSnace = $(rSnace) ;" *
    "Eng.feadd = $(feadd) ;" *
    "Eng.fpylon = $(fpylon) ;" *
    "\n"

    write(NPSS, input_string)
    # write(stdout, input_string)

    out = split(String(readuntil(NPSS.out, "END")), "_") # `readavailable(stream)` is blocking only if no data is available
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

    # neng = parg[igneng    ]

    if (@isdefined mdotf)

        mdotf_tot = mdotf*neng
        parpt[ipt_Ptshaft] = Pshaft_fan #Add this for cost model
        
        
        Wpowertrain = 0.0
        xWpowertrain = 0.0
        parg[igdfan] = Dfan
        
        parg[igWeng  ] = Weng
        

        Wpowertrain = 0.0
        xWpowertrain = 0.0
        
        parg[ igWtesys] = Wpowertrain
        parg[igxWtesys] = xWpowertrain
        
        heatexcess = 0.0

        return NPSS_success, heatexcess, 
        mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, EGT,
        Snace1

    else
        NPSS_success = false
        Ftotal = 0.0
        heatexcess = 0.0
        mdotf_tot = 0.0
        EINOx1  = 0.0
        FAR = 0.0
        Tt3 = 0.0
        OPR = 0.0
        Wc3 = 0.0
        EGT = 0.0
        Snace1 = 0.0

        # println("mdotf not defined")
        return NPSS_success, heatexcess, 
        mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT, Snace1
   
   end

end

""" 
Off-design NPSS for turbofan prop.
"""
function NPSS_TFsysOD(NPSS, alt_in, MN_in, Fn, Tt41,
     first, parg::Array{Float64, 1}, parpt::Array{Union{Float64, Int64},1}, 
     pare, ip::Int64, fanPCT, mofft, Pofft)
    
    #  println("Starting Off des at ip = $ip")
    neng = parg[igneng]
    nfans   = parpt[ipt_nfan]

   if Fn == 0.0
    mode = "222" #NcFan specified
   else
    mode = "333" #Fn specified
   end

   input_string = "$(mode) "*
   "first = $first;" *
   "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
   "Eng.Amb.MN_in  = $MN_in ;" *

   "Eng.CmpH.Fl_B027.Wbld             = $(mofft*2.20462);"*
   "Eng.Gen.ShP             = $(Pofft*0.00134102);"*

   "Fn_target         = $(Fn/neng);"*
   "fanPCT            = $fanPCT;"*
   "Tt41              = $Tt41 ;" *
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

   if (@isdefined mdotf)

    mdotf_tot = mdotf*neng
    Ftotal    = Ftotal*neng

   heatexcess = 0.0

    NPSS_success = true

    return NPSS_success, Ftotal, heatexcess, 
    mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT

   else
    NPSS_success = false
    Ftotal = 0.0
    heatexcess = 0.0
    mdotf_tot = 0.0
    EINOx1  = 0.0
    FAR = 0.0
    Tt3 = 0.0
    OPR = 0.0
    Wc3 = 0.0
    Tt41 = 0.0
    EGT = 0.0

   return NPSS_success, Ftotal, heatexcess, 
    mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT
   
   end

end

"""
Off-design sizing option where Tt41 is specified
"""
function NPSS_TFsysOD2(NPSS, alt_in, MN_in, Fn, Tt41,
    first, parg::Array{Float64, 1}, parpt::Array{Union{Float64, Int64},1}, 
    pare, ip::Int64, fanPCT, mofft, Pofft)
   
  neng = parg[igneng]
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

  "Eng.CmpH.Fl_B027.Wbld             = $(mofft*2.20462);"*
  "Eng.Gen.ShP             = $(Pofft*0.00134102);"*


  "Fn_target         = $(Fn/neng);"*
  "fanPCT            = $fanPCT;"*
  "Tt41              = $Tt41 ;" *   
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

      dummy = parse(Float64, out[14])
      dummy = parse(Float64, out[15])
      dummy = parse(Float64, out[16])
      dummy = parse(Float64, out[17])
      dummy = parse(Float64, out[18])
      dummy = parse(Float64, out[19])
      dummy = parse(Float64, out[20])
      dummy = parse(Float64, out[21])
      dummy = parse(Float64, out[22])
      dummy = parse(Float64, out[23])
   #    pare[iegsNmechF  , ip] = parse(Float64, out[24])

  end

  if (@isdefined mdotf)

   mdotf_tot = mdotf*neng
   Ftotal    = Ftotal*neng

  heatexcess = 0.0

   NPSS_success = true

   return NPSS_success, Ftotal, heatexcess, 
   mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT

  else
   NPSS_success = false
   Ftotal = 0.0
   heatexcess = 0.0
   mdotf_tot = 0.0
   EINOx1  = 0.0
   FAR = 0.0
   Tt3 = 0.0
   OPR = 0.0
   Wc3 = 0.0
   Tt41 = 0.0
   EGT = 0.0

  return NPSS_success, Ftotal
  
  end

end

"""
Off-design sizing option where Tt41 is specified
"""
function NPSS_TFsysOD3(NPSS, alt_in, MN_in, Fn, Tt41,
    first, parg::Array{Float64, 1}, parpt::Array{Union{Float64, Int64},1}, 
    pare, ip::Int64, fanPCT)
   
   #  println("Starting Off des at ip = $ip")
  neng = parg[igneng]
  nfans   = parpt[ipt_nfan]

  mode = "444"

  input_string = "$(mode) "*
  "first = $first;" *
  "Eng.Amb.alt_in = $(alt_in/ft_to_m);" *
  "Eng.Amb.MN_in  = $MN_in ;" *

  "Fn_target         = $(Fn/neng);"*
  "fanPCT            = $fanPCT;"*
  "Tt41              = $Tt41 ;" *
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

  if (@isdefined mdotf)

   mdotf_tot = mdotf*neng
   Ftotal    = Ftotal*neng

  heatexcess = 0.0

   NPSS_success = true

   return NPSS_success, Ftotal, heatexcess, 
   mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT

  else
   NPSS_success = false
   Ftotal = 0.0
   heatexcess = 0.0
   mdotf_tot = 0.0
   EINOx1  = 0.0
   FAR = 0.0
   Tt3 = 0.0
   OPR = 0.0
   Wc3 = 0.0
   Tt41 = 0.0
   EGT = 0.0

  return NPSS_success, Ftotal, heatexcess, 
   mdotf_tot, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT
  
  end


end