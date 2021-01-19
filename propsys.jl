"""
# TurboShaft calculations

## Design pt. method

Inputs: 
- alt_in : Ambient altitude 
- MN_in  : Ambient Mach number
- ShP    : Design ShaftPower in [N]
- π_fan  : Fan pressure ratio

Outputs:


"""
function TurboShaft(alt_in::Float64, MN_in::Float64, 
                    ShP::Float64, 
                    π_LPC::Float64, π_HPC::Float64, Tt41::Float64,
                    cpsi::Float64, w::Float64, lcat::Float64, deNOx_in::Float64; LHV = 120,
                    file_name = "NPSS_Turboshaft/DesScl.int")

    NPSS_TShaft_input(alt_in, MN_in, ShP, 
                        Tt41, π_LPC, π_HPC, 
                        cpsi, w, lcat, deNOx_in; LHV = LHV)

    open("NPSS_Turboshaft/OffDesInputs.inp", "w") do io
        println(io, "//DUMMY since only running design point now")
    end

    NPSS_run("NPSS_Turboshaft/", "TP.bat")

    include("NPSS_Turboshaft/Eng.output")
    # Write all design point scalars to file
    open(file_name, "w") do io
        println(io, "\n// LHV based on fuel")
        println(io, "Eng.FusEng.LHV = ", LHV*429.923, ";") # save to file to avoid later conflicts
        println(io, "Eng.FsEng.W_in = ", mdot*2.205, ";")
        println(io, "Eng.FusEng.Wfuel = ", mdotf*2.205, ";")

        println(io, "\n// Map scalars")
        println(io, "Eng.CmpL.S_map.s_effDes = ", LPCscalars[1], ";")
        println(io, "Eng.CmpL.S_map.s_PRdes  = ", LPCscalars[2], ";")
        println(io, "Eng.CmpL.S_map.s_WcDes  = ", LPCscalars[3], ";")
        println(io, "Eng.CmpL.S_map.s_NcDes  = ", LPCscalars[4], ";")

        println(io, "Eng.CmpH.S_map.s_effDes = ", HPCscalars[1], ";")
        println(io, "Eng.CmpH.S_map.s_PRdes  = ", HPCscalars[2], ";")
        println(io, "Eng.CmpH.S_map.s_WcDes  = ", HPCscalars[3], ";")
        println(io, "Eng.CmpH.S_map.s_NcDes  = ", HPCscalars[4], ";")

        println(io, "Eng.TrbH.S_map.s_effDes = ", HPTscalars[1], ";")
        println(io, "Eng.TrbH.S_map.s_PRdes  = ", HPTscalars[2], ";")
        println(io, "Eng.TrbH.S_map.s_WpDes  = ", HPTscalars[3], ";")
        println(io, "Eng.TrbH.S_map.s_NpDes  = ", HPTscalars[4], ";")
  
        println(io, "Eng.TrbP.S_map.s_effDes = ", PTscalars[1], ";")
        println(io, "Eng.TrbP.S_map.s_PRdes  = ", PTscalars[2], ";")
        println(io, "Eng.TrbP.S_map.s_WpDes  = ", PTscalars[3], ";")
        println(io, "Eng.TrbP.S_map.s_NpDes  = ", PTscalars[4], ";")

        println(io, "\nEng.CmpL.S_map.NcMapDes    = ", LPC[1], ";")
        println(io, "Eng.CmpL.S_map.WcMapDes    = ", LPC[2], ";")
        println(io, "Eng.CmpL.S_map.PRmapDes    = ", LPC[3], ";")
        println(io, "Eng.CmpL.S_map.effMapDes   = ", LPC[4], ";")
        println(io, "Eng.CmpL.S_map.RlineMapDes = ", LPC[5], ";")

        println(io, "Eng.CmpH.S_map.NcMapDes    = ", HPC[1], ";")
        println(io, "Eng.CmpH.S_map.WcMapDes    = ", HPC[2], ";")
        println(io, "Eng.CmpH.S_map.PRmapDes    = ", HPC[3], ";")
        println(io, "Eng.CmpH.S_map.effMapDes   = ", HPC[4], ";")
        println(io, "Eng.CmpH.S_map.RlineMapDes = ", HPC[5], ";")

        println(io, "Eng.TrbH.S_map.PRmapDes = ", HPT[1], ";")
        println(io, "Eng.TrbH.S_map.NpMapDes = ", HPT[2], ";")

        println(io, "Eng.TrbP.S_map.PRmapDes = ", PT[1], ";")
        println(io, "Eng.TrbP.S_map.NpMapDes = ", PT[2], ";")

        
   
        println(io, "\n// Nozzle Area")
        println(io, "Eng.NozPri.AthCold  = ", NozArea, ";")
    
       
        println(io, "\n// PCEC parameters")
        println(io, "Eng.PCEC.l = ", lcat, ";")
        println(io, "Eng.PCEC.w = ", w, ";")
        println(io, "Eng.PCEC.cpsi = ", cpsi, ";")
        println(io, "Eng.PCEC.Af = ", Af, ";")

        println(io, "\n// Bleeds")
        println(io, "Eng.B030.TCLA_NC.fracW = ", TCLA_NC,";")
        println(io, "Eng.B030.TCLA_CH.fracW = ", TCLA_CH,";")

        println(io, "\n// Areas")
        println(io, "Eng.FsEng.Fl_O.Aphy  = ", Areas[ 1], ";" ) 
        println(io, "Eng.InEng.Fl_O.Aphy  = ", Areas[ 2], ";" ) 
        println(io, "Eng.CmpL.Fl_O.Aphy   = ", Areas[ 3], ";" ) 
        println(io, "Eng.D025.Fl_O.Aphy   = ", Areas[ 4], ";" ) 
        println(io, "Eng.CmpH.Fl_O.Aphy   = ", Areas[ 5], ";" ) 
        println(io, "Eng.B030.Fl_O.Aphy   = ", Areas[ 6], ";" ) 
        println(io, "Eng.BrnPri.Fl_O.Aphy = ", Areas[ 7], ";" ) 
        println(io, "Eng.TrbH.Fl_O.Aphy   = ", Areas[ 8], ";" ) 
        println(io, "Eng.D050.Fl_O.Aphy   = ", Areas[ 9], ";" ) 
        println(io, "Eng.TrbP.Fl_O.Aphy   = ", Areas[10], ";" ) 
        println(io, "Eng.D060.Fl_O.Aphy   = ", Areas[11], ";" ) 
        println(io, "Eng.PCEC.Fl_O.Aphy   = ", Areas[12], ";" ) 
        println(io, "Eng.NozPri.Fl_O.Aphy = ", Areas[13], ";" ) 
    end
  
    return eta_thermal, mdotf, BSFC, deNOx, mcat  #, MapScalars, NozArea

end
function TurboShaft(alt_in::Float64, MN_in::Float64, 
                    ShP::Float64)

    NPSS_TShaft_input(alt_in, MN_in, ShP; file_name = "NPSS_Turboshaft/OffDesInputs.inp")

    NPSS_run("NPSS_Turboshaft/", "TP.bat")

    include("NPSS_Turboshaft/Eng.output")


    return eta_thermal, mdotf, BSFC, deNOx, mcat  #, MapScalars, NozArea

end


"""
# Ducted Fan calculations

## Design pt. method

Inputs: 
- alt_in : Ambient altitude 
- MN_in  : Ambient Mach number
- Fn     : Design Thurst reqired in [N]
- π_fan  : Fan pressure ratio

Outputs:
- Dfan        : diameter of fan in [m]
- Fan_power   : Required fan shaft power [kW]
- Torque_fan  : Required fan torque [N*m]
- N_fan       : Fan speed [RPM]
- eta_prop    : propulsive eff
- MapScalars  : Fan map scalars for off-des calcs
- NozArea     : Fan nozzle area for off-des calcs [in^2] in inches so can easily pass back to NPSS

"""
function DuctedFan(alt_in::Float64, MN_in::Float64,  Fn::Float64,
                    Kinl::Float64, Φinl::Float64,
                     π_fan::Float64)
    
    NPSS_Fan_input(alt_in, MN_in, Fn, Kinl, Φinl, π_fan)
    #Run Design model
    NPSS_run("NPSS_Turboshaft/", "FanDes.bat")

    include("NPSS_Turboshaft/Fan.output")
    
    ARfan  = 3   # Blade aspeect ratio
    bladeσ = 0.4 # Blade solidity c/s
    ktech = 0.5 
    Utip  = Dfan/2* (2π * N_fan /60)
    # Sagerser 1971, NASA TM X-2406
    # Note: The term "weight" in Sagerser1971 is actually mass
    mfan = ktech*(135.0 * Dfan^2.7/sqrt(ARfan) * (bladeσ/1.25)^0.3 * (Utip/350.0)^0.3)
    Wfan = mfan*gee

    return Dfan, Fan_power, Torque_fan, N_fan, Mtip,
            eta_prop, eta_DF,
            MapScalars, NozArea, Wfan

end

"""
# Ducted Fan calculations

## Off-design method

Inputs: 
- alt_in       : Ambient altitude 
- MN_in        : Ambient Mach number
- Fn           : Design Thurst reqired 
- Map Scalars  : Fan map scalars for off-des
- NozArea      : Nozzle area for off-des

Outputs:
- Fan_power   : Required fan shaft power [kW]
- Torque_fan  : Required fan torque [N*m]
- N_fan       : Fan speed [RPM]
- eta_prop    : propulsive eff


"""
function DuctedFan(alt_in::Float64, MN_in::Float64,  Fn::Float64,
                    Kinl::Float64, Φinl::Float64)
    
    NPSS_Fan_input(alt_in, MN_in, Fn, Kinl, Φinl)
    #Run Off-Design model
    NPSS_run("NPSS_Turboshaft/", "FanOffDes.bat")

    include("NPSS_Turboshaft/Fan.output")

    return Fan_power, Torque_fan, N_fan, Mtip, eta_prop, eta_DF

end

"""
PowerTrain 

Design method - sizes the powertrain
"""
function PowerTrain(alt_in::Float64, MN_in::Float64, Fn::Float64,
                    parpt::Array{Union{Float64, Int64},1},
                    parmot::Array{Float64, 1},
                    pargen::Array{Float64, 1})
    
    Wpowertrain = 0.0
    Hrej = 0.0
    nfan = parpt[ipt_nfan]
    ngen = parpt[ipt_ngen]
    nTshaft = parpt[ipt_nTshaft]

    Ffan = Fn/nfan
    Kinl, Φinl = 0., 0.
    # Call ducted fan design method
        πfan = parpt[ipt_pifan]

        Dfan, Fan_power, Torque_fan, N_fan, Mtip,
        ηpropul, ηDF,
        FanMapScalars, FanNozArea, Wfan = DuctedFan(alt_in, MN_in, Ffan, Kinl, Φinl, πfan )
        # println("Fan:")
        # println("Fan Ø = ", Dfan, " m")
        # println("Fan Fn = ", Ffan, " N")
        parpt[ipt_Wfan] = Wfan
        Wpowertrain += Wfan*nfan
        Pshaft_mot = -1000. * Fan_power
        parpt[ipt_NdesFan] = N_fan

    # Size motor
        ratAsp   = parpt[ipt_ARmot]
        σAg      = parpt[ipt_sigAgMot]
        ratSplit = parpt[ipt_ratSplitMot]

        Wmot, PreqMot, ηmot, RPMmot,
        PL, PLiron, PLCu, PLwind, SPmot = PMSM(Pshaft_mot, ratAsp, σAg, ratSplit, parmot)
        # println("\nMotor:")
        # println("Losses = ", PL)
        # println("\tIron losses    = ", PLiron, "%")
        # println("\tCopper losses  = ", PLCu)
        # println("\tWindage losses = ", PLwind)
        Hwaste_motor = PreqMot - Pshaft_mot
        Hrej += nfan*Hwaste_motor # Heat rejected from motors
        parpt[ipt_Wmot] = Wmot
        Wpowertrain += Wmot*nfan # Add to total powertrain weight
        parpt[ipt_NdesMot] = RPMmot

    
    # Size Inverter and cables
        ηinv, Winv, SPinv = inverter(PreqMot, RPMmot/60, parmot)
        Hwaste_inv = PreqMot*(1-ηinv)
        Hrej += nfan * Hwaste_inv # Heat rejected from all inverters
        parpt[ipt_Winv] = Winv
        Wpowertrain += Winv*nfan # Add to total powertrain weight
    
        ηcable, Wcable = cable() #TODO cable is dummy right now
        Hwaste_cable = PreqMot*(1-ηcable)
        Hrej += Hwaste_cable  # Heat rejected from all inverters
        parpt[ipt_Wcables] = Wcable
        Wpowertrain += Wcable # Add to total powertrain weight

    # Size generator
        PreqGen = PreqMot * ηinv * ηcable

        ratAsp   = parpt[ipt_ARgen]
        σAg      = parpt[ipt_sigAgGen]
        ratSplit = parpt[ipt_ratSplitGen]

        Wgen, PgenShaft, ηgen, RPMgen,
        PLgen, PLirongen, PLCugen, PLwindgen, SPgen = PMSM(PreqGen*nfan/ngen,
                                                            ratAsp, σAg, ratSplit,
                                                            pargen)
                                                            
        parpt[ipt_NdesGen] = RPMgen
        Hwaste_gen = (PgenShaft - PreqGen*nfan/ngen)
        Hrej += ngen*Hwaste_gen # Heat rejected from motors
        parpt[ipt_Wgen] = Wgen
        Wpowertrain += Wgen*ngen
    
    # Size tubo-shaft and PCEC
        πLPC = parpt[ipt_piLPC]
        πHPC = parpt[ipt_piHPC]
        Tt41 = parpt[ipt_Tt41]
        cpsi = parpt[ipt_cpsi]
        w    = parpt[ipt_wcat]
        lcat = parpt[ipt_lcat]
        deNOx= parpt[ipt_deNOx]
        Ptshaft = PgenShaft*ngen/nTshaft
        ηthermal, mdotf, BSFC, deNOx_out, mcat = TurboShaft(alt_in, MN_in, Ptshaft,
                                            πLPC, πHPC, Tt41,
                                            cpsi, w, lcat, deNOx)

        Wcat = mcat*gee
        parpt[ipt_Wcatalyst] = Wcat
        Wpowertrain += Wcat*nTshaft

        SPtshaft= 10.4e3 # W/kg Based on the RR T406 (4.58 MW power output). The GE38 (~5 MW) has a power density of 11.2 kW/kg
        Wtshaft = gee*(PgenShaft*ngen/nTshaft)/SPtshaft
        parpt[ipt_Wtshaft] = Wtshaft
        Wpowertrain += Wtshaft*nTshaft

    parpt[ipt_Wpttotal] = Wpowertrain

    return [ηpropul, ηmot, ηinv, ηcable, ηgen, ηthermal],
           [Pshaft_mot, PreqMot, PgenShaft, Ptshaft], 
           [Hwaste_motor, Hwaste_inv, Hwaste_cable, Hwaste_gen, Hrej]./1000,
           [Wfan, Wmot, Winv, Wcable, Wgen, Wtshaft, Wcat, Wpowertrain]./gee,
           [SPmot, SPinv, SPgen, SPtshaft], mdotf, BSFC,
           deNOx, FanMapScalars, FanNozArea

end
function PowerTrain(alt_in::Float64, MN_in::Float64, Fn::Float64, FanMapScalars, FanNozArea,
                    parpt::Array{Union{Float64, Int64},1},
                    parmot::Array{Float64, 1},
                    pargen::Array{Float64, 1})
    
    Hrej = 0.0
    nfan = parpt[ipt_nfan]
    ngen = parpt[ipt_ngen]
    nTshaft = parpt[ipt_nTshaft]

    Kinl, Φinl = 0., 0.
    Ffan = Fn/nfan
    # Call ducted fan off-design method
    Fan_power, Torque_fan, N_fan, Mtip, ηpropul, ηDF = DuctedFan(alt_in, MN_in, Ffan, Kinl, Φinl)
    Pshaft_mot = -1000. * Fan_power
    # println("P motor = ", Pshaft_mot)
    
    Nmot = parpt[ipt_NdesMot] * (N_fan/parpt[ipt_NdesFan])
    # Off-des motor call
    PreqMot, ηmot, PL = PMSM(Pshaft_mot, Nmot/60, parmot)
    # println(ηmot)
    Hwaste_motor = PreqMot - Pshaft_mot
    Hrej += nfan*Hwaste_motor # Heat rejected from motors
    
    # Off-des inverter
    ηinv = inverter(PreqMot, parmot)
    Hwaste_inv = PreqMot*(1-ηinv)
    Hrej += nfan * Hwaste_inv # Heat rejected from all inverters
   
    ηcable, Wcable = cable() #TODO cable is dummy right now
    Hwaste_cable = PreqMot*(1-ηcable)
    Hrej += Hwaste_cable  # Heat rejected from all inverters
    # println(ηinv)
    PreqGen = (PreqMot * nfan) * ηinv * ηcable /ngen
    # println("Generator power req = ", PreqGen)
    Ngen = parpt[ipt_NdesGen] 
    PgenShaft, ηgen, PLgen = PMSM(PreqGen, Ngen/60, pargen)
    # println(ηgen)          
    Hwaste_gen = (PgenShaft - PreqGen*nfan/ngen)
    Hrej += ngen*Hwaste_gen # Heat rejected from motors
    # println("Tshaft power req = ", PgenShaft*ngen/nTshaft)
    ηthermal, mdotf, BSFC, deNOx_out, mcat = TurboShaft(alt_in, MN_in, PgenShaft*ngen/nTshaft)

    return [ηmot, ηinv, ηcable, ηgen, ηthermal],
           [Pshaft_mot, PreqMot, PgenShaft], 
           [Hwaste_motor, Hwaste_inv, Hwaste_cable, Hwaste_gen, Hrej]./1000,
           mdotf, BSFC,
           deNOx

end