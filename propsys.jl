"""
# TurboShaft calculations

## Design pt. method

Inputs: 
- alt_in : Ambient altitude [ft]
- MN_in  : Ambient Mach number [-]
- ShP    : Design ShaftPower in [W]
- π_LPC  : LPC pressure ratio [-]
- π_HPC  : HPC pressure ratio [-]
- Tt41   : Design Turbine inlet temp [R]
- cpsi   : Cells per square inch for catalyst [CPSI]
- w      : Catalyst cell wall thickness [mil] 
- lcat   : Reacting channel length [m]
- deNOx_in: Target design deNOx [-]

Outputs:
- η eta_thermal
- ṁf
- BSFC
- deNOx
- mcat

"""
function TurboShaft(alt_in::Float64, MN_in::Float64, 
                    ShP::Float64, 
                    π_LPC::Float64, π_HPC::Float64, Tt41::Float64,
                    cpsi::Float64, w::Float64, lcat::Float64, deNOx_in::Float64; LHV = 120,
                    file_name = "NPSS_Turboshaft/DesScl.int")

    NPSS_TShaft_input(alt_in, MN_in, ShP, 
                        Tt41, π_LPC, π_HPC, 
                        cpsi, w, lcat, deNOx_in; LHV = LHV)

    # open("NPSS_Turboshaft/OffDesInputs.inp", "w") do io
    #     println(io, "//DUMMY since only running design point now")
    # end

    NPSS_run("NPSS_Turboshaft/", "TP.bat")

    include("NPSS_Turboshaft/Eng.output")
    # Write all design point scalars to file
      
    return eta_thermal, mdotf, BSFC, deNOx, mcat  #, MapScalars, NozArea

end

"""
# TurboShaft calculations

## Off-Design method

Inputs: 
- alt_in : Ambient altitude 
- MN_in  : Ambient Mach number
- Tt41   : Operating temp in [R]


Outputs:
- Pshaft
- η eta_thermal
- ṁf
- BSFC
- deNOx
- mcat

"""
function TurboShaft(alt_in::Float64, MN_in::Float64, 
                    Tt41::Float64, Nshaft::Float64)

    NPSS_TShaft_input(alt_in, MN_in,
                    Tt41, Nshaft;
                    file_name = "NPSS_Turboshaft/OffDesInputs.inp")

    NPSS_run("NPSS_Turboshaft/", "TPoffDes.bat")

    include("NPSS_Turboshaft/Eng.output")

    return ShP, eta_thermal, mdotf, BSFC, deNOx #, MapScalars, NozArea

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
function DuctedFan(alt_in::Float64, MN_in::Float64,  Pin::Float64,
                    Kinl::Float64, Φinl::Float64)
    
    NPSS_Fan_input(alt_in, MN_in, Pin, Kinl, Φinl)
    #Run Off-Design model
    NPSS_run("NPSS_Turboshaft/", "FanOffDes.bat")

    include("NPSS_Turboshaft/Fan.output")

    return Fn, Fan_power, Torque_fan, N_fan, Mtip, eta_prop, eta_DF

end

"""
PowerTrain 

Design method - sizes the powertrain
"""
function PowerTrain(alt_in::Float64, MN_in::Float64, Fn::Float64,
                    Kinl::Float64, Φinl::Float64,
                    parpt::Array{Union{Float64, Int64},1},
                    parmot::Array{Float64, 1},
                    pargen::Array{Float64, 1})
    
    # Initialize powertrain weight and waste heat 
        Wpowertrain = 0.0
        Hrej = 0.0
    # Unpack number of powertrain elements
        nfan = parpt[ipt_nfan]
        ngen = parpt[ipt_ngen]
        nTshaft = parpt[ipt_nTshaft]

    # Thrust per fan
        Ffan = Fn/nfan

    # Call ducted fan design method
        πfan = parpt[ipt_pifan]

        Dfan, Fan_power, Torque_fan, N_fan, Mtip,
        ηpropul, ηDF,
        FanMapScalars, FanNozArea, Wfan = DuctedFan(alt_in, MN_in, Ffan, Kinl, Φinl, πfan )
        # println("Fan:")
        # println("Fan Ø = ", Dfan, " m")
        # println("Fan Fn = ", Ffan, " N")
        parpt[ipt_Wfan]    = Wfan # Save fan weight
        parpt[ipt_NdesFan] = N_fan

        Wpowertrain += Wfan*nfan
        Pshaft_mot = -1000. * Fan_power
        
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
        parpt[ipt_Wmot]    = Wmot
        parpt[ipt_NdesMot] = RPMmot

        Hwaste_motor = PreqMot - Pshaft_mot
        Hrej += nfan*Hwaste_motor # Heat rejected from motors
        Wpowertrain += Wmot*nfan  # Add to total powertrain weight
           
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
function PowerTrainOD(alt_in::Float64, MN_in::Float64, Tt41::Float64,
                    Kinl::Float64, Φinl::Float64,
                    parpt::Array{Union{Float64, Int64},1},
                    parmot::Array{Float64, 1},
                    pargen::Array{Float64, 1})
    
    Hrej = 0.0
    nfan = parpt[ipt_nfan]
    ngen = parpt[ipt_ngen]
    nTshaft = parpt[ipt_nTshaft]
    Nshaft = 30000. #RPM

    # Calculate Turboshaft power output
        Pshaft, ηthermal,
        mdotf, BSFC,
        deNOx_out = TurboShaft(alt_in, MN_in, Tt41, Nshaft)

    # Run generator
        Pgen_in = Pshaft * nTshaft/ngen
        Ngen = Nshaft
        Pgen_out, ηgen, PLgen = PMSM(Pgen_in, Ngen/60, pargen)

        Hwaste_gen = (Pgen_in - Pgen_out)
        Hrej += ngen*Hwaste_gen # Heat rejected from motors

    # Off-des cable
        ηcable, Wcable = cable() #TODO cable is dummy right now
        Hwaste_cable = Pgen_out*(1-ηcable)
        Hrej += Hwaste_cable  # Heat rejected from cables
        Pinv_in = Pgen_out*ηcable * ngen/nfan # for each inverter

    # Off-des inverter
        ηinv = inverter(Pinv_in, parmot)
        Hwaste_inv = Pinv_in*(1-ηinv)
        Hrej += nfan * Hwaste_inv # Heat rejected from all inverters
        Pmot_in = Pinv_in * ηinv

    # Off-des motor
        Nmot = parpt[ipt_NdesMot] #* (N_fan/parpt[ipt_NdesFan])
        # Off-des motor call
        Pmot_out, ηmot, PL = PMSM(Pmot_in, Nmot/60, parmot)
        # println(ηmot)
        Hwaste_motor = Pmot_in - Pmot_out
        Hrej += nfan*Hwaste_motor # Heat rejected from motors
        # println("Motor speed = ", Nmot)
        # println("Motor power = ", Pmot_out)
    # Ducted fan
        Pfan_in = Pmot_out
        Fn, Fan_power, Torque_fan, N_fan, Mtip,
        ηpropul, ηDF = DuctedFan(alt_in, MN_in, Pfan_in, Kinl, Φinl)
        # println("Fan speed = ", N_fan) 
        # TODO fan vs motor speed discrepency 

    Ftotal = Fn*nfan

    return Ftotal, [ηmot, ηinv, ηcable, ηgen, ηthermal],
           [Pmot_out, Pmot_in, Pinv_in, Pgen_in, Pshaft], 
           [Hwaste_motor, Hwaste_inv, Hwaste_cable, Hwaste_gen, Hrej]./1000,
           mdotf, BSFC,
           deNOx_out

end