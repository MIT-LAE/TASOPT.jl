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
function TurboShaft(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, 
                    ShP::Float64, 
                    π_LPC::Float64, π_HPC::Float64, Tt41::Float64,
                    cpsi::Float64, w::Float64, lcat::Float64, deNOx_in::Float64, first; LHV = 120,
                    file_name = "NPSS_Turboshaft/DesScl.int")

    NPSS_success, 
    eta_thermal, mdotf, BSFC,
     deNOx, mcat, 
     EINOx1, EINOx2, mdot, Tt3, OPR, Wc3  = NPSS_TShaft_input(NPSS, alt_in, MN_in, ShP, 
                        Tt41, π_LPC, π_HPC, 
                        cpsi, w, lcat, deNOx_in, first; LHV = LHV)

    # include("NPSS_Turboshaft/Eng.output")
      
    return eta_thermal, mdotf, BSFC, deNOx, mcat, EINOx1, EINOx2, mdot, Tt3, OPR, Wc3  #, MapScalars, NozArea

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
function TurboShaft(NPSS_TS, alt_in::Float64, MN_in::Float64, 
                    Tt41::Float64, Nshaft::Float64, first)

    # file_name = "NPSS_Turboshaft/Eng.output"

    NPSS_success, 
    ShP, eta_thermal, mdotf,
    BSFC, deNOx = NPSS_TShaft_run(NPSS_TS, alt_in, MN_in, Tt41, 30_000.0, first)
    
    # include(file_name)

    return ShP, eta_thermal, mdotf, BSFC, deNOx #, MapScalars, NozArea

end
function TurboShaft2(NPSS_TS, alt_in::Float64, MN_in::Float64, 
    ShP_dmd::Float64, Nshaft::Float64, first)

# file_name = "NPSS_Turboshaft/Eng.output"

NPSS_success, 
ShP, eta_thermal, mdotf,
BSFC, deNOx = NPSS_TShaft_run2(NPSS_TS, alt_in, MN_in, ShP_dmd, 30_000.0, first)

# include(file_name)

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
function DuctedFan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64,  Fn::Float64,
                    Kinl::Float64, Φinl::Float64,
                     π_fan::Float64, rSnace, first)
    
    # file_name = "NPSS_Turboshaft/Fan.output"
 
    NPSS_success,
    Dfan, Fan_power, Torque_fan, N_fan, Mtip,
    eta_prop, eta_DF,
    NozArea = runNPSS_Fan(NPSS, alt_in, MN_in, Fn, Kinl, Φinl, π_fan, first)
    
    # Read the data
    # include(file_name)
    
    ARfan  = 3   # Blade aspeect ratio
    bladeσ = 0.4 # Blade solidity c/s
    ktech = 0.5 
    Utip  = Dfan/2* (2π * N_fan /60)
    # Sagerser 1971, NASA TM X-2406
    # Note: The term "weight" in Sagerser1971 is actually mass
    mfan = ktech*(135.0 * Dfan^2.7/sqrt(ARfan) * (bladeσ/1.25)^0.3 * (Utip/350.0)^0.3)

    #Nacelle weight model from NASA CR 151970 but excluding core cowl for this ∵ this is only a ducted fan
    Snace1 = rSnace * 0.25 * π*Dfan^2
    Ainlet = 0.4*Snace1
    Acowl  = 0.2*Snace1
    Aexh   = 0.4*Snace1
    Wnace1 = 4.45*(Ainlet/0.3048^2) * (2.5+0.0238*Dfan/0.0254) +
             4.45*(Acowl /0.3048^2) *  1.9 +
             4.45*(Aexh  /0.3048^2) * (2.5+0.0363*Dfan/0.0254)
  
    Wfan = mfan*gee + Wnace1*0.5 #TODO need more recent weight model for nace - here just using 50% of the DC10 one

    return Dfan, Fan_power, Torque_fan, N_fan, Mtip,
            eta_prop, eta_DF,
            NozArea, Wfan, Snace1

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
function DuctedFan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64,  Pin::Float64,
                    Kinl::Float64, Φinl::Float64, first)
    
    # file_name = "NPSS_Turboshaft/Fan.output"

    NPSS_success, 
    Fn, Fan_power, Torque_fan, N_fan,
     Mtip, eta_prop, eta_DF = runNPSS_Fan(NPSS, alt_in, MN_in, Pin, Kinl, Φinl, first)
    
    # Read the data
    # include(file_name)

    return Fn, Fan_power, Torque_fan, N_fan, Mtip, eta_prop, eta_DF

end

function DuctedFan2(NPSS::Base.Process, alt_in::Float64, MN_in::Float64,  Fn::Float64,
    Kinl::Float64, Φinl::Float64, first)

# file_name = "NPSS_Turboshaft/Fan.output"

NPSS_success, 
Fn, Fan_power, Torque_fan, N_fan,
Mtip, eta_prop, eta_DF = runNPSS_Fan2(NPSS, alt_in, MN_in, Fn, Kinl, Φinl, first)

# Read the data
# include(file_name)

return Fn, Fan_power, Torque_fan, N_fan, Mtip, eta_prop, eta_DF

end
"""
PowerTrain 

Design method - sizes the powertrain

Inputs:
- NPSS_Fan -  NPSS process that is running the fan calculations
- alt_in, MN_in - Flight conditions (alt in m)
- Fn            - Required Thrust [N]
- Kinl, Φinl    - Ingested KE defect and disspation



"""
function PowerTrain(NPSS_TS::Base.Process, NPSS_Fan::Base.Process, NPSS_AftFan::Base.Process,
                    alt_in::Float64, MN_in::Float64, Fn::Float64,
                    Kinl::Float64, Φinl::Float64,
                    parg::Array{Float64, 1},
                    parpt::Array{Union{Float64, Int64},1},
                    parmot::Array{Float64, 1},
                    pargen::Array{Float64, 1}, first)
    
    # Initialize powertrain weight and waste heat 
        Wpowertrain = 0.0
       xWpowertrain = 0.0
        Hrej = 0.0
    # Unpack number of powertrain elements
        nfan    = parpt[ipt_nfan]
        ngen    = parpt[ipt_ngen]
        nTshaft = parpt[ipt_nTshaft]
        FnSplit = parpt[ipt_Fnsplit] # How much of the thrust comes from aft BLI fans

        FnAft = Fn*FnSplit
        rSnace = parg[igrSnace]
        #Assume 2 aft fans:
        DAftFan, AftFan_power, AftFan_torque, AftFan_N, AftFan_Mtip,
        AftFan_ηpropul, AftFan_η,
        AftFanNozArea, WAftfan, AftSnace1 = DuctedFan(NPSS_AftFan, alt_in, MN_in, FnAft/2, Kinl, Φinl, 1.25, rSnace, first )

        parg[igdaftfan] = DAftFan
        
        Pshaft_AftFan = -1000.0 * AftFan_power
        Wpowertrain += WAftfan*2.0
       xWpowertrain += WAftfan*2.0*parg[igxgen]
        # println(DAftFan)
        #Make these 0 for wings
        Kinl = 0.0
        Φinl = 0.0
    # Thrust per fan
        Ffan = (Fn - FnAft)/nfan

    # Call ducted fan design method
        πfan = parpt[ipt_pifan]
        NPSS_time = 0.0
        NPSS_time += @elapsed  Dfan, Fan_power, Torque_fan, N_fan, Mtip,
        ηpropul, ηDF,
        FanNozArea, Wfan, Snace1 = DuctedFan(NPSS_Fan, alt_in, MN_in, Ffan, Kinl, Φinl, πfan, rSnace, first )
        parpt[ipt_calls_NPSS] += 1
        # println("Fan:")
        # println("Fan Ø = ", Dfan, " m")
        # println("Fan Fn = ", Ffan, " N")
        parg[igdfan] = Dfan
        parg[igWfan] = Wfan

        parpt[ipt_Wfan]    = Wfan # Save fan weight
        parpt[ipt_NdesFan] = N_fan

        Wpowertrain += Wfan*nfan
       xWpowertrain += Wfan*nfan*parg[igxfan]

        Pshaft_mot = -1000. * Fan_power # convert to W (minus sign since NPSS output is power into fan which is negative)
        
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
        parg[igWmot] = Wmot

        parpt[ipt_Wmot]    = Wmot
        parpt[ipt_NdesMot] = RPMmot
        parpt[ipt_FanGR] = RPMmot/N_fan

        Hwaste_motor = PreqMot - Pshaft_mot
        Hrej += nfan*Hwaste_motor # Heat rejected from motors
        Wpowertrain += Wmot*nfan  # Add to total powertrain weight
       xWpowertrain += Wmot*nfan*parg[igxmot]           

    # Size Inverter and cables
        ηinv, Winv, SPinv = inverter(PreqMot, RPMmot/60, parmot)
        Hwaste_inv = PreqMot*(1-ηinv)
        Hrej += nfan * Hwaste_inv # Heat rejected from all inverters
        
        parg[igWinv] = Winv
        parpt[ipt_Winv] = Winv
        Wpowertrain += Winv*nfan # Add to total powertrain weight
       xWpowertrain += Winv*nfan*parg[igxinv]    

        ηcable, Wcable = cable() #TODO cable is dummy right now
        Hwaste_cable = PreqMot*(1-ηcable)
        Hrej += Hwaste_cable  # Heat rejected from all inverters
        
        parg[igWcables] = Wcable
        parpt[ipt_Wcables] = Wcable
        Wpowertrain += Wcable # Add to total powertrain weight
    #    xWpowertrain += Wcable*parg[igxcable] #need to add x of cable 
    # Size generator
        PreqGen = PreqMot /( ηinv * ηcable)

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
        
        parg[igWgen] = Wgen
        parpt[ipt_Wgen] = Wgen
        Wpowertrain += Wgen*ngen
       xWpowertrain += Wgen*ngen*parg[igxgen]  

    # Size tubo-shaft and PCEC
        πLPC = parpt[ipt_piLPC]
        πHPC = parpt[ipt_piHPC]
        Tt41 = parpt[ipt_Tt41]
        cpsi = parpt[ipt_cpsi]
        w    = parpt[ipt_wcat]
        lcat = parpt[ipt_lcat]
        deNOx= parpt[ipt_deNOx]

        Ptshaft = PgenShaft*ngen/nTshaft + Pshaft_AftFan # Assume one aft fan per Tshaft
        parpt[ipt_Ptshaft] = Ptshaft
        NPSS_time += @elapsed ηthermal, mdotf, BSFC,
         deNOx_out, mcat, EINOx1, EINOx2, mdot, Tt3, OPR, Wc3 = TurboShaft(NPSS_TS, alt_in, MN_in, Ptshaft,
                                            πLPC, πHPC, Tt41,
                                            cpsi, w, lcat, deNOx, first)
        
        parpt[ipt_calls_NPSS] += 1
        
        FAR = mdotf/mdot # fuel air ratio in each turboshaft engine

        mdotf_tot = mdotf*nTshaft
        Wcat = mcat*gee*1.5 # TODO replace fudge factor 1.5 with ammonia tanks/ pumps etc
        
        parg[igWcat] = Wcat
        parpt[ipt_Wcatalyst] = Wcat
        Wpowertrain += Wcat*nTshaft
       xWpowertrain += Wcat*nTshaft*parg[igxtshaft]

        # SPtshaft = SPadv(Ptshaft/hp_to_W)*hp_to_W*2.205 #returns in hp/lb so convert to W/kg
        SPtshaft= 10.4e3 # Alternate W/kg Based on the RR T406 (4.58 MW power output) power density = 10.4 kW/kg. The GE38 (~5 MW) has a power density of 11.2 kW/kg
        Wtshaft = gee*(Ptshaft)/SPtshaft
        
        parg[igWtshaft] = Wtshaft
        parpt[ipt_Wtshaft] = Wtshaft
        Wpowertrain += Wtshaft*nTshaft
       xWpowertrain += Wtshaft*nTshaft*parg[igxtshaft]

    parpt[ipt_Wpttotal] = Wpowertrain

    parg[ igWtesys] = Wpowertrain
    parg[igxWtesys] = xWpowertrain

    parpt[ipt_time_NPSS] += NPSS_time

    return [ηpropul, ηmot, ηinv, ηcable, ηgen, ηthermal],
           [Pshaft_mot, PreqMot, PgenShaft, Ptshaft], 
           [Hwaste_motor, Hwaste_inv, Hwaste_cable, Hwaste_gen, Hrej]./1000,
           [Wfan, Wmot, Winv, Wcable, Wgen, Wtshaft, Wcat, Wpowertrain]./gee,
           [SPmot, SPinv, SPgen, SPtshaft], mdotf_tot, BSFC,
           deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, FanNozArea, Snace1, AftSnace1

end

function PowerTrainOD(NPSS_TS::Base.Process, NPSS_Fan::Base.Process, NPSS_AftFan::Base.Process,
                    alt_in::Float64, MN_in::Float64, Tt41::Float64,
                    Kinl::Float64, Φinl::Float64,
                    parpt::Array{Union{Float64, Int64},1},
                    parmot::Array{Float64, 1},
                    pargen::Array{Float64, 1}, first, Ldebug)
    
    Hrej = 0.0
    nfan = parpt[ipt_nfan]
    ngen = parpt[ipt_ngen]
    nTshaft = parpt[ipt_nTshaft]
    Nshaft = 30000. #RPM
    NPSS_time = 0.0

    # Calculate Turboshaft power output
        NPSS_time += @elapsed Pshaft, ηthermal,
        mdotf, BSFC,
        deNOx_out = TurboShaft(NPSS_TS, alt_in, MN_in, Tt41, Nshaft, first)
        
        parpt[ipt_calls_NPSS] += 1

        mdotf_tot = mdotf*nTshaft

    #Aft fan
        PAftFan = Pshaft*parpt[ipt_Fnsplit]
        FnAft, AftFan_power, AftFan_torque, AftFan_N, AftFan_Mtip,
        AftFan_ηpropul, AftFan_ηDF = DuctedFan(NPSS_AftFan, alt_in, MN_in, PAftFan, Kinl, Φinl, first)

        Kinl = 0.0
        Φinl = 0.0

    # Run generator
        Pgen_in = Pshaft*(1 - parpt[ipt_Fnsplit]) * nTshaft/ngen # Assume power is split the same as the thrust split at design
        Ngen = Nshaft
        PLgen = PMSM(Pgen_in, Ngen/60, pargen)
        Pgen_out = Pgen_in - PLgen #don't directly use Pgen_out 
        ηgen = Pgen_out/Pgen_in
        Hwaste_gen = PLgen
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

    ## Motor and fan speed need to match through the gearing ratio - currently this is not guaranteed 
    ##  during off-design operation. ∴ we iterate. Usually this should converge in ≈5 iterations. 
    ##  This is def not ideal - it would be more efficient to converge the speeds via Newton simultaneously with the Fan model.

    # Off-des motor
        Nmot = parpt[ipt_NdesMot] #* (N_fan/parpt[ipt_NdesFan])
        # Off-des motor call
        PL = PMSM(Pmot_in, Nmot/60, parmot)
        Pmot_out = Pmot_in - PL
        ηmot = Pmot_out/Pmot_in
        Hwaste_motor = PL
        Hrej += nfan*Hwaste_motor # Heat rejected from motors
        # println("Motor speed = ", Nmot)
        # println("Motor power = ", Pmot_out)

    # Ducted fan
    Pfan_in = Pmot_out

    NPSS_time += @elapsed Fn, Fan_power, Torque_fan, N_fan, Mtip,
        ηpropul, ηDF = DuctedFan(NPSS_Fan, alt_in, MN_in, Pfan_in, Kinl, Φinl, first)

        parpt[ipt_calls_NPSS] += 1

        err = 1
        maxiter = 10
        GR = parpt[ipt_FanGR]
        tol = 1e-5

        if Ldebug
            printstyled("\tMotor speed convergence = \n"; color=:blue)
            printstyled(@sprintf("\t\t%5s  %10s  %10s  %10s  %10s  %10s  %10s  \n",
            "iter", "err", "Nmot/GR", "Nfan", "ηmot", "ηfan", "Mtip"); color =:blue)
        end
        
        #iterate to converge fan and mot speeds
        @inbounds for  i = 1:maxiter

            err = (N_fan*GR - Nmot)/Nmot
            
            Ldebug && printstyled(@sprintf("\t\t%5d  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e\n",
                                            i, abs(err), Nmot/GR, N_fan, ηmot, ηDF, Mtip); color = :blue)

            if abs(err)≤tol
                break;
            end
            Nmot = N_fan*GR
            # Off-des motor call
            PL = PMSM(Pmot_in, Nmot/60, parmot)
            Pmot_out = Pmot_in - PL
            ηmot = Pmot_out/Pmot_in
            Hwaste_motor = PL
            Hrej += nfan*Hwaste_motor

            # Ducted fan
            Pfan_in = Pmot_out

        NPSS_time += @elapsed Fn, Fan_power, Torque_fan, N_fan, Mtip,
                ηpropul, ηDF = DuctedFan(NPSS_Fan, alt_in, MN_in, Pfan_in, Kinl, Φinl, first)
        
        parpt[ipt_calls_NPSS] += 1
        end

        abs(err)≥tol && printstyled("Warning [propsys]: Motor-fan speeds not converged! Error = ",abs(err),"\n"; color=:red)

    Ftotal = Fn*nfan

    parpt[ipt_time_NPSS] += NPSS_time
    return Ftotal, [ηmot, ηinv, ηcable, ηgen, ηthermal],
           [Pmot_out, Pmot_in, Pinv_in, Pgen_in, Pshaft], 
           [Hwaste_motor, Hwaste_inv, Hwaste_cable, Hwaste_gen, Hrej]./1000,
           mdotf_tot, BSFC,
           deNOx_out

end


function PowerTrainOD2(NPSS_TS::Base.Process, NPSS_Fan::Base.Process,
    alt_in::Float64, MN_in::Float64, Fn::Float64,
    Kinl::Float64, Φinl::Float64,
    parpt::Array{Union{Float64, Int64},1},
    parmot::Array{Float64, 1},
    pargen::Array{Float64, 1}, first, Ldebug)

Hrej = 0.0
nfan = parpt[ipt_nfan]
ngen = parpt[ipt_ngen]
nTshaft = parpt[ipt_nTshaft]
Nshaft = 30000. #RPM
NPSS_time = 0.0

# Ducted fan

# Thrust per fan
Ffan = Fn/nfan

NPSS_time += @elapsed Fn, Fan_power, Torque_fan, N_fan, Mtip,
ηpropul, ηDF = DuctedFan2(NPSS_Fan, alt_in, MN_in, Ffan, Kinl, Φinl, first)

parpt[ipt_calls_NPSS] += 1
Pmot_out = -1000.0 * Fan_power

GR = parpt[ipt_FanGR]
Nmot = N_fan*GR

# Off-des motor call
Pmot_in, ηmot, PL = PMSM(Pmot_out, Nmot/60, parmot)

Hwaste_motor =  Pmot_in - Pmot_out
Hrej += nfan*Hwaste_motor # Heat rejected from motors
    
Pinv_out = Pmot_in
# Off-des inverter
ηinv = inverter(Pinv_out, parmot)
Hwaste_inv = Pinv_out*(1-ηinv)
Hrej += nfan * Hwaste_inv # Heat rejected from all inverters
Pinv_in = Pinv_out / ηinv

# Off-des cable
ηcable, Wcable = cable() #TODO cable is dummy right now
Hwaste_cable = Pinv_in*(1-ηcable)
Hrej += Hwaste_cable  # Heat rejected from cables
Pgen_out = Pinv_in/ηcable * ngen/nfan # for each inverter

# Run generator
Ngen = Nshaft
Pgen_in, ηgen, PLgen = PMSM(Pgen_out, Ngen/60, pargen)

Hwaste_gen = (Pgen_in - Pgen_out)
Hrej += ngen*Hwaste_gen # Heat rejected from motors

Pshaft = Pgen_in * nTshaft/ngen
println(Pshaft)
# Calculate Turboshaft power output
NPSS_time += @elapsed Pshaft, ηthermal,
mdotf, BSFC,
deNOx_out = TurboShaft2(NPSS_TS, alt_in, MN_in, Pshaft, Nshaft, first)

parpt[ipt_calls_NPSS] += 1

mdotf_tot = mdotf*nTshaft



parpt[ipt_time_NPSS] += NPSS_time

return [ηmot, ηinv, ηcable, ηgen, ηthermal],
[Pmot_out, Pmot_in, Pinv_in, Pgen_in, Pshaft], 
[Hwaste_motor, Hwaste_inv, Hwaste_cable, Hwaste_gen, Hrej]./1000,
mdotf_tot, BSFC,
deNOx_out

end

"""
Specific power for advanced turboshaft engines
Curve fit for SP of turboshaft engines based on https://ntrs.nasa.gov/api/citations/20190025407/downloads/20190025407.pdf

"""
function SPadv(Php)
    -1.185e-7*Php^2 + 2.008e-3*Php + 1.608
end
