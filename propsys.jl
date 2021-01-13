
include("PMSM.jl")  # Motor/generator functions
include("PMSM.inc") # Motor/generator properties array
include("NPSS_functions.jl") # NPSS functions

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
function DuctedFan(alt_in::Float64, MN_in::Float64,  Fn::Float64, π_fan::Float64,)
    
    NPSS_Fan_input(alt_in, MN_in, Fn, π_fan)
    NPSS_run("NPSS_Turboshaft/", "Fan.bat")

    include("NPSS_Turboshaft/Fan.output")

    return Dfan, Fan_power, Torque_fan, N_fan, eta_prop, MapScalars, NozArea

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
function DuctedFan(alt_in::Float64, MN_in::Float64,  Fn::Float64, MapScalars::Array{Float64, 1}, NozArea::Float64)
    
    NPSS_Fan_input(alt_in, MN_in, Fn, MapScalars, NozArea)
    NPSS_run("NPSS_Turboshaft/", "Fan.bat")

    include("NPSS_Turboshaft/Fan.output")

    return Fan_power, Torque_fan, N_fan, eta_prop

end

"""
PowerTrain 

Design method - sizes the powertrain
"""
function PowerTrain(alt_in::Float64, MN_in::Float64, Fn::Float64, parpt::Array{Union{Float64, Int64},1}, parte::Array{Float64, 1})
    
    Wpowertrain = 0.0
    
    Ffan = Fn/parpt[ipt_neng]
    # Call ducted fan design method
    πfan = parpt[ipt_pifan]
    Dfan, Fan_power, Torque_fan, N_fan, ηpropul, MapScalars, NozArea = DuctedFan(alt_in, MN_in, Ffan, πfan)
    # println("Fan:")
    # println("Fan Ø = ", Dfan, " m")
    # println("Fan Fn = ", Ffan, " N")

    Pshaft_mot = -1000. * Fan_power

    ratAsp   = parpt[ipt_ARmot]
    σAg      = parpt[ipt_sigAgMot]
    ratSplit = parpt[ipt_ratSplitMot]

    Wpmsm, PreqMot, ηmot, RPMmot, PL, PLiron, PLCu, PLwind, SP = PMSM(Pshaft_mot, ratAsp, σAg, ratSplit, parte)
    # println("\nMotor:")
    # println("Losses = ", PL)
    # println("\tIron losses    = ", PLiron, "%")
    # println("\tCopper losses  = ", PLCu)
    # println("\tWindage losses = ", PLwind)

    Wpowertrain += Wpmsm*neng # Add to total powertrain weight
    
    ηinv, Winv = inverter(PreqMot, RPMmot/60, parte)
    Wpowertrain += Winv*neng # Add to total powertrain weight
    
    ηcable, Wcable = cable()
    Wpowertrain += Wcable # Add to total powertrain weight

    PreqGen = PreqMot * ηinv * ηcable

    ratAsp   = parpt[ipt_ARgen]
    σAg      = parpt[ipt_sigAgGen]
    ratSplit = parpt[ipt_ratSplitGen]

    Wgen, PgenShaft, ηgen, RPMgen, PLgen, PLirongen, PLCugen, PLwindgen, SPgen = PMSM(PreqGen*neng/ngen, ratAsp, σAg, ratSplit, parte)
    Wpowertrain += Wgen*ngen

    ηthermal, mdotf, BSFC = TurboShaft(alt_in, MN_in, PgenShaft*ngen/nTshaft, 3.0, 6.0)

    return [ηmot, ηpropul, ηinv, ηcable, ηgen, ηthermal], Pshaft_mot, PreqMot, PgenShaft, Fn, Wpmsm, Winv, Wcable, Wgen, SP, SPgen, mdotf, BSFC 

end