
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
function TurboShaft(alt_in::Float64, MN_in::Float64,  ShP::Float64, π_LPC::Float64, π_HPC::Float64)
    
    ShP_hp = ShP / 746 # Convert shaft power from W to HP
    NPSS_TShaft_input(alt_in, MN_in, ShP_hp, 3200, π_LPC, π_HPC )
    NPSS_run("NPSS_Turboshaft/", "TP.bat")

    include("NPSS_Turboshaft/Eng.output")
    
    return eta_thermal, mdotf, BSFC  #, MapScalars, NozArea

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
function PowerTrain(alt_in::Float64, MN_in::Float64, Fn::Float64, neng::Int64, ngen::Int64, π_fan::Float64, parte::Array{Float64, 1})
    
    Wpowertrain = 0.0
    Ffan = Fn/neng
    # Call ducted fan design method
    Dfan, Fan_power, Torque_fan, N_fan, ηpropul, MapScalars, NozArea = DuctedFan(alt_in, MN_in, Ffan, π_fan)
    # println("Fan:")
    # println("Fan Ø = ", Dfan, " m")
    # println("Fan Fn = ", Ffan, " N")

    Pshaft_mot = -1000. * Fan_power

    ratAsp   = 0.8
    σAg      = 40e3
    ratSplit = 0.7

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

    ratAsp   = 0.6
    σAg      = 48e3
    ratSplit = 0.8

    Wgen, PgenShaft, ηgen, RPMgen, PLgen, PLirongen, PLCugen, PLwindgen, SPgen = PMSM(PreqGen*neng/ngen, ratAsp, σAg, ratSplit, parte)
    Wpowertrain += Wgen*ngen

    ηthermal, mdotf, BSFC = TurboShaft(alt_in, MN_in, PgenShaft, 3.0, 6.0)


    return [ηmot, ηpropul, ηinv, ηcable, ηgen, ηthermal], Pshaft_mot, PreqMot, PgenShaft, Fn, Wpmsm, Winv, Wcable, Wgen, SP, SPgen, mdotf, BSFC 

end