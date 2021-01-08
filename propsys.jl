
include("PMSM.jl")  # Motor/generator functions
include("PMSM.inc") # Motor/generator properties array

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
function PowerTrain(alt_in::Float64, MN_in::Float64, Fn::Float64, π_fan::Float64, parte::Array{Float64, 1})

    # Call ducted fan design method
    Dfan, Fan_power, Torque_fan, N_fan, ηpropul, MapScalars, NozArea = DuctedFan(alt_in, MN_in, Fn, π_fan)

    Pshaft_mot = -1000. * Fan_power

    # ratAsp   = 0.8
    # σAg      = 20e3
    # ratSplit = 0.86
    Wpmsm, Preq, ηmot, RPMmot, _, _, _, _, SP = PMSM(Pshaft_mot, ratAsp, σAg, ratSplit, parte)


    return ηmot, ηpropul, Pshaft_mot, Preq, Fn, Wpmsm, SP

end