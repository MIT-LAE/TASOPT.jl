"""
# Ducted Fan calculations

## Design pt. method

Inputs: 
- alt_in : Ambient altitude 
- MN_in  : Ambient Mach number
- Fn     : Design Thurst reqired 
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