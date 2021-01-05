"""
# Ducted Fan calculations

## Design pt. method

Inputs: 
- alt_in : Ambient altitude 
- MN_in  : Ambient Mach number
- Fn     : Design Thurst reqired 
- π_fan  : Fan pressure ratio

"""
function DuctedFan(alt_in::Float64, MN_in::Float64,  Fn::Float64, π_fan::Float64,)
    Design = true
    NPSS_Fan_input(alt_in, MN_in, Fn, π_fan)
    NPSS_run("NPSS_Turboshaft/", "Fan.bat")

    include("NPSS_Turboshaft/Fan.output")

    return Design, Dfan, Fan_power, Torque_fan, N_fan, MapScalars, NozArea

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

"""
function DuctedFan(alt_in::Float64, MN_in::Float64,  Fn::Float64, MapScalars::Array{Float64, 1}, NozArea::Float64)
    Design = false
    NPSS_Fan_input(alt_in, MN_in, Fn, MapScalars, NozArea)
    NPSS_run("NPSS_Turboshaft/", "Fan.bat")

    include("NPSS_Turboshaft/Fan.output")

    return Design, Dfan, Fan_power, Torque_fan, N_fan

end