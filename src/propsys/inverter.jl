"""
$TYPEDEF

Simplified model of an inverter.

$TYPEDFIELDS
"""
@kwdef mutable struct Inverter
    """Mass [kg]"""
    mass::AbstractFloat = 0.0
    """Design power [W]"""
    P_design::AbstractFloat = 1e3 
    """Design Rotational speed [RPM]"""
    N::AbstractFloat = 1e3
    """Specific power [W/kg]"""
    specific_power::AbstractFloat = 19e3
    """Carrier frequency constant"""
    kcf::AbstractFloat = 20.0
end

"""
    inverter(P::Float64, N::Float64)


Simple inverter model that calculates the efficiency and mass of an inverter or rectifier.
Default specific power is 19 kW/kg per [GE and Uni Illinois using SiC and GaN switches respectively.](https://arc.aiaa.org/doi/pdf/10.2514/6.2017-4701)
"""
function inverter(inverter::Inverter, P::Float64, f::Float64)
    kcf = inverter.kcf
    SPinv  = inverter.specific_power
    P_design = inverter.P_design

    fSwitch = f*kcf
    η100 = -2.5e-7*fSwitch + 0.995
    η20  = η100 - 0.0017
    η10  = η100 - 0.0097

    M = [1.0 0.1 10.0
         1.0 0.2  5.0
         1.0 1.0  1.0]

    k1, k2, k3 = M\[η10; η20; η100]

    η = k1 + k2*(P/P_design) + k3*(P/P_design)^-1 #but here P = Pmax at design

    parte[ite_Pinvdes] = P
    
    Winverter = P/SPinv * gee # Weight in [N]

    return η, Winverter, SPinv

end
"""
Off design method for inverter
"""
function inverter(P::Float64, parte::Array{Float64, 1})

    k1 = parte[ite_k1]
    k2 = parte[ite_k2]
    k3 = parte[ite_k3]

    Pdes = parte[ite_Pinvdes]
    
    η = k1 + k2*(P/Pdes) + k3*(P/Pdes)^-1

    return η

end