"""
$TYPEDEF

Simplified model of an inverter.

$TYPEDFIELDS
"""
@kwdef mutable struct inverter
    """Design power [W]"""
    P::AbstractFloat = 1e3 
    """Number of pole-pairs"""
    N_poles::Integer = 2
    """Design Rotational speed [RPM]"""
    N::AbstractFloat
    
end

"""
    inverter(P::Float64, N::Float64, parte::Array{Float64, 1})


Simple inverter model that calculates the efficiency and mass of an inverter or rectifier.
Default specific power is 19 kW/kg per [GE and Uni Illinois using SiC and GaN switches respectively.](https://arc.aiaa.org/doi/pdf/10.2514/6.2017-4701)
"""
function inverter(P::Float64, N::Float64, parte::Array{Float64, 1})
    kcf = 20.
    SPinv  = 19.0e3 #W/kg 

    p = parte[ite_p]

    f = p * N
    fSwitch = f*kcf
    η100 = -2.5e-7*fSwitch + 0.995
    η20  = η100 - 0.0017
    η10  = η100 - 0.0097

    M = [1.0 0.1 10.0
         1.0 0.2  5.0
         1.0 1.0  1.0]

    k1, k2, k3 = M\[η10; η20; η100]

    # η = k1 + k2*(P/Pmax) + k3*(P/Pmax)^-1 but here P = Pmax at design
    η = k1 + k2*1 + k3/1

    parte[ite_k1] = k1
    parte[ite_k2] = k2
    parte[ite_k3] = k3

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