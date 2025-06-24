"""
$TYPEDEF

Simplified model of an inverter.
Default specific power is 19 kW/kg per [GE and Uni Illinois using SiC and GaN switches respectively.](https://arc.aiaa.org/doi/pdf/10.2514/6.2017-4701)

$TYPEDFIELDS
"""
@kwdef mutable struct Inverter
    """Mass [kg]"""
    mass::AbstractFloat = 0.0
    """Design power [W]"""
    P_design::AbstractFloat = 1e3 
    """Input power [W]"""
    P_input::AbstractFloat = 1e3 
    """Output power [W]"""
    P::AbstractFloat = 1e3 
    """Design Rotational speed [RPM]"""
    N::AbstractFloat = 1e3
    """Specific power [W/kg]"""
    specific_power::AbstractFloat = 19e3
    """Carrier frequency constant"""
    kcf::AbstractFloat = 20.0
end #Inverter

"""
    size_inverter!(inverter::Inverter, P_design::Float64, f_design::Float64)

Simple inverter model that calculates the efficiency and mass of an inverter or rectifier.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `inverter::Inverter`: inverter object
    - `P_design::Float64`: design output electrical power [W]
    - `f_design::Float64`: design output electrical frequency [Hz]
    
    **Outputs:**
    No direct outputs. Inverter object gets modified with the input power and mass.
"""
function size_inverter!(inverter::Inverter, P_design::Float64, f_design::Float64)
    SPinv  = inverter.specific_power
    inverter.P_design = P_design
    inverter.P = P_design

    operate_inverter!(inverter, P_design, f_design)
    
    #Mass is calculated from the specified specific power
    inverter.mass = P_design/SPinv #mass in kg
end #size_inverter!

"""
    operate_inverter!(inverter::Inverter, P::Float64, f::Float64)

Off design method for inverter. The model is based on Faranda et al., 2015 (10.3390/en8064853).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `inverter::Inverter`: inverter object
    - `P::Float64`: output electrical power [W]
    - `f::Float64`: output electrical frequency [Hz]
    
    **Outputs:**
    No direct outputs. Inverter object gets modified with the input power.
"""
function operate_inverter!(inverter::Inverter, P::Float64, f::Float64)
    P_design = inverter.P_design
    inverter.P = P
    kcf = inverter.kcf

    #The linear scaling in efficiency with frequency is from Enders, 2020 (Master's thesis)
    fSwitch = f*kcf
    #TODO: it is unclear where these values come from; the thesis cites a paper but
    #the data in the paper does not quite match these values.
    η100 = -2.5e-7*fSwitch + 0.995
    η20  = η100 - 0.0017
    η10  = η100 - 0.0097

    #Data from Faranda et al.
    M = [1.0 0.1 10.0
         1.0 0.2  5.0
         1.0 1.0  1.0]

    #Find coefficients for the fit
    k1, k2, k3 = M\[η10; η20; η100]

    η = k1 + k2*(P/P_design) + k3*(P/P_design)^-1
    inverter.P_input = P/η #input power in W
end #operate_inverter!