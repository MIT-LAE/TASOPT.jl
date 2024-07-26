"""
Structure with thermodynamic parameters for the vapor portion of a saturated mixture.
"""
mutable struct SaturatedGas
    """Species name"""
    species::String
    """Saturated mixture temperature (K)"""
    T::Float64
    """Saturated mixture pressure (Pa)"""
    p::Float64
    """Gas density (kg/m^3)"""
    ρ::Float64
    """Derivative of density with pressure (kg/m^3/Pa)"""
    ρ_p::Float64
    """Gas specific enthalpy (J/kg)"""
    h::Float64
    """Gas specific internal energy (J/kg)"""
    u::Float64
    """Derivative of internal energy with pressure (J/kg/Pa)"""
    u_p::Float64
end

"""
Structure with thermodynamic parameters for the liquid portion of a saturated mixture.
"""
mutable struct SaturatedLiquid
    """Species name"""
    species::String
    """Saturated mixture temperature (K)"""
    T::Float64
    """Saturated mixture pressure (Pa)"""
    p::Float64
    """Liquid density (kg/m^3)"""
    ρ::Float64
    """Derivative of density with pressure (kg/m^3/Pa)"""
    ρ_p::Float64
    """Gas specific enthalpy (J/kg)"""
    h::Float64
    """Gas specific internal energy (J/kg)"""
    u::Float64
    """Derivative of internal energy with pressure (J/kg/Pa)"""
    u_p::Float64
end

"""
Structure with thermodynamic properties of a saturated mixture.
"""
mutable struct SaturatedMixture
    """Saturated vapor properties"""
    gas::SaturatedGas
    """Saturated liquid properties"""
    liquid::SaturatedLiquid
    """Species name"""
    species::String
    """Liquid fill fraction"""
    β::Float64
    """Vapor quality"""
    x::Float64
    """Saturated mixture temperature (K)"""
    T::Float64
    """Saturated mixture pressure (Pa)"""
    p::Float64
    """Mixture density (kg/m^3)"""
    ρ::Float64
    """Mixture specific enthalpy (J/kg)"""
    h::Float64
    """Mixture specific internal energy (J/kg)"""
    u::Float64
    """Derivative of internal energy with pressure (J/kg/Pa)"""
    u_p::Float64
    """Energy derivative"""
    ϕ::Float64
    """Density ratio, `ρl/(ρl-ρg)`"""
    ρ_star::Float64
    """Enthalpy of vaporization (J/kg)"""
    hvap::Float64
end

"""
    SaturatedGas(species::String, p::Float64)

This function produces a saturated gas at a desired pressure.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: pressure (Pa)
    
    **Outputs:**
    - `gas::SaturatedGas`: saturated gas
"""
function SaturatedGas(species::String, p::Float64)

    Tsat, ρ, ρ_p, h, u, u_p = gas_properties(species, p)

    gas = SaturatedGas(species, Tsat, p, ρ, ρ_p, h, u, u_p)
    return gas
end

"""
    SaturatedGas(species::String, p::Float64)

This function produces a saturated liquid at a desired pressure.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: pressure (Pa)
    
    **Outputs:**
    - `liquid::SaturatedLiquid`: saturated liquid
"""
function SaturatedLiquid(species::String, p::Float64)

    Tsat, ρ, ρ_p, h, u, u_p = liquid_properties(species, p)

    liquid = SaturatedLiquid(species, Tsat, p, ρ, ρ_p, h, u, u_p)
    return liquid
end

"""
    calculate_mixture(gas, liquid, β)

This function calculates the thermodynamic properties of a saturated mixture from those of its gas and liquid phases.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gas::SaturatedGas`: saturated gas
    - `liquid::SaturatedLiquid`: saturated liquid
    - `β::Float64`: liquid fill ratio
    
    **Outputs:**
    - `x::Float64`: vapor quality 
    - `ρ::Float64`: mixture density (kg/m^3)
    - `h::Float64`: mixture specific enthalpy (J/kg)
    - `u::Float64`: mixture specific internal energy (J/kg)
    - `u_p::Float64`: derivative of internal energy with pressure at constant density (J/kg/Pa)
    - `ϕ::Float64`: energy derivative
    - `ρ_star::Float64`: density ratio, `ρl/(ρl-ρg)`
    - `hvap::Float64`: enthalpy of vaporization (J/kg)
"""
function calculate_mixture(gas::SaturatedGas, liquid::SaturatedLiquid, β::Float64)
    x = 1 / (1 + (liquid.ρ / gas.ρ) * (β / (1 - β)))
    ρ = 1 / (x/gas.ρ + (1 - x)/liquid.ρ)
    h = x * gas.h + (1 - x) * liquid.h
    u = x * gas.u + (1 - x) * liquid.u

    #Calculate ∂x/∂p at constant ρ for ∂u/∂p at constant ρ calculation
    ρl = liquid.ρ
    ρl_p = liquid.ρ_p
    ρg = gas.ρ
    ρg_p = gas.ρ_p
    
    x_p = (-x/ρg^2 * ρg_p - (1-x)/ρl^2 * ρl_p) / (1/ρl - 1/ρg)

    #Find ∂u/∂p at constant ρ by chain rule
    u_p = x * gas.u_p + (1 - x) * liquid.u_p + x_p * (gas.u - liquid.u)

    ϕ = 1 / (ρ * u_p)
    ρ_star = gas.ρ / (liquid.ρ - gas.ρ)
    hvap = gas.h - liquid.h
    return x, ρ, h, u, u_p, ϕ, ρ_star, hvap
end

"""
    SaturatedMixture(species::String, p::Float64, β::Float64)

This function produces a saturated mixture at a desired pressure for a given liquid fill volume ratio.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: pressure (Pa)
    - `β::Float64`: liquid fill ratio
    
    **Outputs:**
    - `mixture::SaturatedMixture`: saturated mixture
"""
function SaturatedMixture(species::String, p::Float64, β::Float64)

    gas = SaturatedGas(species, p)
    liquid = SaturatedLiquid(species, p)

    x, ρ, h, u, u_p, ϕ, ρ_star, hvap = calculate_mixture(gas, liquid, β)

    mixture = SaturatedMixture(gas, liquid, species, β, x, gas.T, p, ρ, h, u, u_p, ϕ, ρ_star, hvap)
    return mixture
end

"""
update_pβ!(mixture, p, β)

This function updates a saturated mixture when there is a change in pressure or liquid fill volume ratio.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mixture::SaturatedMixture`: saturated mixture
    - `p::Float64`: pressure (Pa)
    - `β::Float64`: liquid fill ratio
    
    **Outputs:**
    No direct outputs. The `mixture` input is modified.
"""
function update_pβ!(mixture::SaturatedMixture, p::Float64, β::Float64)

    gas = SaturatedGas(mixture.species, p)
    liquid = SaturatedLiquid(mixture.species, p)

    x, ρ, h, u, u_p, ϕ, ρ_star, hvap = calculate_mixture(gas, liquid, β)

    mixture.gas = gas
    mixture.liquid = liquid
    mixture.T = gas.T
    mixture.β = β
    mixture.p = p
    mixture.x = x
    mixture.ρ = ρ
    mixture.h = h
    mixture.u = u
    mixture.u_p = u_p
    mixture.ϕ = ϕ
    mixture.ρ_star = ρ_star
    mixture.hvap = hvap
end

"""
    convert_β_same_ρ(species::String, p::Float64, p0::Float64, β0::Float64)

This function calculates the change in fill fraction for a change of pressure of a mixture at constant mixture density.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: new pressure (Pa)
    - `p0::Float64`: original pressure (Pa)
    - `β0::Float64`: original liquid fill ratio
    
    **Outputs:**
    - `β::Float64`: new liquid fill ratio
"""
function convert_β_same_ρ(species::String, p::Float64, p0::Float64, β0::Float64)

    mix0 = SaturatedMixture(species, p0, β0)

    gas = SaturatedGas(species, p)
    liq = SaturatedLiquid(species, p)

    β = (β0 * mix0.liquid.ρ + (1 - β0) * mix0.gas.ρ - gas.ρ) / (liq.ρ - gas.ρ)
    return β
end