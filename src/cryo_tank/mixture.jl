"""
Structure with thermodynamic properties of a saturated mixture.
"""
mutable struct SaturatedMixture
    """Gas phase properties"""
    gas::SaturatedPhaseProps
    """Liquid phase properties"""
    liquid::SaturatedPhaseProps
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
    calculate_mixture!(mixture)

Update the derived thermodynamic properties of a saturated mixture in-place.
Reads `mixture.gas`, `mixture.liquid`, and `mixture.β`; writes `mixture.x`,
`mixture.ρ`, `mixture.h`, `mixture.u`, `mixture.u_p`, `mixture.ϕ`,
`mixture.ρ_star`, and `mixture.hvap`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mixture::SaturatedMixture`: mixture with `gas`, `liquid`, and `β` already set

    **Outputs:**
    No direct outputs. The derived fields of `mixture` are updated in-place.
"""
function calculate_mixture!(mixture::SaturatedMixture)
    gas    = mixture.gas
    liquid = mixture.liquid
    β      = mixture.β

    x = 1 / (1 + (liquid.ρ / gas.ρ) * (β / (1 - β)))
    ρ = 1 / (x/gas.ρ + (1 - x)/liquid.ρ)
    h = x * gas.h + (1 - x) * liquid.h
    u = x * gas.u + (1 - x) * liquid.u

    #Calculate ∂x/∂p at constant ρ for ∂u/∂p at constant ρ calculation
    ρl   = liquid.ρ
    ρl_p = liquid.ρ_p
    ρg   = gas.ρ
    ρg_p = gas.ρ_p

    x_p = (-x/ρg^2 * ρg_p - (1-x)/ρl^2 * ρl_p) / (1/ρl - 1/ρg)

    #Find ∂u/∂p at constant ρ by chain rule
    u_p = x * gas.u_p + (1 - x) * liquid.u_p + x_p * (gas.u - liquid.u)

    mixture.x      = x
    mixture.ρ      = ρ
    mixture.h      = h
    mixture.u      = u
    mixture.u_p    = u_p
    mixture.ϕ      = 1 / (ρ * u_p)
    mixture.ρ_star = gas.ρ / (liquid.ρ - gas.ρ)
    mixture.hvap   = gas.h - liquid.h
end

"""
    SaturatedMixture(gas, liquid, species, p, β)

Construct a `SaturatedMixture` from its defining inputs. All derived thermodynamic
properties (`x`, `ρ`, `h`, `u`, `u_p`, `ϕ`, `ρ_star`, `hvap`) are computed from
`gas`, `liquid`, and `β`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gas::SaturatedPhaseProps`: saturated gas phase properties
    - `liquid::SaturatedPhaseProps`: saturated liquid phase properties
    - `species::String`: species name
    - `p::Float64`: pressure (Pa)
    - `β::Float64`: liquid fill ratio

    **Outputs:**
    - `mixture::SaturatedMixture`: fully initialised saturated mixture
"""
function SaturatedMixture(gas::SaturatedPhaseProps, liquid::SaturatedPhaseProps,
                           species::String, p::Float64, β::Float64)
    x = 1 / (1 + (liquid.ρ / gas.ρ) * (β / (1 - β)))
    ρ = 1 / (x/gas.ρ + (1 - x)/liquid.ρ)
    h = x * gas.h + (1 - x) * liquid.h
    u = x * gas.u + (1 - x) * liquid.u

    ρl   = liquid.ρ;  ρl_p = liquid.ρ_p
    ρg   = gas.ρ;     ρg_p = gas.ρ_p
    x_p  = (-x/ρg^2 * ρg_p - (1-x)/ρl^2 * ρl_p) / (1/ρl - 1/ρg)
    u_p  = x * gas.u_p + (1 - x) * liquid.u_p + x_p * (gas.u - liquid.u)

    return SaturatedMixture(gas, liquid, species, β, x, gas.Tsat, p,
                             ρ, h, u, u_p, 1/(ρ*u_p), gas.ρ/(liquid.ρ-gas.ρ), gas.h-liquid.h)
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
    return SaturatedMixture(gas_properties(species, p), liquid_properties(species, p),
                             species, p, β)
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
    mixture.gas    = gas_properties(mixture.species, p)
    mixture.liquid = liquid_properties(mixture.species, p)
    mixture.T      = mixture.gas.Tsat
    mixture.β      = β
    mixture.p      = p
    calculate_mixture!(mixture)
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

    gas = gas_properties(species, p)
    liq = liquid_properties(species, p)

    β = (β0 * mix0.liquid.ρ + (1 - β0) * mix0.gas.ρ - gas.ρ) / (liq.ρ - gas.ρ)
    return β
end
