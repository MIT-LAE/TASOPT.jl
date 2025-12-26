module atmosphere
using Roots, Printf
export atmos, find_altitude_from_density, AtmosphericState

include("../utils/constants.jl")

"""
    AtmosphericState{T}

Atmospheric conditions at a point computed from the `atmos()` function.

# Fields
- `alt`: Altitude in meters
- `T`: Static temperature in Kelvin
- `p`: Static pressure in Pascals
- `ρ`: Density in kg/m³
- `a`: Speed of sound in m/s
- `μ`: Dynamic viscosity in kg/(m·s)

# Constructors
- `AtmosphericState(alt)`: Compute from altitude using ISA
- `AtmosphericState(alt, ΔT)`: Compute from altitude with temperature offset
"""
struct AtmosphericState{T<:Real}
    alt::T  # Altitude in m
    T::T    # Temperature in K
    p::T    # Pressure in Pa
    ρ::T    # Density in kg/m^3
    a::T    # Speed of sound in m/s
    μ::T    # Dynamic viscosity in kg/m-s
    ΔT::T   # Temperature offset from ISA in K
    function AtmosphericState(alt::T, ΔT::T=0.0) where T<:Real
        temp, p, ρ, a, μ = atmos(T(alt/1000.0), ΔT)
        new{T}(alt, temp, p, ρ, a, μ, ΔT)
    end
end

function Base.show(io::IO, ::MIME"text/plain", atm::AtmosphericState)
    # Print a clear, formatted header
    println(io, "AtmosphericState @ ", atm.alt, " m [", 
    round(atm.alt/ft_to_m; sigdigits=3), " ft]:")

    # Use @printf for aligned, formatted output
    @printf(io, "  %-22s %-12.2f K\n", "Temperature (T):", atm.T)
    @printf(io, "  %-22s %-12.0f Pa\n", "Pressure (p):", atm.p)
    @printf(io, "  %-22s %-12.4f kg/m³\n", "Density (ρ):", atm.ρ)
    @printf(io, "  %-22s %-12.2f m/s\n", "Speed of Sound (a):", atm.a)
    @printf(io, "  %-22s %-12.3e kg/m-s\n", "Dynamic Viscosity (μ):", atm.μ)
end

function Base.show(io::IO, atm::AtmosphericState)
    print(io, "AtmosphericState(alt=$(round(atm.alt/1000, digits=2)) km, ",
              "T=$(round(atm.T, sigdigits=3)) K, ",
              "p=$(round(atm.p/1000, sigdigits=3)) kPa, ",
              "ρ=$(round(atm.ρ, sigdigits=2)) kg/m³, ",
              "a=$(round(atm.a, sigdigits=3)) m/s, ",
              "μ=$(round(atm.μ*1e6, sigdigits=4))×10⁻⁶ kg/(m·s))")
end

"""
    atmos(h, ΔT)
    
Atmospheric functions ` T(h)`, `ρ(h)` etc
valid to `h`=20km, `p(h)` valid to `h`=70km.

Also calculates viscosity using Sutherland's law. Non-standard sea-level temperatures are allowed
with an ISA + ΔT like model.

Units:
- [h]   = km ASL
- [T]   = Kelvin
- [p]   = Pa
- [ρ]   = kg/m^3
- [a]   = m/s
- [μ]   = kg/m-s 
"""
function atmos(h::Float64, ΔT::Float64 = 0.0) 
#TODO: We probs don't need this anymore with the AtmosphericState struct but retain for now. 

 pSL = 1.0132e5 # Pa
 TSL_std = 288.2   # K, standard sea level temperature
 Tpause = 216.65   # K
 Tblend = 2.0      # K (tropopause blending T range)
 Tlapse = -6.5     # K/km
 
 μSL  = 1.78e-5  # kg/m-s
 Tsuth = 110.0    # K
 
 cp  = 1004.0   # J/kg-K
 ɣ = 1.4 
 
 Tstd = log(1.0+exp((TSL_std+Tlapse*h-Tpause)/Tblend))*Tblend + Tpause
 T = Tstd + ΔT
 
 p = pSL*exp( - 0.11800*h /(1.0 + 0.0020*h) - 0.00198*h^2/(1.0 + 0.0006*h^2) )
 ρ = ɣ*p/((ɣ-1.0)*cp*T)
 a = sqrt(ɣ*p/ρ)
 μ = μSL * sqrt(T/TSL_std)^3 * (TSL_std+Tsuth)/(T+Tsuth)
 
 return T,p,ρ,a,μ

end # atmos

"""
    find_altitude_from_density(ρ::Float64, ΔT::Float64 = 0.0) 
    
Uses a non-linear solver to find the altitude corresponding to a given air density.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ρ::Float64`: air density (kg/m^3)
    - `ΔT::Float64`: temperature difference from standard atmosphere (K)
    
    **Outputs:**
    - `h::Float64`: altitude (km)
"""
function find_altitude_from_density(ρ::Float64, ΔT::Float64 = 0.0) 
    res(x) = atmos(x, ΔT)[3] - ρ #Residual for density
    h = find_zero(res, 0.0)
    return h
end

end