module atmosphere
using Roots
export atmos, find_altitude_from_density, AtmosState

"""
    AtmosState{T<:Real}

Container for atmospheric properties returned by [`atmos`](@ref).
"""
struct AtmosState{T<:Real}
    T::T
    p::T
    ρ::T
    a::T
    μ::T
end

"""
    atmos(h, ΔT) -> AtmosState
    
Atmospheric functions ` T(h)`, `ρ(h)` etc
valid to `h`=20km, `p(h)` valid to `h`=70km.

Also calculates viscosity using Sutherland's law. Non-standard sea-level temperatures are allowed
with an ISA + ΔT like model.

Units:
- [h]   = m ASL
- [T]   = Kelvin
- [p]   = Pa
- [ρ]   = kg/m^3
- [a]   = m/s
- [μ]   = kg/m-s 
"""
function atmos(h::Float64, ΔT::Float64 = 0.0)
 h_km = h / 1e3

 pSL = 1.0132e5 # Pa
 TSL_std    = 288.2    # K, standard sea level temperature
 Tpause = 216.65   # K
 Tblend = 2.0      # K (tropopause blending T range)
 Tlapse = -6.5     # K/km
 
 μSL  = 1.78e-5  # kg/m-s
 Tsuth = 110.0    # K
 
 cp  = 1004.0   # J/kg-K
 ɣ = 1.4 
 
 Tstd = log(1.0+exp((TSL_std+Tlapse*h_km-Tpause)/Tblend))*Tblend + Tpause
 T = Tstd + ΔT
 
 p = pSL*exp( - 0.11800*h_km /(1.0 + 0.0020*h_km) - 0.00198*h_km^2/(1.0 + 0.0006*h_km^2) )
 ρ = ɣ*p/((ɣ-1.0)*cp*T)
 a = sqrt(ɣ*p/ρ)
 μ = μSL * sqrt(T/TSL_std)^3 * (TSL_std+Tsuth)/(T+Tsuth)
 
 return AtmosState(T, p, ρ, a, μ)

end # atmos

"""
    find_altitude_from_density(ρ::Float64, ΔT::Float64 = 0.0) 
    
Uses a non-linear solver to find the altitude corresponding to a given air density.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ρ::Float64`: air density (kg/m^3)
    - `ΔT::Float64`: temperature difference from standard atmosphere (K)
    
    **Outputs:**
    - `h::Float64`: altitude (m)
"""
function find_altitude_from_density(ρ::Float64, ΔT::Float64 = 0.0) 
    res(x) = atmos(x, ΔT).ρ - ρ #Residual for density
    h = find_zero(res, 0.0)
    return h
end

end
