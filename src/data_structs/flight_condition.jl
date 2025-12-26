"""
    FlightCondition

Encapsulates the flight state at a given mission point, replacing scattered
`para`/`pare` array indices with a structured type.

- `atmos::AtmosphericState`: Atmospheric conditions (T, p, ρ, a, μ) at altitude
- `Mach::Float64`: Mach number (was `para[iaMach, ip]`)
- `TAS::Float64`: True airspeed in m/s (was `pare[ieu0, ip]`)
- `Re_unit::Float64`: Unit Reynolds number in 1/m (was `para[iaReunit, ip]`)
- `climb_angle::Float64`: Climb angle γ in radians (was `para[iagamV, ip]`)
```
"""
struct FlightCondition
    atmos::atmosphere.AtmosphericState
    Mach::Float64
    TAS::Float64          # True airspeed [m/s]
    Re_unit::Float64      # Unit Reynolds number [1/m]
    climb_angle::Float64  # γ [rad]
end

"""
    FlightCondition(alt::Real, Mach::Real; ΔT::Real=0.0)

Construct a FlightCondition from altitude and Mach number.
Automatically computes TAS and unit Reynolds number from atmospheric properties.

# Arguments
- `alt`: Altitude in meters
- `Mach`: Mach number
- `ΔT`: Temperature offset from ISA in Kelvin (default: 0.0)
"""
function FlightCondition(alt, Mach; ΔT=0.0, climb_angle=0.0)
    atmos_state = atmosphere.AtmosphericState(alt, ΔT)
    TAS = Mach * atmos_state.a
    Re_unit = atmos_state.ρ * TAS / atmos_state.μ
    return FlightCondition(atmos_state, Mach, TAS, Re_unit, climb_angle)
end

function Base.getproperty(obj::FlightCondition, sym::Symbol)
    if hasfield(FlightCondition, sym)
        return getfield(obj, sym)
    elseif hasfield(typeof(obj.atmos), sym)
        return getproperty(getfield(obj, :atmos), sym)
    else
        return getfield(obj, sym)
    end
end



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Temporary adapters for old array-based code
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

"""
    FlightCondition(para::AbstractArray, pare::AbstractArray, ip::Integer; ΔT::Real=0.0)

Extract a FlightCondition from legacy `para` and `pare` arrays at mission point `ip`.

This adapter enables gradual migration from array-based to struct-based code.
The arrays can be 2D views (indices × flight_points) or 3D (indices × flight_points × missions).

# Arguments
- `ΔT`: Temperature offset from ISA in Kelvin (default: 0.0).
  For full consistency, pass `parm[imDeltaTatm]` from the mission parameters.
"""
function FlightCondition(para::AbstractArray, pare::AbstractArray, ip::Integer; ΔT::Real=0.0)
    alt = _get_value(para, iaalt, ip)
    Mach = _get_value(para, iaMach, ip)
    Re_unit = _get_value(para, iaReunit, ip)
    γ = _get_value(para, iagamV, ip)
    TAS = _get_value(pare, ieu0, ip)

    atmos_state = atmosphere.AtmosphericState(alt, ΔT)
    return FlightCondition(atmos_state, Mach, TAS, Re_unit, γ)
end

# Helper to handle both 2D and 3D arrays
function _get_value(arr::AbstractArray, idx::Integer, ip::Integer)
    if ndims(arr) == 2
        return arr[idx, ip]
    else  # 3D: assume mission index 1
        return arr[idx, ip, 1]
    end
end

"""
    update_arrays!(para, pare, fc::FlightCondition, ip::Integer)

Write a FlightCondition back to legacy `para` and `pare` arrays at mission point `ip`.

This adapter enables gradual migration from array-based to struct-based code.
"""
function update_arrays!(para::AbstractArray, pare::AbstractArray, fc::FlightCondition, ip::Integer)
    _set_value!(para, iaalt, ip, fc.atmos.alt)
    _set_value!(para, iaMach, ip, fc.Mach)
    _set_value!(para, iaReunit, ip, fc.Re_unit)
    _set_value!(para, iagamV, ip, fc.climb_angle)

    _set_value!(pare, ieT0, ip, fc.atmos.T)
    _set_value!(pare, iep0, ip, fc.atmos.p)
    _set_value!(pare, ierho0, ip, fc.atmos.ρ)
    _set_value!(pare, iea0, ip, fc.atmos.a)
    _set_value!(pare, iemu0, ip, fc.atmos.μ)
    _set_value!(pare, ieu0, ip, fc.TAS)
    _set_value!(pare, ieM0, ip, fc.Mach)  # Also stored in pare for engine calculations

    return nothing
end

# Helper to handle both 2D and 3D arrays
function _set_value!(arr::AbstractArray, idx::Integer, ip::Integer, val)
    if ndims(arr) == 2
        arr[idx, ip] = val
    else  # 3D: assume mission index 1
        arr[idx, ip, 1] = val
    end
    return nothing
end

# ----------------------------------------


"""
    altitude(fc::FlightCondition)

Return altitude in meters.
"""
altitude(fc::FlightCondition) = fc.atmos.alt

"""
    altitude_km(fc::FlightCondition)

Return altitude in kilometers.
"""
altitude_km(fc::FlightCondition) = fc.atmos.alt / 1000.0

"""
    dynamic_pressure(fc::FlightCondition)

Compute dynamic pressure q = 0.5 * ρ * V² in Pa.
"""
dynamic_pressure(fc::FlightCondition) = 0.5 * fc.atmos.ρ * fc.TAS^2

"""
    total_temperature(fc::FlightCondition)

Compute total (stagnation) temperature Tt = T * (1 + (γ-1)/2 * M²) in K.
Assumes γ = 1.4 for air.
"""
function total_temperature(fc::FlightCondition)
    γ = 1.4
    return fc.atmos.T * (1.0 + (γ - 1.0) / 2.0 * fc.Mach^2)
end

"""
    total_pressure(fc::FlightCondition)

Compute total (stagnation) pressure pt = p * (1 + (γ-1)/2 * M²)^(γ/(γ-1)) in Pa.
Assumes γ = 1.4 for air.
"""
function total_pressure(fc::FlightCondition)
    γ = 1.4
    return fc.atmos.p * (1.0 + (γ - 1.0) / 2.0 * fc.Mach^2)^(γ / (γ - 1.0))
end

function Base.show(io::IO, fc::FlightCondition)
    print(io, "FlightCondition(alt=$(round(fc.atmos.alt/1000, digits=2)) km, ",
              "M=$(round(fc.Mach, digits=3)), ",
              "TAS=$(round(fc.TAS, digits=1)) m/s, ",
              "γ=$(round(rad2deg(fc.climb_angle), digits=2))°)")
end

function Base.show(io::IO, ::MIME"text/plain", fc::FlightCondition)
    println(io, "FlightCondition:")
    println(io, "  Altitude:           $(round(fc.atmos.alt/1000, sigdigits=3)) km ($(round(fc.atmos.alt, digits=1)) m)")
    println(io, "  Mach number:        $(round(fc.Mach, sigdigits=3))")
    println(io, "  True airspeed:      $(round(fc.TAS, sigdigits=2)) m/s")
    println(io, "  Unit Reynolds:      $(round(fc.Re_unit/1e6, sigdigits=3))×10⁶ 1/m")
    println(io, "  Flight path angle:  $(round(rad2deg(fc.climb_angle), digits=2))°")
    println(io, "  Atmospheric state: ", fc.atmos)
    # println(io, "    Temperature:      $(round(fc.atmos.T, sigdigits=3)) K")
    # println(io, "    Pressure:         $(round(fc.atmos.p/1000, sigdigits=3)) kPa")
    # println(io, "    Density:          $(round(fc.atmos.ρ, sigdigits=4)) kg/m³")
    # println(io, "    Speed of sound:   $(round(fc.atmos.a, sigdigits=3)) m/s")
    # print(io,   "    Viscosity:        $(round(fc.atmos.μ*1e6, sigdigits=4))×10⁻⁶ kg/(m·s)")
end

# Approximate equality for testing
function Base.isapprox(fc1::FlightCondition, fc2::FlightCondition; rtol=1e-6, atol=0.0)
    return isapprox(fc1.atmos.alt, fc2.atmos.alt; rtol, atol) &&
           isapprox(fc1.atmos.T, fc2.atmos.T; rtol, atol) &&
           isapprox(fc1.atmos.p, fc2.atmos.p; rtol, atol) &&
           isapprox(fc1.atmos.ρ, fc2.atmos.ρ; rtol, atol) &&
           isapprox(fc1.atmos.a, fc2.atmos.a; rtol, atol) &&
           isapprox(fc1.atmos.μ, fc2.atmos.μ; rtol, atol) &&
           isapprox(fc1.Mach, fc2.Mach; rtol, atol) &&
           isapprox(fc1.TAS, fc2.TAS; rtol, atol) &&
           isapprox(fc1.Re_unit, fc2.Re_unit; rtol, atol) &&
           isapprox(fc1.climb_angle, fc2.climb_angle; rtol, atol)
end


"""
    initialize_flight_conditions(n_points::Integer, n_missions::Integer)

Create an uninitialized flight conditions matrix of size [n_points × n_missions].
All entries are set to a default sea-level, zero-Mach condition.
Flight conditions are populated during mission setup via `set_flight_condition!`.
"""
function initialize_flight_conditions(n_points::Integer, n_missions::Integer)
    # Default FlightCondition at sea level with zero velocity
    default_fc = FlightCondition(0.0, 0.0)
    return fill(default_fc, n_points, n_missions)
end

