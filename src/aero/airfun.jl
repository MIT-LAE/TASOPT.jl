export airfun
const nAfun::Int = 3
const Ai = @MMatrix zeros(2,2)
const Ai_cl = @MMatrix zeros(2,2)
const Ai_tau = @MMatrix zeros(2,2)
const Ai_cl_tau = @MMatrix zeros(2,2)
const Aij = @MVector zeros(2)
const Aij_tau = @MVector zeros(2)
const Aijk = @MVector zeros(nAfun) # for the three vars cdf, cdp and cm

"""
    airfun(cl, τ, Mach, air::airfoil)

Looks up airfoil performance data at specified conditions, as precomputed and found in `./src/airfoil_data/`.
Outside the `cl` and `τ` limits of the database, a quadratic penalty is added to `cdp`; Mach is extrapolated without a penalty.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**      
      - `cl::Float64`: Airfoil section lift coefficient.
      - `τ::Float64`: Airfoil section thickness-to-chord ratio.
      - `Mach::Float64`: Mach number.
      - `air::TASOPT.aerodynamics.airfoil`: `airfoil` structure with performance data.

      **Outputs:**
      - `cdf::Float64`: Airfoil section skin friction drag.
      - `cdp::Float64`: Airfoil section pressure drag, including the penalty outside the database (see [`airfoil_database_penalty`](@ref)).
      - `cdw::Float64`: Airfoil section wave drag (unused and assumed 0 here; placeholder left for future implementation).
      - `cm::Float64`: Airfoil section pitching moment.
"""
@views function airfun(cl, τ, Mach, air::airfoil)

    return airfun(cl, τ, Mach, 
                    air.cl, air.τ, air.Ma,
                    air.A,
                    air.A_M,
                    air.A_τ,
                    air.A_cl,
                    air.A_M_τ,
                    air.A_M_cl,
                    air.A_cl_τ,
                    air.A_M_cl_τ)

end

@views function airfun(cl, τ, Mach, 
                Acl,
                Aτ,
                AMa,
                A,
                A_M,
                A_τ,
                A_cl,
                A_M_τ,
                A_M_cl,
                A_cl_τ,
                A_M_cl_τ)
    
    io, im, dMa , tMa  = findsegment(Mach, AMa)
    jo, jm, dcl , tcl  = findsegment(cl, Acl)
    ko, km, dtau, ttau = findsegment(τ, Aτ)
    
    @inbounds for l = 1:nAfun
        @inbounds for jd = 1:2
            @inbounds for kd = 1:2
                j = jo + jd-2
                k = ko + kd-2

                Ai[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A, :, j, k, l),
                    dYdX=view(A_M, :, j, k, l))
                
                Ai_cl[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_cl, :, j, k, l),
                    dYdX=view(A_M_cl, :, j, k, l))
                
                Ai_tau[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_τ, :, j, k, l),
                    dYdX=view(A_M_τ, :, j, k, l))
                
                Ai_cl_tau[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_cl_τ, :, j, k, l),
                    dYdX=view(A_M_cl_τ, :, j, k, l))

            end
        end

        @inbounds for kd = 1:2
            Aij[kd] = interpolate(ix1=1, ix2=2, 
                dx=dcl, t=tcl, Y=view(Ai, :, kd), dYdX=view(Ai_cl, :, kd))
            
            Aij_tau[kd] = interpolate(ix1=1, ix2=2,
                dx=dcl, t = tcl, Y=view(Ai_tau, :, kd), dYdX=view(Ai_cl_tau,:,kd))
        end

        Aijk[l] = interpolate(ix1=1, ix2=2,
            dx=dtau, t=ttau, Y=view(Aij, :), dYdX=view(Aij_tau, :))
    end
    cdf = Aijk[1]
    cdp = Aijk[2]
    cm  = Aijk[3]
    cdw = 0.0

    #Add quadratic penalties for exceeding database cl or h/c limits
    cdp = cdp + cl_database_penalty(cl, Acl) + thickness_database_penalty(τ, Aτ)

    return cdf, cdp, cdw, cm

end


"""
    cl_database_penalty(cl, Acl)
"""
function cl_database_penalty(cl, Acl)
    if cl < Acl[1]
        return 1.0*(cl-Acl[1])^2
    elseif cl > Acl[end]
        return 1.0*(cl-Acl[end])^2
    end
    return 0.0
end

"""
    thickness_database_penalty(τ, Aτ)
"""
function thickness_database_penalty(τ, Aτ)
    if τ < Aτ[1]
        return 25.0*(τ - Aτ[1])^2
    elseif τ > Aτ[end]
        return 25.0*(τ - Aτ[end])^2
    end
    return 0.0
end

"""
    airfoil_database_penalty(cl, τ, Acl, Aτ)
    airfoil_database_penalty(cl, τ, air::airfoil)

Pressure drag penalty that [`airfun`](@ref) adds outside the airfoil database. Zero inside the database.
"""
function airfoil_database_penalty(cl, τ, Acl, Aτ)
    return cl_database_penalty(cl, Acl) + thickness_database_penalty(τ, Aτ)
end

airfoil_database_penalty(cl, τ, air::airfoil) = airfoil_database_penalty(cl, τ, air.cl, air.τ)

"""
    airfoil_database_exceedances(cl, τ, Mach, air::airfoil)

Returns the airfoil database bounds exceeded (`:cl`, `:thickness`, `:Mach`) as a `Vector{Symbol}`, empty inside the database.
"""
function airfoil_database_exceedances(cl, τ, Mach, air::airfoil)
    exceeded = Symbol[]
    if cl < air.cl[1] || cl > air.cl[end]
        push!(exceeded, :cl)
    end
    if τ < air.τ[1] || τ > air.τ[end]
        push!(exceeded, :thickness)
    end
    if Mach > air.Ma[end]
        push!(exceeded, :Mach)
    end
    return exceeded
end

"""
    warn_outside_airfoil_database(cl, τ, Mach, air::airfoil)

Logs one warning per Julia session naming the exceeded bounds if the conditions are outside the airfoil database.
"""
function warn_outside_airfoil_database(cl, τ, Mach, air::airfoil)
    exceeded = airfoil_database_exceedances(cl, τ, Mach, air)
    isempty(exceeded) && return nothing

    @warn "Airfoil section outside the airfoil database: t/c = $τ, c_l = $cl, M⟂ = $Mach.\n" *
          "Exceeded bounds: $(join(exceeded, ", ")). " *
          "Database ranges: c_l = [$(air.cl[1]), $(air.cl[end])], t/c = [$(air.τ[1]), $(air.τ[end])], " *
          "M⟂ up to $(air.Ma[end]).\n" *
          "Outside the range of c_l and t/c, the pressure drag includes a quadratic database-extrapolation penalty, " *
          "which is NOT physics. Mach is extrapolated by the spline without a penalty. " *
          "The share of the wing pressure drag due to this penalty is reported as `penalty_fraction` in the " *
          "`WingDragComponents` returned by `wing_profiledrag_direct`.\n" *
          "This warning is shown only once per Julia session (maxlog = 1)." maxlog = 1
    return nothing
end

"""
    AirfoilDatabaseStation

Conditions at one wing station outside the airfoil database: `panel` (`:inboard` or `:outboard`), `station` index,
spanwise position `η`, section `cl`, thickness-to-chord `τ` and perpendicular Mach number `Mach`.
"""
struct AirfoilDatabaseStation
    panel::Symbol
    station::Int
    η::Float64
    cl::Float64
    τ::Float64
    Mach::Float64
end

"""
    record_airfoil_database_station!(outside, panel, station, η, cl, τ, Mach, air::airfoil)

Appends the station to `outside` if it is outside the airfoil database.
"""
function record_airfoil_database_station!(outside, panel, station, η, cl, τ, Mach, air::airfoil)
    isempty(airfoil_database_exceedances(cl, τ, Mach, air)) && return nothing
    push!(outside, AirfoilDatabaseStation(panel, station, η, cl, τ, Mach))
    return nothing
end

"""
    warn_wing_outside_airfoil_database(outside, air::airfoil, penalty_fraction, CDpwing; imission = nothing, ip = nothing)

Reports the stations `outside` the airfoil database (see [`AirfoilDatabaseStation`](@ref)) of one `wing_profiledrag_direct` call.
Does nothing if `outside` is empty. Otherwise logs one `@debug` per station and one `@warn` for the wing,
shown once per Julia session for each flight point and set of exceeded bounds.
"""
function warn_wing_outside_airfoil_database(outside, air::airfoil, penalty_fraction, CDpwing;
                                            imission = nothing, ip = nothing)
    isempty(outside) && return nothing

    flight_point = isnothing(ip) ? "" : " at $(ip_labels[ip])" * (isnothing(imission) ? "" : " (mission $imission)")
    fmt(x) = @sprintf("%.3g", x)
    range_text(low, high) = fmt(low) == fmt(high) ? fmt(low) : "$(fmt(low)) to $(fmt(high))"

    for s in outside
        @debug "Wing $(s.panel) station $(s.station) (η = $(fmt(s.η)))$flight_point is outside the airfoil database: " *
               "c_l = $(fmt(s.cl)), t/c = $(fmt(s.τ)), M⟂ = $(fmt(s.Mach)), " *
               "exceeded: $(join(airfoil_database_exceedances(s.cl, s.τ, s.Mach, air), ", ")), " *
               "penalty cdp = $(fmt(airfoil_database_penalty(s.cl, s.τ, air)))"
    end

    stations = join(("$panel " * join((s.station for s in outside if s.panel == panel), ", ") *
                     " (η = $(range_text(extrema(s.η for s in outside if s.panel == panel)...)))"
                    for panel in (:inboard, :outboard) if any(s.panel == panel for s in outside)), "; ")

    # M⟂ has an upper limit only
    quantities = ((:cl, "c_l", [s.cl for s in outside], air.cl[1], air.cl[end]),
                  (:thickness, "thickness t/c", [s.τ for s in outside], air.τ[1], air.τ[end]),
                  (:Mach, "M⟂", [s.Mach for s in outside], -Inf, air.Ma[end]))
    exceeded = Symbol[]
    bounds = String[]
    for (key, label, values, lower, upper) in quantities
        low, high = extrema(values)
        amounts = String[]
        low < lower && push!(amounts, "$(fmt(lower - low)) below the lower bound $(fmt(lower))")
        high > upper && push!(amounts, "$(fmt(high - upper)) above the upper bound $(fmt(upper))")
        isempty(amounts) && continue
        push!(exceeded, key)
        push!(bounds, "$label = $(range_text(low, high)), " * join(amounts, ", "))
    end

    @warn "Wing outside the airfoil database$flight_point.\n" *
          "Stations outside: $stations.\n" *
          "Exceeded bounds: $(join(bounds, "; ")).\n" *
          "Penalty, not physics: ΔCDp = $(fmt(penalty_fraction * CDpwing)), " *
          "penalty_fraction = $(fmt(penalty_fraction)) of the wing pressure drag CDp. " *
          "Mach is extrapolated by the spline without a penalty.\n" *
          "Shown once per flight point and exceeded bounds per Julia session; " *
          "JULIA_DEBUG=TASOPT lists every station." maxlog = 1 _id = Symbol(:airfoil_database_, ip, :_, join(exceeded, :_))
    return nothing
end

"""
    findsegment(x::T, xarr::Vector{T})

Uses bisection to find the right interval of the array `xarr` where x lies.
Returns im and io s.t. `xarr[im] < x < xarr[io]`.

Additionally returns the interval `dx = xarr[io] - xarr[im]`
"""
function findsegment(x::Float64, xarr)

    if isnan(x)
        error("Oops, you're searching for a NaN! Go fix your bug!")
    end
    
    im::Int = find_bisection(x, xarr)
    io = im+1

    dx = xarr[io] - xarr[im]
    t = (x - xarr[im])/dx

    return io, im, dx, t
end
"""
    find_bisection(x::T, X::AbstractVector{T}) where T

Simple, fast bisection search in a **sorted** array. This does not sort the array for you.
"""
function find_bisection(x::T, X::AbstractVector{T}) where T
    ilow::Int = 1
    i::Int = length(X)
    while (i-ilow > 1)
        imid = (i+ilow)÷2
        if x<X[imid]
            i = imid
        else
            ilow = imid
        end
    end
    return i-1
end  # function find_bisection

"""
    interpolate(;ix1, ix2, dx, t, Y, dYdX)

Convenience function to evaluate the symmetric form of the cubic spline equation, given
the x interval defined by the indices (ix1, ix2),  the interval (dx), normalized x location (t),
the Y array to be interpolated and the slopes (dYdX) at the interval end points.

Note: This is intentionally defined as a kw only argument list to encourage anyone working 
with this deep part of the code knows what they are doing. See [`eval_spline`](@ref).
"""
function interpolate(;ix1, ix2, dx, t, Y, dYdX)

    yim = Y[ix1] 
    yi = Y[ix2]
    Δ = yi - yim
    fxm = dx * dYdX[ix1] - Δ
    fxo = dx * dYdX[ix2] - Δ

    return eval_spline(t = t, yim = yim, yi = yi, fxm = fxm, fxo = fxo)
end

"""
    eval_spline(;t, yim, yi, fxm, fxo)

Given normalized interpolation location t, y₋ and y, calculates the
spline interpolation of the form    
    (1 - t)*y₁ + t*y₂ + t*(1 - t)*((1 - t)fxm + t*fx),
    where t = (x - x₁)/Δx.
"""
function eval_spline(;t, yim, yi, fxm, fxo)
    tmi = 1.0 - t
    return t*yi + tmi*yim + t*tmi*(tmi*fxm - t*fxo)
end