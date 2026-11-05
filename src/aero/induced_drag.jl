"""
      induced_drag!(para, wing, htail, trefftz_config)

Computes the induced drag via the Trefftz plane using a discrete vortex method.
Formerly, `cditrp!()` but now broken up into smaller functions.

If the Trefftz plane configuration has a strut surface, the strut carries the lift `CL_strut` given by
the [`strut_lift_share`](@ref) of `CL - CLhtail`. The wing lift is reduced by the same amount, so that the total lift stays `CL`.
The bound circulation of the strut follows ``(1-\\xi^2)^p`` with ``\\xi = 0`` at the fuselage end and ``\\xi = 1`` at the
wing end of the strut (see [`set_strut_circulation_shape!`](@ref)).

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `para::AbstractArray{Float64}`: Array of `aircraft` model aerodynamic parameters.
      - `wing::TASOPT.Wing`: Wing Structure.
      - `htail::TASOPT.Tail`: Htail Structure.
      - `trefftz_config::TrefftzPlaneConfig`: Trefftz plane analysis configuration (panel discretization and analysis parameters).

      **Outputs:**
      - No explicit outputs. Computed induced drag value, span efficiency and strut lift coefficient (`iaCLstrut`) are saved to `para` of `aircraft` model.

!!! compat "Future Changes"
      In an upcoming revision, an `aircraft` `struct` and auxiliary indices will be passed in lieu of pre-sliced `par` arrays.

"""
function induced_drag!(para, ac, trefftz_config::TrefftzPlaneConfig)
    
    wing = ac.wing
    htail = ac.htail

    CL = para[iaCL]

    if (CL == 0.0)
        para[iaCDi] = 0.
        para[iaspaneff] = 1.0
        para[iaCLstrut] = 0.0
        wing.strut.CL = 0.0
        wing.strut.lift_fraction = 0.0
        return
    end

    CLhtail = para[iaCLh] * htail.layout.S / wing.layout.S #normalize by wing area
    # println("CLhtail: $(para[iaCLh]) $(htail.layout.S) $(parg[igS])")
    bref = wing.layout.span
    Sref = wing.layout.S

    Mach = para[iaMach]

    specifies_CL = true

    #---- number of surfaces (wing, horizontal tail and, for a braced wing, strut)
    nsurf = n_surfaces(trefftz_config)

    gammas = zeros(Float64, nsurf)
    gammat = zeros(Float64, nsurf)
    po = zeros(Float64, nsurf)
    CLsurfsp = zeros(Float64, nsurf)

    #---- wing wake parameters
    fLo = wing.fuse_lift_carryover

    gammas[1] = wing.inboard.λ * para[iarcls]
    gammat[1] = wing.outboard.λ * para[iarclt]
    po[1] = 1.0

    #---- strut lift (zero without a strut surface)
    CL_strut = strut_lift_share(wing, trefftz_config) * (CL - CLhtail)
    para[iaCLstrut] = CL_strut
    wing.strut.CL = CL_strut
    wing.strut.lift_fraction = CL_strut / CL

    CLsurfsp[1] = CL - CLhtail - CL_strut

    #---- velocity-change fractions at wing spanwise locations due to fuselage flow
    fduo = para[iafduo]
    fdus = para[iafdus]
    fdut = para[iafdut]

    #---- horizontal tail wake parameters
    gammas[2] = 1.0
    gammat[2] = htail.outboard.λ
    po[2] = 1.0
    CLsurfsp[2] = CLhtail

    #---- strut: its circulation shape is set below from `spanload_exponent`, so po, gammat and gammas stay zero
    if nsurf == 3
        CLsurfsp[3] = CL_strut
    end

    # Ensure wake_system has current geometry (rebuilds if geometry changed)
    ensure_trefftz_current!(ac, po, gammat, gammas, trefftz_config)

    ifrst = [i_first_wing(trefftz_config), i_first_tail(trefftz_config)]
    ilast = [i_last_wing(trefftz_config), i_last_tail(trefftz_config)]
    if nsurf == 3
        push!(ifrst, i_first_strut(trefftz_config))
        push!(ilast, i_last_strut(trefftz_config))
        set_strut_circulation_shape!(TREFFTZ_GEOM, ifrst[3], ilast[3], wing.strut.spanload_exponent)
    end

    calculate_wake_circulation!(gw, TREFFTZ_GEOM.gc, ifrst, ilast, nsurf) #TODO: next cut avoid these module consts. They're ugly.

    if (specifies_CL)
        # Scale circulations to match specified lift coefficient
        scale_circulation!(TREFFTZ_GEOM.gc, gw, TREFFTZ_GEOM.yp, ifrst, ilast,
            CLsurfsp, bref, Sref, nsurf)
    end

    CLtp, CDtp, sefftp = compute_induced_drag!(vnc, ac.wake_system, TREFFTZ_GEOM.gc, gw, bref, Sref)

    para[iaCDi] = CDtp
    para[iaspaneff] = sefftp

    return
end # induced_drag!

"""
    strut_effective_area(strut)

Lifting area of both struts [m²]: the sloped planform area `strut.S` times `cos(strut.dihedral)^2`. One factor
`cos(dihedral)` gives the area seen from above, and the other accounts for the freestream meeting the sloped strut
at an incidence reduced by `cos(dihedral)`.
"""
strut_effective_area(strut) = strut.S * cos(strut.dihedral)^2

"""
    strut_lift_share(wing)
    strut_lift_share(wing, trefftz_config)

Share ``k`` of the lift carried by the wing and the strut, `CL - CLhtail`, that is carried by both struts: ``k = r / (1 + r)``
with ``r`` = `mean_cl_ratio_to_wing` times the square of the local velocity ratio of the strut times the
[`strut_effective_area`](@ref) divided by the wing area. ``k`` is always below 1. It is zero if the wing has no area or the
strut has no area yet, and, with the `trefftz_config` argument, if the configuration has no strut surface.
"""
function strut_lift_share(wing)
    strut = wing.strut
    wing.layout.S > 0.0 || return 0.0
    r = strut.mean_cl_ratio_to_wing * strut.local_velocity_ratio^2 * strut_effective_area(strut) / wing.layout.S
    return r / (1.0 + r)
end

function strut_lift_share(wing, trefftz_config::TrefftzPlaneConfig)
    n_surfaces(trefftz_config) == 3 || return 0.0
    return strut_lift_share(wing)
end

"""
    strut_lift_centroid_fraction(p)

Centroid ``\\xi_{cp}`` of the strut loading ``(1-\\xi^2)^p`` on ``0 \\le \\xi \\le 1`` with ``\\xi = 0`` at the fuselage end of the strut and
``\\xi = 1`` at the wing end, for the loading exponent `p` (`wing.strut.spanload_exponent`):
``\\xi_{cp} = \\int_0^1 \\xi (1-\\xi^2)^p d\\xi / \\int_0^1 (1-\\xi^2)^p d\\xi
= \\Gamma(p+3/2) / (\\sqrt{\\pi}\\,(p+1)\\,\\Gamma(p+1))``. It is 1/2 for ``p = 0``. This is the fraction of the strut lift that
is carried at the wing end of a strut with pinned ends, and the position of the lift centre of the strut between its ends.
"""
strut_lift_centroid_fraction(p) = exp(log_gamma(p + 1.5) - log_gamma(p + 1.0)) / (sqrt(pi) * (p + 1.0))

"""
    strut_lift_arm(wing)

Streamwise offset of the lift centre of the strut from the wing box [m], with the sign convention of the wing lift arm
`co * CMw1` in [`balance_aircraft!`](@ref): the lift acts at `xwbox - arm`, so a lift centre aft of the wing box has a
negative arm. The lift acts at the spanwise centroid of the strut loading ``(1-\\xi^2)^p`` with ``p`` = `wing.strut.spanload_exponent` (see [`set_strut_circulation_shape!`](@ref)), between the
strut root `wing.strut.y_root` (``\\xi = 0``) and the attachment at the wing break (``\\xi = 1``). The centroid is
``\\xi_{cp}`` of [`strut_lift_centroid_fraction`](@ref).

The strut is assumed to lie under the wing box axis, which is swept back by `wing.layout.sweep`, so a strut lift centre
at the spanwise position `y_cp` is `y_cp tan(sweep)` aft of the wing box root. The chordwise offset of the strut from
the wing box and the pitching moment of the strut section (symmetric section) are neglected.
"""
function strut_lift_arm(wing)
    strut = wing.strut
    ξ_cp = strut_lift_centroid_fraction(strut.spanload_exponent)
    y_attach = 0.5 * wing.layout.break_span
    y_cp = strut.y_root + (y_attach - strut.y_root) * ξ_cp
    return -tand(wing.layout.sweep) * y_cp
end

"""
    strut_blended_lift_arm(co_CMw1, wing, trefftz_config)

Offset of the lift centre of the wing and the strut from the wing box [m], with the sign convention of
[`strut_lift_arm`](@ref): the wing lift arm `co_CMw1` (the wing root chord times `CMw1`) and the
[`strut_lift_arm`](@ref) weighted by the [`strut_lift_share`](@ref). Returns `co_CMw1` unchanged if the strut carries
no lift.
"""
function strut_blended_lift_arm(co_CMw1, wing, trefftz_config::TrefftzPlaneConfig)
    k = strut_lift_share(wing, trefftz_config)
    k == 0.0 && return co_CMw1
    return (1.0 - k) * co_CMw1 + k * strut_lift_arm(wing)
end

"""
    log_gamma(x)

Natural logarithm of the gamma function for `x > 0`. The argument is shifted up to at least 15 with
``\\Gamma(x+1) = x\\Gamma(x)``, where the Stirling series with four terms is accurate to about 1e-14.
"""
function log_gamma(x)
    log_product = 0.0
    while x < 15.0
        log_product += log(x)
        x += 1.0
    end
    return (x - 0.5) * log(x) - x + 0.5 * log(2.0 * pi) +
           1.0 / (12.0 * x) - 1.0 / (360.0 * x^3) + 1.0 / (1260.0 * x^5) - 1.0 / (1680.0 * x^7) - log_product
end

"""
    ensure_trefftz_current!(ac, po, gammat, gammas, trefftz_config)

Checks if the Trefftz plane geometry is up-to-date with the current wing and tail geometry
(spans, root spans, break locations, taper ratios and the heights `layout.z` of the wing and tail points
in the wake) and, if the configuration has a strut surface, the strut geometry and `spanload_exponent`.
If not, rebuilds the geometry. The strut circulation itself is set by [`set_strut_circulation_shape!`](@ref) at every
call of [`induced_drag!`](@ref), because [`scale_circulation!`](@ref) scales it in place.
"""
function ensure_trefftz_current!(ac, po, gammat, gammas, trefftz_config)
    wing = ac.wing
    htail = ac.htail

    geometry = (wing.layout.span, wing.layout.root_span, wing.layout.ηs,
        htail.layout.span, htail.layout.root_span, htail.layout.ηs, wing.layout.z, htail.layout.z,
        gammat..., gammas...)
    if n_surfaces(trefftz_config) == 3
        geometry = (geometry..., wing.strut.y_root, wing.strut.z, wing.strut.spanload_exponent)
    end
    current_hash = hash(geometry)

    # Only rebuild if geometry changed
    if current_hash != TREFFTZ_GEOMETRY_HASH[]
        # println("Rebuilding Trefftz plane geometry...", TREFFTZ_GEOMETRY_HASH[], " --> ", current_hash)
        ac.wake_system = _build_trefftz_geometry!(wing, htail, po, gammat, gammas, trefftz_config)
        TREFFTZ_GEOMETRY_HASH[] = current_hash
    end
end  # function ensure_trefftz_current!

"""
    bunch_transform(t, bunch)

Apply bunching transformation to input `t ∈ \\[0,1\\]`.

Transforms `t` to cluster points near the center (t=0.5) when `bunch > 0`.
The transformation is: `t_bunched = t + bunch * t * (1 - t)`

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: Original normalized angle ∈ \\[0,1\\].
    - `bunch::Float64`: Clustering factor ∈ \\[0,1\\]. Higher values cluster more points near center.

    **Outputs:**
    - `t_bunched::Float64`: Transformed parameter value.

See also [`inv_bunch_transform`](@ref).
"""
@inline function bunch_transform(t, bunch)
    return t + bunch * t * (1.0 - t)
end

"""
    inv_bunch_transform(t_bunched, bunch)

Inverse of the bunching transformation. Recovers original `t` from bunched value.

Solves `t_bunched = t + bunch * t * (1 - t)` for `t` using the quadratic formula.
See also [`bunch_transform`](@ref).
"""
@inline function inv_bunch_transform(t_bunched, bunch)
    return (1.0 + bunch - sqrt((1.0 + bunch)^2 - 4.0 * bunch * t_bunched)) * 0.5 / bunch
end

"""
    calculate_wake_circulation!(gw, gc, ifrst, ilast, nsurf)

Compute wake vortex circulation from bound circulation `γ_wake = -dΓ/dy`.

At each panel edge, the change in bound circulation is shed into the wake:
- Root: `gw = 0` (no upstream circulation at the center of the fuselage)
- Interior: `gw[i] = gc[i-1] - gc[i]` 
- Tip: `gw[i] = gc[i-1]` 
"""
@inline function calculate_wake_circulation!(
    gw::AbstractVector{Float64},
    gc::AbstractVector{Float64},
    i_first::AbstractVector{Int},
    i_last::AbstractVector{Int},
    nsurf::Int
)
    @inbounds for isurf = 1:nsurf
        i = i_first[isurf]
        gw[i] = 0.0 # centerline wake circulations cancel from left and right contributions

        @inbounds for i = (i_first[isurf]+1):(i_last[isurf]-1)
            gw[i] = gc[i-1] - gc[i]
        end

        i = i_last[isurf]
        gw[i] = gc[i-1]
    end

    return nothing
end

"""
    scale_circulation!(gc, gw, yp, ifrst, ilast, CLsurfsp, bref, Sref, nsurf)

Scale circulation distribution to match specified specified lift coefficient for each surface.
"""
@inline function scale_circulation!(
    gc::AbstractVector{Float64},
    gw::AbstractVector{Float64},
    yp::AbstractVector{Float64},
    ifrst::AbstractVector{Int},
    ilast::AbstractVector{Int},
    CLsurfsp::AbstractVector{Float64},
    bref::Float64,
    Sref::Float64,
    nsurf::Int
)
    @inbounds for isurf = 1:nsurf
        # Integrate circulation distribution to find current CL
        cl_test = 0.0
        @inbounds for i = ifrst[isurf]:(ilast[isurf]-1)
            dy = yp[i+1] - yp[i]
            cl_test += gc[i] * dy
        end

        # Convert integrated circulation to lift coefficient
        cl_test = cl_test * 2.0 * bref / (0.5 * Sref)

        # Compute scaling factor. A surface with no circulation (the strut) stays at zero.
        gfac = cl_test == 0.0 ? 0.0 : CLsurfsp[isurf] / cl_test

        # Scale both bound and wake circulations
        @inbounds for i = ifrst[isurf]:ilast[isurf]
            gw[i] *= gfac
        end
        @inbounds for i = ifrst[isurf]:ilast[isurf]-1
            gc[i] *= gfac
        end
    end

    return nothing
end

@enum WAKE_CONTRACTION_TYPE begin
    FUSEWAKE
    WINGWAKE
end

function get_wake_contraction(wake_type::WAKE_CONTRACTION_TYPE, y, yo, yop)
      if wake_type === FUSEWAKE
            # Power law to contract the stream tubes in the fuselage region.
            # Use tolerance for boundary case (y ≈ yo)
            if y > yo * (1.0 + 1e-10)
                  error("y=$y must be within fuselage region (y ≤ yo=$yo)")
            end
            yexp = (yo / yop)^2
            return yop * (y / yo)^yexp
      else
            # Mass conservation based contraction outside fuselage region.
            # Use tolerance for boundary case (y ≈ yo)
            if y < yo * (1.0 - 1e-10)
                  error("y=$y must be outside fuselage region (y ≥ yo=$yo)")
            end
            return sqrt(y^2 - yo^2 + yop^2)
      end
end

"""
    generate_panel_points!(
        i, geom,
        k_start, k_end,
        t_start, t_end,
        e_start, e_end,
        y_start, y_end, z_start, z_end,
        g_start, g_end,
        yo, yop, bunch, ktip
    ) -> i_end

Generate points for a single panel section (image, inboard, or outboard).
This eliminates code duplication between the three panel loops.

Type-stable and efficient for repeated use in different panel sections.
"""
@inline function generate_panel_points!(
    i::Int,
    geom::TrefftzGeometry,
    k_start::Int, k_end::Int,
    t_start::Float64, t_end::Float64,
    e_start::Float64, e_end::Float64,
    y_start::Float64, y_end::Float64,
    z_start::Float64, z_end::Float64,
    g_start::Float64, g_end::Float64,
    yo::Float64, yop::Float64,
    bunch::Float64, ktip::Float64,
    wake_contract::WAKE_CONTRACTION_TYPE)

    @inbounds for k = k_start+1:k_end
        i = i + 1

        # Interpolation fraction
        fk = (k - k_start) / (k_end - k_start)

        # Theta values for field point and control point
        geom.t[i] = t_start * (1.0 - fk) + t_end * fk
        tc = 0.5 * (geom.t[i-1] + geom.t[i])

        # Eta values (spanwise position)
        e  = cos(0.5π * bunch_transform(geom.t[i], bunch))
        ec = cos(0.5π * bunch_transform(tc, bunch))

        # Interpolation fractions in η space
        if e_end - e_start == 0.0
            fi = 1.0
            fc = 0.5
        else
            fi = (e - e_start) / (e_end - e_start)
            fc = (ec - e_start) / (e_end - e_start)
        end

        # Physical coordinates - field points
        geom.y[i] = y_start * (1.0 - fi) + y_end * fi
        geom.z[i] = z_start * (1.0 - fi) + z_end * fi

        # Physical coordinates - control points
        geom.yc[i-1] = y_start * (1.0 - fc) + y_end * fc
        geom.zc[i-1] = z_start * (1.0 - fc) + z_end * fc

        # Circulation with tip rolloff
        geom.gc[i-1] = (g_start * (1.0 - fc) + g_end * fc) * sqrt(1.0 - ec^ktip)

        #Trefftz plane locations (y′ & yc′) of stream lines after wake contraction
        geom.yp[i] = get_wake_contraction(wake_contract, geom.y[i], yo, yop)
        geom.ycp[i-1] = get_wake_contraction(wake_contract, geom.yc[i-1], yo, yop)

        geom.zp[i] = geom.z[i]
        geom.zcp[i-1] = geom.zc[i-1]
    end

    return i
end

"""
    generate_trefftz_points!(
        i_start::Int,
        geom::TrefftzGeometry,
        surface,
        po, gammat, gammas,
        panels::SurfaceDiscretization,
        trefftz_config::TrefftzPlaneConfig,
        root_contraction
    ) -> i_end::Int

Generate field points and control points for a single lifting surface in the Trefftz plane.
Returns the final index after point generation.

This is an incremental refactoring step - extracted from _trefftz_analysis for clarity.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `i_start::Int`: Starting index in the arrays.
    - `geom::TrefftzGeometry`: Geometry arrays struct for point storage.
    - `surface`: Wing or tail structure (TASOPT.Wing or TASOPT.Tail) containing geometry.
    - `po::Float64`: Root circulation scaling factor.
    - `gammat::Float64`: Outer section lift distribution taper ratio.
    - `gammas::Float64`: Inner section lift distribution taper ratio.
    - `panels::SurfaceDiscretization`: Panel discretization for this surface.
    - `trefftz_config::TrefftzPlaneConfig`: Configuration with bunching and tip parameters.
    - `root_contraction::Float64`: Wake contraction factor at fuselage root.

    **Outputs:**
    - `i_end::Int`: Final index after point generation.
"""
@inline function generate_trefftz_points!(
    i_start::Int,
    geom::TrefftzGeometry,
    surface,
    po::Float64, gammat::Float64, gammas::Float64,
    panels::SurfaceDiscretization,
    trefftz_config::TrefftzPlaneConfig,
    root_contraction::Float64
)
    # Extract geometry from surface struct
    b = surface.layout.span
    bs = surface.layout.break_span 
    bo = surface.layout.root_span
    bop = bo * root_contraction
    zcent = surface.layout.z

    # Extract configuration
    bunch = trefftz_config.bunch
    ktip = trefftz_config.k_tip

    # Extract panel counts
    n_img = panels.n_image_panels
    n_inn = panels.n_inner_panels
    n_out = panels.n_outer_panels

    # Key span positions (η = spanwise location normalized by b)
    e0 = 0.0
    eo = bo/b
    es = bs/b
    e1 = 1.0
    eop = bop/b

    # Transform to θ space
    t0 = 1.0
    to = acos(eo)/(0.5*π)
    ts = acos(es)/(0.5*π)
    t1 = 0.0

    # Apply inverse bunching
    if bunch > 0.0
        to = inv_bunch_transform(to, bunch)
        ts = inv_bunch_transform(ts, bunch)
    end

    # Physical coordinates at key stations
    y0 = 0.0
    yo = 0.5*bo
    ys = 0.5*bs
    y1 = 0.5*b
    yop = 0.5*bop

    z0 = zcent
    zo = zcent
    zs = zcent
    z1 = zcent

    # Circulation at key stations
    g0 = po
    go = po
    gs = po*gammas
    g1 = po*gammat

    # Panel indices
    k0 = 1
    ko = 1 + n_img
    ks = 1 + n_img + n_inn
    k1 = 1 + n_img + n_inn + n_out

    # Start with center point
    i = i_start
    geom.t[i] = t0
    geom.y[i] = y0
    geom.z[i] = z0
    geom.yp[i] = y0
    geom.zp[i] = z0

    # Image panels (fuselage region)
    i = generate_panel_points!(
        i, geom,
        k0, ko, t0, to, e0, eo,
        y0, yo, z0, zo, g0, go,
        yo, yop, bunch, ktip, FUSEWAKE)

    # Inboard panels
    i = generate_panel_points!(
        i, geom,
        ko, ks, to, ts, eo, es,
        yo, ys, zo, zs, go, gs,
        yo, yop, bunch, ktip, WINGWAKE)

    # Outboard panels
    i = generate_panel_points!(
        i, geom,
        ks, k1, ts, t1, es, e1,
        ys, y1, zs, z1, gs, g1,
        yo, yop, bunch, ktip, WINGWAKE)

    return i
end

"""
    generate_trefftz_points_inclined!(i_start::Int, geom::TrefftzGeometry,
    y_start, z_start, y_end, z_end, n_panels::Int) -> i_end::Int

Generate field points and control points for a zero-lift line in the Trefftz plane, evenly spaced along the
straight segment from (`y_start`, `z_start`) to (`y_end`, `z_end`). Returns the final index after point generation.

The `n_panels + 1` field points start at index `i_start`. Each control point is at the midpoint of two consecutive
field points. The circulation `gc` is zero (see [`set_strut_circulation_shape!`](@ref) for the strut) and the wake is
not contracted (`yp = y`, `zp = z`).
"""
function generate_trefftz_points_inclined!(
    i_start::Int,
    geom::TrefftzGeometry,
    y_start, z_start,
    y_end, z_end,
    n_panels::Int
)
    i = i_start
    geom.y[i] = y_start
    geom.yp[i] = y_start
    geom.z[i] = z_start
    geom.zp[i] = z_start

    @inbounds for k = 1:n_panels
        i = i + 1
        geom.y[i] = y_start + (y_end - y_start) * k / n_panels
        geom.yp[i] = geom.y[i]
        geom.z[i] = z_start + (z_end - z_start) * k / n_panels
        geom.zp[i] = geom.z[i]

        geom.yc[i-1] = 0.5 * (geom.y[i-1] + geom.y[i])
        geom.ycp[i-1] = geom.yc[i-1]
        geom.zc[i-1] = 0.5 * (geom.z[i-1] + geom.z[i])
        geom.zcp[i-1] = geom.zc[i-1]
        geom.gc[i-1] = 0.0
    end

    return i
end

"""
    generate_strut_points!(i_start::Int, geom::TrefftzGeometry, wing, panels::SurfaceDiscretization) -> i_end::Int

Generate field points and control points for the strut in the Trefftz plane. Returns the final index
after point generation.

The points are evenly spaced along the straight line (see [`generate_trefftz_points_inclined!`](@ref)) from the
strut root at the fuselage, at lateral position `wing.strut.y_root` and a height `wing.strut.z` below the wing box
(`wing.layout.z - wing.strut.z`), to the strut attachment at the wing break (`ηs b/2`), at the height of the
wing box `wing.layout.z`. The circulation `gc` is zero here and set by [`set_strut_circulation_shape!`](@ref) for
each [`induced_drag!`](@ref) call. The wake is not contracted (`yp = y`).
"""
function generate_strut_points!(
    i_start::Int,
    geom::TrefftzGeometry,
    wing,
    panels::SurfaceDiscretization
)
    return generate_trefftz_points_inclined!(
        i_start, geom,
        wing.strut.y_root, wing.layout.z - wing.strut.z,
        0.5 * wing.layout.break_span, wing.layout.z,
        n_total(panels))
end

"""
    set_strut_circulation_shape!(geom::TrefftzGeometry, i_first::Int, i_last::Int, spanload_exponent)

Set the bound circulation `gc` of the strut control points `i_first:i_last-1` to the shape ``(1-\\xi^2)^p`` with
``p`` = `spanload_exponent` and ``\\xi = (k-1/2)/n`` for the control point `k` of the `n` strut panels. Each strut is
treated as half a wing mirrored about the fuselage centreline: ``\\xi = 0`` at the fuselage end, where the loading is
largest like at a wing root, and ``\\xi = 1`` at the wing attachment. The last panel (at the wing attachment) has zero
circulation, because the strut is joined to the wing and sheds no tip vortex there. The first panel keeps its shape
value; the wake circulation at the fuselage centreline cancels between the left and right strut, see
[`calculate_wake_circulation!`](@ref). [`scale_circulation!`](@ref) then scales the shape to the strut lift.
"""
function set_strut_circulation_shape!(geom::TrefftzGeometry, i_first::Int, i_last::Int, spanload_exponent)
    n = i_last - i_first
    @inbounds for k = 1:n
        ξ = (k - 0.5) / n
        geom.gc[i_first+k-1] = (1.0 - ξ^2)^spanload_exponent
    end
    geom.gc[i_last-1] = 0.0
    return nothing
end

"""
    _build_trefftz_geometry!(wing, htail, po, gammat, gammas,
    trefftz_config::TrefftzPlaneConfig)

Construct the Trefftz plane geometry for the given wing and horizontal tail and, if `trefftz_config`
has strut panels, the strut of the wing (see [`generate_strut_points!`](@ref)).
"""
function _build_trefftz_geometry!(wing, htail, po, gammat, gammas,
    trefftz_config::TrefftzPlaneConfig)

    i::Int64 = 0

    nsurf = n_surfaces(trefftz_config)
    @inbounds for isurf = 1:nsurf
        i = i + 1
        if isurf == 3
            # Generate points for the strut
            i = generate_strut_points!(i, TREFFTZ_GEOM, wing, trefftz_config.strut_panels)
        else
            # Select the right surface and panel config
            surface = isurf == 1 ? wing : htail
            panels = isurf == 1 ? trefftz_config.wing_panels : trefftz_config.tail_panels
            root_contraction = isurf == 1 ? trefftz_config.wing_root_contraction : trefftz_config.tail_root_contraction

            # Handle T-tail special case (root_span = 0)
            if isurf == 2 && surface.layout.root_span == 0.0
                panels = SurfaceDiscretization(panels.n_outer_panels, panels.n_inner_panels, 0)
            end

            # Generate points for this surface
            i = generate_trefftz_points!(
                i, TREFFTZ_GEOM,
                surface,
                po[isurf], gammat[isurf], gammas[isurf],
                panels,
                trefftz_config,
                root_contraction
            )
        end

        #---- dummy control point between surfaces
        TREFFTZ_GEOM.yc[i] = 0.0
        TREFFTZ_GEOM.zc[i] = 0.0
        TREFFTZ_GEOM.gc[i] = 0.0
    end # nsurf loop
    # ilast = i_last_strut or i_last_tail (trefftz_config) # could use this instead of accumulating i...
    WS = WakeSystem(
        TREFFTZ_GEOM.yp[1:i],
        TREFFTZ_GEOM.zp[1:i],
        TREFFTZ_GEOM.ycp[1:i-1],
        TREFFTZ_GEOM.zcp[1:i-1])
    return WS
end


"""
    compute_induced_drag!(vnc, WS::WakeSystem{Npts, Nelems}, gc, gw, bref, Sref) where {Npts, Nelems}

Compute the induced drag and span efficiency from the wake system.
"""
function compute_induced_drag!(vnc, WS::WakeSystem{Npts, Nelems}, gc, gw, bref, Sref) where {Npts, Nelems}

    @views mul!(vnc[1:Nelems], WS.influence_matrix[1:Nelems, 1:Npts], gw[1:Npts])

    # Normalize in-place
    factor = bref / (2.0π)
    @inbounds for i in 1:Nelems
        vnc[i] *= factor
    end

    # Compute CL and CD
    # Need to extract element geometry on-the-fly
    CL = 0.0
    CD = 0.0

    @inbounds for i in 1:Nelems
        # Extract element geometry from WakeSystem
        elem = WS.elements[i]
        dy = elem.Δy
        ds = elem.length

        CD -= gc[i] * vnc[i] * ds
        CL += 2.0 * gc[i] * dy
    end

    # Normalize
    CD = CD * bref / (0.5 * Sref)
    CL = CL * bref / (0.5 * Sref)

    AR = bref^2 / Sref
    CD_elliptical = CL^2 / (π * AR) #CD for elliptical loading here. 
    # For general CD = CL²/(π*AR*spaneff) ∴spaneff calculated as:
    span_eff = CD_elliptical / CD

    return CL, CD, span_eff
end