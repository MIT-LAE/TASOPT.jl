"""
      induced_drag!(para, wing, htail, trefftz_config)

Computes the induced drag via the Trefftz plane using a discrete vortex method.
Formerly, `cditrp!()` but now broken up into smaller functions.

If the Trefftz plane configuration has a strut surface, the strut carries the lift `CL_strut`, the [`strut_lift_share`](@ref) of `CL - CLhtail`,
and the wing lift is reduced by the same amount.

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
        return
    end

    CLhtail = para[iaCLh] * htail.layout.S / wing.layout.S #normalize by wing area
    # println("CLhtail: $(para[iaCLh]) $(htail.layout.S) $(parg[igS])")
    bref = wing.layout.span
    Sref = wing.layout.S

    Mach = para[iaMach]

    specifies_CL = true

    #---- number of surfaces (wing, horizontal tail and, if braced, strut)
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

    #---- strut lift
    CL_strut = strut_lift_share(wing, trefftz_config) * (CL - CLhtail)
    para[iaCLstrut] = CL_strut

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

    #---- strut circulation shape is set below
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
    if nsurf == 3
        # Unlike the wing, the strut sheds its root circulation at the fuselage end
        gw[ifrst[3]] = -TREFFTZ_GEOM.gc[ifrst[3]]
    end

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
    strut_lift_share(wing)
    strut_lift_share(wing, trefftz_config)

Share of the wing and strut lift, `CL - CLhtail`, carried by both struts: `strut.lift_fraction_strut_to_total`.
Zero if the wing or strut has no area or, with `trefftz_config`, if there is no strut surface.
"""
function strut_lift_share(wing)
    (wing.layout.S > 0.0 && wing.strut.S > 0.0) || return 0.0
    return wing.strut.lift_fraction_strut_to_total
end

function strut_lift_share(wing, trefftz_config::TrefftzPlaneConfig)
    n_surfaces(trefftz_config) == 3 || return 0.0
    return strut_lift_share(wing)
end

"""
    strut_lift_centroid_fraction(p)

Centroid fraction ``\\xi_{cp}`` of the strut loading ``(4\\xi(1-\\xi))^p`` on ``0 \\le \\xi \\le 1`` (``\\xi = 0`` at the fuselage end).
The loading is symmetric, so ``\\xi_{cp} = 1/2`` for every exponent `p`; `p` is not used.
"""
strut_lift_centroid_fraction(p) = 0.5

"""
    strut_lift_arm(wing)

Streamwise offset of the strut lift centre from the wing box [m], with the sign convention of the wing lift arm
`co * CMw1` in [`balance_aircraft!`](@ref): the lift acts at `xwbox - arm`.
The lift acts at the centroid of the strut loading ([`strut_lift_centroid_fraction`](@ref)) between the strut root and the wing break.
The box axis is unswept across the centre box and sweeps from the wing root station, so the sweep offset is measured
from the root station (as in `wing_CM`). The chordwise term is the quarter-chord offset of the box axis,
`(spar_box_x_c - 0.25) c_cp cos²Λ`, with the strut chord `c_cp` at the lift centroid (linear taper from `root_chord` to `chord`).
The strut is assumed to lie under the swept wing box axis; the strut section pitching moment is zero (symmetric section).
"""
function strut_lift_arm(wing)
    strut = wing.strut
    ξ_cp = strut_lift_centroid_fraction(strut.spanload_exponent)
    y_root = 0.5 * wing.layout.root_span
    y_attach = 0.5 * wing.layout.break_span
    c_cp = strut.root_chord + (strut.chord - strut.root_chord) * ξ_cp
    sweep_arm = -tand(wing.layout.sweep) * ξ_cp * (y_attach - y_root)
    chord_arm = (wing.layout.spar_box_x_c - 0.25) * c_cp * strut.cos_lambda^2
    return sweep_arm + chord_arm
end

"""
    strut_blended_lift_arm(co_CMw1, wing, trefftz_config)

Offset of the wing and strut lift centre from the wing box [m]: the wing lift arm `co_CMw1` and the
[`strut_lift_arm`](@ref) weighted by the [`strut_lift_share`](@ref). Returns `co_CMw1` if the strut carries no lift.
"""
function strut_blended_lift_arm(co_CMw1, wing, trefftz_config::TrefftzPlaneConfig)
    k = strut_lift_share(wing, trefftz_config)
    k == 0.0 && return co_CMw1
    return (1.0 - k) * co_CMw1 + k * strut_lift_arm(wing)
end

"""
    ensure_trefftz_current!(ac, po, gammat, gammas, trefftz_config)

Checks if the Trefftz plane geometry is up-to-date with the current wing and tail geometry
(and strut geometry, if present). If not, rebuilds the geometry.
"""
function ensure_trefftz_current!(ac, po, gammat, gammas, trefftz_config)
    wing = ac.wing
    htail = ac.htail

    geometry = (wing.layout.span, wing.layout.root_span, wing.layout.ηs,
        htail.layout.span, htail.layout.root_span, htail.layout.ηs, wing.layout.z, htail.layout.z,
        gammat..., gammas...)
    if n_surfaces(trefftz_config) == 3
        geometry = (geometry..., wing.strut.z, wing.strut.spanload_exponent)
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

        # Compute scaling factor
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

Generate field points and control points in the Trefftz plane, evenly spaced along the
straight segment from (`y_start`, `z_start`) to (`y_end`, `z_end`). Returns the final index after point generation.
The circulation `gc` is zero and the wake is not contracted.
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
    generate_strut_points!(i_start::Int, geom::TrefftzGeometry, wing, panels::SurfaceDiscretization, wing_root_contraction::Float64) -> i_end::Int

Generate field points and control points for the strut in the Trefftz plane. Returns the final index
after point generation.

The points run along a straight line (see [`generate_trefftz_points_inclined!`](@ref)) from the strut root at the fuselage
(the wing root station `wing.layout.root_span / 2`, `wing.strut.z` below the wing box) to the attachment at the wing break
(`wing.layout.break_span / 2`, wing box height). The field and control points are contracted with the wing wake law
([`get_wake_contraction`](@ref)) from the wing root, so the strut tip meets the wing break in the Trefftz plane.
The circulation `gc` is set by [`set_strut_circulation_shape!`](@ref).
"""
function generate_strut_points!(
    i_start::Int,
    geom::TrefftzGeometry,
    wing,
    panels::SurfaceDiscretization,
    wing_root_contraction::Float64
)
    yo = 0.5 * wing.layout.root_span
    yop = yo * wing_root_contraction
    i_end = generate_trefftz_points_inclined!(
        i_start, geom,
        yo, wing.layout.z - wing.strut.z,
        0.5 * wing.layout.break_span, wing.layout.z,
        n_total(panels))

    @inbounds for i = i_start:i_end
        geom.yp[i] = get_wake_contraction(WINGWAKE, geom.y[i], yo, yop)
    end
    @inbounds for i = i_start:(i_end-1)
        geom.ycp[i] = get_wake_contraction(WINGWAKE, geom.yc[i], yo, yop)
    end
    return i_end
end

"""
    set_strut_circulation_shape!(geom::TrefftzGeometry, i_first::Int, i_last::Int, spanload_exponent)

Set the bound circulation `gc` of the strut control points `i_first:i_last-1` to ``(4\\xi(1-\\xi))^p`` with
``p`` = `spanload_exponent` and ``\\xi = (k-1/2)/n`` (0 at the fuselage, 1 at the wing). The loading is nonzero on every
panel. [`scale_circulation!`](@ref) then scales the shape to the strut lift.
"""
function set_strut_circulation_shape!(geom::TrefftzGeometry, i_first::Int, i_last::Int, spanload_exponent)
    n = i_last - i_first
    @inbounds for k = 1:n
        ξ = (k - 0.5) / n
        geom.gc[i_first+k-1] = (4.0 * ξ * (1.0 - ξ))^spanload_exponent
    end
    return nothing
end

"""
    _build_trefftz_geometry!(wing, htail, po, gammat, gammas,
    trefftz_config::TrefftzPlaneConfig)

Construct the Trefftz plane geometry for the given wing and horizontal tail and, if present, the strut.
"""
function _build_trefftz_geometry!(wing, htail, po, gammat, gammas,
    trefftz_config::TrefftzPlaneConfig)

    i::Int64 = 0

    nsurf = n_surfaces(trefftz_config)
    @inbounds for isurf = 1:nsurf
        i = i + 1
        if isurf == 3
            # Generate points for the strut
            i = generate_strut_points!(i, TREFFTZ_GEOM, wing, trefftz_config.strut_panels, trefftz_config.wing_root_contraction)
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