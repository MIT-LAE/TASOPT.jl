"""
    WingSectionDimensions{T<:Real}

Structural dimensions of a sized wing spar-box section.
Fields are ordered as: thicknesses first (web, cap), then areas (web, cap).

!!! details "Fields"
    - `tbweb::T`: Normalized Web thickness
    - `tbcap::T`: Normalized Cap thickness
    - `Abweb::T`: Normalized Web cross-sectional area (both webs combined)
    - `Abcap::T`: Normalized Cap cross-sectional area (both caps combined)
"""
struct WingSectionDimensions{T<:Real}
    tbweb::T
    tbcap::T
    Abweb::T
    Abcap::T
end

"""
    size_wing_section!(section, sweep, sigfac; t_cap_min = 0.0, t_web_min = 0.0, P_axial = 0.0)

Calculates Loads and thicknesses for wing sections

An axial force `P_axial` in the section is carried by the caps together with the bending moment. With the
caps-only idealisation (two caps, section modulus ``Z = A_{cap} h_{rms} c_\\perp / 2``) the cap stress is
``\\sigma = M / Z + P / A_{cap}``, which is the stress of the positive bending moment
``M_{eff} = M + |P|\\, h_{rms}\\, c_\\perp / 2`` with ``c_\\perp = c_s \\cos\\Lambda``. The section is sized for ``M_{eff}``.

The section's absolute minimum gauges `min_cap_thickness` and `min_web_thickness` [m] are divided by ``c_\\perp`` and,
if larger than `t_cap_min` or `t_web_min`, replace them as the minimum normalized thicknesses.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `section::TASOPT.structures.Wing.WingSection`: Wing Section to be sized.
    - `sweep::Float64`: Wing sweep.
    - `sigfac::Float64`: Stress factor
    - `t_cap_min::Float64`: (keyword, default 0.0) Minimum cap thickness
    - `t_web_min::Float64`: (keyword, default 0.0) Minimum web thickness
    - `P_axial::Float64`: (keyword, default 0.0) Axial force in the section [N]; tension or compression add to the cap stress

    **Outputs:** [`WingSectionDimensions`](@ref)
"""
function size_wing_section!(section, sweep, sigfac; t_cap_min = 0.0, t_web_min = 0.0, P_axial = 0.0)
    shear_load = section.max_shear_load
    moment = section.moment

    cross_section = section.cross_section

    Eweb = section.webs.material.E
    Ecap = section.caps.material.E
    Gcap = section.caps.material.G
    Gweb = section.webs.material.G
    sigcap = section.caps.material.σmax * sigfac
    tauweb = section.webs.material.τmax * sigfac

    cs = section.co
    cosL = cosd(sweep)
    cp = cs*cosL

    h_avg, h_rms = get_average_sparbox_heights(section.cross_section)

    web_height = cross_section.web_to_box_height * cross_section.thickness_to_chord

    t_web_floor = max(t_web_min, section.min_web_thickness / cp)
    t_cap_floor = max(t_cap_min, section.min_cap_thickness / cp)

    tbweb, Abweb = size_web(tauweb, shear_load, cp, web_height; t_min = t_web_floor)
    
    moment_effective = moment + abs(P_axial) * h_rms * cp / 2
    tbcap, Abcap = size_cap(sigcap, moment_effective, cross_section.thickness_to_chord,
        cross_section.width_to_chord, h_rms, cs, cosL; t_min = t_cap_floor)

    # EI_xx
    section.EI[1] = Ecap * cp^4 * (h_rms^3 - (h_rms - 2.0 * tbcap)^3) * cross_section.width_to_chord / 12.0 +
                                 Eweb * cp^4 * tbweb * web_height^3 / 6.0
    # EI_yy 
    section.EI[4] = Ecap * cp^4 * tbcap * cross_section.width_to_chord^3 / 6.0 +
                                Eweb * cp^4 * tbweb * web_height * 0.5 * cross_section.width_to_chord^2
    # println("GJ SURFW: cp = $cp, width_to_chord = $(cross_section.width_to_chord), tbweb = $tbweb, h_avg = $h_avg, tbcap = $tbcap, web_to_box_height = $(cross_section.web_to_box_height), thickness_to_chord = $(cross_section.thickness_to_chord), Gweb = $Gweb, Gcap = $Gcap")
    section.GJ = cp^4 * 2.0*((cross_section.width_to_chord-tbweb)*(h_avg-tbcap))^2 /
        (  (cross_section.web_to_box_height*section.cross_section.thickness_to_chord-tbcap)/(Gweb*tbweb) +
        (   cross_section.width_to_chord -tbweb)/(Gcap*tbcap) )

    return WingSectionDimensions(tbweb, tbcap, Abweb, Abcap)
end

"""
    Cantilever

Marker type for a wing without a strut. Its inboard loads follow from the outboard loads and the
inboard panel weights; see [`wing_section_loads`](@ref).
"""
struct Cantilever end

"""
    Braced

Marker type for a wing braced by a strut. Its inboard loads follow from the outboard loads, the inboard
panel weights and the strut force; see [`wing_case_loads`](@ref).
"""
struct Braced end

"""
    outboard_section_loads!(wing, po, gammat, gammas, Nload, We, neout, dyeout)

Sets the outboard section shear and bending moment at the wing break (strut attachment point), from the
outer-wing lift distribution, the tip roll-off lift loss, the outboard panel weight and any outboard engines.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure.
    - `po::Float64`: Wing root loading magnitude.
    - `gammat::Float64`: Wing outer section taper ratio.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `Nload::Int`: Max vertical load factor for wing bending loads.
    - `We::Float64`: Weight of one engine [N].
    - `neout::Int`: Number of engines per wing mounted outboard of the wing break.
    - `dyeout::Float64`: Average spanwise moment arm of the outboard engines about the wing break station ηs [m].

    **Outputs:**
    - `wing.outboard.max_shear_load` and `wing.outboard.moment` are updated in place.
"""
function outboard_section_loads!(wing, po, gammat, gammas, Nload, We, neout, dyeout)
    etas = wing.layout.ηs

    # Tip roll off Lift (modeled as a point load) and it's moment about ηs
    dLt = wing.tip_lift_loss * po * wing.layout.root_chord * gammat * wing.outboard.λ
    dMt = dLt * 0.5 * wing.layout.span * (1.0 - etas)

    # Outboard section:
    #---- strut-attach shear,moment from outer-wing loading. 
    #     Note added term to account for any outboard engines.
    #     If neout = 0 this simplifies to Drela's version which assumes engine
    #     fixed at ηs locations only.
    wing.outboard.max_shear_load = (po * wing.layout.span / 4.0) * (gammas + gammat) * (1.0 - etas) +
                                   dLt - Nload * wing.outboard.weight - Nload * neout * We
    wing.outboard.moment = (po * wing.layout.span^2 / 24.0) * (gammas + 2.0 * gammat) * (1.0 - etas)^2 +
                           dMt - Nload * wing.outboard.dyW - Nload * neout * We * dyeout
    return nothing
end # outboard_section_loads!

"""
    wing_engine_arms(wing)

Returns the number of engines per wing and their spanwise moment arms, in the form that the wing load
functions take, for the one engine per wing at `wing.η_engine`. An engine at or inboard of the wing break
`wing.layout.ηs` counts as an inboard engine and an engine outboard of the wing break as an outboard engine.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the engine location `η_engine`, the span and the wing break
      and wing root stations in `wing.layout`.

    **Outputs:**
    - `(neinn, dyeinn, neout, dyeout)`:
      - `neinn::Int`: Number of engines per wing mounted at or inboard of the wing break.
      - `dyeinn::Float64`: Spanwise moment arm of the inboard engine about the wing root station ηo [m], zero for an engine inboard of ηo.
      - `neout::Int`: Number of engines per wing mounted outboard of the wing break.
      - `dyeout::Float64`: Spanwise moment arm of the outboard engine about the wing break station ηs [m].
"""
function wing_engine_arms(wing)
    half_span = 0.5 * wing.layout.span
    if wing.η_engine <= wing.layout.ηs
        return (neinn = 1, dyeinn = max(half_span * (wing.η_engine - wing.layout.ηo), 0.0), neout = 0, dyeout = 0.0)
    else
        return (neinn = 0, dyeinn = 0.0, neout = 1, dyeout = half_span * (wing.η_engine - wing.layout.ηs))
    end
end # wing_engine_arms

"""
    wing_section_loads(::Cantilever, wing, po, gammas, Nload, We, neinn, dyeinn)

Sets the inboard section shear and bending moment of a wing without a strut, from the outboard
section loads (`wing.outboard.max_shear_load`, `wing.outboard.moment`), the inboard lift distribution
and the inboard panel weight, including bending relief from inboard engines.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the outboard section loads already set.
    - `po::Float64`: Wing root loading magnitude.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `Nload::Int`: Max vertical load factor for wing bending loads.
    - `We::Float64`: Weight of one engine [N].
    - `neinn::Int`: Number of engines per wing mounted at or inboard of the wing break.
    - `dyeinn::Float64`: Average spanwise moment arm of the inboard engines about the wing root station ηo [m].

    **Outputs:**
    - `wing.inboard.max_shear_load` and `wing.inboard.moment` are updated in place.
"""
function wing_section_loads(::Cantilever, wing, po, gammas, Nload, We, neinn, dyeinn)
    etao = wing.layout.ηo
    etas = wing.layout.ηs

    # Modifed to account for bending relief from multiple engines.
    # dyeinn allows engine to be in locations other than ηs
    So = wing.outboard.max_shear_load - Nload*neinn*We +
        0.25*po*wing.layout.span*(1.0+gammas)*(etas-etao) -
        Nload*wing.inboard.weight
    Mo = wing.outboard.moment + wing.outboard.max_shear_load*0.5*wing.layout.span*(etas-etao) +
        (1.0/24.0)*po*wing.layout.span^2*(1.0+2.0*gammas)*(etas-etao)^2 -
        Nload*wing.inboard.dyW - Nload*neinn*We*dyeinn

    #----- limit So,Mo to Ss,Ms, which might be needed with heavy outboard engine
    #-      (rules out negatively-tapered structure, deemed not feasible for downloads)
    wing.inboard.max_shear_load = max(So ,wing.outboard.max_shear_load)
    wing.inboard.moment = max(Mo ,wing.outboard.moment)
    return nothing
end # wing_section_loads

"""
    wing_case_loads(::Cantilever, wing, po, gammat, gammas, Nload, We, neout, dyeout, neinn, dyeinn)

Returns the root and break shear and bending moment of a wing without a strut for one load case, as a
[`WingCaseLoads`](@ref). The outboard loads follow from [`outboard_section_loads!`](@ref) and the inboard
loads from [`wing_section_loads`](@ref).

This function sizes nothing, but it mutates the section load fields of the `wing` passed in
(`wing.outboard.max_shear_load`, `wing.outboard.moment`, `wing.inboard.max_shear_load` and
`wing.inboard.moment`); pass a copy to keep the wing loads unchanged.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the panel weights of the load case.
    - `po::Float64`: Wing root loading magnitude.
    - `gammat::Float64`: Wing outer section taper ratio.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `Nload::Int`: Max vertical load factor for wing bending loads.
    - `We::Float64`: Weight of one engine [N].
    - `neout::Int`: Number of engines per wing mounted outboard of the wing break.
    - `dyeout::Float64`: Average spanwise moment arm of the outboard engines about the wing break station ηs [m].
    - `neinn::Int`: Number of engines per wing mounted at or inboard of the wing break.
    - `dyeinn::Float64`: Average spanwise moment arm of the inboard engines about the wing root station ηo [m].

    **Outputs:**
    - `WingCaseLoads`: Shear and moment at the wing root and at the wing break.
"""
function wing_case_loads(::Cantilever, wing, po, gammat, gammas, Nload, We, neout, dyeout, neinn, dyeinn)
    outboard_section_loads!(wing, po, gammat, gammas, Nload, We, neout, dyeout)
    wing_section_loads(Cantilever(), wing, po, gammas, Nload, We, neinn, dyeinn)
    return WingCaseLoads(S_root = wing.inboard.max_shear_load, M_root = wing.inboard.moment,
                         S_break = wing.outboard.max_shear_load, M_break = wing.outboard.moment)
end # wing_case_loads

"""
    strut_force_compatibility(M0, EI, L, f_strut; npts = 21)

Returns the vertical reaction `R` of a strut on a wing bay, from the compatibility of the wing and strut
deflections at the strut attachment. The bay is a beam clamped at the wing root (``x = 0``) and supported by
the strut at ``x = L``. With ``M_0(x)`` the bending moment of the unbraced bay, the moment of the braced bay is
``M(x) = M_0(x) - R\\,(L - x)`` and the strut reaction follows from a zero vertical displacement of the wing
minus the strut stretch at the attachment:
```math
R = \\frac{\\delta_0}{f_{ww} + f_s}, \\qquad
\\delta_0 = \\int_0^L \\frac{M_0\\,(L - x)}{EI}\\,dx, \\qquad
f_{ww} = \\int_0^L \\frac{(L - x)^2}{EI}\\,dx
```
where ``\\delta_0`` is the deflection of the unbraced bay at the attachment, ``f_{ww}`` the flexibility of the
bay there and ``f_s`` the strut flexibility `f_strut`. Both integrals use the trapezoid rule with `npts` points.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `M0::Function`: Bending moment of the unbraced bay [N m] at the distance `x` from the root, `M0(x)`.
    - `EI::Function`: Bending stiffness of the bay [N m²] at the distance `x` from the root, `EI(x)`.
    - `L::Float64`: Length of the bay, from the root to the strut attachment [m].
    - `f_strut::Float64`: Axial flexibility of the strut, referred to the vertical displacement at the attachment [m/N].
    - `npts::Int`: (keyword, default 21) Number of points of the trapezoid rule.

    **Outputs:**
    - `R::Float64`: Vertical strut reaction on the wing, against the wing lift [N], with the sign of the load factor.
"""
function strut_force_compatibility(M0, EI, L, f_strut; npts = 21)
    dx = L / (npts - 1)
    d0 = 0.0
    f_ww = 0.0
    for i in 1:npts
        x = (i - 1) * dx
        weight = (i == 1 || i == npts) ? 0.5 * dx : dx # trapezoid weights
        d0 += weight * M0(x) * (L - x) / EI(x)
        f_ww += weight * (L - x)^2 / EI(x)
    end
    return d0 / (f_ww + f_strut)
end # strut_force_compatibility

"""
    euler_buckling_load(EI, K, length)

Returns the Euler buckling load ``P_{cr} = \\pi^2 EI / (K l)^2`` [N] of a member of bending stiffness `EI` [N m²],
effective-length factor `K` [-] and length `length` [m].
"""
function euler_buckling_load(EI, K, length)
    return pi^2 * EI / (K * length)^2
end # euler_buckling_load

"""
    size_member_box(caps, webs; chord, height_to_chord, width_to_chord, web_to_box_height,
        tension, compression, buckling_length, K, shear_load, bending_moment = 0.0, sigfac, min_cap_thickness, min_web_thickness)

Sizes the caps and webs of a rectangular box member (strut or jury strut). The caps carry the whole axial force, in tension,
in compression and against Euler buckling, and the webs carry the shear load. The caps also carry the bending moment of a
transverse load (the own weight of the strut), for yield only. Lengths are normalised by the member chord ``c``,
as in [`WingCrossSection`](@ref). With the box height ``h`` (`height_to_chord`), width ``w`` (`width_to_chord`) and web
height ``h_w`` = `web_to_box_height` ``h``, the caps (top and bottom, total) and the webs (two, total) have the areas
```math
A_{cap} = 2\\,\\bar t_{cap}\\,w\\,c^2, \\qquad A_{web} = 2\\,\\bar t_{web}\\,h_w\\,c^2
```
and the bending stiffnesses about the box axes are
```math
EI_{xx} = c^4\\left[E_{cap}\\,\\frac{w\\left(h^3 - (h - 2\\bar t_{cap})^3\\right)}{12} + E_{web}\\,\\frac{\\bar t_{web} h_w^3}{6}\\right], \\qquad
EI_{yy} = c^4\\left[E_{cap}\\,\\frac{\\bar t_{cap} w^3}{6} + E_{web}\\,\\frac{\\bar t_{web} h_w w^2}{2}\\right]
```
with ``EI_{min} = \\min(EI_{xx}, EI_{yy})``. The sweep of the member is ignored.

The cap thickness ``\\bar t_{cap}`` is the largest of
- the yield thickness ``F / (2 w c^2 \\sigma_{cap})``, from the larger ``F`` of the tension and the compression and
  ``\\sigma_{cap}`` = `caps.material.σmax * sigfac`. With a `bending_moment` ``M > 0`` it is the larger of this thickness and
  the thickness of [`calc_cap_thickness`](@ref) for the effective moment ``M + F h c / 2`` (the moment ``M`` plus the axial force
  ``F`` on the lever arm of the caps, as for the wing caps in [`size_wing_section!`](@ref)) of a box of height ``h`` and width ``w``,
- the thickness for which both ``EI_{xx}`` (flapwise bending, box height direction) and ``EI_{yy}`` (chordwise bending,
  box width direction) reach ``f_u |P_c| (K l)^2 / \\pi^2``, and
- the minimum gauge `min_cap_thickness / chord`.

The web thickness is the larger of [`calc_web_thickness`](@ref) for the shear load and the shear stress
`webs.material.τmax * sigfac`, and the minimum gauge `min_web_thickness / chord`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `caps::StructuralMember`: Cap material.
    - `webs::StructuralMember`: Web material.
    - `chord::Float64`: Chord of the member [m].
    - `height_to_chord::Float64`: Box height to chord ratio, the thickness-to-chord ratio [-].
    - `width_to_chord::Float64`: Box width to chord ratio [-].
    - `web_to_box_height::Float64`: Web height to box height ratio [-].
    - `tension::Float64`: Axial force of the most tensile case, zero or positive [N].
    - `compression::Float64`: Axial force of the most compressive case, zero or negative [N].
    - `buckling_length::Float64`: Member length for Euler buckling [m].
    - `K::Float64`: Effective-length factor [-].
    - `shear_load::Float64`: Shear load of the webs [N].
    - `bending_moment::Float64`: (keyword, default `0.0`) Bending moment of the transverse load carried by the caps [N m].
    - `sigfac::Float64`: Stress factor applied to the material strengths [-].
    - `min_cap_thickness::Float64`: Minimum cap gauge [m].
    - `min_web_thickness::Float64`: Minimum web gauge [m].

    **Outputs:**
    - `tbcap::Float64`: Cap thickness to chord ratio [-].
    - `tbweb::Float64`: Web thickness to chord ratio [-].
    - `area_caps::Float64`: Material area of the caps [m²].
    - `area_webs::Float64`: Material area of the webs [m²].
    - `area::Float64`: Material area of the caps and webs [m²].
    - `EI_min::Float64`: Smaller of the bending stiffnesses ``EI_{xx}`` and ``EI_{yy}`` of the final box [N m²].
    - `buckling_load::Float64`: Euler buckling load of the box [N].
    - `criterion::Symbol`: Requirement that sets the cap thickness, `:tension` (yield), `:buckling` or `:gauge`.
    - `buckling_axis::Symbol`: Bending axis whose buckling requirement needs the larger cap thickness, `:xx` (flapwise),
    `:yy` (chordwise) or `:none` without compression.
"""
function size_member_box(caps, webs; chord, height_to_chord, width_to_chord, web_to_box_height,
        tension, compression, buckling_length, K, shear_load, bending_moment = 0.0, sigfac, min_cap_thickness, min_web_thickness)
    c, h, w = chord, height_to_chord, width_to_chord
    hw = web_to_box_height * h
    E_cap, E_web = caps.material.E, webs.material.E

    # Webs: shear, with the minimum gauge
    tbweb = max(calc_web_thickness(webs.material.τmax * sigfac, shear_load, c, hw), min_web_thickness / c)

    # Caps: yield, Euler buckling and minimum gauge
    F = max(tension, abs(compression))
    σ_cap = caps.material.σmax * sigfac
    tb_yield = F / (2 * w * c^2 * σ_cap)
    if bending_moment > 0.0
        # the axial force acts on the lever arm h c / 2 of the caps, as for the wing caps
        try
            tb_yield = max(tb_yield, calc_cap_thickness(σ_cap, bending_moment + F * h * c / 2, h, w, h, c, 1.0))
        catch err
            err isa ErrorException || rethrow()
            error("Strut cap sizing infeasible: the caps cannot carry the bending moment within the box height. ",
                "Increase the chord or the thickness-to-chord ratio of the member.")
        end
    end
    tb_gauge = min_cap_thickness / c

    # Euler buckling: the caps must give both EI_xx (flapwise, box height direction) and EI_yy (chordwise, box width
    # direction) at least the required stiffness; the more constraining axis sets the cap thickness
    tb_buckling = 0.0
    buckling_axis = :none
    if compression < 0.0
        EI_req = caps.material.safety_factor * abs(compression) * (K * buckling_length)^2 / pi^2

        # Stiffness the caps must add to the web contribution of each axis, floored at zero if the webs alone suffice
        EI_caps_xx = max(EI_req - E_web * c^4 * tbweb * hw^3 / 6, 0.0)
        EI_caps_yy = max(EI_req - E_web * c^4 * tbweb * hw * w^2 / 2, 0.0)

        # EI_xx cap term: E_cap c^4 w (h^3 - (h - 2 tb)^3) / 12 = EI_caps_xx
        X = 12 * EI_caps_xx / (E_cap * c^4 * w)
        X < h^3 || error("Strut cap sizing infeasible: EI_xx (flapwise bending, box height direction) cannot reach ",
            "the buckling stiffness within the box height. Increase the chord or the thickness-to-chord ratio of the member.")
        tb_xx = 0.5 * (h - cbrt(h^3 - X))

        # EI_yy cap term: E_cap c^4 tb w^3 / 6 = EI_caps_yy
        tb_yy = 6 * EI_caps_yy / (E_cap * c^4 * w^3)
        tb_yy <= h / 2 || error("Strut cap sizing infeasible: EI_yy (chordwise bending, box width direction) cannot reach ",
            "the buckling stiffness within the box height. Increase the chord or the width-to-chord ratio of the member.")

        tb_buckling = max(tb_xx, tb_yy)
        buckling_axis = tb_xx >= tb_yy ? :xx : :yy
    end

    tbcap = max(tb_yield, tb_buckling, tb_gauge)
    criterion = tbcap == tb_yield ? :tension : (tbcap == tb_buckling ? :buckling : :gauge)

    EI_xx = E_cap * c^4 * w * (h^3 - (h - 2 * tbcap)^3) / 12 + E_web * c^4 * tbweb * hw^3 / 6
    EI_yy = E_cap * c^4 * tbcap * w^3 / 6 + E_web * c^4 * tbweb * hw * w^2 / 2
    EI_min = min(EI_xx, EI_yy)

    area_caps = 2 * tbcap * w * c^2
    area_webs = 2 * tbweb * hw * c^2
    return (tbcap = tbcap, tbweb = tbweb, area_caps = area_caps, area_webs = area_webs, area = area_caps + area_webs,
        EI_min = EI_min, buckling_load = euler_buckling_load(EI_min, K, buckling_length), criterion = criterion,
        buckling_axis = buckling_axis)
end # size_member_box

"""
    strut_buckling_length(strut_length, eta_strut)
    strut_buckling_length(strut_length, ::Nothing)

Returns the length of the main strut that is used for its Euler buckling. A vertical jury strut at the fraction
`eta_strut` of the strut length, from the fuselage end, holds the strut laterally, so the strut buckles over the longer of
its two parts, ``\\max(\\eta_{strut}, 1 - \\eta_{strut})\\,l``. Without a jury strut (`nothing`), it is the full strut length.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `strut_length::Float64`: Strut length [m].
    - `eta_strut::Union{Nothing,Float64}`: Position of the jury strut along the strut, as a fraction of the strut length, or `nothing`.

    **Outputs:**
    - `l_buckling::Float64`: Strut length for buckling [m].
"""
strut_buckling_length(strut_length, eta_strut) = max(eta_strut, 1.0 - eta_strut) * strut_length
strut_buckling_length(strut_length, ::Nothing) = strut_length

"""
Jury strut axial load as a fraction of the compression of the main strut [-].
A conventional bracing-force fraction for a member that only restrains the main strut laterally; an assumption, not an input.
"""
const jury_load_fraction = 0.02

"""
    strut_chord_at(strut, fraction)

Chord of the main strut [m] at the `fraction` of its length from the fuselage end (0) to the wing end (1).
The chord varies linearly from `strut.root_chord` to `strut.chord`.
"""
strut_chord_at(strut, fraction) = strut.root_chord + (strut.chord - strut.root_chord) * fraction

"""
    size_jury!(wing, etao, etas, sigfac, sinL, cosL)

Sizes the vertical jury strut `wing.jury` of a truss-braced wing and sets its geometry, axial force, area, buckling load,
weight and weight moment. The jury is a lateral support only, so its load is the fraction `jury_load_fraction` of the compression
of the main strut (`wing.strut.axial_force_compression`), taken as compression. Its length is ``z (1 - \\eta_{strut})`` with the
strut height ``z`` (`wing.strut.z`) and its wing attachment is ``\\eta_o + \\eta_{strut} (\\eta_s - \\eta_o)``.
The jury is a pinned (``K = 1``) box of caps and webs, sized by [`size_member_box`](@ref) for this compression over its length,
without tension or shear (the webs are at their minimum gauge). Its chord is `jury.chord_to_strut_chord` times the chord of the strut at the jury attachment ([`strut_chord_at`](@ref)), and the
thickness-to-chord, width and web height ratios are those of the main strut. `jury.area` is the material area of caps and webs.
The weight is the member weight of caps and webs of n_wings times ``1 + f_{added}`` (`jury.weight_frac_added`), with the bare areas. Its moment `jury.dxW` is the weight times the
streamwise offset of the wing attachment station `jury.eta_wing` from the wing box root, ``0.5\\,b\\,(\\eta_{wing} - \\eta_o) \\tan\\Lambda``,
with the wing span ``b`` and sweep ``\\Lambda``.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the main strut already sized.
    - `etao::Float64`: Non-dimensional span coordinate of the wing root station.
    - `etas::Float64`: Non-dimensional span coordinate of the strut attachment.
    - `sigfac::Float64`: Stress Factor.
    - `sinL::Float64`: Sine of the wing sweep.
    - `cosL::Float64`: Cosine of the wing sweep.
"""
function size_jury!(wing, etao, etas, sigfac, sinL, cosL; n_wings = 2.0)
    jury, strut = wing.jury, wing.strut

    jury.eta_wing = etao + jury.eta_strut * (etas - etao)
    jury.length = strut.z * (1.0 - jury.eta_strut) # only correct if the wing and truss/strut have the same sweep and are in the same plane. 3D effects on length missing
    jury.axial_force = -jury_load_fraction * abs(strut.axial_force_compression)

    jury.chord = jury.chord_to_strut_chord * strut_chord_at(strut, jury.eta_strut)
    box = size_member_box(jury.caps, jury.webs;
        chord = jury.chord, height_to_chord = strut.thickness_to_chord, width_to_chord = strut.box_width_to_chord,
        web_to_box_height = strut.web_to_box_height, tension = 0.0, compression = jury.axial_force,
        buckling_length = jury.length, K = 1.0, shear_load = 0.0, sigfac = sigfac,
        min_cap_thickness = jury.min_cap_thickness, min_web_thickness = jury.min_web_thickness)
    jury.area = box.area
    jury.buckling_load = box.buckling_load

    jury.weight = n_wings * gee * jury.length * (jury.caps.material.ρ * box.area_caps + jury.webs.material.ρ * box.area_webs) *
                  (1.0 + jury.weight_frac_added)
    jury.dxW = jury.weight * 0.5 * wing.layout.span * (jury.eta_wing - etao) * sinL / cosL
    return nothing
end # size_jury!

"""
    unbraced_inboard_loads(wing, po, gammas, Nload, We, neinn, dyeinn, x; F_wing = 0.0)

Returns the shear and bending moment at the distance `x` from the wing root station ηo, in the inboard bay
of the wing clamped at the root without a strut. The bay is loaded by the section loads at the wing break
(`wing.outboard.max_shear_load`, `wing.outboard.moment`), a linear spanload from `po` at the root to `po*gammas`
at the break, the inboard panel weight, acting at the moment arm `wing.inboard.dyW / wing.inboard.weight`
from the root, and the inboard engines `neinn`, acting at `dyeinn`. The point load `F_wing` at the strut attachment,
upward like the lift of the wing (the part of the strut lift that the strut passes to the wing), loads the bay
too. Only the loads outboard of `x` count, so `F_wing` counts for `x` inboard of the attachment and adds `F_wing (L - x)`
to the moment and `F_wing` to the shear. At `x = 0` they are the root loads of [`wing_section_loads`](@ref), before its limiting to the break loads.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the outboard section loads already set.
    - `po::Float64`: Wing root loading magnitude.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `Nload::Int`: Max vertical load factor for wing bending loads.
    - `We::Float64`: Weight of one engine [N].
    - `neinn::Int`: Number of engines per wing mounted at or inboard of the wing break.
    - `dyeinn::Float64`: Average spanwise moment arm of the inboard engines about the wing root station ηo [m].
    - `x::Float64`: Distance from the wing root station ηo, between 0 and the bay length [m].
    - `F_wing::Float64`: (keyword, default 0.0) Upward point load at the strut attachment, at the end of the bay [N].

    **Outputs:**
    - `(S, M)`: Shear [N] and bending moment [N m], positive for upward lift loading.
"""
function unbraced_inboard_loads(wing, po, gammas, Nload, We, neinn, dyeinn, x; F_wing = 0.0)
    L = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
    u = L - x # length of the bay outboard of x

    # Spanload at x and its linear decrease towards the break
    p_x = po * (1.0 - (1.0 - gammas) * x / L)
    dp = po * (1.0 - gammas) / L

    S = wing.outboard.max_shear_load + p_x * u - dp * u^2 / 2
    M = wing.outboard.moment + wing.outboard.max_shear_load * u + p_x * u^2 / 2 - dp * u^3 / 3

    if x < L
        S += F_wing
        M += F_wing * u
    end
    if wing.inboard.weight != 0.0 && x <= wing.inboard.dyW / wing.inboard.weight
        S -= Nload * wing.inboard.weight
        M -= Nload * (wing.inboard.dyW - wing.inboard.weight * x)
    end
    if x <= dyeinn
        S -= Nload * neinn * We
        M -= Nload * neinn * We * (dyeinn - x)
    end
    return (S = S, M = M)
end # unbraced_inboard_loads

"""
    strut_reaction(wing, po, gammas, M0, L, ls, F_wing)

Returns the vertical strut reaction on the inboard bay of a braced wing, positive upward, for the strut
force model `wing.strut.opt_strut_force`:
- `Compatibility`: [`strut_force_compatibility`](@ref) with the unbraced bay moment `M0(x)`, a bending stiffness
  that varies linearly from `wing.inboard.EI[1]` at the root to `wing.outboard.EI[1]` at the strut attachment
  and the axial flexibility of the strut referred to the vertical displacement at the attachment
  ```math
  f_s = \\frac{l_s}{E_s\\,A_s\\,\\sin^2\\theta}, \\qquad \\sin\\theta = \\frac{z}{l_s}
  ```
  with the strut modulus ``E_s`` (`wing.strut.caps.material.E`), the strut area ``A_s`` (`wing.strut.area`) and the
  strut length ``l_s`` and height ``z`` (`wing.strut.z`) in the plane of the strut.
  The stiffness and the strut area are those of the previous sizing iteration.
- `Prescribed`: the rule `R = (po b / 12) (ηs - ηo) (1 + 2 γs) + S_break + F_wing`, which makes the bending moment at the
  root equal to the one at the break; the point load `F_wing` is part of the shear that the strut supports.

The strut force, in the direction of the strut, is `R ls / z`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the outboard section loads already set.
    - `po::Float64`: Wing root loading magnitude.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `M0::Function`: Bending moment of the unbraced inboard bay [N m] at the distance `x` from the root, `M0(x)`.
    - `L::Float64`: Length of the inboard bay, from the root to the strut attachment [m].
    - `ls::Float64`: Length of the strut in the plane of the strut [m].
    - `F_wing::Float64`: Upward point load at the strut attachment from the lift of the strut [N], see [`strut_wing_end_lift`](@ref).

    **Outputs:**
    - `R::Float64`: Vertical strut reaction on the wing [N].
"""
function strut_reaction(wing, po, gammas, M0, L, ls, F_wing)
    EI_root, EI_attach, A_strut = wing.inboard.EI[1], wing.outboard.EI[1], wing.strut.area

    # TEMPORARY FALLBACK: the stiffness and the strut area are zero in the first sizing iteration, so the
    # compatibility force cannot be computed and the prescribed force is used instead. This is a temporary
    # measure and is to be replaced later.
    stiffness_known = EI_root > 0.0 && EI_attach > 0.0 && A_strut > 0.0

    if wing.strut.opt_strut_force == StrutForceModel.Compatibility && stiffness_known
        EI(x) = EI_root + (EI_attach - EI_root) * x / L
        f_strut = ls / (wing.strut.caps.material.E * A_strut * (wing.strut.z / ls)^2)
        return strut_force_compatibility(M0, EI, L, f_strut)
    else
        return (po * wing.layout.span / 12.0) * (wing.layout.ηs - wing.layout.ηo) * (1.0 + 2.0 * gammas) +
               wing.outboard.max_shear_load + F_wing
    end
end # strut_reaction

"""
    strut_lift_on_strut(wing, po, gammat, gammas)

Returns the lift on the strut of one wing for the wing root loading `po` of a wing that is relieved by the strut lift. The wing
carries the share ``1 - k`` of the lift of the wing and the strut, with the strut lift share ``k`` of
[`strut_lift_share`](@ref) (``k = 0`` if the strut has no lift ratio, no area, or the wing has no area yet). The wing lift
of one wing is ``p_o K_p b / 2`` with the span factor ``K_p`` of [`wing_lift_span_factor`](@ref), so the lift on the strut
of one wing is
```math
F_{strut} = \\frac{k}{1 - k}\\, \\frac{p_o K_p b}{2}
```
This is the share ``k`` of the net lift ``(N W - L_h)/2`` of one wing and strut. It has the sign of `po`, which is the sign of the
load factor.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the strut data.
    - `po::Float64`: Wing root loading magnitude, already reduced by the factor ``1 - k``.
    - `gammat::Float64`: Wing outer section taper ratio of the spanload of `po`, `wing.outboard.λ * rclt`.
    - `gammas::Float64`: Wing inner section taper ratio of the spanload of `po`, `wing.inboard.λ * rcls`.

    **Outputs:**
    - `F_strut::Float64`: Lift on the strut of one wing [N].
"""
function strut_lift_on_strut(wing, po, gammat, gammas)
    k = strut_lift_share(wing)
    return k / (1.0 - k) * po * wing_lift_span_factor(wing, gammat, gammas) * wing.layout.span / 2.0
end # strut_lift_on_strut

"""
    strut_wing_end_lift(wing, F_strut)

Returns the part of the lift `F_strut` on the strut of one wing that loads the wing, as a point load at the strut attachment.
The strut has pinned ends and its loading ``(1-\\xi^2)^p`` has the centroid fraction ``\\xi_{cp}`` of
[`strut_lift_centroid_fraction`](@ref) from the fuselage end, so the part at the wing end is ``F_{wing} = \\xi_{cp} F_{strut}``
and the rest, ``(1 - \\xi_{cp}) F_{strut}``, goes to the fuselage.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with `wing.strut.spanload_exponent`.
    - `F_strut::Float64`: Lift on the strut of one wing [N].

    **Outputs:**
    - `F_wing::Float64`: Part of the strut lift at the wing end of the strut [N], with the sign of `F_strut`.
"""
strut_wing_end_lift(wing, F_strut) = strut_lift_centroid_fraction(wing.strut.spanload_exponent) * F_strut

"""
    strut_fuselage_end_force(wing, R, ls, F_strut)

Returns the axial force of the strut at its fuselage end, positive in tension, for the strut reaction `R` on the wing, the strut
length `ls` in the plane of the strut and the lift `F_strut` on the strut of one wing (see [`wing_case_loads`](@ref) for the
equilibrium of the strut):
```math
T_f = R\\,\\frac{l_s}{z} + (1 - \\xi_{cp})\\,F_{strut}\\,\\frac{z}{l_s}
```
The second term is the along-axis component of the part of the strut lift at the fuselage end. `R` and `F_strut` have the same
sign, so ``|T_f|`` is the largest axial force of the strut. Without strut lift it is ``R l_s / z``.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the strut height `wing.strut.z` and `wing.strut.spanload_exponent`.
    - `R::Float64`: Strut reaction on the wing [N], see [`strut_reaction`](@ref).
    - `ls::Float64`: Length of the strut in the plane of the strut [m].
    - `F_strut::Float64`: Lift on the strut of one wing [N], with the sign of the load factor.

    **Outputs:**
    - `T_f::Float64`: Axial force of the strut at the fuselage end [N].
"""
strut_fuselage_end_force(wing, R, ls, F_strut) =
    R * ls / wing.strut.z + (F_strut - strut_wing_end_lift(wing, F_strut)) * wing.strut.z / ls

"""
    wing_case_loads(::Braced, wing, po, gammat, gammas, Nload, We, neout, dyeout, neinn, dyeinn; strut_lift = 0.0)

Returns the root and break shear and bending moment, the strut force and the axial force of a wing braced
by a strut for one load case, as a [`WingCaseLoads`](@ref). The outboard loads follow from
[`outboard_section_loads!`](@ref). The inboard bay is the bay of [`unbraced_inboard_loads`](@ref), clamped at
the root and supported by the strut at the break, where the strut reaction `R` follows from
[`strut_reaction`](@ref). The bay length is ``L = (b/2)(η_s - η_o)`` and the strut height is `wing.strut.z`.

The strut carries the share ``k`` of the lift (see [`strut_lift_share`](@ref)), and `po` is the root loading of the
wing that carries the rest, ``(1-k)`` times the loading of the wing and strut together (see [`wing_case_root_loading`](@ref)).

Equilibrium of the strut. The strut has pinned ends, the fuselage end at the origin and the wing end at ``(L, z)``, and the
vertical lift ``F_{strut}`` on it with the centroid ``\\xi_{cp} L`` from the fuselage end. The strut exerts on the wing the
horizontal force ``-R L / z`` and the vertical force ``-R`` of the reaction, together with ``F_{wing} = \\xi_{cp} F_{strut}``
of the lift (the point load of the wing). 
```math
T_w = R\\,\\frac{l_s}{z} - F_{wing} \\sin\\Gamma \\text{ at the wing end}, \\qquad
T_f = T_w + F_{strut} \\sin\\Gamma = R\\,\\frac{l_s}{z} + (1-\\xi_{cp})\\,F_{strut}\\,\\frac{z}{l_s} \\text{ at the fuselage end}
```
``R`` and ``F_{strut}`` have the same sign, so ``|T_f| \\ge |T_w|``: the fuselage end governs, in tension (positive load factor) and in
compression. The axial force of the strut and the axial force in the inboard wing are
```math
T = T_f, \\qquad P = -R\\,\\frac{l_s}{z}\\,\\frac{L}{l_s} = -R\\,\\frac{L}{z}, \\qquad l_s = \\sqrt{z^2 + L^2}
```
The strut force `T` is positive in tension and negative in compression (negative load factors). The axial
force `P` in the inboard wing is positive in tension, so that an under-wing strut in a pull-up compresses the wing; it
follows from the horizontal force ``R L / z`` of the strut on the wing, which the strut lift does not change.

The returned loads are
- `S_root`, `M_root`: shear and bending moment at the root, ``S_0(0) - R`` and ``M_0(0) - R L`` (the loads of
  the cantilever without its limiting to the break loads, with the strut reaction applied);
- `S_break`, `M_break`: the outboard section loads;
- `inboard_stations`, `inboard_shear`, `inboard_moment`: the distance ``x`` from the root of 11 stations
  between the root and the strut attachment, and the shear ``S_0(x) - R`` and bending moment
  ``M_0(x) - R (L - x)`` there (the shear jumps by `R` at the strut attachment, the moment is continuous; the last
  station is the attachment itself, where the shear is that of the outboard side and does not include the point load `F_wing`);
- `M_inboard_extreme`: the bending moment of the largest magnitude of `inboard_moment`;
- `strut_reaction`: the strut reaction ``R``, which acts on the wing against the wing lift;
- `strut_force` and `P_axial_inboard`: ``T = T_f`` (see [`strut_fuselage_end_force`](@ref)) and ``P``;
- `strut_lift`: the lift on the strut of one wing ``F_{strut}``, with the sign of the load factor.

Limitations: the offset torque of the strut joint on the wing box is not modelled, and the panel weights of
the wing are those of the design cruise mass state in every load case. The strut force and the strut load on the fuselage are those of
the strut and fuselage ends in the plane of the strut: the exact fuselage pin force has the vertical part ``R + (1-\\xi_{cp}) F_{strut}``
and the horizontal part ``R L / z``, whereas `wing.strut.force_vertical` and `wing.strut.force_horizontal` of [`wing_weights!`](@ref)
are the components of ``T_f`` along the strut axis.

This function sizes nothing, but it mutates the outboard section load fields of the `wing` passed in
(`wing.outboard.max_shear_load` and `wing.outboard.moment`); pass a copy to keep the wing loads unchanged.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the panel weights of the load case, the strut
      geometry and the stiffness and strut area of the previous sizing iteration.
    - `po::Float64`: Wing root loading magnitude, reduced by the strut lift share (see above).
    - `gammat::Float64`: Wing outer section taper ratio.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `Nload::Int`: Max vertical load factor for wing bending loads.
    - `We::Float64`: Weight of one engine [N].
    - `neout::Int`: Number of engines per wing mounted outboard of the wing break.
    - `dyeout::Float64`: Average spanwise moment arm of the outboard engines about the wing break station ηs [m].
    - `neinn::Int`: Number of engines per wing mounted at or inboard of the wing break.
    - `dyeinn::Float64`: Average spanwise moment arm of the inboard engines about the wing root station ηo [m].
    - `strut_lift::Float64`: (keyword, default 0.0) Lift on the strut of one wing [N], with the sign of the load factor.

    **Outputs:**
    - `WingCaseLoads`: Shear and moment at the wing root and at the wing break, the shear and moment along the
      inboard wing, the extreme inboard moment, the strut reaction and force, the lift on the strut and the axial force
      in the inboard wing.
"""
function wing_case_loads(::Braced, wing, po, gammat, gammas, Nload, We, neout, dyeout, neinn, dyeinn; strut_lift = 0.0)
    outboard_section_loads!(wing, po, gammat, gammas, Nload, We, neout, dyeout)

    L = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
    z = wing.strut.z
    ls = sqrt(z^2 + L^2)

    F_wing = strut_wing_end_lift(wing, strut_lift)
    inboard_loads(x) = unbraced_inboard_loads(wing, po, gammas, Nload, We, neinn, dyeinn, x; F_wing)
    M0(x) = inboard_loads(x).M
    R = strut_reaction(wing, po, gammas, M0, L, ls, F_wing)

    stations = range(0.0, L; length = 11)
    M_stations = [M0(x) - R * (L - x) for x in stations]
    S_stations = [inboard_loads(x).S - R for x in stations]
    strut_force = strut_fuselage_end_force(wing, R, ls, strut_lift)

    root = inboard_loads(0.0)
    return WingCaseLoads(S_root = root.S - R, M_root = root.M - R * L,
                         S_break = wing.outboard.max_shear_load, M_break = wing.outboard.moment,
                         M_inboard_extreme = M_stations[argmax(abs.(M_stations))],
                         P_axial_inboard = -(R * ls / z) * L / ls, strut_force = strut_force,
                         inboard_stations = collect(stations), inboard_shear = S_stations,
                         inboard_moment = M_stations, strut_reaction = R, strut_lift = strut_lift)
end # wing_case_loads

"""
    far_manoeuvre_load_factor(W_MTO)

Returns the positive manoeuvre load factor required by FAR/CS 25.337 (n+) for an aircraft of maximum takeoff
weight `W_MTO`:
```math
n = 2.1 + \\frac{24000}{W_{lb} + 10000}, \\quad 2.5 \\le n \\le 3.8
```
where ``W_{lb}`` is `W_MTO` in pounds-force. The load factor is limited to 2.5 for heavy aircraft and to 3.8 for light ones.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `W_MTO::Float64`: Maximum takeoff weight [N].

    **Outputs:**
    - `n::Float64`: Positive manoeuvre load factor [-].
"""
far_manoeuvre_load_factor(W_MTO) = clamp(2.1 + 24000 / (W_MTO / 4.4482216152605 + 10000), 2.5, 3.8)

"""
    wing_case_load_factor(name, case, W_MTO, far_load_factor; wing = nothing, W_lift = NaN, Mach_cruise = NaN, V_cruise_EAS = 0.0)

Returns the load factor that a wing sizing load case is evaluated at. The `case` is not modified.
- `Gust` cases: the load factor is ``n = 1 + \\mathrm{sign}(n_{case})\\,\\Delta n``, where `n_case` is `case.load_factor` (only its
  sign is used) and ``\\Delta n`` is the Pratt–Walker increment of [`gust_load_increment`](@ref). 
- Other cases: the load factor is the `case.load_factor` entered by the user, except for the positive manoeuvre cases
  `:pullup` and `:pullup_light` if `far_load_factor` is `true`: their load factor is then
  [`far_manoeuvre_load_factor`](@ref) of `W_MTO`. The gust keywords are ignored and not validated.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `name::Symbol`: Name of the case in `WingSizingLoads`.
    - `case::WingLoadCase`: Wing sizing load case.
    - `W_MTO::Float64`: Maximum takeoff weight of the current sizing iteration [N].
    - `far_load_factor::Bool`: Whether the positive manoeuvre cases use the FAR/CS 25.337 load factor (`Options.opt_far_load_factor`).
    - `wing::TASOPT.structures.Wing`: (keyword, `Gust` cases only) Wing structure, with the layout area `S`, span `span`, `AR` and `sweep`.
    - `W_lift::Float64`: (keyword, `Gust` cases only) Aircraft weight carried by the wing lift in the mass state of the case [N].
    - `Mach_cruise::Float64`: (keyword, `Gust` cases only) Cruise Mach number, for the lift curve slope.
    - `V_cruise_EAS::Float64`: (keyword, `Gust` cases only, default 0.0) Cruise equivalent airspeed [m/s], used if the case has no speed of its own.

    **Outputs:**
    - `n::Float64`: Load factor of the case [-].
"""
function wing_case_load_factor(name, case, W_MTO, far_load_factor;
        wing = nothing, W_lift = NaN, Mach_cruise = NaN, V_cruise_EAS = 0.0)
    if case.kind == WingLoadKind.Gust
        V_EAS = case.gust_speed_EAS > 0 ? case.gust_speed_EAS : V_cruise_EAS
        if !(V_EAS > 0) || !(0 <= Mach_cruise < 1)
            throw(ArgumentError("gust sizing cases need a positive speed (case.gust_speed_EAS or V_cruise_EAS) and a subsonic Mach_cruise; got V_EAS = $V_EAS, Mach_cruise = $Mach_cruise"))
        end
        if wing === nothing || isnan(W_lift)
            throw(ArgumentError("gust sizing cases need the keywords wing and W_lift"))
        end
        S, span = wing.layout.S, wing.layout.span
        CL_alpha = wing_lift_curve_slope(wing.layout.AR, wing.layout.sweep, Mach_cruise)
        Δn = gust_load_increment(W_lift / S, S / span, CL_alpha, case.gust_velocity_EAS, V_EAS)
        return 1 + sign(case.load_factor) * Δn
    end
    if far_load_factor && (name === :pullup || name === :pullup_light)
        return far_manoeuvre_load_factor(W_MTO)
    end
    return case.load_factor
end # wing_case_load_factor

"""
    wing_case_root_loading(case, wing, rclt, rcls, W_lift, Lhtail, n_pullup; load_factor = case.load_factor, strut_share = 0.0)

Returns the load factor and the wing root loading magnitude of one wing sizing load case.

A `Maneuver` case uses the load factor of the case (`load_factor`, by default `case.load_factor`) and a tail load of
`case.tail_load_scale * (load factor / n_pullup) * Lhtail`, and its root loading follows from
[`wing_loading`](@ref). A `Gust` case is treated in the same way, with its computed load factor 1 ± Δn passed as
`load_factor` (see [`wing_case_load_factor`](@ref)). A `Ground` case has the load factor of the case and no wing lift.
If the struts carry the share `strut_share` of the lift (see [`strut_lift_share`](@ref)), the wing carries the rest and
its root loading is `(1 - strut_share)` times the one of the wing and the struts together.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `case::WingLoadCase`: Wing sizing load case.
    - `wing::TASOPT.structures.Wing`: Wing structure.
    - `rclt::Float64`: Tip/root cl ratio (clt/clo).
    - `rcls::Float64`: Break/root cl ratio (cls/clo).
    - `W_lift::Float64`: Aircraft weight carried by the wing lift in the case [N].
    - `Lhtail::Float64`: Worst-case (most negative) tail lift at the pull-up load factor [N].
    - `n_pullup::Float64`: Pull-up load factor that `Lhtail` belongs to.
    - `load_factor::Float64`: (keyword, default `case.load_factor`) Load factor of the case, see [`wing_case_load_factor`](@ref).
    - `strut_share::Float64`: (keyword, default 0.0) Share of the lift carried by the struts, from 0 (cantilever wing) to below 1.

    **Outputs:**
    - `(Nload, po)`: Load factor for wing bending loads and wing root loading magnitude of the wing, reduced by `strut_share`.
"""
function wing_case_root_loading(case, wing, rclt, rcls, W_lift, Lhtail, n_pullup; load_factor = case.load_factor, strut_share = 0.0)
    if case.kind == WingLoadKind.Ground
        return (Nload = load_factor, po = 0.0)
    end
    # Maneuver and Gust cases
    n = load_factor
    Lh = case.tail_load_scale * (n / n_pullup) * Lhtail
    return (Nload = n, po = (1.0 - strut_share) * wing_loading(wing, rclt, rcls, n, W_lift, Lh))
end # wing_case_root_loading

"""
    wing_lift_curve_slope(AR, sweep, Mach)

Returns the lift curve slope of a finite wing from the Helmbold (DATCOM) formula
```math
C_{L_\\alpha} = \\frac{2\\pi AR}{2 + \\sqrt{4 + (AR\\,\\beta)^2 \\left(1 + \\tan^2\\Lambda / \\beta^2\\right)}}, \\qquad \\beta = \\sqrt{1 - M^2}
```
The formula is valid for subsonic Mach numbers. DATCOM uses the mid-chord sweep ``\\Lambda``; the
quarter-chord sweep of the wing layout is used as an approximation of it.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `AR::Float64`: Wing aspect ratio.
    - `sweep::Float64`: Wing sweep angle [deg].
    - `Mach::Float64`: Flight Mach number.

    **Outputs:**
    - `CL_alpha::Float64`: Wing lift curve slope [1/rad].
"""
function wing_lift_curve_slope(AR, sweep, Mach) # standard range of the formulation is up to Mach <= 0.80
    β = sqrt(1 - Mach^2)
    return 2π * AR / (2 + sqrt(4 + (AR * β)^2 * (1 + tand(sweep)^2 / β^2))) # Polhamus formulation
end # wing_lift_curve_slope

"""
    gust_load_increment(W_over_S, c_mean, CL_alpha, U_de, V_EAS)

Returns the incremental load factor of an aircraft in a vertical gust from the Pratt–Walker gust alleviation
formula (https://ntrs.nasa.gov/api/citations/19930090988/downloads/19930090988.pdf)
```math
\\mu = \\frac{2\\,(W/S)}{\\rho_0\\, \\bar{c}\\, C_{L_\\alpha}\\, g}, \\qquad K_g = \\frac{0.88\\,\\mu}{5.3 + \\mu}, \\qquad
\\Delta n = \\frac{\\rho_0\\, U_{de}\\, V_{EAS}\\, C_{L_\\alpha}\\, K_g}{2\\,(W/S)}
```
where ``\\rho_0`` is the sea-level air density. The gust velocity and the airspeed are equivalent airspeeds (EAS).
The wing sizing cases use this function for their gust cases, see [`wing_case_load_factor`](@ref).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `W_over_S::Float64`: Wing loading [N/m²].
    - `c_mean::Float64`: Mean geometric chord, wing area over span [m].
    - `CL_alpha::Float64`: Wing lift curve slope [1/rad].
    - `U_de::Float64`: Derived gust velocity [m/s EAS].
    - `V_EAS::Float64`: Equivalent airspeed [m/s].

    **Outputs:**
    - `Δn::Float64`: Magnitude of the incremental load factor.
"""
function gust_load_increment(W_over_S, c_mean, CL_alpha, U_de, V_EAS)
    ρ0 = atmos(0.0).ρ
    μ = 2 * W_over_S / (ρ0 * c_mean * CL_alpha * gee)
    K_g = 0.88 * μ / (5.3 + μ)
    return ρ0 * U_de * V_EAS * CL_alpha * K_g / (2 * W_over_S)
end # gust_load_increment

"""
    wing_sizing_masses(; W_cruise, W_MTO, W_pay, W_paymax, W_fuel, freserve, W_fuel_cap, has_wing_fuel)

Returns the weights that [`evaluate_wing_sizing_cases!`](@ref) passes to [`wing_mass_state`](@ref), from the
aircraft weights of a sizing iteration. The reserve, empty and fuel-at-MTO weights are simple estimates:
- `W_reserve = W_fuel * freserve / (1 + freserve)`: `freserve` is the reserve fuel over the burned fuel, and
  `W_fuel` includes the reserve.
- `W_empty = W_MTO - W_pay - W_fuel`: the operating empty weight, without payload and fuel.
- `W_fuel_at_MTO = W_pay - W_paymax + W_fuel`: the fuel that fits at the maximum takeoff weight with the
  maximum payload, capped at `W_fuel_cap` if the fuel is stored in the wing.
- `r_fuel_design`: the same fuel weight over `W_fuel_cap` (zero if the fuel is not stored in the wing). It is
  not capped at one.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `W_cruise::Float64`: Aircraft weight at the design cruise point [N].
    - `W_MTO::Float64`: Maximum takeoff weight [N].
    - `W_pay::Float64`: Design payload weight [N].
    - `W_paymax::Float64`: Maximum payload weight [N].
    - `W_fuel::Float64`: Design mission fuel weight, including the reserve [N].
    - `freserve::Float64`: Reserve fuel as a fraction of the burned fuel.
    - `W_fuel_cap::Float64`: Fuel weight that fills the wing tanks [N].
    - `has_wing_fuel::Bool`: Whether fuel is stored in the wing.

    **Outputs:**
    - `masses::NamedTuple`: `W_cruise`, `W_MTO`, `W_fuel`, `W_reserve`, `W_empty`, `W_fuel_cap`, `W_fuel_at_MTO`
      and `r_fuel_design`, as defined above [N, except `r_fuel_design`].
"""
function wing_sizing_masses(; W_cruise, W_MTO, W_pay, W_paymax, W_fuel, freserve, W_fuel_cap, has_wing_fuel)
    W_fuel_design_MTO = W_pay - W_paymax + W_fuel
    return (W_cruise = W_cruise,
            W_MTO = W_MTO,
            W_fuel = W_fuel,
            W_reserve = W_fuel * freserve / (1 + freserve),
            W_empty = W_MTO - W_pay - W_fuel,
            W_fuel_cap = W_fuel_cap,
            W_fuel_at_MTO = has_wing_fuel ? min(W_fuel_cap, W_fuel_design_MTO) : W_fuel_design_MTO,
            r_fuel_design = (has_wing_fuel && W_fuel_cap > 0) ? W_fuel_design_MTO / W_fuel_cap : 0.0)
end # wing_sizing_masses

"""
    evaluate_wing_sizing_cases!(wing, masses, has_wing_fuel, rclt, rcls, Lhtail, gammat, gammas, We, dyeinn;
        has_strut = false, neinn = 1, neout = 0, dyeout = 0.0, far_load_factor = false,
        Mach_cruise = NaN, V_cruise_EAS = 0.0)

Computes the root and break loads of every active wing sizing load case and returns their envelope. The
wing is a cantilever wing, or a wing braced by a strut if `has_strut` is `true`; the loads of each case follow
from [`wing_case_loads`](@ref) for [`Cantilever`](@ref) or [`Braced`](@ref). The loads of each case are stored in
`wing.sizing_loads.results`, the load factor of the case in its `n_computed`, and the cases with the largest root
and break bending moment in `wing.sizing_loads.governing_root` and `wing.sizing_loads.governing_break`.

For a braced wing, the case with the largest strut tension is stored in `wing.sizing_loads.governing_strut`,
and the envelope of the inboard section is the largest over the cases of the root and interior bending moment
and of the root and break shear, instead of the root loads only. The case with the most negative strut force,
if the strut is compressed in any case, is stored in `wing.sizing_loads.governing_strut_compression`.

Active gust cases are evaluated like the manoeuvre cases, with the load factor ``1 \\pm \\Delta n`` of
[`wing_case_load_factor`](@ref). By default the wing has one engine per wing at the inboard side; the keywords
`neinn`, `neout` and `dyeout` set other engine locations (see [`wing_engine_arms`](@ref)). The section loads
of the `wing` are not modified.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure, with the panel weights of the design cruise mass state.
    - `masses::NamedTuple`: Weights for [`wing_mass_state`](@ref): `W_cruise`, `W_MTO`, `W_fuel`, `W_reserve`,
      `W_empty`, `W_fuel_cap`, `W_fuel_at_MTO` and `r_fuel_design`.
    - `has_wing_fuel::Bool`: Whether fuel is stored in the wing.
    - `rclt::Float64`: Tip/root cl ratio (clt/clo).
    - `rcls::Float64`: Break/root cl ratio (cls/clo).
    - `Lhtail::Float64`: Worst-case (most negative) tail lift at the pull-up load factor [N].
    - `gammat::Float64`: Wing outer section taper ratio.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `We::Float64`: Weight of one engine [N].
    - `dyeinn::Float64`: Average spanwise moment arm of the inboard engine about the wing root station ηo [m].
    - `has_strut::Bool`: (keyword, default `false`) Whether the wing is braced by a strut.
    - `neinn::Int`: (keyword, default 1) Number of engines per wing mounted at or inboard of the wing break.
    - `neout::Int`: (keyword, default 0) Number of engines per wing mounted outboard of the wing break.
    - `dyeout::Float64`: (keyword, default 0.0) Average spanwise moment arm of the outboard engines about the wing break station ηs [m].
    - `far_load_factor::Bool`: (keyword, default `false`) Whether the positive manoeuvre cases use the FAR/CS 25.337 load factor.
    - `Mach_cruise::Float64`: (keyword, default `NaN`) Cruise Mach number, for the lift curve slope of the gust cases.
    - `V_cruise_EAS::Float64`: (keyword, default 0.0) Cruise equivalent airspeed [m/s] of the gust cases that have no `gust_speed_EAS` of their own.

    **Outputs:**
    - Cantilever wing: `(S_root, M_root, S_break, M_break)`: Largest magnitude of the shear and bending moment
      at the wing root and at the wing break, each from its governing case, or `nothing` if no case was evaluated.
    - Braced wing: `(S_root, M_root, S_break, M_break, T_strut, T_strut_compression, P_axial)`, where `S_root`
      and `M_root` are the loads for the inboard section sized as described above, `T_strut` is the largest strut
      tension [N] over the cases (0 if the strut is compressed in all cases), `T_strut_compression` the most
      negative strut force [N] over the cases (0 if the strut is not compressed in any case) and `P_axial` the
      axial force of largest magnitude in the inboard wing [N] (negative in compression), or `nothing` if no
      case was evaluated.
"""
function evaluate_wing_sizing_cases!(wing, masses, has_wing_fuel, rclt, rcls, Lhtail, gammat, gammas, We, dyeinn;
        has_strut = false, neinn = 1, neout = 0, dyeout = 0.0, far_load_factor = false,
        Mach_cruise = NaN, V_cruise_EAS = 0.0)
    sl = wing.sizing_loads
    n_pullup = wing_case_load_factor(:pullup, sl.pullup, masses.W_MTO, far_load_factor)
    config = has_strut ? Braced() : Cantilever()
    strut_share = has_strut ? strut_lift_share(wing) : 0.0
    empty!(sl.results)
    sl.governing_strut = :none
    sl.governing_strut_compression = :none

    # Limitation: every case uses the panel weights and inertia stored in the wing (design cruise mass
    # state). Only the lift weight, the load factor and the tail load follow the mass state of the case.
    # Per-case panel inertia, from the structural and fuel weights of each mass state, is not implemented yet.
    for name in fieldnames(WingSizingLoads)
        case = getfield(sl, name)
        case isa WingLoadCase || continue
        case.active || continue

        ms = wing_mass_state(case.mass_state, has_wing_fuel; masses...)
        n_case = wing_case_load_factor(name, case, masses.W_MTO, far_load_factor;
                                       wing, W_lift = ms.W_lift, Mach_cruise, V_cruise_EAS)
        (Nload, po) = wing_case_root_loading(case, wing, rclt, rcls, ms.W_lift, Lhtail, n_pullup;
                                             load_factor = n_case, strut_share)
        case.n_computed = Nload
        if has_strut
            # the spanload of po is the one of wing_loading, with rclt and rcls
            F_strut = strut_lift_on_strut(wing, po, wing.outboard.λ * rclt, wing.inboard.λ * rcls)
            sl.results[name] = wing_case_loads(config, deepcopy(wing), po, gammat, gammas, Nload, We, neout, dyeout, neinn, dyeinn;
                                               strut_lift = F_strut)
        else
            sl.results[name] = wing_case_loads(config, deepcopy(wing), po, gammat, gammas, Nload, We, neout, dyeout, neinn, dyeinn)
        end
    end
    isempty(sl.results) && return nothing

    # first case in declaration order wins ties
    names = [name for name in fieldnames(WingSizingLoads) if haskey(sl.results, name)]
    sl.governing_root = names[argmax([abs(sl.results[name].M_root) for name in names])]
    sl.governing_break = names[argmax([abs(sl.results[name].M_break) for name in names])]

    root, brk = sl.results[sl.governing_root], sl.results[sl.governing_break]
    if has_strut
        strut_forces = [sl.results[name].strut_force for name in names]
        sl.governing_strut = names[argmax(strut_forces)]
        strut_compression, i_compression = findmin(strut_forces)
        strut_compression < 0 && (sl.governing_strut_compression = names[i_compression])
        cases = [sl.results[name] for name in names]
        return (S_root = maximum(max(abs(c.S_root), abs(c.S_break)) for c in cases),
                M_root = maximum(max(abs(c.M_root), abs(c.M_inboard_extreme)) for c in cases),
                S_break = abs(brk.S_break), M_break = abs(brk.M_break),
                T_strut = max(0.0, maximum(strut_forces)),
                T_strut_compression = min(0.0, strut_compression),
                P_axial = cases[argmax([abs(c.P_axial_inboard) for c in cases])].P_axial_inboard)
    end
    return (S_root = abs(root.S_root), M_root = abs(root.M_root),
            S_break = abs(brk.S_break), M_break = abs(brk.M_break))
end # evaluate_wing_sizing_cases!

"""
    wing_mass_state(mass_state, has_wing_fuel; W_cruise, W_MTO, W_fuel, W_reserve, W_empty,
        W_fuel_cap, W_fuel_at_MTO, r_fuel_design)

Returns the lift weight and the fuel split of a wing mass state: the weight the wing lifts, the fuel
carried in the wing as a fraction of the wing tank capacity, and the fuel carried in the fuselage.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mass_state::WingMassState.T`: Mass state (`DesignCruise`, `MTOW_MaxFuel`, `MZFW` or `Light`).
    - `has_wing_fuel::Bool`: Whether fuel is stored in the wing. If `false`, the wing fuel fraction is zero
      and `W_fuel_cap` is not used; if `true`, all fuel is in the wing and the fuselage fuel weight is zero.
    - `W_cruise::Float64`: Aircraft weight at the design cruise point [N].
    - `W_MTO::Float64`: Maximum takeoff weight [N].
    - `W_fuel::Float64`: Design mission fuel weight [N].
    - `W_reserve::Float64`: Reserve fuel weight [N].
    - `W_empty::Float64`: Operating empty weight, without fuel [N].
    - `W_fuel_cap::Float64`: Fuel weight that fills the wing tanks [N].
    - `W_fuel_at_MTO::Float64`: Fuel weight on board at maximum takeoff weight with the maximum payload [N].
    - `r_fuel_design::Float64`: Wing fuel fraction (fuel weight over `W_fuel_cap`) at the design cruise point.

    **Outputs:**
    - `W_lift::Float64`: Weight carried by the wing lift [N].
    - `r_fuel_wing::Float64`: Fuel weight in the wing as a fraction of `W_fuel_cap`; at most one at `MTOW_MaxFuel`.
    - `W_fuel_fuse::Float64`: Fuel weight in the fuselage [N].
"""
function wing_mass_state(mass_state, has_wing_fuel;
        W_cruise, W_MTO, W_fuel, W_reserve, W_empty, W_fuel_cap, W_fuel_at_MTO, r_fuel_design)

    if mass_state == WingMassState.DesignCruise
        W_lift, W_fuel_state = W_cruise, W_fuel
        r_fuel_wing = has_wing_fuel ? r_fuel_design : 0.0
    elseif mass_state == WingMassState.MTOW_MaxFuel
        W_lift, W_fuel_state = W_MTO, W_fuel_at_MTO
        r_fuel_wing = has_wing_fuel ? min(1.0, W_fuel_at_MTO / W_fuel_cap) : 0.0
    elseif mass_state == WingMassState.MZFW
        W_lift, W_fuel_state = W_MTO - W_fuel + W_reserve, W_reserve
        r_fuel_wing = has_wing_fuel ? W_reserve / W_fuel_cap : 0.0
    else # WingMassState.Light
        W_lift, W_fuel_state = W_empty + W_reserve, W_reserve
        r_fuel_wing = has_wing_fuel ? W_reserve / W_fuel_cap : 0.0
    end

    W_fuel_fuse = has_wing_fuel ? 0.0 : W_fuel_state
    return (W_lift = W_lift, r_fuel_wing = r_fuel_wing, W_fuel_fuse = W_fuel_fuse)
end # wing_mass_state

"""
    wing_weights!(wing, po, gammat, gammas, 
       Nload, We, neout, dyeout, neinn, dyeinn, sigfac, rhofuel; n_wings=2.0, has_strut=false, has_jury=false,
       design_wing_fuel=nothing, has_wing_fuel=false, has_centerbox_fuel=false, section_loads_envelope=nothing)

Calculates Wing or Tail loads, stresses, weights of individual wing sections.
Also returns the material gauges, torsional and bending stiffness. Formerly, `get_wing_weights!()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure.
    - `po::Float64`: Wing root loading magnitude.
    - `gammat::Float64`: Wing outer section taper ratio.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `Nload::Int`: Max vertical load factor for wing bending loads.
    - `We::Float64`: Weight of one engine [N].
    - `neout::Int`:  Number of engines per wing mounted outboard of the wing break (strut attachment point).
    - `dyeout::Float64`: Average spanwise moment arm of the outboard engines about the wing break station ηs, i.e. ``y_e - y_s`` [m].
    - `neinn::Int`: Number of engines per wing mounted at or inboard of the wing break (strut attachment point).
    - `dyeinn::Float64`: Average spanwise moment arm of the inboard engines about the wing root station ηo, i.e. ``y_e - y_o`` [m]. For an engine at the break, ``(b/2)(η_s - η_o)``.
    - `sigfac::Float64`: Stress Factor.
    - `rhofuel::Float64`: Density of the fuel.
    - `n_wings::Int64`: Number of total wings (1 for Vtail).
    - `has_strut::Bool`: Whether the wing is braced by a strut (`false` for a cantilever wing).
    - `has_jury::Bool`: Whether the strut has a vertical jury strut (truss-braced wing), only used with `has_strut`. The jury strut is a
        lateral support: it shortens the buckling length of the strut, see `strut_buckling_length`, and is sized by `size_jury!`.
        It does not change the wing loads or the strut force.
    - `design_wing_fuel::Union{Nothing,Float64}`: Fuel weight carried at the design point [N]. If given (not `nothing`), the panel
        weights `wing.inboard.weight`, `wing.outboard.weight` and moments `wing.inboard.dyW`, `wing.outboard.dyW` are updated with the
        panel structure weight (including added weight) plus the wing fuel weight, scaled by the fuel fraction `design_wing_fuel / Wfmax`.
        If `nothing` (default), these panel fields are not modified.
    - `has_wing_fuel::Bool`: Whether fuel is stored in the wing (the fuel fraction is zero otherwise). Only used with `design_wing_fuel`.
    - `has_centerbox_fuel::Bool`: Whether fuel is stored in the wing center box. Only used with `design_wing_fuel`.
    - `section_loads_envelope::Union{Nothing,NamedTuple}`: Envelope of the section loads `(S_root, M_root, S_break, M_break)` over the wing sizing cases.
        If given (not `nothing`), it replaces the computed shear and moment of the inboard (root) and outboard (break) sections before they are sized.

See [Geometry](@ref geometry),  [Wing/Tail Structures](@ref wingtail), and Section 2.7  of the [TASOPT Technical Description](@ref dreladocs). 
"""
function wing_weights!(wing, po, gammat, gammas,
       Nload, We, neout, dyeout, neinn, dyeinn, sigfac, rhofuel; n_wings=2.0, has_strut=false, has_jury=false,
       design_wing_fuel=nothing, has_wing_fuel=false, has_centerbox_fuel=false, section_loads_envelope=nothing)

    tauweb = wing.inboard.webs.material.τmax * sigfac

    cosL = cosd(wing.layout.sweep)
    sinL = sind(wing.layout.sweep)

    # Calculate non-dim span coordinate at span break and root (ηs and ηo resp.)
    etao = wing.layout.ηo
    etas = wing.layout.ηs

    h_avgo, h_rmso = get_average_sparbox_heights(wing.inboard.cross_section)
    h_avgs, h_rmss = get_average_sparbox_heights(wing.outboard.cross_section)

    outboard_section_loads!(wing, po, gammat, gammas, Nload, We, neout, dyeout)
    has_braced_envelope = has_strut && !isnothing(section_loads_envelope) &&
                          wing.strut.opt_strut_force == StrutForceModel.Compatibility
    if (!has_strut || has_braced_envelope) && !isnothing(section_loads_envelope)
        wing.outboard.max_shear_load = section_loads_envelope.S_break
        wing.outboard.moment = section_loads_envelope.M_break
    end

    #---- size strut-attach station at etas
    cs = wing.layout.root_chord*wing.inboard.λ

    outboard_dims = size_wing_section!(wing.outboard, wing.sweep, sigfac)
    # Inboard Section:
    if !has_strut
        #----- no strut, with or without engine at etas
        ls = 0.
        Tstrut = 0.
        Rstrut = 0.
        Pstrut = 0.

        wing_section_loads(Cantilever(), wing, po, gammas, Nload, We, neinn, dyeinn)
        if !isnothing(section_loads_envelope)
            wing.inboard.max_shear_load = section_loads_envelope.S_root
            wing.inboard.moment = section_loads_envelope.M_root
        end

        inboard_dims = size_wing_section!(wing.inboard, wing.sweep, sigfac)

        lsp = 0.
        Tstrutp = 0.
        if wing isa Wing # Tail has no strut
            wing.strut.cos_lambda = 1.0
            wing.strut.dihedral = 0.0
            wing.strut.axial_force_compression = 0.0
            wing.strut.force_vertical = 0.0
            wing.strut.force_horizontal = 0.0
            wing.strut.reaction_share = 0.0
        end

    else
        #----- strut present
        ls = sqrt(wing.strut.z^2 + (0.5*wing.layout.span*(etas-etao))^2)
        Tstrut_compression = 0.
        if has_braced_envelope
            #----- strut tension, and inboard shear, moment and axial force from the sizing load cases
            Tstrut = section_loads_envelope.T_strut
            Tstrut_compression = get(section_loads_envelope, :T_strut_compression, 0.0)
            wing.inboard.max_shear_load = section_loads_envelope.S_root
            wing.inboard.moment = section_loads_envelope.M_root

            #----- size inboard station at etao
            inboard_dims = size_wing_section!(wing.inboard, wing.sweep, sigfac;
                                              P_axial = section_loads_envelope.P_axial)
        else
            F_strut = strut_lift_on_strut(wing, po, gammat, gammas)
            F_wing = strut_wing_end_lift(wing, F_strut)
            Rstrut = (po*wing.layout.span/12.0)*(etas-etao)*(1.0+2.0*gammas) + wing.outboard.max_shear_load + F_wing
            Tstrut = strut_fuselage_end_force(wing, Rstrut, ls, F_strut)
            #c     Pstrut = Rstrut*0.5*wing.layout.span*(etas-etao)/zs

            #----- inboard shear,moment used for sparbox sizing
            wing.inboard.max_shear_load = wing.outboard.max_shear_load
            wing.inboard.moment = wing.outboard.moment
            #
            #----- size inboard station at etao
            inboard_dims = size_wing_section!(wing.inboard, wing.sweep, sigfac)
        end

        #----- total strut length, tension
        lsp = sqrt(wing.strut.z^2 + (0.5*wing.layout.span*(etas-etao)/cosL)^2)
        Tstrutp = Tstrut*lsp/ls
        wing.strut.cos_lambda = ls/lsp
        wing.strut.dihedral = atan(wing.strut.z/(0.5*wing.layout.span*(etas-etao))) # TODO: Update formulation # in-plane strut angle to the horizontal
        wing.strut.axial_force = Tstrutp
        # The compression is scaled by lsp/ls like the tension: it is the force along the total (swept, 3-D)
        # strut length lsp, whereas the strut force of the load cases is in the plane, along the length ls.
        wing.strut.axial_force_compression = Tstrut_compression*lsp/ls
        wing.strut.length = lsp

        # Force of the strut on the fuselage for the fuselage sizing: the largest in-plane force (tension or compression),
        # and the share of the vertical wing support that the strut takes in the case with the largest strut force.
        L_bay = 0.5*wing.layout.span*(etas-etao)
        T_plane = max(Tstrut, abs(Tstrut_compression))
        wing.strut.force_vertical = n_wings*T_plane*wing.strut.z/ls
        wing.strut.force_horizontal = T_plane*L_bay/ls
        results = wing.sizing_loads.results
        wing.strut.reaction_share = 0.0
        if !isempty(results)
            worst = results[argmax(name -> abs(results[name].strut_force), collect(keys(results)))]
            support = worst.S_root + worst.strut_reaction
            iszero(support) || (wing.strut.reaction_share = clamp(worst.strut_reaction/support, 0.0, 1.0))
        end
        # The strut is a box of caps and webs of the chord `chord_to_wing_chord` times the wing chord at the strut attachment.
        # With a tapered strut, this wing-end chord is the smaller chord, so the box is sized for the thinnest section.
        # A vertical jury strut shortens the buckling length of the strut and loads the webs in shear (previous iteration).
        strut = wing.strut
        strut.chord = strut.chord_to_wing_chord * cs
        strut.root_chord = strut.chord / strut.taper_ratio
        buckling_length = strut_buckling_length(lsp, has_jury ? wing.jury.eta_strut : nothing)

        # The load factor is the largest of the active sizing cases, or Nload without the load case envelope.
        n_max = Nload
        if has_braced_envelope
            sl = wing.sizing_loads
            n_max = maximum(abs(getfield(sl, name).n_computed) for name in fieldnames(WingSizingLoads)
                            if getfield(sl, name) isa WingLoadCase && getfield(sl, name).active; init = 0.0)
        end
        W_strut_one_wing = strut.weight / n_wings
        perpendicular_load(F) = abs(F) / lsp * sqrt(lsp^2 - strut.z^2) / lsp
        if any(r.strut_lift != 0 for r in values(results))
            sl = wing.sizing_loads
            q_perp = maximum(perpendicular_load(r.strut_lift - getfield(sl, name).n_computed * W_strut_one_wing)
                             for (name, r) in results)
        else
            q_perp = perpendicular_load(n_max * W_strut_one_wing)
        end
        V_transverse = q_perp * buckling_length / 2
        strut.bending_moment = q_perp * buckling_length^2 / 8
        strut.shear_load = (has_jury ? abs(wing.jury.axial_force) : 0.0) + V_transverse

        box = size_member_box(strut.caps, strut.webs;
            chord = strut.chord, height_to_chord = strut.thickness_to_chord, width_to_chord = strut.box_width_to_chord,
            web_to_box_height = strut.web_to_box_height, tension = Tstrutp, compression = strut.axial_force_compression,
            buckling_length = buckling_length, K = strut.end_fixity_K,
            shear_load = strut.shear_load, bending_moment = strut.bending_moment, sigfac = sigfac,
            min_cap_thickness = strut.min_cap_thickness, min_web_thickness = strut.min_web_thickness)
        strut.area = box.area
        strut.EI_min = box.EI_min
        strut.buckling_load = box.buckling_load
        strut.sizing_criterion = box.criterion
        # The strut weight is the bare member weight of both wings times the added weight fraction (fairings and fittings);
        # the area stays bare, and the weight moment uses the total weight.
        # The box area scales with the chord, so the weight is that of the mean chord over the wing-end chord of the box.
        mean_to_wing_end_chord = (1.0 + 1.0/strut.taper_ratio) / 2
        strut.weight = n_wings*gee*lsp*(strut.caps.material.ρ*box.area_caps + strut.webs.material.ρ*box.area_webs) *
                       mean_to_wing_end_chord * (1.0 + strut.weight_frac_added)
        wing.strut.dxW = wing.strut.weight * 0.25*wing.layout.span*(etas-etao) * sinL/cosL
        has_jury && size_jury!(wing, etao, etas, sigfac, sinL, cosL; n_wings = n_wings) 
    end

    Abfuels = calc_sparbox_internal_area(wing.inboard.cross_section.width_to_chord,
                                    h_avgs, outboard_dims.tbcap, outboard_dims.tbweb)
    Abfuelo = calc_sparbox_internal_area(wing.inboard.cross_section.width_to_chord,
                                    h_avgo, inboard_dims.tbcap, inboard_dims.tbweb)
    
    

    Vcen = wing.layout.root_chord^2*wing.layout.span*  etao / 2.0

    Vinn = wing.layout.root_chord^2*wing.layout.span* (etas-etao) *
        (1.0 + wing.inboard.λ + wing.inboard.λ^2)/6.0 *
        cosL
    Vout = wing.layout.root_chord^2*wing.layout.span* (1.0 -etas) *
        (wing.inboard.λ^2 + wing.inboard.λ*wing.outboard.λ + wing.outboard.λ^2)/6.0 *
        cosL

    dxVinn = wing.layout.root_chord^2*wing.layout.span^2 * (etas-etao)^2 *
        (1.0 + 2.0*wing.inboard.λ + 3.0*wing.inboard.λ^2)/48.0 *
        sinL
    dxVout = wing.layout.root_chord^2*wing.layout.span^2 * (1.0 -etas)^2 *
        (wing.inboard.λ^2 + 2.0*wing.inboard.λ*wing.outboard.λ + 3.0*wing.outboard.λ^2)/48.0 *
        sinL +
        wing.layout.root_chord^2*wing.layout.span^2 * (etas-etao)*(1.0 -etas) *
        (wing.inboard.λ^2 + wing.inboard.λ*wing.outboard.λ + wing.outboard.λ^2)/12.0 *
        sinL

    dyVinn = wing.layout.root_chord^2*wing.layout.span^2 * (etas-etao)^2 *
        (1.0 + 2.0*wing.inboard.λ + 3.0*wing.inboard.λ^2)/48.0 *
        cosL
    dyVout = wing.layout.root_chord^2*wing.layout.span^2 * (1.0 -etas)^2 *
        (wing.inboard.λ^2 + 2.0*wing.inboard.λ*wing.outboard.λ + 3.0*wing.outboard.λ^2)/48.0 *
        cosL

    #---- set chord^2 weighted average areas for inner panel
    Abcapi = (inboard_dims.Abcap + outboard_dims.Abcap*wing.inboard.λ^2)/(1.0+wing.inboard.λ^2)
    Abwebi = (inboard_dims.Abweb + outboard_dims.Abweb*wing.inboard.λ^2)/(1.0+wing.inboard.λ^2)

    Wscen   = (wing.inboard.caps.ρ*inboard_dims.Abcap  + wing.inboard.webs.ρ*inboard_dims.Abweb )*gee*Vcen
    Wsinn   = (wing.inboard.caps.ρ*Abcapi              + wing.inboard.webs.ρ*Abwebi             )*gee*Vinn
    Wsout   = (wing.inboard.caps.ρ*outboard_dims.Abcap + wing.inboard.webs.ρ*outboard_dims.Abweb)*gee*Vout

    dxWsinn = (wing.inboard.caps.ρ*Abcapi              + wing.inboard.webs.ρ*Abwebi             )*gee*dxVinn
    dxWsout = (wing.inboard.caps.ρ*outboard_dims.Abcap + wing.inboard.webs.ρ*outboard_dims.Abweb)*gee*dxVout

    dyWsinn = (wing.inboard.caps.ρ*Abcapi              + wing.inboard.webs.ρ*Abwebi             )*gee*dyVinn
    dyWsout = (wing.inboard.caps.ρ*outboard_dims.Abcap + wing.inboard.webs.ρ*outboard_dims.Abweb)*gee*dyVout


    Abfueli= (Abfuelo + Abfuels*wing.inboard.λ^2)/(1.0+wing.inboard.λ^2)

    Wfcen   = rhofuel*Abfuelo *gee*Vcen
    Wfinn   = rhofuel*Abfueli *gee*Vinn
    Wfout   = rhofuel*Abfuels *gee*Vout

    dxWfinn = rhofuel*Abfueli *gee*dxVinn
    dxWfout = rhofuel*Abfuels *gee*dxVout

    dyWfinn = rhofuel*Abfueli *gee*dyVinn
    dyWfout = rhofuel*Abfuels *gee*dyVout

    wing.inboard.caps.weight.W   = 2.0*wing.inboard.caps.ρ*gee*(inboard_dims.Abcap*Vcen + Abcapi*Vinn + outboard_dims.Abcap*Vout)
    wing.inboard.webs.weight.W   = 2.0*wing.inboard.webs.ρ*gee*(inboard_dims.Abweb*Vcen + Abwebi*Vinn + outboard_dims.Abweb*Vout)

    dxWcap = 2.0*wing.inboard.caps.ρ*gee*( Abcapi*dxVinn + outboard_dims.Abcap*dxVout )
    dxWweb = 2.0*wing.inboard.webs.ρ*gee*( Abwebi*dxVinn + outboard_dims.Abweb*dxVout )

    Vout = wing.layout.root_chord^2*wing.layout.span* (1.0 -etas) *
    (wing.inboard.λ^2 + wing.inboard.λ*wing.outboard.λ + wing.outboard.λ^2)/6.0 *
    cosL

    wing.inboard.caps.thickness = inboard_dims.tbcap
    wing.inboard.webs.thickness = inboard_dims.tbweb
    wing.outboard.caps.thickness = outboard_dims.tbcap
    wing.outboard.webs.thickness = outboard_dims.tbweb

    fwadd = wing_additional_weight(wing)
    Wwing = n_wings * (Wscen + Wsinn + Wsout) * (1.0 + fwadd)
    wing.dxW = n_wings * (dxWsinn + dxWsout) * (1.0 + fwadd)

    if !isnothing(design_wing_fuel)
        Wfmax = has_wing_fuel ? 2.0 * ((has_centerbox_fuel ? Wfcen : 0.0) + Wfinn + Wfout) : 0.0
        rfmax = has_wing_fuel ? design_wing_fuel / Wfmax : 0.0

        wing.inboard.weight = Wsinn * (1.0 + fwadd) + rfmax * Wfinn
        wing.outboard.weight = Wsout * (1.0 + fwadd) + rfmax * Wfout
        wing.inboard.dyW = dyWsinn * (1.0 + fwadd) + rfmax * dyWfinn
        wing.outboard.dyW = dyWsout * (1.0 + fwadd) + rfmax * dyWfout
    end

    return Wwing, Wsinn, Wsout, dyWsinn, dyWsout, Wfcen, Wfinn, Wfout,
    dxWfinn, dxWfout, dyWfinn, dyWfout, lsp

end # wing_weights


"""
    size_cap(σmax, moment, h̄, w̄, h_rms, c, cosΛ; t_min = 0.0)

Calculates area and thickness of wing caps based on maximum stress and geometric parameters.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `σmax::Float64`: Maximum allowable stress
    - `moment::Float64`: Bending moment
    - `h̄::Float64`: Height-to-chord ratio
    - `w̄::Float64`: Width-to-chord ratio
    - `h_rms::Float64`: RMS height of sparbox
    - `c::Float64`: Chord length
    - `cosΛ::Float64`: Cosine of sweep angle
    - `t_min::Float64`: (keyword, default 0.0) Minimum cap thickness; `t_cap = max(required, t_min)`
    
    **Outputs:**
    - `t_cap::Float64`: Thickness of cap
    - `Ab_cap::Float64`: Cross-sectional area of cap
"""
function size_cap(σmax, moment, h̄, w̄, h_rms, c, cosΛ; t_min = 0.0)
    t_cap = max(calc_cap_thickness(σmax, moment, h̄, w̄, h_rms, c, cosΛ), t_min)
    Ab_cap = 2 * t_cap * w̄
    return t_cap, Ab_cap
end  # function size_cap

"""
    size_web(τmax, shear, c_perp, web_height; t_min = 0.0)

Calculates area and thickness of wing webs based on maximum shear stress and geometric parameters.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `τmax::Float64`: Maximum allowable shear stress
    - `shear::Float64`: Shear load
    - `c_perp::Float64`: Perpendicular chord length
    - `web_height::Float64`: Height of web
    - `t_min::Float64`: (keyword, default 0.0) Minimum web thickness; `t_web = max(required, t_min)`
    
    **Outputs:**
    - `t_web::Float64`: Thickness of web
    - `Ab_web::Float64`: Cross-sectional area of web
"""
function size_web(τmax, shear, c_perp, web_height; t_min = 0.0)
    t_web = max(calc_web_thickness(τmax, shear, c_perp, web_height), t_min)
    Ab_web = 2 * t_web * web_height
    return t_web, Ab_web
end  # function size_web

"""
    calc_web_thickness(τmax, shear, c_perp, web_height)

Calculates the required web thickness based on maximum shear stress and loading conditions.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `τmax::Float64`: Maximum allowable shear stress
    - `shear::Float64`: Shear load
    - `c_perp::Float64`: Perpendicular chord length
    - `web_height::Float64`: Height of web
    
    **Outputs:**
    - `t_web::Float64`: Required thickness of web
"""
function calc_web_thickness(τmax, shear, c_perp, web_height)
    t_web = shear / (c_perp^2 * 2 * web_height * τmax)
    return t_web
end  # function calc_web_thickness

"""
    calc_cap_thickness(σmax, moment, h̄, w̄, h_rms, c, cosΛ)

Calculates the required cap thickness based on maximum stress and geometric parameters.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `σmax::Float64`: Maximum allowable stress
    - `moment::Float64`: Bending moment
    - `h̄::Float64`: Height-to-chord ratio
    - `w̄::Float64`: Width-to-chord ratio
    - `h_rms::Float64`: RMS height of sparbox
    - `c::Float64`: Chord length
    - `cosΛ::Float64`: Cosine of sweep angle
    
    **Outputs:**
    - `t_cap::Float64`: Required thickness of cap

Throws an `ErrorException` if the cap sizing is infeasible (`con >= h_rms^3`).
"""
function calc_cap_thickness(σmax, moment, h̄, w̄, h_rms, c, cosΛ)
    con = moment * 6h̄ / w̄ * 1 / (c^3 * σmax * cosΛ^4)
    if con >= h_rms^3
        error("Wing cap sizing infeasible: bending load parameter con = $con >= h_rms^3 = $(h_rms^3). " *
              "The cap cannot carry the moment within the sparbox height.")
    end
    t_cap = 0.5 * (h_rms - ∛(h_rms^3 - con))
    return t_cap
end  # function calc_cap_thickness

"""
    calc_sparbox_internal_area(width, height, t_cap, t_web)

Calculates the internal aera of the sparbox, accounting for the thickness of
the webs and caps.
A = (w - 2tweb)×(h - 2tcap)
"""
function calc_sparbox_internal_area(width, height, t_cap, t_web)
    return (width - 2*t_web)*(height - 2*t_cap)
end  # function calc_internal_area