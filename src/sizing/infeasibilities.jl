"""
    record_infeasibility!(ac::aircraft, flag::Symbol)

Records that the design violates the constraint `flag`, without interrupting the sizing.
Each flag is stored once. Flags are reported by [`report_infeasibilities`](@ref).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object whose `infeasibilities` are updated.
    - `flag::Symbol`: name of the violated constraint.

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function record_infeasibility!(ac::aircraft, flag::Symbol)
    flag in ac.infeasibilities || push!(ac.infeasibilities, flag)
    return ac
end

"""
    INFEASIBILITY_DESCRIPTIONS

Description of each infeasibility flag, shown by [`report_infeasibilities`](@ref).
"""
const INFEASIBILITY_DESCRIPTIONS = Dict{Symbol,String}(
    :wing_cap_min_gauge => "minimum gauge sets the wing cap thickness",
    :airfoil_database_exceeded => "wing outside the airfoil database, drag includes the extrapolation penalty",
    :strut_airfoil_database_exceeded => "strut or jury strut outside the airfoil database",
    :strut_attachment_outboard_limit => "strut attached outboard of the limit",
    :wing_beam_column_margin => "inboard wing compression too close to its Euler load (placeholder metric)",
    :wing_fuel_volume_exceeded => "fuel weight exceeds the capacity of the wing tanks",
    :wing_span_exceeded => "wing span exceeds the span limit",
)

"""
    report_infeasibilities(ac::aircraft)

Issues a single warning naming every infeasibility flag recorded in `ac`, if any.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object whose `infeasibilities` are reported.

    **Outputs:**
    - `nothing`
"""
function report_infeasibilities(ac::aircraft)
    isempty(ac.infeasibilities) && return nothing
    describe(flag) = haskey(INFEASIBILITY_DESCRIPTIONS, flag) ? "$flag ($(INFEASIBILITY_DESCRIPTIONS[flag]))" : string(flag)
    @warn "Infeasible design: constraints violated in the final sizing iteration: " *
          join(describe.(ac.infeasibilities), ", ")
    return nothing
end

"""
    check_wing_span!(ac::aircraft)

Records `:wing_span_exceeded` if `wing.layout.span` exceeds `wing.layout.max_span` (zero or less means no limit).
Must be called after the wing geometry has been set. It only reports; it does not change the sizing.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object with the wing geometry set.

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function check_wing_span!(ac::aircraft)
    layout = ac.wing.layout
    if layout.max_span > 0.0 && layout.span > layout.max_span
        record_infeasibility!(ac, :wing_span_exceeded)
    end
    return ac
end

"""
    check_wing_cap_min_gauge!(ac::aircraft)

Records `:wing_cap_min_gauge` if the minimum gauge `min_cap_thickness` sets the cap thickness of the
inboard or outboard wing section. Must be called after the wing has been sized.
The floor is `min_cap_thickness / (co cos(Λ))`, as in [`size_wing_section!`](@ref).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object with a sized wing.

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function check_wing_cap_min_gauge!(ac::aircraft)
    wing = ac.wing
    for section in (wing.inboard, wing.outboard)
        section.min_cap_thickness > 0.0 || continue
        t_cap_floor = section.min_cap_thickness / (section.co * cosd(wing.sweep))
        section.caps.thickness <= t_cap_floor && record_infeasibility!(ac, :wing_cap_min_gauge)
    end
    return ac
end

"""
    STRUT_ETA_S_MAX

Placeholder limit of the strut attachment station `ηs`, flagged by [`check_strut_attachment!`](@ref).
"""
const STRUT_ETA_S_MAX = 0.9

"""
    check_strut_attachment!(ac::aircraft)

Records `:strut_attachment_outboard_limit` if the wing has a strut and `wing.layout.ηs` exceeds [`STRUT_ETA_S_MAX`](@ref).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object with the wing layout set.

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function check_strut_attachment!(ac::aircraft)
    if ac.options.has_strut && ac.wing.layout.ηs > STRUT_ETA_S_MAX
        record_infeasibility!(ac, :strut_attachment_outboard_limit)
    end
    return ac
end

"""
    WING_BEAM_COLUMN_RATIO_MAX

Placeholder limit of the inboard wing compression to Euler load ratio, see [`check_wing_beam_column!`](@ref).
"""
const WING_BEAM_COLUMN_RATIO_MAX = 0.5

"""
    WING_INBOARD_PANEL_EFFECTIVE_LENGTH_FACTOR

Effective-length factor ``K`` of the inboard wing panel (clamped at the fuselage, pinned at the strut attachment),
see [`wing_beam_column_ratio`](@ref). Ignores the outboard panel restraint, which would lower ``K`` toward 0.5.
"""
const WING_INBOARD_PANEL_EFFECTIVE_LENGTH_FACTOR = 0.7

"""
    wing_beam_column_ratio(wing::Wing)

Returns ``P_c / P_E``, the largest inboard wing compression over the sizing load cases over the Euler load
``P_E = \\pi^2 EI / (K L)^2``, with ``K =`` [`WING_INBOARD_PANEL_EFFECTIVE_LENGTH_FACTOR`](@ref) and
``L = (\\eta_s - \\eta_o) b / 2``. ``EI`` is the mean of the flapwise stiffnesses at the root and at the strut
attachment. Returns 0 if either stiffness is not sized.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::Wing`: sized wing.

    **Outputs:**
    - `ratio::Float64`: compression over Euler load [-].
"""
function wing_beam_column_ratio(wing::Wing)
    EI_root = wing.inboard.EI[1]
    EI_strut = wing.outboard.EI[1]
    (EI_root > 0.0 && EI_strut > 0.0) || return 0.0
    EI = 0.5 * (EI_root + EI_strut)

    L = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
    P_c = maximum(max(-case.P_axial_inboard, 0.0) for case in values(wing.sizing_loads.results); init = 0.0) # tension positive
    P_E = π^2 * EI / (WING_INBOARD_PANEL_EFFECTIVE_LENGTH_FACTOR * L)^2
    return P_c / P_E
end

"""
    check_wing_beam_column!(ac::aircraft)

Records `:wing_beam_column_margin` if the compression in the inboard wing of a braced wing exceeds
[`WING_BEAM_COLUMN_RATIO_MAX`](@ref) times the Euler load, see [`wing_beam_column_ratio`](@ref).
Must be called after the wing has been sized.

!!! warning "Placeholder metric"
    Chordwise (in-plane) buckling and the beam-column moment amplification in the box sizing are not included.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object with a sized wing.

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function check_wing_beam_column!(ac::aircraft)
    ac.options.has_strut || return ac

    if wing_beam_column_ratio(ac.wing) > WING_BEAM_COLUMN_RATIO_MAX
        record_infeasibility!(ac, :wing_beam_column_margin)
    end
    return ac
end

"""
    check_wing_fuel_volume!(ac::aircraft, W_fuel, W_fuel_capacity)

Records `:wing_fuel_volume_exceeded` if the wing stores fuel (`options.has_wing_fuel`) and `W_fuel` exceeds the
capacity `W_fuel_capacity` of the wing tanks (and centre-box tank if `options.has_centerbox_fuel`).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object.
    - `W_fuel::Float64`: Design mission fuel to be stored in the wing [N].
    - `W_fuel_capacity::Float64`: Weight of fuel that fills the wing tanks [N].

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function check_wing_fuel_volume!(ac::aircraft, W_fuel, W_fuel_capacity)
    if ac.options.has_wing_fuel && W_fuel > W_fuel_capacity
        record_infeasibility!(ac, :wing_fuel_volume_exceeded)
    end
    return ac
end
