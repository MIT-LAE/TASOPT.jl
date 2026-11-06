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
    :strut_attachment_outboard_limit => "strut attached outboard of the limit",
    :wing_beam_column_margin => "inboard wing axial load too close to its Euler load (placeholder metric)",
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

Placeholder limit of the inboard wing axial force to Euler load ratio, see [`check_wing_beam_column!`](@ref).
"""
const WING_BEAM_COLUMN_RATIO_MAX = 0.5

"""
    check_wing_beam_column!(ac::aircraft)

Records `:wing_beam_column_margin` if the axial force in the inboard wing of a braced wing exceeds
[`WING_BEAM_COLUMN_RATIO_MAX`](@ref) times the Euler load ``\\pi^2 EI / (4 L^2)`` of a clamped-free column.
Must be called after the wing has been sized.

!!! warning "Placeholder metric"
    Uses the largest axial force of the sizing load cases, the root stiffness `inboard.EI[1]` and the inboard panel length `(ηs - ηo) b / 2`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object with a sized wing.

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function check_wing_beam_column!(ac::aircraft)
    wing = ac.wing
    ac.options.has_strut || return ac

    EI = wing.inboard.EI[1]
    EI > 0.0 || return ac

    L = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
    P_axial = maximum(abs(case.P_axial_inboard) for case in values(wing.sizing_loads.results); init = 0.0)
    if abs(P_axial) / (π^2 * EI / (4 * L^2)) > WING_BEAM_COLUMN_RATIO_MAX # axial force over Euler load
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
