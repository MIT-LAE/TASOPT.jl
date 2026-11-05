"""
    record_infeasibility!(ac::aircraft, flag::Symbol)

Records that the current design violates a constraint, identified by `flag`, without
interrupting the sizing. Each flag is stored once, in order of first occurrence.
The flags are reported by [`report_infeasibilities`](@ref) at the end of sizing.

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
    clear_infeasibilities!(ac::aircraft)

Removes all recorded infeasibility flags from `ac`. Called in each iteration of the weight sizing loop
that sizes the wing, after the convergence test, so the flags describe the final sizing pass only.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object whose `infeasibilities` are emptied.

    **Outputs:**
    - `ac::aircraft`: the same aircraft object.
"""
function clear_infeasibilities!(ac::aircraft)
    empty!(ac.infeasibilities)
    return ac
end

"""
    INFEASIBILITY_DESCRIPTIONS

Short description of each infeasibility flag recorded by the sizing, shown by [`report_infeasibilities`](@ref).
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

Issues a single warning naming every infeasibility flag recorded in `ac`.
Does nothing if no flags are recorded. The flags are left in place.

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

Records the flag `:wing_span_exceeded` if the wing span `wing.layout.span` strictly exceeds the span limit
`wing.layout.max_span` (`maxSpan` in the input TOML). A limit of zero or less means no limit, so the wing is never flagged.
Must be called after the wing geometry has been set.

The span limit is normally enforced by the optimiser outside `size_aircraft!` (see `example/example_opt.jl`).
This check only reports a violation, it does not change the sizing.

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

Records the flag `:wing_cap_min_gauge` if the minimum gauge `min_cap_thickness` sets the cap thickness of the
inboard or outboard wing section, i.e. the strength-required cap thickness is thinner than the gauge.
Must be called after the wing has been sized. The section floor is the same normalized thickness that
[`size_wing_section!`](@ref) uses, `min_cap_thickness / (co cos(Λ))`. A section without a minimum gauge is never flagged.

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

Placeholder limit of the spanwise strut attachment station `ηs` (fraction of the half-span) of a braced wing.
Attaching the strut further outboard is flagged by [`check_strut_attachment!`](@ref).
"""
const STRUT_ETA_S_MAX = 0.9

"""
    check_strut_attachment!(ac::aircraft)

Records the flag `:strut_attachment_outboard_limit` if the wing is braced by a strut and its attachment
station `wing.layout.ηs` is outboard of [`STRUT_ETA_S_MAX`](@ref). A cantilever wing is never flagged.

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

Placeholder limit of the ratio of the axial force in the inboard wing of a braced wing to its cantilever Euler
load, see [`wing_beam_column_ratio`](@ref).
"""
const WING_BEAM_COLUMN_RATIO_MAX = 0.5

"""
    wing_beam_column_ratio(P_axial, EI, L)

Ratio of an axial force to the Euler buckling load of a column of bending stiffness `EI` and length `L`
that is clamped at one end and free at the other, ``|P| / (\\pi^2 EI / (4 L^2))``.

!!! warning "Placeholder metric"
    This is a temporary, simple estimate, not a validated beam-column model. It treats the inboard wing panel as a
    prismatic cantilever column and ignores the strut support, the bending moment and the taper.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `P_axial::Float64`: Axial force [N]; tension or compression give the same ratio.
    - `EI::Float64`: Bending stiffness [N m²].
    - `L::Float64`: Column length [m].

    **Outputs:**
    - `ratio::Float64`: Axial force over Euler load.
"""
wing_beam_column_ratio(P_axial, EI, L) = abs(P_axial) / (π^2 * EI / (4 * L^2))

"""
    check_wing_beam_column!(ac::aircraft)

Records the flag `:wing_beam_column_margin` if the axial force in the inboard wing of a braced wing exceeds
[`WING_BEAM_COLUMN_RATIO_MAX`](@ref) times its Euler load from [`wing_beam_column_ratio`](@ref). Must be called after
the wing has been sized. A cantilever wing is never flagged, nor is a wing whose inboard stiffness is not positive.

!!! warning "Placeholder metric"
    This check is temporary, to be refined. It uses the largest axial force magnitude of the wing sizing load cases,
    the flapwise stiffness `inboard.EI[1]` of the inboard section (at the wing root station `ηo`, where it is
    largest) and the length of the inboard panel `(ηs - ηo) b / 2` between the root station and the strut attachment.

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
    if wing_beam_column_ratio(P_axial, EI, L) > WING_BEAM_COLUMN_RATIO_MAX
        record_infeasibility!(ac, :wing_beam_column_margin)
    end
    return ac
end

"""
    check_wing_fuel_volume!(ac::aircraft, W_fuel, W_fuel_capacity)

Records the flag `:wing_fuel_volume_exceeded` if fuel is stored in the wing (`options.has_wing_fuel`) and the
fuel weight `W_fuel` strictly exceeds the weight `W_fuel_capacity` that the geometric wing tanks
(plus the centre-box tank if `options.has_centerbox_fuel`) can hold. The limit is 100% of the geometric capacity,
without margin. Does nothing if there is no fuel in the wing.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft object.
    - `W_fuel::Float64`: Weight of the design mission fuel to be stored in the wing [N]. This is not the fuel at
      maximum payload that is used for the structural sizing.
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
