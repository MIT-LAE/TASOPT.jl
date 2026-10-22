"""
    OperatingPointMaximumReport

Values underlying a maximum-over-operating-points constraint. `values` contains
one value for every design-mission point in `point_indices`; `maximum` is the
aggregate value used by the compact constraint representation. A non-`nothing`
`fixed_sizing_point` identifies the point at which a fixed-limit sizing model
is applied rather than merely evaluated.
"""
struct OperatingPointMaximumReport
    name::Symbol
    point_indices::Vector{Int}
    point_labels::Vector{Symbol}
    values::Vector{Float64}
    limit::Float64
    maximum::Float64
    maximizing_point::Int
    fixed_sizing_point::Union{Nothing, Int}

    function OperatingPointMaximumReport(name::Symbol, point_indices,
        point_labels, values, limit::Real; fixed_sizing_point=nothing)
        indices = Int.(collect(point_indices))
        labels = Symbol.(collect(point_labels))
        point_values = Float64.(collect(values))
        length(indices) == length(labels) == length(point_values) ||
            throw(DimensionMismatch("operating-point indices, labels, and values must align"))
        isempty(indices) && throw(ArgumentError("operating-point report cannot be empty"))
        isfinite(limit) && limit > 0 ||
            throw(ArgumentError("operating-point limit must be finite and positive"))
        maximum_value, local_index = findmax(point_values)
        fixed = isnothing(fixed_sizing_point) ? nothing : Int(fixed_sizing_point)
        isnothing(fixed) || fixed in indices || throw(ArgumentError(
            "fixed sizing point must be present in the operating-point report"))
        return new(name, indices, labels, point_values, Float64(limit),
            maximum_value, indices[local_index], fixed)
    end
end

"""Raw physical residuals, scaled residuals, physical values, limits, and operating-point details."""
struct ConstraintReport{R,S,V,L,O}
    raw::R
    scaled::S
    values::V
    limits::L
    operating_points::O
end

ConstraintReport(raw, scaled, values, limits) =
    ConstraintReport(raw, scaled, values, limits, NamedTuple())

"""Immutable context supplied to every per-mission constraint policy."""
struct MissionConstraintContext
    mission_index::Int
    design_mtow::Float64
    is_off_design::Bool

    function MissionConstraintContext(mission_index::Integer,
        design_mtow::Real, is_off_design::Bool)
        mission_index > 0 || throw(ArgumentError(
            "mission index must be a positive integer"))
        isfinite(design_mtow) && design_mtow > 0 || throw(ArgumentError(
            "design MTOW must be finite and positive"))
        return new(Int(mission_index), Float64(design_mtow), is_off_design)
    end
end

function _operating_point_maximum(name, point_values, limit;
    fixed_sizing_point=nothing)
    points = collect(ipstatic:ipdescentn)
    labels = Symbol.(ip_labels[points])
    return OperatingPointMaximumReport(name, points, labels, point_values,
        limit; fixed_sizing_point)
end

"""
    aircraft_constraint_report(ac)

Report aircraft constraints that depend on the shared aircraft design rather than
on an individual mission: wing span, fan diameter, and fan/core-LPC pressure
ratio ordering. These constraints are intended to be evaluated once after
sizing in a multi-mission analysis.
"""
function aircraft_constraint_report(ac)
    values = (wing_span=ac.wing.layout.span,
        fan_diameter=ac.parg[igdfan],
        fan_pressure_ratio=ac.pare[iepif, ipcruise1, 1])
    limits = (wing_span=ac.wing.layout.max_span,
        fan_diameter=ac.engine.max_fan_diameter,
        fan_pressure_ratio=ac.pare[iepilc, ipcruise1, 1])
    raw = (wing_span=values.wing_span-limits.wing_span,
        fan_diameter=values.fan_diameter-limits.fan_diameter,
        fan_to_lpc_pressure_ratio=
            values.fan_pressure_ratio-limits.fan_pressure_ratio)
    scaled = (wing_span=values.wing_span/limits.wing_span-1.0,
        fan_diameter=values.fan_diameter/limits.fan_diameter-1.0,
        fan_to_lpc_pressure_ratio=
            values.fan_pressure_ratio/limits.fan_pressure_ratio-1.0)
    return ConstraintReport(raw, scaled, values, limits)
end

"""
    mission_constraint_report(ac, mission_index, context)

Report the aircraft constraints that must hold for every mission: fuel capacity,
balanced field length, top-of-climb gradient, compressor-exit temperature,
and turbine-metal temperature. The takeoff weight row compares every mission
against the MTOW fixed by design-mission sizing; it is an identity for mission
1 and the required fixed-MTOW constraint for off-design missions.
"""
function mission_constraint_report(ac, mission_index::Integer,
    context)
    1 <= mission_index <= size(ac.parm, 2) || throw(ArgumentError(
        "mission index is outside the aircraft mission set"))
    points = ipstatic:ipdescentn
    t3_report = _operating_point_maximum(:compressor_exit_temperature,
        @view(ac.pare[ieTt3, points, mission_index]), ac.engine.max_compressor_exit_temperature)
    metal_by_point = vec(maximum(@view(ac.pare[
        ieTmet1:ieTmet1+ncrowx-1, points, mission_index]); dims=1))
    fixed_point = mission_index == 1 ? iprotate : nothing
    tmetal_report = _operating_point_maximum(:turbine_metal_temperature,
        metal_by_point, ac.parg[igTmetal]; fixed_sizing_point=fixed_point)
    design_mtow = Float64(context.design_mtow)
    isfinite(design_mtow) && design_mtow > 0 || throw(ArgumentError(
        "design MTOW must be finite and positive"))
    values = (fuel_weight=ac.parm[imWfuel, mission_index],
        balanced_field_length=ac.parm[imlBF, mission_index],
        top_of_climb_gradient=ac.para[iagamV, ipclimbn, mission_index],
        compressor_exit_temperature=t3_report.maximum,
        turbine_metal_temperature=tmetal_report.maximum,
        takeoff_weight=ac.parm[imWTO, mission_index])
    limits = (fuel_weight=ac.parg[igWfmax],
        balanced_field_length=ac.parg[iglBFmax],
        top_of_climb_gradient=ac.parg[iggtocmin],
        compressor_exit_temperature=ac.engine.max_compressor_exit_temperature,
        turbine_metal_temperature=ac.parg[igTmetal],
        takeoff_weight=design_mtow)
    raw = (fuel_capacity=values.fuel_weight-limits.fuel_weight,
        balanced_field_length=
            values.balanced_field_length-limits.balanced_field_length,
        top_of_climb_gradient=
            limits.top_of_climb_gradient-values.top_of_climb_gradient,
        compressor_exit_temperature=
            values.compressor_exit_temperature-limits.compressor_exit_temperature,
        turbine_metal_temperature=
            values.turbine_metal_temperature-limits.turbine_metal_temperature,
        design_mtow=values.takeoff_weight-limits.takeoff_weight)
    scaled = (fuel_capacity=values.fuel_weight/limits.fuel_weight-1.0,
        balanced_field_length=
            values.balanced_field_length/limits.balanced_field_length-1.0,
        top_of_climb_gradient=
            1.0-values.top_of_climb_gradient/limits.top_of_climb_gradient,
        compressor_exit_temperature=
            values.compressor_exit_temperature/limits.compressor_exit_temperature-1.0,
        turbine_metal_temperature=
            values.turbine_metal_temperature/limits.turbine_metal_temperature-1.0,
        design_mtow=values.takeoff_weight/limits.takeoff_weight-1.0)
    operating_points = (compressor_exit_temperature=t3_report,
        turbine_metal_temperature=tmetal_report)
    return ConstraintReport(raw, scaled, values, limits, operating_points)
end

"""Return whether every scaled constraint residual is within `tolerance`."""
function is_feasible(report::ConstraintReport, tolerance::Real=0.0)
    isfinite(tolerance) && tolerance >= 0 ||
        throw(ArgumentError("feasibility tolerance must be finite and nonnegative"))
    residuals = report.scaled
    return all(value -> isfinite(value) && value <= tolerance, values(residuals))
end
