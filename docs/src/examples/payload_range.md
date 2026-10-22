# Payload-range diagram

After sizing a feasible aircraft, plot its payload-range capability:

```julia
using TASOPT
ac = read_aircraft_model("aircraft.toml")
result = size_aircraft!(ac)
@assert result.success
TASOPT.PayloadRange(ac; filename="payload_range.png")
```

The single-panel diagram connects four characteristic points:

0. Maximum payload at zero plotted range, verified with a 500 nmi mission.
1. Maximum feasible range at maximum payload.
2. Maximum feasible range at the nominal full-fuel payload: MTOW minus empty
   weight minus fuel capacity.
3. Maximum feasible range at zero payload (ferry range).

The original design mission is marked separately with a star, colored orange and
labeled when outside the constraint tolerance. Point 0 is an axis intercept, not
a zero-distance cruise mission.
If the full-fuel payload is outside the interval from zero to maximum payload,
it is clamped to that interval and the coincident corner is omitted.

## Constraints and altitude selection

Every range trial searches starting cruise altitudes **28,000:1,000:50,000 ft**
and selects the lowest fuel burn **among feasible, converged candidates**.
Design Mach, reserves, and other mission settings are preserved. TASOPT retains
its cruise-climb profile; altitude bounds apply to the start of cruise.

The default policies reuse `aircraft_constraint_report` and
`mission_constraint_report`:

- Shared aircraft: wing span, model-defined fan diameter, and fan/LPC pressure-ratio ordering.
- Each mission: fuel capacity, sized MTOW, balanced field length, minimum top-of-climb
  gradient, model-defined compressor-exit temperature, and turbine-metal temperature.
- Temperature maxima include all mission operating points; trajectories must
  also have positive climb gradients and nonnegative cruise distance.

The default scaled constraint tolerance is `1e-6` (one part per million of each
limit). Constraint limits come from the aircraft model. Other aircraft can provide
`aircraft_constraints(ac)` and `mission_constraints(ac, mission_index, context)`
callbacks returning `ConstraintReport` objects with their applicable limits.

Range refinement stops at **any** active constraint. BFL, climb, or temperature
limits may truncate the classical weight/fuel envelope, so point 2 need not
actually fill the tanks. Straight connecting segments are illustrative; only
the solved points and the short-mission origin check have been verified.
A shared aircraft violation, infeasible verification mission, or failure to
obtain any converged trial produces an error rather than a fabricated point.

```julia
pr = TASOPT.payload_range_data(ac;
    altitude_ft=28000:1000:50000,
    range_tol=100.0,           # metres
    short_range=500 * 1852.0, # verification mission, metres
    constraint_tol=1e-6,
    itermax=35)
```

Returned arrays contain range (m), payload (N), fuel including reserves (N),
takeoff weight (N), PFEI, and starting altitude (m). `constraints` contains each
point's mission report; entry 1 corresponds to `origin_check` at 500 nmi.
`aircraft_report`, `design_report`, and `design_feasible` provide the shared and
sizing-mission checks. The zero-range intercept has NaN fuel, takeoff weight,
PFEI, and altitude. Ferry PFEI is also NaN. PFEI is not plotted.

## Axis and sampling options

```julia
# Four characteristic points, plus the design marker (default).
TASOPT.PayloadRange(ac; plots_OEW=false)

# OEW + payload, with every solved boundary sample marked.
TASOPT.PayloadRange(ac; plots_OEW=true, all_points=true, Rpts=10, Ppts=8)
```

`plots_OEW=false` shows payload alone; `true` adds empty weight to both the
boundary and design marker. The input aircraft is unchanged.

With `all_points=true`, `Rpts` is the number of equally spaced range samples
from the short verification mission to corner 1, including both endpoints.
The zero-range intercept is additional; shorter missions are not solved.
`Ppts` is the number of equally spaced payload levels on **each** descending
segment, including its corners. Shared corners appear only once. Both counts
default to 10 and must be at least 2; they have no effect in four-point mode.
For four distinct corners, the dense plot has `Rpts + 2Ppts - 1` boundary markers.

Every additional mission uses the same altitude search and constraint checks;
intermediate payload levels get their own maximum-range solve. Thus the denser
boundary can curve away from the four-point straight segments and takes longer
to compute. The design marker is always added separately.

## Engine limits in the model TOML

```toml
[Propulsion]
max_fan_diameter = "2 m"

[Propulsion.temperature_limits]
compressor_exit = "900 K"
turbine_metal = "1280 K"
```

These are maximum values, checked across all mission operating points for the
engine temperatures. Numbers without units use metres and kelvin. Limits must
be finite and positive. The turbine-metal value also controls cooling sizing;
it takes precedence over the legacy `Propulsion.T_max_metal` key.

For older TOMLs, omitted fan-diameter and compressor-exit limits default to
2 m and 900 K (or values provided by the input template). Omitted turbine-metal
limits use `T_max_metal`. Set appropriate limits explicitly for each validated
model. `save_aircraft_model` writes all three limits and keeps `T_max_metal`
consistent for legacy readers.

In memory, fan diameter and compressor-exit limits are stored in
`ac.engine.max_fan_diameter` and `ac.engine.max_compressor_exit_temperature`;
the existing `ac.parg[igTmetal]` remains the turbine-metal sizing/checking limit.
