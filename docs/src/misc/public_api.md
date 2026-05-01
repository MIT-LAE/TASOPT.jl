# Public API and SemVer Scope

This note defines the public API and Semantic Versioning policy for `TASOPT.jl` starting from v3.0.

## Stable public API

Only the symbols listed below are semver-protected.

### Core model lifecycle

These are the key functions that are used while reading an aircraft model definition `.toml` file and sizing the aircraft. [`fly_mission!`](@ref TASOPT.fly_mission!) and [`balance_aircraft!`](@ref TASOPT.balance_aircraft!) let users evaluate ad-hoc scenarios outside the typical sizing routine.

- [`read_aircraft_model`](@ref)
- `load_default_model`
- [`save_aircraft_model`](@ref)
- [`size_aircraft!`](@ref TASOPT.size_aircraft!)
- [`fly_mission!`](@ref TASOPT.fly_mission!)
- [`balance_aircraft!`](@ref TASOPT.balance_aircraft!)

### Primary model types

- [`aircraft`](@ref)

### Data and reporting helpers

- [`output_csv`](@ref)
- [`plot_airf`](@ref)
- [`aeroperf_sweep`](@ref)
- [`stickfig`](@ref)
- [`plot_details`](@ref)
- [`plot_drag_breakdown`](@ref)
- [`PayloadRange`](@ref TASOPT.PayloadRange)
- [`DragPolar`](@ref)

### Public configuration types and enums

- [`StructuralAlloy`](@ref TASOPT.materials.StructuralAlloy), [`Conductor`](@ref TASOPT.materials.Conductor), [`Insulator`](@ref TASOPT.materials.Insulator), [`ThermalInsulator`](@ref TASOPT.materials.ThermalInsulator)
- [`EngineLocation`](@ref), [`PropSysArch`](@ref), [`WingMove`](@ref), [`FuelType`](@ref), [`TrimVar`](@ref), [`TailSizing`](@ref)

### Unit conversion and constants

- `convertMass`, `convertForce`, `convertDist`, `convertSpeed`, `convertPower`, `convertAngle`
- `gee`, `gamSL`, `cpSL`, `μAir`, `pref`, `Tref`

This API still supports the main user workflow:

1. Build or load a model ([`read_aircraft_model`](@ref) / `load_default_model`).
2. Run sizing and mission analysis ([`size_aircraft!`](@ref TASOPT.size_aircraft!), [`fly_mission!`](@ref TASOPT.fly_mission!), [`balance_aircraft!`](@ref TASOPT.balance_aircraft!)).
3. Save and export results ([`save_aircraft_model`](@ref), [`output_csv`](@ref), plotting helpers).

Note: Users can still call any submodule functions directly, but those calls are treated as advanced/internal use and are not semver-stable and can have updated syntax or behavior in minor releases.

## Explicitly non-public API

The following are **not** semver-protected and may change in minor/patch releases:

- `__TASOPTindices__` and symbols brought in via `include(__TASOPTindices__)`. Refactoring these global arrays is a priority on the LAE development plan for `TASOPT.jl`.
- Submodule exports and internals (for example, `TASOPT.engine.*`, `TASOPT.aerodynamics.*`, `TASOPT.CryoTank.*`).
- Non-listed top-level exports.
- `quicksave_aircraft` and `quickload_aircraft` (JLD2 snapshot helpers).
- CSV index helper exports (`default_output_indices`, `output_indices_all`, `output_indices_wGeom`, `output_indices_wEngine`).
- Internal struct field layouts unless explicitly documented in the stable list above.
- Internal solver, gas model, and mission iteration details.


## SemVer policy that will be followed

### MAJOR version bump required

- Breaking changes to any symbol in the stable public API list above.

### MINOR version bump

- Adding new stable API symbols without breaking existing stable public API.
- Behavioral additions that are backward compatible for stable symbols.
- Changes to non-public internals and submodule exports.

### PATCH version bump

- Bug fixes and internal refactors that do not break current submodule exports.

## Release checklist

Before release, if yes to any item below, use a MAJOR bump:

1. Did any stable-listed symbol in the public API list above disappear, rename, or break behavior?
2. Did any stable-listed function reject previously valid calls?
