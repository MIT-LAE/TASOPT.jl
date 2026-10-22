# Verification and Validation against real world aircraft

Public information, including airport planning documents (APD) and the ICAO Engine Emissions Databank (EDB), provides a basis for TASOPT models of different aircraft classes.

The primary comparison is mission performance, represented by payload-range diagrams. Aircraft weights such as maximum takeoff mass (MTOM) and operating empty weight (OEW) are secondary comparisons; improving payload-range agreement can also improve these estimates. Each aircraft report identifies its source data and compares geometry, payload range, aircraft weights, and available EDB thrust and fuel-flow evidence.

Each model was optimized to minimize mission-weighted fuel burn. Five off-design ranges supplement the sizing mission, with weights based on their frequency in annual schedules. The exact weights retained from each study are listed in its report; the schedule dataset is not distributed here. These are local optimization results, not proofs of global optimality.

| Class | Reference aircraft | Model | V&V analysis |
|---|---|---|---|
| Narrow body | Boeing 737-800 | [default_narrow.toml](../defaults/default_narrow.toml) | [Narrow body](narrow.md) |
| Regional jet | Embraer E195 | [default_regional.toml](../defaults/default_regional.toml) | [Regional jet](regional.md) |
| Wide body | Boeing 777-300ER | [default_wide.toml](../defaults/default_wide.toml) | [Wide body](wide.md) |

## Evidence and reproducibility

[Source evidence hashes](data/source_evidence.json) identify the archived APDs, EDB workbook and source reports. Numerical CSV/TOML receipts preserve their original units and source-study paths; those paths are provenance, not runtime dependencies.

The reports preserve the selected 21 October 2026 studies and their original analysis methods. Their main plots and numerical tables are archived comparisons. Captured plot captions retain the source study’s terminology; “current runtime” inside those images refers to that archived runtime. Separately labeled current plotting-function checks use this branch and the script below. Calibration against an APD curve is not independent validation. Outstanding checks and operating-limit violations remain visible in each report; these models are not presented as formally validated releases.

The source models were expanded using their previous TASOPT default templates before replacing those templates, so omitted inputs retain their original values. Model names/descriptions were updated and model-defined fan/temperature limits added. [Model provenance](data/model_provenance.toml) records source, template and imported SHA-256 hashes. `default_input.toml` remains unchanged.

The current `PayloadRange` function searches starting cruise altitude from 28,000 to 50,000 ft in 1,000-ft steps and enforces its aircraft and mission constraints with scaled tolerance `1e-6`. That is a different calculation from some archived studies below. This altitude interval bounds the starting altitude, not the whole cruise-climb trajectory. Range brackets terminate at 100 m (about 0.054 nmi); the zero-range intercept is backed by a 500-nmi maximum-payload check. Straight connecting lines do not certify every intermediate state.

To size the imported models and generate current four-point diagrams separately from the archived comparisons, run from the repository root:

```sh
julia --project example/validation/generate_plots.jl
```

The script writes to `example/validation/generated/`, reports constraint failures without relaxing them, and retains the design point. It does not rerun the original optimizers. Percent error is `100 × (model/reference − 1)`; zero-reference errors are undefined. Empty-weight and MZFW proxies, thrust conventions and digitization limitations are identified in the aircraft reports.

## Import verification

On Julia 1.11.2 with the plotting implementation at TASOPT commit `da6f4f50`, all three imported models size successfully and reproduce the archived model takeoff masses within the regression tolerance of 0.01%. The full suite passes **7634/7634 tests**, including the five off-design missions for each model; the documentation build also passes.

Fresh four-point plots are included at the end of each aircraft report. Every solved boundary point passes the current constraint checks. The E195 design point passes; the B738 design point exceeds BFL/climb tolerance, and the 777 design point exceeds BFL tolerance. Exact residuals are preserved in each `data/<class>/current_runtime.toml`. Sizing convergence alone does not establish operating feasibility or independent validation.
