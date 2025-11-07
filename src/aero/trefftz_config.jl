# This file has a set of configuration structs to neatly hold parameters
# related to the Trefftz plane drag estimation method.

"""
Struct to hold the discretization parameters for surface panels
in the far field drag estimates in the Trefftz Plane.
"""
struct SurfaceDiscretization
    n_outer_panels::Int # Number of panels in the outboard wing section
    n_inner_panels::Int # Number of panels in the inboard wing section
    n_image_panels::Int # Number of panels "inside" the fuselage

    function SurfaceDiscretization(n_outer_panels, n_inner_panels, n_image_panels)
        all(>=(0), (n_outer_panels, n_inner_panels, n_image_panels)) ||
            throw(ArgumentError("Number of panels must be non-negative integers."))
        new(n_outer_panels, n_inner_panels, n_image_panels)
    end
end

"""
    TrefftzPlaneConfig

Configuration for Trefftz plane induced drag analysis.
Contains both panel discretization and analysis constants/ parameters.

# Example
```julia
config = TrefftzPlaneConfig(
    wing_panels = SurfaceDiscretization(40, 12, 6),
    tail_panels = SurfaceDiscretization(20, 0, 4),
    k_tip = 16.0,
    bunch = 0.5,
    root_contraction = 0.2
)
```
"""
struct TrefftzPlaneConfig
    wing_panels::SurfaceDiscretization
    tail_panels::SurfaceDiscretization
    k_tip::Float64
    bunch::Float64
    root_contraction::Float64

    function TrefftzPlaneConfig(wing_panels::SurfaceDiscretization,
                                tail_panels::SurfaceDiscretization;
                                k_tip::Float64=16.0,
                                bunch::Float64=0.5,
                                root_contraction::Float64=0.2)
        0.0 <= bunch <= 1.0 || throw(ArgumentError("bunch must be in [0,1]"))
        k_tip > 0.0 || throw(ArgumentError("k_tip must be positive"))
        root_contraction > 0.0 || throw(ArgumentError("root_contraction must be positive"))
        
        new(wing_panels, tail_panels, k_tip, bunch, root_contraction)
    end
end

"""
    get_trefftz_config(quality::String; k_tip=16.0, bunch=0.5, root_contraction=0.2)

Create a [`TrefftzPlaneConfig`](@ref).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `quality::String`: Discretization quality level ("COARSE", "MEDIUM", or "FINE").
    - `k_tip::Float64=16.0`: Tip loading exponent (controls circulation decay at wing tips).
    - `bunch::Float64=0.5`: Center clustering factor ∈ \\[0,1\\] (controls spanwise panel spacing near root).
    - `root_contraction::Float64=0.2`: Root streamline contraction factor (stream tube contraction near fuselage).

    **Outputs:**
    - `TrefftzPlaneConfig`: Complete configuration for Trefftz plane induced drag analysis.

!!! info "Quality Levels"
    - **"FINE"** (328 panels): Reference accuracy, slowest
      - Wing: 160 outer, 48 inner, 24 image
      - Tail: 80 outer, 0 inner, 16 image
    - **"MEDIUM"** (84 panels): ~0.6% higher CDi vs FINE, recomended default
      - Wing: 40 outer, 12 inner, 6 image
      - Tail: 20 outer, 0 inner, 4 image
    - **"COARSE"** (43 panels): ~1.8% higher CDi vs FINE, fastest
      - Wing: 20 outer, 6 inner, 3 image
      - Tail: 10 outer, 0 inner, 2 image
"""
function get_trefftz_config(quality::String;
                           k_tip::Float64=16.0,
                           bunch::Float64=0.5,
                           root_contraction::Float64=0.2)
    quality_upper = uppercase(quality)

    if quality_upper == "FINE" # 328 panels. Reference case.
        wing_panels = SurfaceDiscretization(160, 48, 24)
        tail_panels = SurfaceDiscretization(80, 0, 16)
    elseif quality_upper == "MEDIUM" # ~0.60% higher CDi relative to FINE
        wing_panels = SurfaceDiscretization(40, 12, 6)
        tail_panels = SurfaceDiscretization(20, 0, 4)
    elseif quality_upper == "COARSE" # ~1.8% higher CDi relative to FINE
        wing_panels = SurfaceDiscretization(20, 6, 3)
        tail_panels = SurfaceDiscretization(10, 0, 2)
    else
        throw(ArgumentError("Invalid discretization resolution: $quality. Must be 'COARSE', 'MEDIUM', or 'FINE'."))
    end

    return TrefftzPlaneConfig(wing_panels, tail_panels;
                             k_tip=k_tip,
                             bunch=bunch,
                             root_contraction=root_contraction)
end

const DEFAULT_TREFFTZ_CONFIG = get_trefftz_config("MEDIUM")
