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

n_total(sd::SurfaceDiscretization) = sd.n_outer_panels + sd.n_inner_panels + sd.n_image_panels

"""
    TrefftzPlaneConfig

Configuration for Trefftz plane induced drag analysis.
Contains both panel discretization and analysis constants/ parameters.
The surfaces are ordered wing, tail, strut. The strut is only part of the analysis when `strut_panels`
has panels; by default it has none and the analysis has two surfaces. [`get_trefftz_config`](@ref) gives
the strut halfway between the panel counts of the wing and the tail.

# Example
```julia
config = TrefftzPlaneConfig(
    wing_panels = SurfaceDiscretization(40, 12, 6),
    tail_panels = SurfaceDiscretization(20, 0, 4),
    k_tip = 16.0,
    bunch = 0.5,
    wing_root_contraction = 0.2,
    tail_root_contraction = 1.0,
    strut_panels = SurfaceDiscretization(24, 0, 0) # optional: adds a strut as a third surface
)
```
"""
struct TrefftzPlaneConfig
    wing_panels::SurfaceDiscretization
    tail_panels::SurfaceDiscretization
    k_tip::Float64
    bunch::Float64
    wing_root_contraction::Float64
    tail_root_contraction::Float64
    strut_panels::SurfaceDiscretization

    function TrefftzPlaneConfig(wing_panels::SurfaceDiscretization,
                                tail_panels::SurfaceDiscretization;
                                k_tip::Float64=16.0,
                                bunch::Float64=0.5,
                                wing_root_contraction::Float64=0.2,
                                tail_root_contraction::Float64=1.0,
                                strut_panels::SurfaceDiscretization=SurfaceDiscretization(0, 0, 0))

        0.0 <= bunch <= 1.0 || throw(ArgumentError("bunch must be in [0,1]"))
        k_tip > 0.0 || throw(ArgumentError("k_tip must be positive"))
        0.0 < wing_root_contraction ≤ 1.0 || throw(ArgumentError("wing_root_contraction must be [0..1]"))
        0.0 < tail_root_contraction ≤ 1.0 || throw(ArgumentError("tail_root_contraction must be [0..1]"))

        i_sum = wing_panels.n_outer_panels +
                wing_panels.n_inner_panels +
                wing_panels.n_image_panels + 1 +
                tail_panels.n_outer_panels +
                tail_panels.n_inner_panels +
                tail_panels.n_image_panels + 1 +
                n_strut_points(strut_panels)
        geom_size = length(TREFFTZ_GEOM.y)
        if i_sum > geom_size
            error("TREFFTZ: Geometry array overflow. Required $i_sum panels but TrefftzGeometry{$geom_size} only has $geom_size slots. " *
                  "Check panel discretization settings.")
        end
        new(wing_panels, tail_panels, k_tip, bunch, wing_root_contraction, tail_root_contraction, strut_panels)
    end
end

"""
    n_points_used(config::TrefftzPlaneConfig)

Returns the total number of points used in the Trefftz plane calculations
based on the provided configuration.
"""
function n_points_used(trefftz_config::TrefftzPlaneConfig)
    n_wings = trefftz_config.wing_panels.n_outer_panels +
              trefftz_config.wing_panels.n_inner_panels +
              trefftz_config.wing_panels.n_image_panels + 1 
              #add 1 for dummy point between surfaces
    n_tails = trefftz_config.tail_panels.n_outer_panels +
             trefftz_config.tail_panels.n_inner_panels +
             trefftz_config.tail_panels.n_image_panels + 1 
             #add 1 for dummy point between surfaces
    return n_wings + n_tails + n_strut_points(trefftz_config.strut_panels)
end

"""
    n_strut_points(strut_panels::SurfaceDiscretization)

Number of points in the strut block of the Trefftz plane calculations: zero if the strut has no panels,
otherwise its panels plus 1 for the dummy point after the surface.
"""
n_strut_points(strut_panels::SurfaceDiscretization) =
    n_total(strut_panels) == 0 ? 0 : n_total(strut_panels) + 1

"""
    n_surfaces(config::TrefftzPlaneConfig)

Number of lifting surfaces in the Trefftz plane analysis: 2 (wing, tail), or 3 (wing, tail, strut)
if the configuration has strut panels.
"""
n_surfaces(config::TrefftzPlaneConfig) = n_total(config.strut_panels) == 0 ? 2 : 3

i_first_wing(config::TrefftzPlaneConfig) = 1
i_last_wing(config::TrefftzPlaneConfig) = i_first_wing(config) +
                                          n_total(config.wing_panels)
i_first_tail(config::TrefftzPlaneConfig) = i_last_wing(config) + 1
i_last_tail(config::TrefftzPlaneConfig) = i_first_tail(config) +
                                          n_total(config.tail_panels)
i_first_strut(config::TrefftzPlaneConfig) = i_last_tail(config) + 1
i_last_strut(config::TrefftzPlaneConfig) = i_first_strut(config) +
                                           n_total(config.strut_panels)

"""
    get_trefftz_config(quality::String; braced=false, k_tip=16.0, bunch=0.5, wing_root_contraction=0.2, tail_root_contraction=1.0)

Create a [`TrefftzPlaneConfig`](@ref).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `quality::String`: Discretization quality level ("COARSE", "MEDIUM", or "FINE").
    - `braced::Bool=false`: If `true`, the strut is included as a third surface with halfway between the total panel
      counts of the wing and the tail, rounded down (20, 41 or 164 for "COARSE", "MEDIUM" or "FINE"), all laid out
      along the strut. If `false`, the analysis has the wing and tail only.
    - `k_tip::Float64=16.0`: Tip loading exponent (controls circulation decay at wing tips).
    - `bunch::Float64=0.5`: Center clustering factor ∈ \\[0,1\\] (controls spanwise panel spacing near root).
    - `wing_root_contraction::Float64=0.2`: Wing root streamline contraction factor (stream tube contraction near fuselage).
    - `tail_root_contraction::Float64=1.0`: Tail root streamline contraction factor (typically 1.0 = no contraction).

    **Outputs:**
    - `TrefftzPlaneConfig`: Complete configuration for Trefftz plane induced drag analysis.

!!! info "Quality Levels"
    - **"FINE"** (328 panels): Reference accuracy, slowest
      - Wing: 160 outer, 48 inner, 24 image
      - Tail: 80 outer, 0 inner, 16 image
    - **"MEDIUM"** (84 panels): ~0.6% higher CDi vs FINE,
      - Wing: 40 outer, 12 inner, 6 image
      - Tail: 20 outer, 0 inner, 4 image
    - **"COARSE"** (43 panels): ~1.8% higher CDi vs FINE, fastest
      - Wing: 20 outer, 6 inner, 3 image
      - Tail: 10 outer, 0 inner, 2 image
"""
function get_trefftz_config(quality::String;
                           braced::Bool=false,
                           k_tip::Float64=16.0,
                           bunch::Float64=0.5,
                           wing_root_contraction::Float64=0.2,
                           tail_root_contraction::Float64=1.0)
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

    # The strut has no fuselage image region, so its panels are all ordinary panels
    n_strut_panels = div(n_total(wing_panels) + n_total(tail_panels), 2)
    strut_panels = braced ? SurfaceDiscretization(n_strut_panels, 0, 0) : SurfaceDiscretization(0, 0, 0)

    return TrefftzPlaneConfig(wing_panels, tail_panels;
                             k_tip=k_tip,
                             bunch=bunch,
                             wing_root_contraction=wing_root_contraction,
                             tail_root_contraction=tail_root_contraction,
                             strut_panels=strut_panels)
end

"""
    TrefftzGeometry{N}

Geometry work arrays for Trefftz calculations.

# Fields
- `t, y, yp, z, zp`: Vortex point coordinates _p indicates far downstream wake points.
- `yc, ycp, zc, zcp`: Control point coordinates
- `gc`: Circulation strength at the control points.

"""
struct TrefftzGeometry{N}
    t::MVector{N, Float64}
    y::MVector{N, Float64}
    yp::MVector{N, Float64}
    z::MVector{N, Float64}
    zp::MVector{N, Float64}
    yc::MVector{N, Float64}
    ycp::MVector{N, Float64}
    zc::MVector{N, Float64}
    zcp::MVector{N, Float64}
    gc::MVector{N, Float64}
end

# These functions are to create properly sized TrefftzGeometry structs but for now
# we use the fixed-size module arrays of aero.jl to avoid changing too much code at once.
"""
    calculate_array_size(config::TrefftzPlaneConfig, htail)

Calculate required array size for geometry arrays based on panel configuration.
Accounts for wing, tail and (if present) strut panels, handling T-tail case where root_span == 0.
"""
function calculate_array_size(config::TrefftzPlaneConfig, htail)
    # Wing contribution
    n_panels = config.wing_panels.n_outer_panels +
               config.wing_panels.n_inner_panels +
               config.wing_panels.n_image_panels + 1

    # Tail contribution (handle T-tail case where root_span == 0)
    tail_image_panels = (htail.layout.root_span == 0.0) ? 0 : config.tail_panels.n_image_panels
    n_panels += config.tail_panels.n_outer_panels +
                config.tail_panels.n_inner_panels +
                tail_image_panels + 1

    # Strut contribution (none for a cantilever wing)
    n_panels += n_strut_points(config.strut_panels)

    return n_panels
end

"""
    TrefftzGeometry{N}() where {N}

Create zero-initialized geometry arrays of size N.
"""
function TrefftzGeometry{N}() where {N}
    TrefftzGeometry{N}(
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64}),
        zeros(MVector{N, Float64})
    )
end

"""
    TrefftzGeometry(config::TrefftzPlaneConfig, htail)

Create geometry arrays sized for the given configuration.
Uses StaticArrays for stack allocation, providing excellent performance.
"""
function TrefftzGeometry(config::TrefftzPlaneConfig, htail)
    N = calculate_array_size(config, htail)
    return TrefftzGeometry{N}()
end
