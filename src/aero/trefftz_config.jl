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

function get_trefftz_discretization(s::String)
    quality = uppercase(s)
    if quality == "FINE" #328 panels. Reference case. 
        return (wing_panels = SurfaceDiscretization(160, 48, 24),
                tail_panels = SurfaceDiscretization(80, 0, 16))
    elseif quality == "MEDIUM" #gives about 0.60% higher CDi relative to the FINE
        return (wing_panels = SurfaceDiscretization(40, 12, 6),
                tail_panels = SurfaceDiscretization(20, 0, 4))
    elseif quality == "COARSE" # COARSE #gives about 1.8% higher CDi relative to the FINE
        return (wing_panels = SurfaceDiscretization(20, 6, 3),
                tail_panels = SurfaceDiscretization(10, 0, 2))
    else
        throw(ArgumentError("Invalid discretization resolution: $s. Must be 'COARSE', 'MEDIUM', or 'FINE'."))
    end
end

const DEFAULT_TREFFTZ_DISCRETIZATION = get_trefftz_discretization("MEDIUM")

"""
Trefftz plane analysis configuration parameters.
Controls point clustering and tip correction behavior.
"""
struct TrefftzConfig
    k_tip::Float64          # Tip loading exponent (default: 16.0)
    bunch::Float64          # Center clustering factor ∈ [0,1] (default: 0.5)
    root_contraction::Float64  # Root streamline contraction (default: 0.2)

    function TrefftzConfig(; k_tip=16.0, bunch=0.5, root_contraction=0.2)
        0.0 <= bunch <= 1.0 || throw(ArgumentError("bunch must be in [0,1]"))
        k_tip > 0.0 || throw(ArgumentError("k_tip must be positive"))
        root_contraction > 0.0 || throw(ArgumentError("root_contraction must be positive"))
        new(k_tip, bunch, root_contraction)
    end
end