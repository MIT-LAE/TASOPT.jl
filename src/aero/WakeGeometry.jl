
using StaticArrays
using LinearAlgebra
using DocStringExtensions

const Point2D = SVector{2, Float64} # [y, z] really just a shorthand for SVector

struct WakeElement
    "Starting point of the wake element"
    p1::Point2D
    "Ending point of the wake element"
    p2::Point2D
    "Control point (midpoint) of the element"
    control_point::Point2D
    "Unit normal vector to the element"
    unit_normal::Point2D
    
    # Inner constructor to help in calculating the control point and the normal
    function WakeElement(p1::Point2D, p2::Point2D)
        Δs = p2 - p1
        # Use LinearAlgebra's norm for potential SIMD optimization
        length = norm(Δs)
        # Construct unit normal ensuring right-hand convention
        unit_normal = Point2D(Δs[2] / length, -Δs[1] / length)
        # Midpoint calculation
        control_point = (p1 + p2) * 0.5
        new(p1, p2, control_point, unit_normal)
    end
end

@views function WakeElement(wp::SVector{N, Point2D}) where N
    SVector{N-1,WakeElement}([WakeElement(a,b) for (a,b) in zip(wp[begin:end-1], wp[begin+1:end])])
end
