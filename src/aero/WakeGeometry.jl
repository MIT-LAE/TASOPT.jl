
using StaticArrays
using LinearAlgebra
using DocStringExtensions

"""
    $(TYPEDEF)

Two-dimensional point type using StaticArrays for stack allocation.
First component is y-coordinate, second is z-coordinate.
"""
const Point2D = SVector{2, Float64} # [y, z] really just a shorthand for SVector

"""
    $(TYPEDEF)

Represents a single wake element defined by two points in the y-z plane.

# Fields
$(FIELDS)

It store a few more things inside WakeElement to help with calculations later.
The unit normal is computed as [Δz, -Δy]/|Δs| where Δs is the vector connecting the end points.
"""
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

"""
    WakeElement(wp::SVector{N, Point2D}) where N

Returns an SVector of WakeElements
"""
@views function WakeElement(wp::SVector{N, Point2D}) where N
    SVector{N-1,WakeElement}([WakeElement(a,b) for (a,b) in zip(wp[begin:end-1], wp[begin+1:end])])
end

struct WakeSystem{NP, NE}
    points::SVector{NP, Point2D}
    elements::SVector{NE, WakeElement}
    influence_matrix::SMatrix{NE, NP, Float64}
    
    function WakeSystem(points::SVector{NP, Point2D}, 
        elements::SVector{NE, WakeElement},
        influence_matrix::SMatrix{NE, NP, Float64}) where {NP, NE}
        if NE != NP - 1
            throw(ArgumentError("Number of elements (NE) must be exactly one less than the number of points (NP)."))
        end
        new{NP, NE}(points, elements, influence_matrix)
    end
end

function WakeSystem(points::SVector{NP, Point2D}) where NP
    elements = WakeElement(points)
    NE = NP - 1
    influence_matrix = @MMatrix zeros(NE, NP)
    ws = WakeSystem(points,elements,influence_matrix)
    calculate_influence_matrix!(ws)
    return ws
end  # function WakeSystem
