
using StaticArrays
using LinearAlgebra
using DocStringExtensions

"""
    $(TYPEDEF)

Two-dimensional point type using StaticArrays for stack allocation.
First component is y-coordinate, second is z-coordinate.
# Examples
```julia-repl
julia> p = Point2D(1.0, 2.0) # Creates point at y=1.0, z=2.0
2-element SVector{2, Float64} with indices SOneTo(2):
 1.0
 2.0
```
"""
const Point2D = SVector{2, Float64} # [y, z] really just a shorthand for SVector


"""
    $(TYPEDEF)

Represents a single wake element defined by two points in the y-z plane and a control
point. Note the control point does *not* have to be the midpoint, but it defaults
to that if nothing is provided.

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
    "Length of element"
    length::Float64
    "y-extent of element"
    Δy::Float64
    "z-extent of element"
    Δz::Float64
    "Unit normal vector to the element"
    unit_normal::Point2D
    
    # Inner constructor to help in calculating the control point and the normal
    function WakeElement(p1::Point2D, p2::Point2D;
        control_point::Union{Point2D, Nothing}=nothing)
        Δs = p2 - p1
        length = norm(Δs)
        # Normal points "up" in s-n-l coordinate system along the sheet where 
        # s is along wake sheet (left to right), l is out of the page like x̂ and so n is up.
        unit_normal = Point2D(-Δs[2] / length, Δs[1] / length)
        # Midpoint calculation as control point if not provided
        control_point = isnothing(control_point) ? (p1 + p2) * 0.5 : control_point
        new(p1, p2, control_point, length, Δs[1], Δs[2], unit_normal)
    end
end


"""
    generate_wake_elements(points::SVector{N, Point2D};
    control_points::Union{SVector{M, Point2D}, Nothing}=nothing) where {N, M}

Returns an SVector of WakeElements
"""
@inline function generate_wake_elements(points::SVector{N, Point2D};
    control_points::Union{SVector{M, Point2D}, Nothing}=nothing) where {N, M}
    if isnothing(control_points)
        # If no control points are provided, calculate midpoints and set that as the control point
        return SVector{N-1,WakeElement}(WakeElement(points[i], points[i+1]) for i in 1:N-1)
    else
        if M != N - 1
            throw(ArgumentError("Number of control points must be exactly one less than the number of points."))
        end
    return SVector{N-1,WakeElement}(WakeElement(points[i], points[i+1]; 
                                        control_point=control_points[i]) for i in 1:N-1)
    end
end

function element_lengths(wake_elements::SVector{N, WakeElement}) where N
    SVector{N, Float64}(wake_elements[i].length for i in 1:N)
end

"""
"""
function element_dys(wake_elements::SVector{N, WakeElement}) where N
    SVector{N, Float64}(wake_elements[i].Δy for i in 1:N)
end  # function element_dys

function element_dzs(wake_elements::SVector{N, WakeElement}) where N
    SVector{N, Float64}(wake_elements[i].Δz for i in 1:N)
end  # function element_dys