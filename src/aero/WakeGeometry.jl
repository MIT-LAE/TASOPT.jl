
using StaticArrays
using LinearAlgebra
using DocStringExtensions

const MIN_DISTANCE_SQUARED = eps()

"""
    $(TYPEDEF)

Two-dimensional point type using StaticArrays for stack allocation.
First component is y-coordinate (along wing direction), second is z-coordinate (up/down).
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
        if (length == 0.0)
            throw(ArgumentError("p1 and p2 must be distinct"))
        end
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

#Convenience functions to get element lengths, Δy's and Δz's
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


"""
    $(TYPEDEF)

A complete wake system with N wake points forming N-1 wake elements.

# Type Parameters
$(TYPEDFIELDS)

The `influence_matrix` is purely geometric! So compute once and reuse it for different
loading conditions. For each element-point pair, the coefficient is:
```math
\\left(\\hat{x} \\times (\\mathbf{r}_{cp} - \\mathbf{r}_{wp})\\right) \\cdot \\hat{n}
```
where:
- ``\\hat{x}`` is the unit vector in x-direction
- ``\\mathbf{r}_{cp}`` is the control point position vector
- ``\\mathbf{r}_{wp}`` is the wake point position vector
- ``\\hat{n}`` is the unit normal to the wake element
"""
struct WakeSystem{NP, NE}
    points::SVector{NP, Point2D}
    elements::SVector{NE, WakeElement}
    influence_matrix::MMatrix{NE, NP, Float64}
    
    function WakeSystem(points::SVector{NP, Point2D}, 
        elements::SVector{NE, WakeElement},
        influence_matrix::MMatrix{NE, NP, Float64}) where {NP, NE}
        if NE != NP - 1
            throw(ArgumentError("Number of elements (NE) must be exactly one less than the number of points (NP)."))
        end
        new{NP, NE}(points, elements, influence_matrix)
    end
end

element_lengths(WS::WakeSystem) = element_lengths(WS.elements)
element_dys(WS::WakeSystem) = element_dys(WS.elements)
element_dzs(WS::WakeSystem) = element_dzs(WS.elements)

"""
"""
function WakeSystem(points::SVector{NP, Point2D}; 
    control_points::Union{SVector{NC, Point2D}, Nothing}=nothing) where {NP, NC}

    elements = generate_wake_elements(points; control_points = control_points)
    NE = NP - 1
    influence_matrix = @MMatrix zeros(NE, NP)
    ws = WakeSystem(points,elements,influence_matrix)
    calculate_influence_matrix!(ws)
    return ws
end  # function WakeSystem

"""
    calculate_influence_coefficient(r_vec, normal)

Returns the 2D vortex influence coefficient induced at a control point by a unit vortex
at the wake point. 

`r_vec` is the displacement vector `[Δy, Δz]` from the vortex to the evaluation point.
⟹`  n⋅(x̂ × r)/|r|²`,
  `= (-ny*z + nz*y)/|r|²  ⟸ x̂ × r = [-z, y]`

If `|r|²` is below a threshold (`MIN_DISTANCE_SQUARED`), the value is clipped to zero to avoid singularities.
"""
@inline function calculate_influence_coefficient(r_vec::Point2D, normal::Point2D)
    r_squared = dot(r_vec, r_vec)
    # Effectively checking if the control point is too close to a wake point. 
    # Could be improved with some model of a finite core vortex?
    if r_squared > MIN_DISTANCE_SQUARED
        # Cross product x̂ × r gives [-z, y], dot with normal [ny, nz] gives = y*nz - z*ny
        return (r_vec[1] * normal[2] - r_vec[2] * normal[1]) / r_squared
    end
    return 0.0
end

"""
    mirror_point(p1::Point2D)

Mirrors point about the x-z plane by negating the y coordinate
p1 = [y, z] --> [-y, z]
"""
function mirror_point(p1::Point2D)
    return Point2D(-p1[1], p1[2])
end 

function calculate_influence_matrix!(ws::WakeSystem{NP,NE}) where {NP,NE}
    @inbounds for j in 1:NP
        wake_point = ws.points[j]
        mirrored_point = mirror_point(wake_point)
        @inbounds for i in 1:NE
            element = ws.elements[i]
            # Get contribution from j'th wake point
            r_vec = element.control_point - wake_point
            influence = calculate_influence_coefficient(r_vec, element.unit_normal)
            # Now do the same for the mirrored point:
            r_vec_mirror = element.control_point - mirrored_point
            # Note the subtraction below is because the contribution of the "image" leg of the 
            # horse shoe vortex is going to be in the opposite orientation so the negation.
            influence = influence - calculate_influence_coefficient(r_vec_mirror, element.unit_normal)
            
            ws.influence_matrix[i,j] = influence
        end
    end

    return nothing
end

"""
    generate_wake_system(y::AbstractVector{Float64}, z::AbstractVector{Float64}, 
                            ycp::AbstractVector{Float64}, zcp::AbstractVector{Float64})

Constructs a WakeSystem from vectors of y and z coordinates for wake points and control points.
!!! details "🔃 Inputs and Outputs"
**Inputs:**
    - `yp`: Vector of y-coordinates for wake points
    - `zp`: Vector of z-coordinates for wake points
    - `ycp`: Vector of y-coordinates for control points
    - `zcp`: Vector of z-coordinates for control points
    (the _p suffix is to remember that these are locations far downstream in the wake (y′,z′))
"""
function generate_wake_system(yp::AbstractVector{Float64}, zp::AbstractVector{Float64}, 
                            ycp::AbstractVector{Float64}, zcp::AbstractVector{Float64})
    # Check dimensions
    np = length(yp)
    @assert length(zp) == np "yp and zp must have same length"
    @assert length(ycp) == np-1 "Control points must have one less control point than wake points"
    @assert length(zcp) == np-1 "Control points must have one less control point than wake points"
    
    # Create arrays of wake points and control points 
    points = SVector{np, Point2D}(Point2D(yp[i], zp[i]) for i in 1:np)
    control_points = SVector{np-1, Point2D}(Point2D(ycp[i], zcp[i]) for i in 1:np-1)
    WakeSystem(points, control_points = control_points)
end