
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
