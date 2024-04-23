using StaticArrays, DocStringExtensions
using LinearAlgebra
abstract type AbstractLoad end

# Direction Vectors
next_id = Ref(0)
struct frame
    """some way of id'ing this frame"""
    id::Int64
    """name of frame"""
    Name::String
    """Origin of this frame"""
    origin::SVector{3, Float64}

    function frame(name::String, origin::AbstractVector)
        instance = new(next_id.x, name, origin)
        next_id.x = next_id.x + 1
        return instance
    end

    """
        frame(name::String)

    Creates a frame with default origin at global origin.
    """
    function frame(name::String)
        instance = new(next_id.x, name, SA[0.0, 0.0, 0.0])
        next_id.x += 1
        return instance
    end
end

const WORLD = frame("World Frame")

"""
$TYPEDEF

Point load.

$TYPEDFIELDS
"""
@kwdef struct PointLoad <: AbstractLoad
    """Point force F̲ = {Fx, Fy, Fz} [N]"""
    force::AbstractVector{Float64} = [0, 0, 0]
    """Origin of reference frame where load is applied [m]"""
    origin::AbstractVector{Float64} = [0,0,0]
    """Centroid (position) of applied force {x, y, z} [m]"""
    r::AbstractVector{Float64}
end

"""
$TYPEDEF

Uniformly distributed load.

$TYPEDFIELDS
"""
@kwdef struct UniformLoad <: AbstractLoad
    """Total equivalent force at centroid F̲ = {Fx, Fy, Fz} [N]"""
    force::AbstractVector{Float64} = [0, 0, 0]
    """Start location of load [m]"""
    x1::AbstractVector{Float64} = [0,0,0]
    """End location of load [m]"""
    x2::AbstractVector{Float64} = [1,0,0]
    """Origin of reference frame where load is applied [m]"""
    origin::AbstractVector{Float64} = [0,0,0]
end

function Base.getproperty(L::UniformLoad, sym::Symbol)
    if sym === :r
        return 0.5*(L.x1 + L.x2)
    elseif sym === :p
        return L.force ./ norm(L.x2 - L.x1)
    else 
        return getfield(L, sym)
    end
end


function Base.show(io::IO, L::AbstractLoad)
    fx,fy,fz = L.force
    x,y,z = L.r
    println(typeof(L))
    print("F̲ = ", fx, "î + ",fy, "ĵ + ", fz,"k̂ N\n")
    print("r̲ = ", x, "î + ",y, "ĵ + ", z,"k̂ m")
end 

"""
    moment(L::AbstractLoad)

Calculates the moment of a given load `L` about the `origin` of that load.
"""
function moment(L::AbstractLoad)
    return L.r × L.force #Cross product to get moment
end  # function moment

"""
    moment_arm(L::AbstractLoad)
    
Calculates the moment arm of a given load `L` about the `origin` of that load.
"""
function moment_arm(L::AbstractLoad)
    return L.r - L.origin
end  # function moment_arm

"""
    moment_about(L::AbstractLoad, point::AbstractVector{Float64})

Calculates the moment of a given load `L` about a specified `point`
"""
function moment_about(L::AbstractLoad, point::AbstractVector{Float64})
    Δ = L.origin - point # position vector of load origin relative to given point
    moment_arm = L.r + Δ # translate load position to relative to given point
    return moment_arm × L.force
end  # function moment_about 



## Sketches:

# import Base.+
# """
# """
# function +(L1::T, L2::T;
#      origin::AbstractVector{Float64}=[0.0,0.0,0.0]) where T<:AbstractLoad
    
#     fx,fy,fz = L1.force + L2.force
#     # Calculate moments
#     m1 = moment_about(L1, origin) 
#     m2 = moment_about(L2, origin) 
#     mx,my,mz = m1 + m2


# end  
# """
# """
# function moment(fuse::Fuselage)
#     moment(internal loads) + moment( structures)
# end  # function moment
