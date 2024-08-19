using StaticArrays, DocStringExtensions
using LinearAlgebra
abstract type AbstractLoad end

# Direction Vectors
const î = SA[1.0, 0.0, 0.0]
const ĵ = SA[0.0, 1.0, 0.0]
const k̂ = SA[0.0, 0.0, 1.0]

next_id = Ref(0)
"""
$(TYPEDEF)

Represents a coordinate frame at `origin` w.r.t. the WORLD frame.

$(TYPEDFIELDS)
"""
struct Frame
    """some way of id'ing this frame"""
    id::Int64
    """Origin of this frame"""
    origin::SVector{3, Float64}

    function Frame(origin::AbstractVector)
        instance = new(next_id.x, origin)
        next_id.x = next_id.x + 1
        return instance
    end

    """
        frame(name::String)

    Creates a frame with default origin at global origin.
    """
    function Frame()
        instance = new(next_id.x, SA[0.0, 0.0, 0.0])
        next_id.x += 1
        return instance
    end
end

#WORLD is a privileged global frame.
const WORLD = Frame()

"""
    toWORLD(L::AbstractLoad)

Returns a new equivalent load in the WORLD frame.
"""
function copytoWORLD(L::AbstractLoad)
    O′= L.frame.origin
    O = WORLD.origin
    Δ = O′ - O
    L2 = deepcopy(L)
    L2.frame = WORLD
    L2.r = L.r + Δ
    return L2
end  # function toWORLD

"""
    movetoWORLD!(L::AbstractLoad)

Mutates the load by moving it to the WORLD frame
"""
function movetoWORLD!(L::AbstractLoad)
    O′= L.frame.origin
    O = WORLD.origin
    Δ = O′ - O
    L.frame = WORLD
    L.r = L.r + Δ
    return nothing
end  # function movetoWORLD!

"""
$(TYPEDEF) 

Represents a weight at a given location `r` with respect to a 
specified coordinate `frame`

$(TYPEDFIELDS)
"""
mutable struct Weight <: AbstractLoad
    """Weight [N]"""
    W::Float64
    """Location {x,y,z} [m]"""
    r::SVector{3, Float64}
    """Coordinate Frame"""
    frame::Frame
end

"""
$(TYPEDSIGNATURES)
"""
function Weight(W::Float64, r::AbstractVector)
    Weight(W, r, WORLD)
end  # function Weight

"""
    Weight(;W::Float64, x::Float64=0.0, y::Float64=0.0, z::Float64=0.0, frame::Frame=WORLD)

$(TYPEDSIGNATURES)

"""
function Weight(;W::Float64=0.0, x::Float64=0.0, y::Float64=0.0, z::Float64=0.0, frame::Frame=WORLD)
    Weight(W, SA[x,y,z], frame)
end

import Base.+, Base.*

function +(W1::T, W2::T) where T<:Weight
    W1.frame == W2.frame || error("Cannot add weights in different frames")
    total_W = W1.W + W2.W
    Weight(total_W, (W1.r*W1.W + W2.r*W2.W)/total_W)
end  # function +

function *(W::Weight, fac::Float64)
    return Weight(W.W*fac, W.r)
end

function scale!(W::Weight, fac::Float64)
    W.W = W.W*fac
    return nothing
end

"""
    center_of_weight(W_array::AbstractArray{Weight})

Calculates the coordinates of the center of mass/weight and returns a `Weight`
type of the equivalent weight and at the center of mass.

Note: This is notably faster than doing `sum(W_array)`.
"""
function center_of_weight(W_array::AbstractArray{Weight}, frame::Frame = WORLD)
    total_weight = 0.0

    ## The following three SVector definitions are all aliases:
    # r̄ = SA[0.0, 0.0, 0.0]
    # r̄ = SVector{3}(0.0, 0.0, 0.0)
    r̄ = SArray{Tuple{3}}(0.0, 0.0, 0.0)

    for weight in W_array
        weight.frame == frame || error("Weights not in same frame $weight")
        total_weight += weight.W
        r̄ = r̄ + weight.W * weight.r
    end
    return Weight(total_weight, r̄./total_weight)
end

function Base.getproperty(obj::Weight, sym::Symbol)
    if sym === :x
        return obj.r[1]
    elseif sym === :y
        return obj.r[2]
    elseif sym === :z
        return obj.r[3]
    else
        return getfield(obj, sym)
    end
end

"""
    moment(weight::Weight)

Specialized `moment` function for `Weight`s
    ┌─  ─┐   ┌─  ─┐ ┌─  ─┐   ┌─  ─┐ ┌─  ─┐   ┌─     ─┐   ┌─      ─┐  
    │ Mx │   │ rx │ │ Fx │   │ rx │ │ 0. │   │ ry*Fz │   │ -ry|W| │  
    │ My │ = │ ry │x│ Fy │ = │ ry │x│ 0. │ = │-rx*Fz │ = │ +rx|W| │   
    │ Mz │   │ rz │ │ Fz │   │ rz │ │ Fz │   │  0.0  │   │   0.0  │  
    └─  ─┘   └─  ─┘ └─  ─┘   └─  ─┘ └─  ─┘   └─     ─┘   └─      ─┘                            

"""
function moment(weight::Weight)
    #Weight always points in -z direction
    SA[x_moment(weight), y_moment(weight), z_moment(weight)]
end  # function moment

# Key assumption here is that weights are always forces in the -z direction.
x_moment(weight::Weight) = - weight.W * weight.r[2]
y_moment(weight::Weight) = + weight.W * weight.r[1] # + comes from - * -
z_moment(weight::Weight) = 0.0

function Base.show(io::IO, W::Weight)
    x,y,z = W.r
    xo,yo,zo = W.frame.origin
    println(typeof(W))
    print("F̲ = ", -W.W,"k̂ N; ")
    print("r̲ = [", x, "î, ",y, "ĵ, ", z,"k̂] m; ")
    print("in Frame[",W.frame.id,"], o̲ = ", xo, "î + ",yo, "ĵ + ", zo,"k̂ m\n")
end 

"""
$TYPEDEF

Point load.

$TYPEDFIELDS
"""
@kwdef mutable struct PointLoad <: AbstractLoad
    """Point force F̲ = {Fx, Fy, Fz} [N]"""
    force::AbstractVector{Float64} = [0, 0, 0]
    """Reference frame"""
    frame::Frame = WORLD
    """Centroid (position) of applied force {x, y, z} [m]"""
    r::AbstractVector{Float64} = [0,0,0]
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
    """Reference frame"""
    frame::Frame = WORLD
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
    xo,yo,zo = L.frame.origin
    println(typeof(L))
    print("o̲ = ", xo, "î + ",yo, "ĵ + ", zo,"k̂ m\n")
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
    return L.r - L.frame.origin
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
