using DocStringExtensions
"""
$TYPEDEF

StructuralMember:
Contains structural and material properties of structural members 

$TYPEDFIELDS
"""
@kwdef mutable struct StructuralMember
    """Material: Automatically sets stress and density of StructuralMember [StructuralAlloy]"""
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    """Weight [N]"""
    weight::Float64 = 0 #OR Union{Float64, Nothing} ::Weight
    """Stress [Pa]"""
    σ::Float64 = material.σmax
    """Horizontal/Bending Stiffness [N m^2]""" 
    EIh::Float64 = 0
    """Vertical/Normal Stiffness [N m^2]"""
    EIv::Float64 = 0
    """Torsional Rigidity [N m^2]"""
    GJ::Float64 = 0
    """Thickness [m]"""
    thickness::Float64 = 0
    """Density [kg/m^3]"""
    ρ::Float64 = material.ρ
    """Position [m]"""
    x::Float64 = 0
    """Weight Lateral Distribution"""
    dxW::Float64 = 0
    #Material = StructuralAlloy
end

function Base.show(io::IO, x::StructuralMember)
    print(io, "StructuralMember(", x.material,", ", x.weight, ", ", x.x, ")")
end

"""
$TYPEDEF

InternalMember:
Contains structural and material properties of internal members 

$TYPEDFIELDS
"""
@kwdef mutable struct InternalMember
    """Material: Automatically sets stress and density of StructuralMember [StructuralAlloy]"""
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    """Weight [N]"""
    weight::Float64 = 0 #OR Union{Float64, Nothing} ::Weight
    """Thickness [m]"""
    thickness::Float64 = 0
    """Weight per length [N/m]"""
    W_per_length::Float64 = 0
    """Weight per area [N/m^2]"""
    W_per_area::Float64 = 0
   
end

function Base.show(io::IO, x::InternalMember)
    print(io, "InternalMember(", x.material,", ", x.weight, ", ", x.x, ")")
end