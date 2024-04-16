using DocStringExtensions
"""
$TYPEDEF

StructuralMember:
Contains structural and material properties of internal members 

$TYPEDFIELDS
"""
@kwdef mutable struct StructuralMember
    """Material: Automatically sets stress and density of StructuralMember [StructuralAlloy]"""
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    """Weight [N]"""
    weight::Float64 = 0 #OR Union{Float64, Nothing}
    """Stress [Pa]"""
    σ::Float64 = material.σmax
    """Horizontal Stiffnesss [N m^2]""" 
    EIh::Float64 = 0
    """Vertical Stiffnesss [N m^2]"""
    EIv::Float64 = 0
    """Torsional Rigidity [N m^2]"""
    GJ::Float64 = 0
    """Thickness [m]"""
    thickness::Float64 = 0
    """Density [kg/m^3]"""
    ρ::Float64 = material.ρ
    """Position [m]"""
    x::Float64 = 0
    # xv::Float64 = 0
    #Material = StructuralAlloy
end