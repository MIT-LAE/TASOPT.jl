abstract type AbstractWingSection end

"""
$TYPEDEF

Wing Section

$TYPEDFIELDS
"""
@kwdef mutable struct WingSection <: AbstractWingSection
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    """Wing Section layout"""
    layout::WingSectionLayout = WingSectionLayout()
    """Wing Section webs"""
    webs::StructuralMember = StructuralMember(material=material)
    """Wing Section caps"""
    caps::StructuralMember = StructuralMember(material=material)
    """Bending Stiffness EI matrix [N m^2]"""
    EI::Matrix{Float64} = zeros(Float64, 2, 2)
    """Torsional Stiffness GJ [N m^2]"""
    GJ::Float64 = 0
    """Lift Rolloff"""
    lift_rolloff::Float64 = 0
    """Max shear load [N]"""
    max_shear_load::Float64 = 0
    """Moment [N m]"""
    moment::Float64 = 0
    """Weight [N]"""
    weight::Float64 = 0
    """Wing root moment contribution from wing weight section of engine [N m]"""
    dyW::Float64 = 0
end

@kwdef mutable struct TailSection <: AbstractWingSection
    """Tail Section layout"""
    layout::WingSectionLayout = WingSectionLayout()
    """Bending Stiffness EI matrix [N m^2]"""
    EI::Matrix{Float64} = zeros(Float64, 2, 2)
    """Torsional Stiffness GJ [N m^2]"""
    GJ::Float64 = 0
    """Caps Thickness [m]"""
    thickness_cap::Float64 = 0 
    """Webs Thickness [m]"""
    thickness_web::Float64 = 0
    dxW::Float64 = 0
end