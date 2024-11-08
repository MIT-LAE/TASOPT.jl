abstract type AbstractWingSection end

"""
$TYPEDEF

Wing Section

$TYPEDFIELDS
"""
@kwdef mutable struct WingSection <: AbstractWingSection
    """Wing Section layout"""
    cross_section::WingCrossSection = WingCrossSection()
    """Starting section chord [m]"""
    co::Float64 = 0.0
    """Taper ratio λ = c_end/c_start"""
    λ::Float64 = 0.0
    """Wing section sweep (assumed to be constant across sections)"""
    sweep::Float64 = 0.0
    """Wing Section webs"""
    webs::StructuralMember = StructuralMember(material=StructuralAlloy("TASOPT-Al"))
    """Wing Section caps"""
    caps::StructuralMember = StructuralMember(material=StructuralAlloy("TASOPT-Al"))
    """Bending Stiffness EI matrix [N m^2]"""
    EI::Matrix{Float64} = zeros(Float64, 2, 2)
    """Torsional Stiffness GJ [N m^2]"""
    GJ::Float64 = 0
    """Max shear load [N]"""
    max_shear_load::Float64 = 0
    """Moment [N m]"""
    moment::Float64 = 0
    """Weight [N]"""
    weight::Float64 = 0 # make it a weight
    """Wing root moment contribution from wing weight section of engine [N m]"""
    dyW::Float64 = 0
end

@kwdef mutable struct TailSection <: AbstractWingSection
    """Tail Section layout"""
    cross_section::WingCrossSection = WingCrossSection()
    """Taper ratio λ = c_end/c_start"""
    λ::Float64 = 0.0
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