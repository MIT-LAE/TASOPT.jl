abstract type AbstractWingSection end

"""
$TYPEDEF

Strut Webs and Caps

$TYPEDFIELDS
"""
@kwdef mutable struct Web_Cap_Strut
    """Bending Stiffness [N m^2] """
    EI_bending::Float64 = 0
    """Normal Stiffness [N m^2] """
    EI_normal::Float64 = 0
    """Torsional Rigidity [N m^2]"""
    GJ::Float64 = 0
end

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
    """Wing Strut web and caps"""
    web_cap::Web_Cap_Strut = Web_Cap_Strut() # Do we need this and the webs + caps??
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
    """Bending Stiffness [N m^2] """
    EI_bending::Float64 = 0
    """Normal Stiffness [N m^2] """
    EI_normal::Float64 = 0
    """Torsional Rigidity [N m^2]"""
    GJ::Float64 = 0
    """Caps Thickness [m]"""
    thickness_cap::Float64 = 0 
    """Webs Thickness [m]"""
    thickness_web::Float64 = 0
    dxW::Float64 = 0
end