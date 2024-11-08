using DocStringExtensions
using ..aerodynamics
"""
$TYPEDEF

Wing Strut

$TYPEDFIELDS
"""
@kwdef mutable struct Strut
    """Strut Material """
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    """Strut Area [m^2] """
    S::Float64 = 0
    """Strut Length [m] """
    length::Float64 = 0
    """Strut z position [m]"""
    z::Float64 = 0
    """Strut thickness to chord"""
    thickness_to_chord::Float64 = 0
    """Strut local velocity ratio"""
    local_velocity_ratio::Float64 = 0
    """Cosine of strut sweep angle"""
    cos_lambda::Float64 = 0 
    """Strut weight [N]"""
    weight::Float64 = 0
    """Strut axial force [N]"""
    axial_force::Float64 = 0
    """Strut chord [m]"""
    chord::Float64 = 0
    """Aircraft pitching moment contribution from the weight distribution of the strut [N m]"""
    dxW::Float64 = 0
end

"""
$TYPEDEF

Wing Structure:
    Divided into 6 modules
    1. General Properties
    2. Wing Layout
    3. Material
    4. Wing Sections
    5. Strut
    6. Weight Fractions

$TYPEDFIELDS
"""
@kwdef mutable struct Wing
    """Wing Weight [N] """
    weight::Float64 = 0
    """Aircraft pitching moment contribution from the weight distribution of the wing [Nm]"""
    dxW::Float64 = 0
    """Wing Layout """
    layout::WingLayout = WingLayout()
    """Wing Planform (0: wing catilever, plain; 1: wing cantilever with engine"""
    planform::Int64 = 0 # 0: Wing cantilever, plain
                      # 1: Wing cantilever with engine
    """Wing Material """
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")

    """Airfoil data"""
    airsection::aerodynamics.airfoil = aerodynamics.airtable(joinpath(__TASOPTroot__,"airfoil_data/C.air"))

    """Inboard Wing Section (at wing root)"""
    inboard::WingSection = WingSection() # at wing root 
    """Outboard Wing Section (at strut attachment point)"""
    outboard::WingSection = WingSection() # at strut attachment point
    """Span fraction of inner wing break ("snag")"""
    ηs::Float64 = 0 

    """Fuselage lift carryover factor"""
    fuse_lift_carryover::Float64 = 0.0
    """Tip lift roll-off factor"""
    tip_lift_loss::Float64 = 0.0

    """Mean Aerodynamic Chord"""
    mean_aero_chord::Float64 = 0

    """Wing Strut"""
    has_strut::Bool = false
    strut::Strut = Strut()
    
    """Wing flap weight fraction"""
    weight_frac_flap::Float64 = 0
    """Wing slats weight fraction"""
    weight_frac_slat::Float64 = 0
    """Wing ailerons weight fraction"""
    weight_frac_ailerons::Float64 = 0
    """Wing leading_trailing_edge weight fraction"""
    weight_frac_leading_trailing_edge::Float64 = 0
    """Wing ribs weight fraction"""
    weight_frac_ribs::Float64 = 0
    """Wing spoilers weight fraction"""
    weight_frac_spoilers::Float64 = 0
    """Wing attachments weight fraction"""
    weight_frac_attachments::Float64 = 0

end

function wing_additional_weight(wing::Wing)
    return wing.weight_frac_flap + wing.weight_frac_slat + wing.weight_frac_ailerons + 
            wing.weight_frac_leading_trailing_edge + wing.weight_frac_ribs +
            wing.weight_frac_spoilers + wing.weight_frac_attachments
end

"""
"""
function Base.getproperty(obj::WingLayout, sym::Symbol)
    if sym === :ηo
        getfield(obj, :root_span)/getfield(obj, :span)
    elseif sym === :break_span
        getfield(obj, :ηs)*getfield(obj, :span)
    else
        getfield(obj, sym)
    end
    
end  # function Base.get_property