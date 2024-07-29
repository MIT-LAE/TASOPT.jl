# using DocStringExtensions
# """
# $TYPEDEF



# $TYPEDFIELDS
# """
@kwdef mutable struct Web_Cap_Strut
    EI_bending::Float64 = 0
    EI_normal::Float64 = 0
    GJ::Float64 = 0
end

@kwdef mutable struct WingSection
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    webs::StructuralMember = StructuralMember(material=material)
    caps::StructuralMember = StructuralMember(material=material)
    web_cap::Web_Cap_Strut = Web_Cap_Strut()
    shear_load::Float64 = 0
    moment::Float64 = 0
    cos_sweep::Float64 = 0
    weight::Float64 = 0
    dxW::Float64 = 0
    dyW::Float64 = 0
end

@kwdef mutable struct Strut
    # Strut
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    S::Float64 = 0
    length::Float64 = 0
    z::Float64 = 0
    toc::Float64 = 0
    local_velocity_ratio::Float64 = 0
    weight::Float64 = 0
    axial_force::Float64 = 0
    chord::Float64 = 0
    dxW::Float64 = 0
end

@kwdef mutable struct Wing
    weight::Float64 = 0
    dxW::Float64 = 0
    # Layout
    layout::WingLayout = WingLayout()
    planform::Int64 = 0 # 0: Wing cantilever, plain
                      # 1: Wing cantilever with engine

    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")

    # Wing Sections
    inboard::WingSection = WingSection(material=material) # at wing root 
    outboard::WingSection = WingSection(material=material) # at strut attachment point
    weight_center::Float64 = 0

    # Strut
    strut::Strut = Strut()

    # Weight fractions
    weight_frac_flap::Float64 = 0
    weight_frac_slat::Float64 = 0
    weight_frac_ailerons::Float64 = 0
    weight_frac_leading_trailing_edge::Float64 = 0
    weight_frac_ribs::Float64 = 0
    weight_frac_spoilers::Float64 = 0
    weight_frac_attachments::Float64 = 0

end

function wing_additional_weight(wing::Wing)
    return wing.weight_frac_flap + wing.weight_frac_slat + wing.weight_frac_ailerons + 
            wing.weight_frac_leading_trailing_edge + wing.weight_frac_ribs +
            wing.weight_frac_spoilers + wing.weight_frac_attachments
end