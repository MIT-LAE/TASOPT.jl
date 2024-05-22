#Wing class
mutable struct WingSection
    webs::StructuralMember
    caps::StructuralMember
    web_cap::StructuralMember
end

function WingSection(material)
    return WingSection(StructuralMember(material=material), StructuralMember(material=material), StructuralMember(material=material))
end

@kwdef mutable struct Strut
    # Strut
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    area::Float64 = 0
    length::Float64 = 0
    z::Float64 = 0
    toc::Float64 = 0
    local_velocity_ratio::Float64 = 0
    weight::Float64 = 0
end

@kwdef mutable struct Wing
    weight::Float64 = 0
    # Layout
    layout::WingLayout = WingLayout()
    planform::Int64 = 0 # 0: Wing cantilever, plain
                      # 1: Wing cantilever with engine

    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")

    # Wing Sections
    inboard::WingSection = WingSection(material) # at wing root 
    outboard::WingSection = WingSection(material) # at strut attachment point

    # Strut
    strut::Strut = Strut()

    # Weight fractions
    flap_weight_frac::Float64 = 0
    slat_weight_frac::Float64 = 0
    aileron_weight_frac::Float64 = 0
    leading_trailing_edge_weight_frac::Float64 = 0
    ribs_weight_frac::Float64 = 0
    spoilers_weight_frac::Float64 = 0
    attachments_weight_frac::Float64 = 0

end

function wing_additional_weight(wing::Wing)
    return wing.flap_weight_frac + wing.slat_weight_frac + wing.aileron_weight_frac + 
            wing.leading_trailing_edge_weight_frac + wing.ribs_weight_frac +
            wing.spoilers_weight_frac + wing.attachments_weight_frac
end