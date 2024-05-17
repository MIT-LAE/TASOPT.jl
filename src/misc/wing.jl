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
end

