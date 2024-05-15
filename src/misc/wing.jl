#Wing class
mutable struct WingSection
    webs::StructuralMember
    caps::StructuralMember
    web_cap::StructuralMember
end

function WingSection(material)
    return WingSection(StructuralMember(material=material), StructuralMember(material=material), StructuralMember(material=material))
end

@kwdef mutable struct Wing
    # Layout
    layout::WingLayout = WingLayout()
    planform::Int64 = 0 # 0: Wing cantilever, plain
                      # 1: Wing cantilever with engine

    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")

    # Wing Sections
    inboard::WingSection = WingSection(material)
    outboard::WingSection = WingSection(material)

    # Strut
    strut::StructuralMember = StructuralMember(material=material)
    strut_area::Float64 = 0
    strut_length::Float64 = 0
    strut_z::Float64 = 0
    strut_toc::Float64 = 0
    strut_local_velocity_ratio::Float64 = 0
end

