#Wing class
@kwdef mutable struct Wing
    # Layout
    layout::WingLayout = WingLayout()

    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")

    # Wing Sections
    inboard::WingSection = WingSection(material)
    outboard::WingSection = WingSection(material)

    # Strut
    strut::StructuralMember = StructuralMember(material=material)
    strut_area::Float64 = 0
    strut_length::Float64 = 0
end
mutable struct WingSection
    webs::StructuralMember = StructuralMember(material=material)
    caps::StructuralMember = StructuralMember(material=material)
    web_cap::StructuralMember = StructuralMember(material=material)
end

function WingSection(material)
    return WingSection(StructuralMember(material=material), StructuralMember(material=material), StructuralMember(material=material))
end
