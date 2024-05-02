#Wing class
@kwdef mutable struct Wing
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    webs::StructuralMember = StructuralMember(material=material)
    caps::StructuralMember = StructuralMember(material=material)
    strut::StructuralMember = StructuralMember(material=material)
    web_cap_strut::StructuralMember = StructuralMember(material=material)
    web_cap_root::StructuralMember = StructuralMember(material=material)

    strut_area::Float64 = 0
    strut_length::Float64 = 0
    web_thickness_strut_attach_point::Float64 = 0
    cap_thickness_strut_attach_point::Float64 = 0
end