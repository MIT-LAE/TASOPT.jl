using DocStringExtensions
"""
$TYPEDEF

Fuselage Structure:
    Divided into 5 modules
    1. General Properties
    2. Internal Structure
    3. External Loads
    4. Fuselage Layout
    5. Misc Properties

$TYPEDFIELDS
"""
@kwdef mutable struct Fuselage
    # General Properties
    """Fuselage Weight [N] """
    weight::Float64 = 0.0
    """Fuselage Volume [m^3] """
    volume::Float64 = 0.0
    """Fuselage Weight [Nm^3] """
    moment::Float64 = 0.0

    """Fuselage Layout"""
    layout::FuselageLayout = FuselageLayout()

    # Internal Structure
    """Fuselage Material"""
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")

    """Structural Members"""
    skin::StructuralMember = StructuralMember(material=material)
    shell::StructuralMember = StructuralMember(material=material) # IS just Skin + Additional
    cone::StructuralMember = StructuralMember()

    """Internal Members"""
    floor::StructuralMember = StructuralMember()
    insulation::structures.Weight = structures.Weight()
    window::structures.Weight = structures.Weight()
    floor_W_per_area::Float64 = 0.0
    insulation_W_per_area::Float64 = 0.0
    window_W_per_length::Float64 = 0.0

    """Bending Material"""
    bendingmaterial_h::StructuralMember = StructuralMember(material=material)
    bendingmaterial_v::StructuralMember = StructuralMember(material=material)

    """External Weights"""
    APU::structures.Weight = structures.Weight()
    seat::structures.Weight = structures.Weight()
    added_payload::structures.Weight = structures.Weight()
    HPE_sys::structures.Weight = structures.Weight()
    fixed::structures.Weight = structures.Weight()

    # Misc properties
    """Number of decks in fuselage"""
    n_decks::Float64 = 0
    """Fuselage Weight fraction of stringers """
    weight_frac_stringers::Float64 = 0
    """Fuselage Weight fraction of frame """
    weight_frac_frame::Float64 = 0
    """Fuselage Weight fraction of additional weights on skin """
    weight_frac_skin_addl::Float64 = 0
    """Fuselage Shell Modulus Ratio Ebend/Eskin"""
    ratio_young_mod_fuse_bending::Float64 = 0
end

function dx_cabin(fuse::Fuselage)
    return fuse.layout.x_end_cylinder - fuse.layout.x_start_cylinder
end

# fuselage = Fuselage()


# tfweb, cabVol = TASOPT.fusew!(fuselage,13344.666000000001, 219964.5779, 76987.602265, 21996.45779, 7698.7602265000005, 0.0, 0.0, 0.0, 0.0, 0.0, "", 54911.323281976234, 435.0, 22.0, 60.0, 8607.309570000001, 8607.309570000001, 0.4, 0.7, 439929.1558, 439929.15580000007, 5.1591026057637395, 0.3, 1.0, 34.8996, 33.528, 20.165065369407014, 17.3736, 0.0, 2.1336, 36.576, 15.8496, 0.0)

# fuselage_autodiff(fuselage,13344.666000000001, 219964.5779, 76987.602265, 21996.45779, 7698.7602265000005, 0.0, 0.0, 0.0, 0.0, 0.0, "", 54911.323281976234, 435.0, 22.0, 60.0, 8607.309570000001, 8607.309570000001, 0.4, 0.7, 439929.1558, 439929.15580000007, 5.1591026057637395, 0.3, 1.0, 34.8996, 33.528, 20.165065369407014, 17.3736, 0.0, 2.1336, 36.576, 15.8496, 0.0)

# using Zygote
# function fusewAD(fuselage,i)
#     global ac
#     # Create an array of outputs
#     ac.fuselage = fuselage
#     size_aircraft!(ac)
#     oPAr = [ac.fuselage.weight]

#     return oPAr[i]
# end

# function fuselage_autodiff()
#     global ac
#     outputParams = [ac.fuselage.weight]
#     inp = ac.fuselage
#     for index = 1:length(outputParams)
#         fuse_sens = gradient((inp,i)-> fusewAD(inp,i),inp,index)
#         println(index," - ",fuse_sens)
#     end
# end

# ac = load_default_model()
# size_aircraft!(ac)
# fuselage_autodiff(ac)

# ac = load_default_model()
# size_aircraft!(ac)
# fuselage_autodiff(ac)