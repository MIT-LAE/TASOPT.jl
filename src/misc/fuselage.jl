#Fuselage class
@kwdef mutable struct Fuselage
    # General Fuselage Properties
    weight::Float64 = 0.0
    volume::Float64 = 0.0
    # sigma_factor::Float64 = 1.0
    moment::Float64 = 0.0

    # Layout
    layout::FuselageLayout = FuselageLayout()

    # Structures
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    #to do: Add material-> Structural members
    skin::StructuralMember = StructuralMember(material=material)
    shell::StructuralMember = StructuralMember(material=material) # IS just Skin + Additional
    cone::StructuralMember = StructuralMember()
    floor::StructuralMember = StructuralMember()
    insulation::StructuralMember = StructuralMember()
    window::StructuralMember = StructuralMember()
    bending_h::StructuralMember = StructuralMember(material=material)
    bending_v::StructuralMember = StructuralMember(material=material)

    # Loads
    
    # Misc properties
    # Nland::Float64 = 6.0
    n_decks::Float64 = 0
    weight_frac_string::Float64 = 0
    weight_frac_frame::Float64 = 0
    # ffadd::Float64 = 0.2
    # nftanks::Int64 = 1
    rEshell::Float64 = 0
    # tank_placement
end

function dx_cabin(fuse::Fuselage)
    return fuse.layout.x_end_cylinder - fuse.layout.x_start_cylinder
end

# fuselage = Fuselage()


# tfweb, cabVol = TASOPT.fusew!(fuselage,13344.666000000001, 219964.5779, 76987.602265, 21996.45779, 7698.7602265000005, 0.0, 0.0, 0.0, 0.0, 0.0, "", 54911.323281976234, 435.0, 22.0, 60.0, 8607.309570000001, 8607.309570000001, 0.4, 0.7, 439929.1558, 439929.15580000007, 5.1591026057637395, 0.3, 1.0, 34.8996, 33.528, 20.165065369407014, 17.3736, 0.0, 2.1336, 36.576, 15.8496, 0.0)

# fuselage_autodiff(fuselage,13344.666000000001, 219964.5779, 76987.602265, 21996.45779, 7698.7602265000005, 0.0, 0.0, 0.0, 0.0, 0.0, "", 54911.323281976234, 435.0, 22.0, 60.0, 8607.309570000001, 8607.309570000001, 0.4, 0.7, 439929.1558, 439929.15580000007, 5.1591026057637395, 0.3, 1.0, 34.8996, 33.528, 20.165065369407014, 17.3736, 0.0, 2.1336, 36.576, 15.8496, 0.0)


# function fusewAD(fuselage,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
#     Waftfuel, Wftank, ltank, xftankaft, tank_placement,
#     deltap,
#     Wpwindow,Wppinsul,Wppfloor, 
#     Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
#     bv,lambdav,nvtail,xhtail,xvtail,
#     xwing,xwbox,cbox,
#     xfix,xapu,xeng,xfuel,i)
#     # Create an array of outputs
#     tfweb, cabVol = TASOPT.fusew!(fuselage,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
#     Waftfuel, Wftank, ltank, xftankaft, tank_placement,
#     deltap,
#     Wpwindow,Wppinsul,Wppfloor, 
#     Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
#     bv,lambdav,nvtail,xhtail,xvtail,
#     xwing,xwbox,cbox,
#     xfix,xapu,xeng,xfuel)
#     oPAr = [fuselage.weight]

#     return oPAr[i]
# end

# function fuselage_autodiff(fuselage,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
#     Waftfuel, Wftank, ltank, xftankaft, tank_placement,
#     deltap,
#     Wpwindow,Wppinsul,Wppfloor, 
#     Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
#     bv,lambdav,nvtail,xhtail,xvtail,
#     xwing,xwbox,cbox,
#     xfix,xapu,xeng,xfuel)
#     outputParams = [fuselage.weight]
#     for index = 1:length(outputParams)
#         fuse_sens = gradient((fuselage,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
#         Waftfuel, Wftank, ltank, xftankaft, tank_placement,
#         deltap,
#         Wpwindow,Wppinsul,Wppfloor, 
#         Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
#         bv,lambdav,nvtail,xhtail,xvtail,
#         xwing,xwbox,cbox,
#         xfix,xapu,xeng,xfuel,i)-> fusewAD(fuselage,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
#         Waftfuel, Wftank, ltank, xftankaft, tank_placement,
#         deltap,
#         Wpwindow,Wppinsul,Wppfloor, 
#         Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
#         bv,lambdav,nvtail,xhtail,xvtail,
#         xwing,xwbox,cbox,
#         xfix,xapu,xeng,xfuel,i),fuselage,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
#         Waftfuel, Wftank, ltank, xftankaft, tank_placement,
#         deltap,
#         Wpwindow,Wppinsul,Wppfloor, 
#         Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
#         bv,lambdav,nvtail,xhtail,xvtail,
#         xwing,xwbox,cbox,
#         xfix,xapu,xeng,xfuel,index)
#         println(index," - ",fuse_sens)
#     end
# end


    