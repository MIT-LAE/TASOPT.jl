using EnumX
using ..aerodynamics
abstract type AbstractWing end

"""
$TYPEDEF

Wing Strut

$TYPEDFIELDS
"""
@kwdef mutable struct Strut
    """Strut box caps"""
    caps::StructuralMember = StructuralMember(material = StructuralAlloy("TASOPT-Al"))
    """Strut box webs"""
    webs::StructuralMember = StructuralMember(material = StructuralAlloy("TASOPT-Al"))

    # Geometry inputs
    """Planform area of both struts [m^2]"""
    S::Float64 = 0
    """Strut z position [m]"""
    z::Float64 = 0
    """Strut y position at fuselage [m]"""
    y_root::Float64 = 0
    """Strut thickness to chord"""
    thickness_to_chord::Float64 = 0
    """Strut chord to wing chord at the strut attachment (placeholder default)"""
    chord_to_wing_chord::Float64 = 0.5
    """Strut chord taper ratio, wing end to fuselage end (0 < taper_ratio <= 1)"""
    taper_ratio::Float64 = 1.0
    """Strut box width to chord ratio"""
    box_width_to_chord::Float64 = 0.5
    """Strut box web height to box height ratio"""
    web_to_box_height::Float64 = 0.75
    """Euler effective-length factor"""
    end_fixity_K::Float64 = 1.0
    """Minimum strut cap thickness [m]"""
    min_cap_thickness::Float64 = 1.5e-3
    """Minimum strut web thickness [m]"""
    min_web_thickness::Float64 = 1.0e-3
    """Strut local velocity ratio"""
    local_velocity_ratio::Float64 = 0
    """Strut force model"""
    opt_strut_force::StrutForceModel.T = StrutForceModel.Compatibility
    """Strut mean section cl over wing mean section cl (0 = no strut lift)"""
    mean_cl_ratio_to_wing::Float64 = 0.33
    """Strut spanload exponent ``p``, ``\\Gamma \\propto (1-\\xi^2)^p`` with ``\\xi`` = 0 at fuselage, 1 at wing (``p = 0.5`` is elliptic)"""
    spanload_exponent::Float64 = 0.0
    """Fuselage x of the strut attach [m]. NaN follows the wing box, `wing.layout.box_x`"""
    x_strut_attach::Float64 = NaN

    # Computed geometry
    """Strut Length [m] """
    length::Float64 = 0
    """Strut chord at the wing end [m]"""
    chord::Float64 = 0
    """Strut chord at the fuselage end [m]"""
    root_chord::Float64 = 0
    """Strut caps and webs material area [m^2]"""
    area::Float64 = 0
    """Cosine of strut sweep angle"""
    cos_lambda::Float64 = 0 
    """Strut dihedral angle [rad]"""
    dihedral::Float64 = 0

    # Loads
    """Maximum strut tension force [N]"""
    axial_force::Float64 = 0
    """Strut axial force in the most compressive case [N] (<= 0)"""
    axial_force_compression::Float64 = 0
    """Strut web shear load, jury load plus strut own weight [N]"""
    shear_load::Float64 = 0
    """Strut midspan bending moment from own weight [N m]"""
    bending_moment::Float64 = 0
    """Strut Euler buckling load P_cr [N]"""
    buckling_load::Float64 = 0
    """Strut minor-axis bending stiffness [N m^2]"""
    EI_min::Float64 = 0
    """Governing strut sizing criterion: `:tension`, `:buckling` or `:gauge`"""
    sizing_criterion::Symbol = :none

    """Vertical force of both struts on the fuselage [N]"""
    force_vertical::Float64 = 0
    """Horizontal force of one strut on the fuselage [N]"""
    force_horizontal::Float64 = 0
    """Share of vertical wing support taken by the strut, in the largest-strut-force case"""
    reaction_share::Float64 = 0

    # Aerodynamics
    """Lift coefficient of both struts, referenced to wing area"""
    CL::Float64 = 0
    """Strut lift over total lift"""
    lift_fraction::Float64 = 0

    # Weights
    """Strut weight [N]"""
    weight::Float64 = 0
    """Strut added weight fraction (fairings, edges, fittings) on the bare box weight. Placeholder from Boeing SUGAR Phase IV"""
    weight_frac_added::Float64 = 0.60
    """Strut attachment weight fraction (twice the wing default, as it attaches at both ends)"""
    weight_frac_attachments::Float64 = 0.06
    """Aircraft pitching moment contribution from the weight distribution of the strut [N m]"""
    dxW::Float64 = 0
end

"""
$TYPEDEF

Vertical jury strut of a truss-braced wing, connecting the wing to the main strut.
It shortens the buckling length of the main strut but does not change the wing or strut loads.

$TYPEDFIELDS
"""
@kwdef mutable struct JuryStrut
    """Jury strut box caps"""
    caps::StructuralMember = StructuralMember(material = StructuralAlloy("TASOPT-Al"))
    """Jury strut box webs"""
    webs::StructuralMember = StructuralMember(material = StructuralAlloy("TASOPT-Al"))
    # Input
    """Jury position as a fraction of the main strut length from the fuselage end"""
    eta_strut::Float64 = 0.5
    """Jury strut chord to main strut chord"""
    chord_to_strut_chord::Float64 = 0.5
    """Minimum jury strut cap thickness [m]"""
    min_cap_thickness::Float64 = 1.5e-3
    """Minimum jury strut web thickness [m]"""
    min_web_thickness::Float64 = 1.0e-3

    # Computed geometry
    """Wing-side attachment as a fraction of the semispan: `ηo + eta_strut (ηs - ηo)`"""
    eta_wing::Float64 = 0
    """Jury strut length [m]"""
    length::Float64 = 0
    """Jury strut chord [m]"""
    chord::Float64 = 0
    """Jury strut caps and webs material area [m^2]"""
    area::Float64 = 0

    # Loads
    """Jury strut axial force [N] (compression is negative)"""
    axial_force::Float64 = 0
    """Jury strut buckling load [N]"""
    buckling_load::Float64 = 0

    # Weights
    """Jury strut weight [N], both wings"""
    weight::Float64 = 0
    """Jury strut added weight fraction (fairings, fittings) on the bare weight. Placeholder, same as the strut"""
    weight_frac_added::Float64 = 0.60
    """Aircraft pitching moment contribution from the weight distribution of the jury strut [N m]"""
    dxW::Float64 = 0
end

"""
Mass state of a wing sizing load case.

- `DesignCruise`: start-of-cruise weight
- `MTOW_MaxFuel`: maximum takeoff weight, wing tanks as full as possible
- `MZFW`: maximum zero-fuel weight, reserve wing fuel only
- `Light`: operating empty weight plus reserve fuel
"""
@enumx WingMassState DesignCruise MTOW_MaxFuel MZFW Light

"""
How the load factor of a wing sizing load case is obtained.

- `Maneuver`: prescribed limit load factor in flight
- `Gust`: from a discrete gust; `load_factor` gives the sign only
- `Ground`: prescribed load factor on the ground, no wing lift
"""
@enumx WingLoadKind Maneuver Gust Ground

"""
$TYPEDEF

One wing sizing load case.

$TYPEDFIELDS
"""
@kwdef mutable struct WingLoadCase
    """How the load factor is obtained"""
    kind::WingLoadKind.T = WingLoadKind.Maneuver
    """Limit load factor (sign only for `Gust`)"""
    load_factor::Float64 = 0.0
    """Design gust velocity U_de [m/s EAS] (56 ft/s, FAR 25.341(a)(5))"""
    gust_velocity_EAS::Float64 = 17.0688
    """Gust flight speed [m/s EAS]; 0 uses the mission cruise speed"""
    gust_speed_EAS::Float64 = 0.0
    """Weight and wing fuel state"""
    mass_state::WingMassState.T
    """Tail load scale: tail_load_scale · (load factor / pull-up load factor) · maximum tail load"""
    tail_load_scale::Float64 = 1.0
    """Whether the case is used in wing sizing"""
    active::Bool = true
    """Load factor used (gust: 1 ± Δn)"""
    n_computed::Float64 = 0.0
end

"""
$TYPEDEF

Wing root and break loads of one sizing load case.

$TYPEDFIELDS
"""
@kwdef mutable struct WingCaseLoads
    """Shear load at the wing root [N]"""
    S_root::Float64 = 0
    """Bending moment at the wing root [N m]"""
    M_root::Float64 = 0
    """Shear load at the wing break [N]"""
    S_break::Float64 = 0
    """Bending moment at the wing break [N m]"""
    M_break::Float64 = 0
    """Largest bending moment along the inboard bay of a braced wing [N m]"""
    M_inboard_extreme::Float64 = 0
    """Inboard wing axial force from the strut, tension positive [N]"""
    P_axial_inboard::Float64 = 0
    """Strut axial force at the fuselage end, tension positive [N]"""
    strut_force::Float64 = 0
    """Distance from the root of the 11 stations along the inboard bay [m]; empty without a strut"""
    inboard_stations::Vector{Float64} = Float64[]
    """Shear load at `inboard_stations`, with the strut reaction applied [N]"""
    inboard_shear::Vector{Float64} = Float64[]
    """Bending moment at `inboard_stations`, with the strut reaction applied [N m]"""
    inboard_moment::Vector{Float64} = Float64[]
    """Vertical strut reaction on the wing, against the wing lift [N]"""
    strut_reaction::Float64 = 0
    """Lift on one strut [N]"""
    strut_lift::Float64 = 0
end

"""
$TYPEDEF

The fixed set of wing sizing load cases, with their results.
Only the `pullup` case is active by default.

$TYPEDFIELDS
"""
@kwdef mutable struct WingSizingLoads
    """Pull-up, design cruise mass state"""
    pullup::WingLoadCase = WingLoadCase(load_factor = 3.0, mass_state = WingMassState.DesignCruise)
    """Push-over, MTOW with full wing fuel"""
    pushover::WingLoadCase = WingLoadCase(load_factor = -1.0, mass_state = WingMassState.MTOW_MaxFuel, active = false)
    """Pull-up, light mass state"""
    pullup_light::WingLoadCase = WingLoadCase(load_factor = 3.0, mass_state = WingMassState.Light, active = false)
    """Push-over, light mass state"""
    pushover_light::WingLoadCase = WingLoadCase(load_factor = -1.0, mass_state = WingMassState.Light, active = false)
    """Taxi, MTOW with full wing fuel"""
    taxi::WingLoadCase = WingLoadCase(kind = WingLoadKind.Ground, load_factor = -2.0, mass_state = WingMassState.MTOW_MaxFuel, active = false)
    """Upward discrete gust, light mass state"""
    gust_up_light::WingLoadCase = WingLoadCase(kind = WingLoadKind.Gust, load_factor = 1.0, mass_state = WingMassState.Light, active = false)
    """Downward discrete gust, light mass state"""
    gust_down_light::WingLoadCase = WingLoadCase(kind = WingLoadKind.Gust, load_factor = -1.0, mass_state = WingMassState.Light, active = false)

    """Loads of each case, keyed by the case field name"""
    results::Dict{Symbol,WingCaseLoads} = Dict{Symbol,WingCaseLoads}()
    """Case that sets the wing root loads"""
    governing_root::Symbol = :pullup
    """Case that sets the wing break loads"""
    governing_break::Symbol = :pullup
    """Case with the largest strut tension (`:none` without a strut)"""
    governing_strut::Symbol = :none
    """Case with the most negative strut force (`:none` if none)"""
    governing_strut_compression::Symbol = :none
end

"""
$TYPEDEF

The Wing structure is composed of 6 sub-structures as follow and are visualized [here](../assets/wing_struct.png).

1. General Properties
2. Wing Layout
3. Material
4. Wing Sections
5. Strut
6. Weight Fractions

See [`structures.WingSection`](@ref), [`structures.WingLayout`](@ref) or [`structures.WingCrossSection`](@ref) for more detail.

$TYPEDFIELDS
"""
@kwdef mutable struct Wing{air<:aerodynamics.airfoil} <: AbstractWing
    """Wing Weight [N] """
    weight::Float64 = 0.0
    """Aircraft pitching moment contribution from the weight distribution of the wing [Nm]"""
    dxW::Float64 = 0.0
    """Wing Layout """
    layout::WingLayout = WingLayout()
    """Wing Material """
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")

    """Airfoil data"""
    airsection::air = aerodynamics.airtable(joinpath(__TASOPTroot__,"airfoil_data/C.air"))

    """Inboard Wing Section (at wing root)"""
    inboard::WingSection = WingSection() # at wing root 
    """Outboard Wing Section (at strut attachment point)"""
    outboard::WingSection = WingSection() # at strut attachment point
    """Engine spanwise location, ``\\eta = y/(b/2)``"""
    η_engine::Float64 = 0.0

    """Span fraction of inner wing break ("snag")"""
    # ηs::Float64 = 0 

    """Fuselage lift carryover factor"""
    fuse_lift_carryover::Float64 = 0.0
    """Tip lift roll-off factor"""
    tip_lift_loss::Float64 = 0.0

    """Mean Aerodynamic Chord"""
    mean_aero_chord::Float64 = 0.0

    """Wing Strut"""
    strut::Strut = Strut()

    """Vertical jury strut of a truss-braced wing"""
    jury::JuryStrut = JuryStrut()

    """Wing sizing load cases and their results"""
    sizing_loads::WingSizingLoads = WingSizingLoads()

    """Wing flap weight fraction"""
    weight_frac_flap::Float64 = 0.0
    """Wing slats weight fraction"""
    weight_frac_slat::Float64 = 0.0
    """Wing ailerons weight fraction"""
    weight_frac_ailerons::Float64 = 0.0
    """Wing leading_trailing_edge weight fraction"""
    weight_frac_leading_trailing_edge::Float64 = 0.0
    """Wing ribs weight fraction"""
    weight_frac_ribs::Float64 = 0.0
    """Wing spoilers weight fraction"""
    weight_frac_spoilers::Float64 = 0.0
    """Wing attachments weight fraction"""
    weight_frac_attachments::Float64 = 0.0

end

"""
$TYPEDEF

Tail

$TYPEDFIELDS
"""
@kwdef mutable struct Tail <: AbstractWing
    """Tail Layout """
    layout::WingLayout = WingLayout()
    """Tail Sections """
    outboard::WingSection= WingSection()
    inboard::WingSection = WingSection()
    """Tip lift roll-off factor"""
    tip_lift_loss::Float64 = 0.0
    """Aircraft pitching moment contribution from the weight distribution of the strut [N m]"""
    dxW::Float64 = 0
    """Tail Weight [N] """
    weight::Float64 = 0
    """Tail Added Weight Fraction"""
    weight_fraction_added::Float64 = 0
    """Tail Max CL """
    CL_max::Float64 = 0
    """Tail Volume [m^3] """
    volume::Float64 = 0
    """Tail Sizing assumption selection - different for HTail vs VTail """
    opt_sizing::TailSizing.T = TailSizing.FixedVh
    """Tail Downwash factor dε/dα """
    downwash_factor::Float64 = 0
    """Tail max fwd CG (only used if opt_sizing == "CLmax_fwdCG" for HTail) """
    CL_max_fwd_CG::Float64 = 0
    """Tail Minimum static margin"""
    SM_min::Float64 = 0
    """Max Tail down load. Tail download param at max load case"""
    CL_CLmax::Float64 = 0
    """Number of Tails"""
    ntails::Float64 = 0
end

function wing_additional_weight(wing::AbstractWing)
    if wing isa Wing
        return wing.weight_frac_flap + wing.weight_frac_slat + wing.weight_frac_ailerons + 
        wing.weight_frac_leading_trailing_edge + wing.weight_frac_ribs +
        wing.weight_frac_spoilers + wing.weight_frac_attachments
    elseif wing isa Tail
        return wing.weight_fraction_added
    end
    
end


"""
"""
function Base.getproperty(obj::AbstractWing, sym::Symbol)
    concrete_type = typeof(obj)  # Gets the actual type (Wing or Tail)
    if hasfield(concrete_type, sym)
        return getfield(obj, sym)
    elseif hasfield(WingLayout, sym)
        return getfield(obj.layout, sym)
    else
        throw(KeyError("Property $sym not found in $(concrete_type) or WingLayout"))
    end
end