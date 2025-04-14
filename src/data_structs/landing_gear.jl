using DocStringExtensions

"""
$TYPEDEF
Individual landing gear (nose or main): this object stores the geometric and mass properties of
a landing gear.
$TYPEDFIELDS
"""
@kwdef mutable struct IndividualLandingGear
    """Weight and location"""
    weight::Weight = Weight()
    """Landing gear length (m)"""
    length::Float64 = 0.0
    """Number of shock struts"""
    number_struts::Int64 = 0
    """Number of wheels per strut"""
    wheels_per_strut::Int64 = 0
    """CG-to-landing gear distance, for main gear (m)"""
    distance_CG_to_landing_gear::Float64 = 0.0
    """y-offset as fraction of halfspan"""
    y_offset_halfspan_fraction::Float64 = 0.0
    """Overall mass fraction of MTOW"""
    overall_mass_fraction::Float64 = 0.0
end

function Base.getproperty(obj::IndividualLandingGear, sym::Symbol)
    if sym === :moment
        return structures.y_moment(getfield(obj, :weight))

    else
        # Use `getfield` to directly access fields of `IndividualLandingGear`
        return getfield(obj, sym)
    end
end

"""
$TYPEDEF
Landing Gear structure:
    Divided into 2 modules
    1. Main and nose gear objects
    2. Misc Properties
$TYPEDFIELDS
"""
@kwdef mutable struct LandingGear
    model::String = ""

    """Front and main landing gears"""
    nose_gear::IndividualLandingGear = IndividualLandingGear()
    main_gear::IndividualLandingGear = IndividualLandingGear()

    #Miscellaneous
    """Tailstrike angle (rad)"""
    tailstrike_angle::Float64 = 0.0
    """Wing dihedral angle (rad)"""
    wing_dihedral_angle::Float64 = 0.0 #TODO this should be a wing param, but it does not appear to be used elsewhere
    """Engine ground clearance (m)"""
    engine_ground_clearance::Float64 = 0.0
end