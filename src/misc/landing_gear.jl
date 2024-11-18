@kwdef mutable struct IndividualLandingGear
    """Weight and location"""
    weight::structures.Weight = structures.Weight()
    """Landing gear length (m)"""
    length::Float64 = 0.0
    """Number of wheels"""
    number_wheels::Int64 = 0
    """Number of shock struts"""
    number_struts::Int64 = 0
    """Wing-to-landing gear distance, for main gear (m)"""
    distance_wing_to_landing_gear::Float64 = 0.0
end

function Base.getproperty(obj::IndividualLandingGear, sym::Symbol)
    if sym === :moment
        # Use `getfield` to directly access fields of `Weight`
        weight = getfield(obj, :weight)
        return getfield(weight, :W) * getfield(weight, :r)[1]

    else
        # Use `getfield` to directly access fields of `IndividualLandingGear`
        return getfield(obj, sym)
    end
end

@kwdef mutable struct LandingGear
   
    """Front and main landing gears"""
    nose_gear::IndividualLandingGear = IndividualLandingGear()
    main_gear::IndividualLandingGear = IndividualLandingGear()

    # Misc properties
    """Design load factor"""
    load_factor::Float64 = 0.0
end