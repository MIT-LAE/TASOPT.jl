@kwdef mutable struct IndividualLandingGear
    """Weight and location"""
    weight::structures.Weight = structures.Weight()
    """Landing gear length (m)"""
    length::Float64 = 0.0
    """Number of wheels"""
    number_wheels::Int64 = 0
    """Number of shock struts"""
    number_struts::Int64 = 0
end

@kwdef mutable struct LandingGear
   
    """Front and main landing gears"""
    nose_gear::IndividualLandingGear = IndividualLandingGear()
    main_gear::IndividualLandingGear = IndividualLandingGear()

    # Misc properties
    """Design load factor"""
    load_factor::Float64 = 0.0
end