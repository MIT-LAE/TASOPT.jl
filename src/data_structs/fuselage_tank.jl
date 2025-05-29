"""
$TYPEDEF
Fuselage tank component. Usually for Hydrogen aircraft
$TYPEDFIELDS
"""
@kwdef mutable struct fuselage_tank
    """Fuel type name"""
    fueltype::String = ""
    """Fuel tank count"""
    tank_count::Int64 = 0
    """Fuel tank location"""
    placement::String = ""
    """Flag for insulation sizing"""
    sizes_insulation::Bool = false
    """Weight of fuel in one tank (N)"""
    Wfuelintank::Float64 = 0.0

    clearance_fuse::Float64 = 0.0

    """Vector with insulation layer thickness (m)"""
    t_insul::Vector{Float64} = Float64[]
    """Vector with insulation materials"""
    material_insul::Vector{ThermalInsulator} = Float64[]
    """Vector with insulation layer design indices"""
    iinsuldes::Vector{Int64} = Float64[]
    """Length of cylindrical portion of tank (m)"""
    l_cyl_inner::Float64 = 0.0
    """Length of inner tank (m)"""
    l_inner::Float64 = 0.0
    """Inner tank radius (m)"""
    Rinnertank::Float64 = 0.0
    """Vector with surface areas of insulation tank heads (m^2)"""
    Shead_insul::Vector{Float64} = Float64[]

    """Inner vessel material"""
    inner_material::StructuralAlloy = StructuralAlloy("Al-2219-T87")
    """Outer vessel material"""
    outer_material::StructuralAlloy = StructuralAlloy("Al-2219-T87")
    """Tank head aspect ratio"""
    ARtank::Float64 = 0.0
    """Angular location of inner vessel stiffeners"""
    theta_inner::Float64 = 0.0
    """Vector with angular location of outer vessel stiffeners"""
    theta_outer::Vector{Float64} = Float64[]
    """Number of intermediate stiffeners in outer vessel"""
    Ninterm::Float64 = 1.0
    
    """Venting pressure (Pa)"""
    pvent::Float64 = 0.0
    """Fill pressure (Pa)"""
    pinitial::Float64 = 0.0
    """Minimum allowable tank pressure (Pa)"""
    pmin::Float64 = 0.0
    """Departure hold time (s)"""
    t_hold_orig::Float64 = 0.0
    """Arrival hold time (s)"""
    t_hold_dest::Float64 = 0.0
    """Sea-level temperature for tank design (K)"""
    TSLtank::Vector{Float64} = []

    """Liquid fuel density (kg/m^3)"""
    rhofuel::Float64 = 0.0
    """Liquid fuel temperature in tank (K)"""
    Tfuel::Float64 = 0.0
    """Gas fuel density (kg/m^3)"""
    rhofuelgas::Float64 = 0.0
    """Fuel specific enthalpy of vaporization (J/kg)"""
    hvap::Float64 = 0.0
    """Percentage tank boiloff rate at start of cruise (%/h)"""
    boiloff_rate::Float64 = 0.0

    """Vessel additional mass fraction"""
    ftankadd::Float64 = 0.0
    """Vessel weld efficiency"""
    ew::Float64 = 0.0
    """Minimum ullage fraction"""
    ullage_frac::Float64 = 0.0
    """Heat leakage factor"""
    qfac::Float64 = 0.0
    """Pressure rise factor"""
    pfac::Float64 = 0.0
end