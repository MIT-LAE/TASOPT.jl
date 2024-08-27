using DocStringExtensions
abstract type AbstractComponent end
abstract type AbstractTurbomachine end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a conventional
single-bubble type fuselage.

$TYPEDFIELDS
"""
@kwdef mutable struct EngineElement <: AbstractComponent
    """Flag for inclusion"""
    included::Bool = false
    """Pressure ratio"""
    pressure_ratio::Float64 = 0.0
end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a conventional
single-bubble type fuselage.

$TYPEDFIELDS
"""
@kwdef mutable struct Turbomachine <: AbstractTurbomachine
    """Flag for inclusion"""
    included::Bool = false
    """Pressure ratio"""
    pressure_ratio::Vector{Float64} = []
    """Polytropic efficiency"""
    polytropic_efficiency::Vector{Float64} = []
    """Corrected mass flow rate"""
    corrected_mass_flow_rate::Vector{Float64} = []
    """Design pressure ratio"""
    design_pressure_ratio::Float64 = 0.0
    """Design polytropic efficiency"""
    design_polytropic_efficiency::Float64 = 0.0
    """Deisgn corrected mass flow rate"""
    design_corrected_mass_flow_rate::Float64 = 0.0
end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a conventional
single-bubble type fuselage.

$TYPEDFIELDS
"""
@kwdef mutable struct EngineComponents
    """Diffuser"""
    diffuser::AbstractComponent = EngineElement()
    """Fan"""
    fan::AbstractTurbomachine = Turbomachine()
    """Low-pressure compressor"""
    LPC::AbstractTurbomachine = Turbomachine()
    """High-pressure compressor"""
    HPC::AbstractTurbomachine = Turbomachine()
    """Burner"""
    burner::AbstractComponent = EngineElement()
    """High-pressure turbine"""
    HPT::AbstractTurbomachine = Turbomachine()
    """Low-pressure turbine"""
    LPT::AbstractTurbomachine = Turbomachine()
    """Fan nozzle"""
    fan_nozzle::AbstractComponent = EngineElement()
    """Turbine nozzle"""
    turbine_nozzle::AbstractComponent = EngineElement()
end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a conventional
single-bubble type fuselage.

$TYPEDFIELDS
"""
@kwdef mutable struct FlowStation
    """Total temperature [K]"""
    Tt::Float64 = 0.0
    """Total specific enthalpy [J/kg]"""
    ht::Float64 = 0.0
    """Total pressure [Pa]"""
    pt::Float64 = 0.0
    """Total specific heat at constant pressure [J/(kg K)]"""
    cpt::Float64 = 0.0
    """Gas constant [J/(kg K)]"""
    R::Float64 = 0.0
    """Mach number"""
    M::Float64 = 0.0
end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a conventional
single-bubble type fuselage.

$TYPEDFIELDS
"""
@kwdef mutable struct EngineStation
    """Flag for inclusion"""
    included::Bool = false
    """Cross-sectional area [m^2]"""
    area::Float64 = 0.0
    """Flow properties"""
    flow::Vector{FlowStation} = []
end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a conventional
single-bubble type fuselage.

$TYPEDFIELDS
"""
@kwdef mutable struct EngineStations
    """Freestream"""
    S0::EngineStation = EngineStation()
    """Fan face outside of casing BLs"""
    S18::EngineStation = EngineStation()
    """Fan face over LPC portion"""
    S19::EngineStation = EngineStation()
    """Fan face over fan portion"""
    S2::EngineStation = EngineStation()
    """Fan exit in bypass stream"""
    S21::EngineStation = EngineStation()
    """Fan precooler outlet, fan+LPC inlet"""
    S19c::EngineStation = EngineStation()
    """LPC exit, intercooler inlet """
    S25::EngineStation = EngineStation()
    """Intercooler exit, HPC inlet"""
    S25c::EngineStation = EngineStation()
    """Compressor exit"""
    S3::EngineStation = EngineStation()
    """Combustor exit before cooling air addition"""
    S4::EngineStation = EngineStation()
    """Turbine  inlet after  cooling air addition"""
    S41::EngineStation = EngineStation()
    """HPT exit, LPT inlet"""
    S45::EngineStation = EngineStation()
    """LPT exit, regenerative cooler inlet"""
    S49::EngineStation = EngineStation()
    """Regenerative cooler outlet"""
    S49c::EngineStation = EngineStation()
    """Core nozzle"""
    S5::EngineStation = EngineStation()
    """Core flow downstream"""
    S6::EngineStation = EngineStation()
    """Fan nozzle"""
    S7::EngineStation = EngineStation()
    """Fan flow downstream"""
    S8::EngineStation = EngineStation()
end

"""
$TYPEDEF

Engine Structure:
    Divided into 5 modules
    1. General Properties
    2. Internal Structure
    3. External Loads
    4. Fuselage Layout
    5. Misc Properties

$TYPEDFIELDS
"""
@kwdef mutable struct Engine
    # General Properties
    """Engine Weight [N] """
    weight::Float64 = 0.0
    """Engine Volume [m^3] """
    volume::Float64 = 0.0
    """Engine Weight-moment [Nm] """
    moment::Float64 = 0.0

    # """Fuel properties"""
    # fuel::Fuel = Fuel()

    """Stations"""
    station::EngineStations = EngineStations()

    """Components"""
    component::EngineComponents = EngineComponents()

    # """Performance"""
    # performance::EnginePerformance
end