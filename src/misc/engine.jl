using DocStringExtensions
abstract type AbstractComponent end
abstract type AbstractTurbomachine end
iptotal = 17
"""
$TYPEDEF

Contains properties of engine elements that deliver no heat or work to the flow
but produce a pressure drop.

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

Contains properties of the burner (combustor).

$TYPEDFIELDS
"""
@kwdef mutable struct Burner <: AbstractComponent
    """Flag for inclusion"""
    included::Bool = false
    """Pressure ratio"""
    pressure_ratio::Float64 = 0.0
    """Combustion efficiency"""
    combustion_efficiency::Float64 = 0.0
end

"""
$TYPEDEF

Contains the heat and pressure drop properties of a heat exchanger.

$TYPEDFIELDS
"""
@kwdef mutable struct HeatExchanger <: AbstractComponent
    """Flag for inclusion"""
    included::Bool = false
    """Pressure drop [Pa]"""
    pressure_drop::Vector{Float64} = zeros(iptotal)
    """Specific enthalpy rise [J/kg]"""
    enthalpy_rise::Vector{Float64} =  zeros(iptotal)
end

"""
$TYPEDEF

Contains the physical properties of turbomachines doing or exctracting work from the flow. 
It can be used for fans, compressors, or turbines.

$TYPEDFIELDS
"""
@kwdef mutable struct Turbomachine <: AbstractComponent
    """Flag for inclusion"""
    included::Bool = false
    """Pressure ratio"""
    pressure_ratio::Vector{Float64} = zeros(iptotal)
    """Polytropic efficiency"""
    polytropic_efficiency::Vector{Float64} = zeros(iptotal)
    """Corrected mass flow rate"""
    corrected_mass_flow_rate::Vector{Float64} = zeros(iptotal)
    """Design pressure ratio"""
    corrected_speed::Vector{Float64} = zeros(iptotal)
    """Isentropic efficiency"""
    isentropic_efficiency::Vector{Float64} =zeros(iptotal)
    """Design pressure ratio"""
    design_pressure_ratio::Float64 = 0.0
    """Design polytropic efficiency"""
    design_polytropic_efficiency::Float64 = 0.0
    """Deisgn corrected mass flow rate"""
    design_corrected_mass_flow_rate::Float64 = 0.0
    """Deisgn corrected speed"""
    design_corrected_speed::Float64 = 0.0
end

"""
$TYPEDEF

Contains the properties of all the components in an engine.

$TYPEDFIELDS
"""
@kwdef mutable struct EngineComponents
    """Diffuser"""
    diffuser::AbstractComponent = EngineElement()
    """Fan"""
    fan::AbstractComponent = Turbomachine()
    """Fan precooler"""
    fan_precooler::AbstractComponent = HeatExchanger()
    """Low-pressure compressor"""
    LPC::AbstractComponent = Turbomachine()
    """Compressor intercooler"""
    compressor_intercooler::AbstractComponent = HeatExchanger()
    """High-pressure compressor"""
    HPC::AbstractComponent = Turbomachine()
    """Burner"""
    burner::AbstractComponent = Burner()
    """High-pressure turbine"""
    HPT::AbstractComponent = Turbomachine()
    """Low-pressure turbine"""
    LPT::AbstractComponent = Turbomachine()
    """Regenerative cooler"""
    regenerative_cooler::AbstractComponent = HeatExchanger()
    """Fan nozzle"""
    fan_nozzle::AbstractComponent = EngineElement()
    """Turbine nozzle"""
    turbine_nozzle::AbstractComponent = EngineElement()
end

"""
$TYPEDEF

Contains the physical properties of the flow at a given station.

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

Contains all the properties of a given engine station, including area and flow conditions.

$TYPEDFIELDS
"""
@kwdef mutable struct EngineStation
    """Flag for inclusion"""
    included::Bool = false
    """Cross-sectional area [m^2]"""
    area::Float64 = 0.0
    """Flow properties"""
    flow::Vector{FlowStation} = fill(FlowStation(), iptotal)
end

"""
$TYPEDEF

Contains all the engine stations.

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
    """Fan precooler outlet, core fan+LPC inlet"""
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

Contains the parameters that describe the engine performance, such as thrust.

$TYPEDFIELDS
"""
@kwdef mutable struct EnginePerformance
    """Thrust [N]"""
    thrust::Float64 = 0.0
    """Specific thrust"""
    specific_thrust::Float64 = 0.0
    """Thrust-specific fuel consumption [kg/(N s)]"""
    TSFC::Float64 = 0.0
    """Thrust-specific energy consumption [W/N]"""
    TSEC::Float64 = 0.0
    """Core mass flow rate [kg/s]"""
    core_mass_flow_rate::Float64 = 0.0
    """Fan mass flow rate [kg/s]"""
    fan_mass_flow_rate::Float64 = 0.0
    """Fan power [W]"""
    fan_power::Float64 = 0.0
end


"""
$TYPEDEF

Contains the fuel-related parameters, such as type and physical properties.

$TYPEDFIELDS
"""
@kwdef mutable struct Fuel
    """Fuel type"""
    type::String = ""
    """Fuel index"""
    ifuel::Int64 = 0
    """Fuel properties"""
    flow::Vector{FlowStation} = fill(FlowStation(), iptotal)
end

"""
$TYPEDEF

Engine:
    Divided into 6 modules
    1. General Properties
    2. Fuel
    3. Stations
    4. Components
    5. Performance
    6. Misc Properties

$TYPEDFIELDS
"""
@kwdef mutable struct Engine
    # General Properties
    """Engine Weight [N] """
    weight::Float64 = 0.0
    """Engine Weight-moment [Nm] """
    moment::Float64 = 0.0
    """Number of engines"""
    number_engines::Int64 = 0

    """Fuel properties"""
    fuel::Fuel = Fuel()

    """Stations"""
    station::EngineStations = EngineStations()

    """Components"""
    component::EngineComponents = EngineComponents()

    """Performance"""
    performance::Vector{EnginePerformance} = fill(EnginePerformance(), iptotal)

    #Miscellaneous properties
    """Design bypass ratio"""
    design_bypass_ratio::Float64 = 0.0
    """Fan gear ratio"""
    fan_gear_ratio::Float64 = 0.0
end