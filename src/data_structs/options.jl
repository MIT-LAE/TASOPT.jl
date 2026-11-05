"""
    $TYPEDEF

Field of an [`aircraft`](@ref) containing configuration-level design choices.

$TYPEDFIELDS
"""
@kwdef mutable struct Options
    #fuel options
    """Fuel type"""
    opt_fuel::FuelType.T
    """Fuel option index (non-driving; determined and used by gas calcs)"""
    ifuel::Int
    """Indicates presence of centerbox fuel tank, can only be true if has_wing_fuel is true"""
    has_centerbox_fuel::Bool
    """Indicates presence of wing fuel tanks """
    has_wing_fuel::Bool
    """Indicates presence of fuselage fuel tanks (non-driving; set by `fuse_tank` inputs)"""
    has_fuselage_fuel::Bool 
      #TODO: consider making ^ a driving parameter, rather than a reflection of fuse_tank parameters
      #Note: right now fuel can only be stored in the wings or the fuselage, not both
    
    #engine options
    """Engine location"""
    opt_engine_location::EngineLocation.T
    """Propulsion system architecture, performance and weight models set in ac.Engine"""
    opt_prop_sys_arch::PropSysArch.T
    """Calculate takeoff length and engine performance"""
    calculate_takeoff::Bool
    
    #fuselage/cabin options

    """Indicates if the aircraft has a double-decker fuselage configuration"""
    is_doubledecker::Bool

    #wing/stability options
    """Wing position strategy for longitudinal stability analysis: `Fixed` = static wing position, `FixedCLh` = move wing to achieve `CLh = CLhspec` in cruise, `MinStaticMargin` = move wing to achieve minimum static margin = `SMmin`"""
    opt_move_wing::WingMove.T
    """Wing bracing configuration: `Cantilever` = unbraced, `StrutBraced` = braced by a strut, `TrussBraced` = braced by a truss"""
    opt_wing_config::WingConfig.T = WingConfig.Cantilever
    """If `true`, the positive manoeuvre wing sizing cases (`pullup` and `pullup_light`) use the FAR/CS 25.337 load factor of the current maximum takeoff weight, instead of their `load_factor` input"""
    opt_far_load_factor::Bool = false

    #Trefftz plane options
    """Trefftz plane induced drag analysis configuration (discretization, k_tip, bunch, root_contraction)"""
    trefftz_config::aerodynamics.TrefftzPlaneConfig
end

# `options.has_strut` is derived from `opt_wing_config` (true for `StrutBraced` and `TrussBraced`,
# false for `Cantilever`). It is read-only: there is no such field, so assigning to it throws.
function Base.getproperty(opt::Options, name::Symbol)
  name === :has_strut && return getfield(opt, :opt_wing_config) != WingConfig.Cantilever
  return getfield(opt, name)
end
Base.propertynames(opt::Options, private::Bool=false) = (fieldnames(Options)..., :has_strut)

function Base.summary(opt::Options)
  println("\n-------- Options Summary --------")
  println("Fuel Type: ", string(opt.opt_fuel))
  println("Fuel stored in: "*(opt.has_wing_fuel ? "wing "*(opt.has_centerbox_fuel ? " wingbox " : "(none in wingbox)") : "")*(opt.has_fuselage_fuel ? "fuselage " : ""))
  println("Propulsion Architecture: ", string(opt.opt_prop_sys_arch))
  println("Engine Location: ", string(opt.opt_engine_location))
  println("Cabin decks: ", opt.is_doubledecker ? "double" : "single")
  println("Wing Configuration: ", string(opt.opt_wing_config))
  println("FAR Load Factor: ", opt.opt_far_load_factor)
  println("---------------------------------")
end
function Base.show(io::IO, opt::Options)
  print(io, "Options(Fuel: $(string(opt.opt_fuel)); Fuel Storage: " * (opt.has_wing_fuel ? "wing " *
      (opt.has_centerbox_fuel ? "wingbox " : "(none in wingbox) ") : "") * (opt.has_fuselage_fuel ? "fuselage " : "") * 
      "; Propulsion Arch.: $(string(opt.opt_prop_sys_arch)); Engine Location: $(string(opt.opt_engine_location)); Cabin decks: $(opt.is_doubledecker ? "double" : "single"); Wing Configuration: $(string(opt.opt_wing_config)); FAR Load Factor: $(opt.opt_far_load_factor))")
end
