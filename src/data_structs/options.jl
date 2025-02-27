"""
    $TYPEDEF

Field of an [`aircraft`](@ref) containing configuration-level design choices.

$TYPEDFIELDS
"""
@kwdef mutable struct Options
    #fuel options
    """Fuel type (e.g., Jet-A, LH2)"""
    opt_fuel::String
    """Fuel option index (non-driving; determined and used by gas calcs)"""
    ifuel::Integer
    """Indicates presence of centerbox fuel tank, can only be true if has_wing_fuel is true"""
    has_centerbox_fuel::Bool
    """Indicates presence of wing fuel tanks """
    has_wing_fuel::Bool
    """Indicates presence of fuselage fuel tanks (non-driving; set by `fuse_tank` inputs)"""
    has_fuselage_fuel::Bool 
      #TODO: consider making ^ a driving parameter, rather than a reflection of fuse_tank parameters
      #Note: right now fuel can only be stored in the wings or the fuselage, not both
    
    #engine options
    """Engine location ("wing", "fuselage")"""
    opt_engine_location::String
    """Propulsion system architecture (e.g., "tf" for turbofan, "te" for turboelectric), performance and weight models set in ac.Engine"""
    opt_prop_sys_arch::String
    
    #fuselage/cabin options

    """Indicates if the aircraft has a double-decker fuselage configuration"""
    is_doubledecker::Bool

    #wing/stability options
    """Move wingbox selection for longitudinal stability analysis. "fixed" = static wing position ,"fixed_CLh" move wing to get CLh="CLhspec" in cruise, "min_static_margin" = move wing to get min static margin = "SMmin"  """
    opt_move_wing::String
end

function Base.summary(opt::Options)
  println("\n-------- Options Summary --------")
  println("Fuel Type: ", opt.opt_fuel)
  println("Fuel stored in: "*(opt.has_wing_fuel ? "wing "*(opt.has_centerbox_fuel ? " wingbox " : "(none in wingbox)") : "")*(opt.has_fuselage_fuel ? "fuselage " : ""))
  println("Propulsion Architecture: ", opt.opt_prop_sys_arch)
  println("Engine Location: ", opt.opt_engine_location)
  println("Cabin decks: ", opt.is_doubledecker ? "double" : "single")
  println("---------------------------------")
end
function Base.show(io::IO, opt::Options)
  print(io, "Options(Fuel: $(opt.opt_fuel); Fuel Storage: " * (opt.has_wing_fuel ? "wing " * 
      (opt.has_centerbox_fuel ? "wingbox " : "(none in wingbox) ") : "") * (opt.has_fuselage_fuel ? "fuselage " : "") * 
      "; Propulsion Arch.: $(opt.opt_prop_sys_arch); Engine Location: $(opt.opt_engine_location); Cabin decks: $(opt.is_doubledecker ? "double" : "single"))")
end
