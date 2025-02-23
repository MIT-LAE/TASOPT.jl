"""
    $TYPEDEF

Configuration options for an [`aircraft`](@ref).

$TYPEDFIELDS
"""
@kwdef mutable struct Options
    #fuel options
    """Fuel type (e.g., Jet-A, LH2)"""
    opt_fuel::String
    """Fuel option index (non-driving; determined and used by gas calcs)"""
    ifuel::Integer
    """Indicates presence of centerbox fuel tank [true/false], can only be true if has_wing_fuel is true"""
    has_centerbox_fuel::Bool
    """Indicates presence of wing fuel tanks [true/false]"""
    has_wing_fuel::Bool
    """Indicates presence of fuselage fuel tanks [true/false] (non-driving; set by `fuse_tank`` inputs)"""
    has_fuselage_fuel::Bool 
      #TODO: consider making ^ a driving parameter, rather than a reflection of fuse_tank parameters
      #Note: right now fuel can only be stored in the wings or the fuselage, not both
    
    #engine options
    """Engine location ("wing", "fuselage")"""
    opt_engine_location::String
    """Propulsion system architecture (e.g., "tf" for turbofan, "te" for turboelectric), performance and weight models set in ac.Engine"""
    opt_prop_sys_arch::String
    
    #fuselage/cabin options
    is_doubledecker::Bool
end
