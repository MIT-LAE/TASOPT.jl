"""
    $TYPEDEF

Configuration options for an [`aircraft`](@ref).

$TYPEDFIELDS
"""
@kwdef mutable struct Options
    #fuel options
    """Fuel type_________________""" #TODO: finish these descriptions
    opt_fuel::String
    ifuel::Integer
    has_centerbox_fuel::Bool
    has_wing_fuel::Bool
    has_fuselage_fuel::Bool 
        #TODO: consider making ^ a driving parameter, rather than a reflection of fuse_tank parameters
    
    #engine options
    opt_engine_location::String
    opt_prop_sys_arch::String
    opt_engine_model::String
    opt_engine_weight_model::String
    
    #fuselage/cabin options
    is_doubledecker::Bool
end
