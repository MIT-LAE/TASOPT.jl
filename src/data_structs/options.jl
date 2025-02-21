"""
    $TYPEDEF

Configuration options for an [`aircraft`](@ref).

$TYPEDFIELDS
"""
@kwdef mutable struct Options
    #fuel options
    """Fuel type_________________""" #TODO: finish these descriptions
    opt_fuel::String
    has_centerbox_fuel::Bool
    has_wing_fuel::Bool
    # has_fuselage_fuel::Bool       #TODO: put this in and...
    fuselage_fueltank_count::Int    #... move this field to the fuselage class.
    
    #engine options
    opt_engine_location::String
    opt_prop_sys_arch::String
    opt_engine_model::String
    opt_engine_weight_model::String
    has_BLI_cores::Bool
    #TODO: has_BLI_cores should be in the engine class^
    
    #fuselage/cabin options
    is_doubledecker::Bool
end
