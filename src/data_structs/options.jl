"""
    $TYPEDEF

Configuration options for an [`aircraft`](@ref).

$TYPEDFIELDS
"""
@kwdef mutable struct Options
    #fuel options
    opt_fuel::String
    has_centerbox_fuel::Bool
    has_wing_fuel::Bool
    # has_fuselage_fuel::Bool       #TODO: put this in and...
    fuselage_fueltank_count::Int    #... move this field to the fuselage class.

    #wing options
    opt_wing_type::String
    moves_wingbox_forbalance::Bool
    
    #engine options
    opt_engine_location::String
    opt_engine_type::String
    opt_engine_model::String
    opt_engine_weight_model::String
    has_BLI_cores::Bool
    #TODO: has_BLI_cores should be in the engine class^
    
    #fuselage options
    opt_fuselage_taper::String
    is_doubledecker::Bool
    #TODO: these should be in the fuselage class^
end
