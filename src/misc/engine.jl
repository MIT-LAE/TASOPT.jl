using DocStringExtensions
"""
$TYPEDEF

Engine and models

$TYPEDFIELDS
"""
@kwdef mutable struct Engine
    model_name::String = ""
    model::Function = (x->x)
    weight_model_name::String = ""
    weight_model::Function = (x->x)
end

function store_engine_model!(engine)
    if engine.model_name == "turbofan_md"
        model = tfwrap!
    end

    if engine.weight_model_name == "turbofan"
        weight_model = tfweightwrap!
    end

    engine.model = model
    engine.weight_model = weight_model
end

function extract_engine_model(engine)
    enginecalc! = engine.model
    engineweight! = engine.weight_model
    return enginecalc!, engineweight!
end