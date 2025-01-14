using DocStringExtensions
"""
$TYPEDEF

Engine and models

$TYPEDFIELDS
"""
@kwdef mutable struct Engine
    """Engine performance model identifier"""
    model_name::String = ""
    """Engine function to be used by TASOPT"""
    model::Function = (x->x)
    """Weight model identifier"""
    weight_model_name::String = ""
    """Weight model to be used by TASOPT"""
    weight_model::Function = (x->x)
end

"""
    store_engine_model!(engine)
Function to produce and store an engine model based on identifiers.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `engine::Engine`: engine object

    **Outputs:**
    No direct outputs. The `engine` object is modified with the functions.
"""
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

"""
    extract_engine_model(engine)
Function that returns handles corresponding to the engine models.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `engine::Engine`: engine object

    **Outputs:**
    - `enginecalc!::Function`: engine performance function
    - `engineweight!::Function`: engine weight function
"""
function extract_engine_model(engine)
    enginecalc! = engine.model
    engineweight! = engine.weight_model
    return enginecalc!, engineweight!
end