using DocStringExtensions
abstract type AbstractModel end
"""
$TYPEDEF

Engine and models

$TYPEDFIELDS
"""
@kwdef mutable struct Engine
    model::AbstractModel

    """Heat exchanger parameters and data"""
    heat_exchangers::Vector{Any} = []
end

struct TurbofanModel <: AbstractModel
    """Engine performance model identifier"""
    model_name::String
    """Engine function to be used by TASOPT"""
    enginecalc!::Function
    """Weight model identifier"""
    weight_model_name::String 
    """Weight model to be used by TASOPT"""
    engineweight!::Function
end

# Override Engine getproperty to return default values
function Base.getproperty(obj::Engine, sym::Symbol)
    if sym === :type
        if typeof(obj.model) == TurbofanModel
            return "turbofan"
        end
    elseif sym === :enginecalc! || sym ===:engineweight! #Access model directly from engine
        return getfield(obj.model, sym)

    else
        return getfield(obj, sym)
    end
end  # function Base.getproperty