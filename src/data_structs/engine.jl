using DocStringExtensions
abstract type AbstractModel end
"""
$TYPEDEF

Engine and models

$TYPEDFIELDS
"""
mutable struct Engine{M<:AbstractModel}
    model::M

    #TODO add more engine elements are pare gets undone
    """Heat exchanger parameters and data"""
    heat_exchangers::Vector{HX_struct}
end

"""
$TYPEDEF

Turbofan model

$TYPEDFIELDS
"""
struct TurbofanModel{F1, F2} <: AbstractModel
    """Engine performance model identifier"""
    model_name::String
    """Engine function to be used by TASOPT"""
    enginecalc!::F1
    """Weight model identifier"""
    weight_model_name::String 
    """Weight model to be used by TASOPT"""
    engineweight!::F2
    """Flag if engine core ingests upstream BL. `false` for clean flow, `true` if ingests KE defect """
    has_BLI_cores::Bool
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