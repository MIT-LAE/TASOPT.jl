using DocStringExtensions
abstract type AbstractModel end
abstract type AbstractData end
"""
$TYPEDEF

Engine and models

$TYPEDFIELDS
"""
mutable struct Engine{M<:AbstractModel}
    model::M

    #TODO add more engine elements are pare gets undone
    """Engine data storage"""
    data::AbstractData
    """Heat exchanger parameters and data"""
    heat_exchangers::Vector{HX_struct}

end

struct EmptyData <: AbstractData 
end
#----------------------------
# Turbofan model and data
#----------------------------
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
end

#----------------------------
# Fuel cell with ducted fan model and data
#----------------------------
struct FuelCellDuctedFan{F1, F2} <: AbstractModel
    """Engine performance model identifier"""
    model_name::String
    """Engine function to be used by TASOPT"""
    enginecalc!::F1
    """Weight model identifier"""
    weight_model_name::String 
    """Weight model to be used by TASOPT"""
    engineweight!::F2
end

mutable struct FuelCellDuctedFanData <: AbstractData
    """Fuel cell type"""
    type::String
    """Number of cells"""
    number_cells::Float64
    """Cell area (m^2)"""
    area_cell::Float64
    """Membrane thickness (m)"""
    thickness_membrane::Float64
    """Anode thickness (m)"""
    thickness_anode::Float64
    """Cathode thickness (m)"""
    thickness_cathode::Float64
    """Design_stack_power (W)"""
    design_power::Float64
    """Fuel cell total power (W)"""
    fuel_cell_power::Array{Float64, 2}
    """Fuel cell heat rate (W)"""
    fuel_cell_heat::Array{Float64, 2}
    """Stack voltage (V)"""
    stack_voltage::Array{Float64, 2}
    """Fuel cell temperature (K)"""
    FC_temperature::Array{Float64, 2}
    """Fuel cell pressure (Pa)""" 
    FC_pressure::Array{Float64, 2}
    """Current density (A/m^2)"""
    current_density::Array{Float64, 2}
    """Water flux ratio"""
    α_water::Array{Float64, 2}
    """Water concentration at anode""" 
    water_concentration_anode::Array{Float64, 2}
    """Water concentration at cathode""" 
    water_concentration_cathode::Array{Float64, 2}
    """Stoichiometric ratio of hydrogen"""
    λ_H2::Array{Float64, 2}
    """Stoichiometric ratio of oxygen"""
    λ_O2::Array{Float64, 2}
   
    function FuelCellDuctedFanData(nmis::Int;
        type::String="", 
        number_cells::Float64=0.0, 
        area_cell::Float64=0.00, 
        thickness_membrane::Float64=0.0, 
        thickness_anode::Float64=0.0, 
        thickness_cathode::Float64=0.0,
        design_power::Float64=0.0)
        # Initialize matrix fields
        matrices = ntuple(_ -> zeros(iptotal, nmis), 11)  # Adjust count for number of matrix fields

        return new(type, number_cells, area_cell, thickness_membrane, 
        thickness_anode, thickness_cathode, design_power,matrices...)
    end
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