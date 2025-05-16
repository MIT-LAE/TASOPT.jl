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
    """Flag if engine core ingests upstream BL. `false` for clean flow, `true` if ingests KE defect """
    has_BLI_cores::Bool
end

# Data object for a ducted fan powered by a fuel cell
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
    """Design stack voltage (V)"""
    design_voltage::Float64
    """Design stack_power (W)"""
    design_power::Float64
    """Fuel cell total power (W)"""
    FC_power::Array{Float64, 2}
    """Fuel cell heat rate (W)"""
    FC_heat::Array{Float64, 2}
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
        design_voltage::Float64=0.0,
        design_power::Float64=0.0)
        # Initialize matrix fields
        matrices = ntuple(_ -> zeros(iptotal, nmis), 11)  # Adjust count for number of matrix fields

        return new(type, number_cells, area_cell, thickness_membrane, 
        thickness_anode, thickness_cathode, design_voltage, design_power,matrices...)
    end
end