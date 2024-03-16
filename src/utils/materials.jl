abstract type AbstractMaterials end
abstract type Metals <: AbstractMaterials end
abstract type Dielectrics <: AbstractMaterials end

MatProp = TOML.parsefile("src/material_data/MaterialProperties.toml")

"""
$TYPEDEF

Material properties covering a wide range of disciplines:

$TYPEDFIELDS
"""
@kwdef struct MaterialProperties
    # ---------------------
    # Mass Properties
    # ---------------------
    """Density [kg/m³]"""
    ρ::Float64
    # ---------------------
    # Structural Properties
    # ---------------------
    """Young's Modulus [Pa]"""
    E::Float64
    """Shear Modulus [Pa]"""
    G::Float64
    """Poisson's Ratio [-]"""
    ν::Float64
    """Maximum Stress [Pa] (Yield or Ultimate Strength)"""
    σmax::Float64
    """Maximum Shear [Pa]"""
    τmax::Float64
    # ---------------------
    # Electric Properties
    # ---------------------
    """Resistivity [Ω⋅m]"""
    resistivity::Float64
    """Thermal coefficient of resitivity [K⁻¹]"""
    α::Float64
    """Temperature at base resistivity [K]"""
    T0::Float64 = 293.15 # 20°C
    """Dielectric strength [V/m]"""
    Emax::Float64

end

"""
$TYPEDEF

Generic conductor.

$TYPEDFIELDS
"""
@kwdef struct conductor <: Metals
    """Density [kg/m³]"""
    ρ::Float64
    """Resistivity [Ω⋅m]"""
    resistivity::Float64
    """Thermal coefficient of resitivity [K⁻¹]"""
    α::Float64
    """Temperature at base resistivity [K]"""
    T0::Float64 = 293.15 # 20°C
end


function conductor(material::String)
    local dict, ρ, resistivity, α, T0
    try
        dict = MatProp[material]
    catch
        error("Cannot find $material in Material Properties database")
    end

    try
        ρ = dict["density"]
        resistivity = dict["resistivity"]
        α = dict["alpha"]
        T0 = dict["T0"]
    catch 
        error("Insufficient data in database for $material to build a conductor")
    end
    conductor(ρ, resistivity, α, T0)
end

"""
$TYPEDEF

Generic insulator.

$TYPEDFIELDS
"""
@kwdef struct insulator <: Dielectrics
    """Density [kg/m³]"""
    ρ::Float64
    """Dielectric strength [V/m]"""
    Emax::Float64
end

# struct SturctAlloy <: Metals
#     UTS
#     E
#    < and so on... >
# end

"""
    resxden(cond::conductor)

Returns the resisitivity-density product in kg⋅Ω/m²
"""
function resxden(cond::conductor)
    return cond.resistivity*cond.ρ
end

"""
    resistivity(cond::conductor, T::Float64)

Returns the resistivity of the conductor at the given temperature. Defaults to 
T = 293.15 K = 20°C
"""
function resistivity(cond::conductor, T::Float64=293.15)
    ΔT = T - cond.T0
    return cond.resistivity*(1 + cond.α*ΔT)
end  # function resistivity

#-------------------------------------------------------------
# Standard materials taken from here: https://en.wikipedia.org/wiki/Electrical_resistivity_and_conductivity#Resistivity_and_conductivity_of_various_materials
const Al = aluminium = conductor(ρ = 2700.0, resistivity = 2.65e-8, α = 3.90e-3)
const Cu = copper = conductor(ρ =  8960.0, resistivity = 1.68e-8, α = 4.04e-3)
const Ag = silver = conductor(ρ = 10490.0, resistivity = 1.59e-8, α = 3.80e-3)
const Au = gold   = conductor(ρ = 19300.0, resistivity = 2.44e-8, α = 3.40e-3)

# Insulators
const PTFE = insulator(ρ = 2150.0, Emax = 19.7e6)
const PEEK = insulator(ρ = 1320.0, Emax = 23e6)
const polyimide = insulator(ρ = 1700, Emax = 10e6) # Dowdle et al https://doi.org/10.2514/6.2018-5026



"""
    create_dict(material::AbstractMaterials)

Creates a dictionary from a given AbstractMaterials subtype
"""
function create_dict(material::AbstractMaterials)
    fn = fieldnames(typeof(material))
    dict = Dict("Material" => Dict(fn .=> getproperty.([material], fn)))
end

"""
    save_material_toml(filename::String, dict::AbstractDict)

Takes a filename and dict and saves a TOML file
"""
function save_material_toml(filename::String, dict::AbstractDict)
    if splitext(filename)[end] == ""
        filename = filename*".toml"
    end
    open(filename, "w") do io
        TOML.print(io, dict)
    end
end  # function save_material_toml

save_material_toml(filename::String, material::AbstractMaterials) = 
save_material_toml(filename, create_dict(material))