abstract type AbstractMaterials end

MaterialProperties = TOML.parsefile("src/material_data/MaterialProperties.toml")


"""
$TYPEDEF

Generic conductor.

$TYPEDFIELDS
"""
@kwdef struct conductor <: AbstractMaterials
    """Density [kg/m³]"""
    ρ::Float64
    """Resistivity [Ω⋅m]"""
    resistivity::Float64
    """Thermal coefficient of resitivity [K⁻¹]"""
    α::Float64
    """Temperature at base resistivity [K]"""
    T0::Float64 = 293.15 # 20°C
end


"""
    conductor(material::String)

Outer constructor for `conductor` types. 
Material specified needs to have the following data in the database:
- ρ (density): Density [kg/m³]
- resistivity: Resistivity [Ω⋅m]
- α (alpha: Thermal coefficient of resisitivity [K⁻¹]
- T0: Temperature at base resistivity [K]
"""
function conductor(material::String)
    local MatProp, ρ, resistivity, α, T0
    try
        MatProp = MaterialProperties[material]
    catch
        error("Cannot find $material in Material Properties database")
    else
        try
            ρ = MatProp["density"]
            resistivity = MatProp["resistivity"]
            α = MatProp["alpha"]
            T0 = MatProp["T0"]
        catch 
            error("Insufficient data in database for $material to build a conductor")
        else
            conductor(ρ, resistivity, α, T0)
        end
    end

end

"""
$TYPEDEF

Generic insulator.

$TYPEDFIELDS
"""
@kwdef struct insulator <: AbstractMaterials
    """Density [kg/m³]"""
    ρ::Float64
    """Dielectric strength [V/m]"""
    Emax::Float64
end

"""
    insulator(material::String)

Outer constructor for `insulator` types. 
Material specified needs to have the following data in the database:
- ρ (density): Density [kg/m³]
- Emax (dielectric strength): Dielectric strength [V/m]
"""
function insulator(material::String)
    local MatProp, ρ, Emax
    try
        MatProp = MaterialProperties[material]
    catch
        error("Cannot find $material in Material Properties database")
    else
        try
            ρ = MatProp["density"]
            Emax = MatProp["dielectric_strength"]
        catch 
            error("Insufficient data in database for $material to build an insulator")
        else
            insulator(ρ, Emax)
        end
    end

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