"""
`materials` is a module that provides basic functionality to represent
various materials such as `StructuralAlloy`s,`Conductor`s, and `Insulator`s. 
"""
module materials

using TOML, DocStringExtensions
import ..TASOPT: __TASOPTroot__

export StructuralAlloy, Conductor, Insulator, ElectricSteel, 
ThermalInsulator, thermal_conductivity
export resistivity, resxden

MaterialProperties = TOML.parsefile(joinpath(__TASOPTroot__,"material_data/MaterialProperties.toml"))

"""
$TYPEDEF

Generic structural alloy.

$TYPEDFIELDS
"""
@kwdef mutable struct StructuralAlloy
    """Name"""
    name::String = ""
    """Density [kg/m³]"""
    ρ::Float64
    """Young's Modulus [Pa]"""
    E::Float64
    """Shear Modulus [Pa]"""
    G::Float64
    """Poisson's Ratio [-]"""
    ν::Float64
    """Maximum Stress (Yield or Ultimate Strength) [Pa]"""
    σmax::Float64
    """Maximum Shear [Pa]"""
    τmax::Float64
    """Yield tensile strength [Pa]"""
    YTS::Float64
    """Ultimate tensile strength [Pa]"""
    UTS::Float64
    """Ultimate shear strength [Pa]"""
    USS::Float64
    """Thermal conductivity [W/(m K)]"""
    k::Float64
end

"""
    StructuralAlloy(material::String; max_avg_stress = 1.1, safety_factor = 1.0)

Outer constructor for `StructuralAlloy` types. 
Material specified needs to have the following data in the database:
- ρ   : Density [kg/m³]
- E   : Young's Modulus [Pa]
- G   : Shear Modulus [Pa]
- ν   : Poisson's Ratio [-]
- σmax: Maximum Stress (Yield or Ultimate Strength) [Pa]
- τmax: Maximum Shear [Pa]
"""
function StructuralAlloy(material::String; max_avg_stress = 1.1, safety_factor = 1.5)
    local MatProp, ρ, E, G, ν, σmax, τmax, YTS, UTS, USS, k
    try
        MatProp = MaterialProperties[material]
    catch
        error("Cannot find $material in Material Properties database")
    else
        try
            ρ    = MatProp["density"]
            E    = MatProp["youngs_modulus"]
            G    = MatProp["shear_modulus"]
            ν    = MatProp["poissons_ratio"]
            YTS  = MatProp["YTS"]
            UTS  = MatProp["UTS"]
            USS  = MatProp["shear_strength"]
            σmax = YTS/max_avg_stress/safety_factor
            τmax = USS/max_avg_stress/safety_factor
            k = MatProp["thermal_conductivity"]
        catch 
            error("Insufficient data in database for $material to build a StructuralAlloy")
        else
            StructuralAlloy(material, ρ, E, G, ν, σmax, τmax, YTS, UTS, USS, k)
        end
    end

end

"""
$TYPEDEF

Generic conductor.

$TYPEDFIELDS
"""
@kwdef struct Conductor 
    """Name"""
    name::String = ""
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
    Conductor(material::String)

Outer constructor for `conductor` types. 
Material specified needs to have the following data in the database:
- ρ (density): Density [kg/m³]
- resistivity: Resistivity [Ω⋅m]
- α (alpha: Thermal coefficient of resisitivity [K⁻¹]
- T0: Temperature at base resistivity [K]
"""
function Conductor(material::String)
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
            error("Insufficient data in database for $material to build a Conductor")
        else
            Conductor(material, ρ, resistivity, α, T0)
        end
    end

end

"""
$TYPEDEF

Generic insulator.

$TYPEDFIELDS
"""
@kwdef struct Insulator
    """Name"""
    name::String = ""
    """Density [kg/m³]"""
    ρ::Float64
    """Dielectric strength [V/m]"""
    Emax::Float64
end

"""
    Insulator(material::String)

Outer constructor for `insulator` types. 
Material specified needs to have the following data in the database:
- ρ (density): Density [kg/m³]
- Emax (dielectric strength): Dielectric strength [V/m]
"""
function Insulator(material::String)
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
            error("Insufficient data in database for $material to build an Insulator")
        else
            Insulator(material, ρ, Emax)
        end
    end

end
"""
$TYPEDEF

ElectricSteel.

$TYPEDFIELDS
"""
@kwdef struct ElectricSteel 
    """Name"""
    name::String = ""
    """Density [kg/m³]"""
    ρ::Float64
    """Eddy current loss coefficient [W/lbm/Hz²/T²]"""
    kₑ::Float64
    """Hysteresis loss coefficient [W/lbm/Hz]"""
    kₕ::Float64
    """Exponential fit coefficient for hysteresis loss"""
    α::Float64
end
"""
    ElectricSteel(material::String)

Outer constructor for `ElectricSteel` types. 
Material specified needs to have the following data in the database:
- ρ (density): Density [kg/m³]
- ke
- kh
- α 
"""
function ElectricSteel(material::String)
    local MatProp, ρ, ke, kh, α
    try
        MatProp = MaterialProperties[material]
    catch
        error("Cannot find $material in Material Properties database")
    else
        try
            ρ = MatProp["density"]
            ke = MatProp["ke"]
            kh = MatProp["kh"]
            α = MatProp["alpha"]
        catch 
            error("Insufficient data in database for $material to build a Conductor")
        else
            ElectricSteel(material, ρ, ke, kh, α)
        end
    end

end


"""
    resxden(cond::conductor)

Returns the resisitivity-density product in kg⋅Ω/m²
"""
function resxden(cond::Conductor)
    return cond.resistivity*cond.ρ
end

"""
    resistivity(cond::conductor, T::Float64)

Returns the resistivity of the conductor at the given temperature. Defaults to 
T = 293.15 K = 20°C
"""
function resistivity(cond::Conductor, T::Float64=293.15)
    ΔT = T - cond.T0
    return cond.resistivity*(1 + cond.α*ΔT)
end  # function resistivity

"""
$TYPEDEF

Generic thermal insulator.

$TYPEDFIELDS
"""
@kwdef struct ThermalInsulator 
    """Name"""
    name::String = ""
    """Density [kg/m³]"""
    ρ::Float64 = 0.0
    """Coefficients for thermal conductivity as a function of temperature [W/(m⋅K)]"""
    conductivity_coeffs::Vector{Float64} = ""
end

"""
    ThermalInsulator(material::String)

Outer constructor for `ThermalInsulator` types. 
Material specified needs to have the following data in the database:
- ρ (density): Density [kg/m³]
- conductivity (thermal conductivity): a string with the thermal conductivity as a function of `T` [W/(m⋅K)]
"""
function ThermalInsulator(material::String)
    local MatProp, ρ, cond_coeffs
    try
        MatProp = MaterialProperties[material]
    catch
        error("Cannot find $material in Material Properties database")
    else
        try
            ρ = MatProp["density"]
            cond_coeffs = MatProp["conductivity_coeffs"]

        catch 
            error("Insufficient data in database for $material to build a ThermalInsulator")
        else
            ThermalInsulator(material, ρ, cond_coeffs)
        end
    end

end

"""
    thermal_conductivity(material::ThermalInsulator, T::Float64)

This function evaluates a thermal conductivity polynomial.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `material::ThermalInsulator`: material object`
      - `T::Float64`: temperature (K)

      **Outputs:**
      - `k::Float64`: thermal conductivity ([W/(m⋅K)]).
"""
function thermal_conductivity(material::ThermalInsulator, T::Float64)
    coeffs = material.conductivity_coeffs
    k = 0.0
    for i = 1:length(coeffs)
        k += coeffs[i] * T^(i-1)
    end
    return k
end

"""
    create_dict(material)

Creates a dictionary from a given material type
"""
function create_material_dict(material)
    fn = fieldnames(typeof(material))
    isempty(fn) ? error("Not an appropriate struct") :
    Dict("Material" => Dict(fn .=> getproperty.([material], fn)))
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

save_material_toml(filename::String, material) = 
save_material_toml(filename, create_material_dict(material))

function Base.show(io::IO, alloy::StructuralAlloy)
    print(io, "StructuralAlloy(", alloy.name, ")")
end

function Base.show(io::IO, ::MIME"text/plain", alloy::StructuralAlloy)
    print("StructuralAlloy(",alloy.name,"):")
    print("\n ρ    = ",alloy.ρ   ," kg/m³")
    print("\n E    = ",alloy.E   ," Pa")
    print("\n G    = ",alloy.G   ," Pa")
    print("\n ν    = ",alloy.ν)
    print("\n YTS  = ",alloy.YTS ," Pa")
    print("\n UTS  = ",alloy.UTS ," Pa")
    print("\n USS  = ",alloy.USS ," Pa")
    print("\n σmax = ",round(alloy.σmax,sigdigits=3)," Pa")
    print("\n τmax = ",round(alloy.τmax,sigdigits=3)," Pa")
end

end