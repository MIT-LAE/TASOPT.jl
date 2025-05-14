"""
    aircraft

A type representing a TASOPT aircraft model including, geometric,
aerodynamic, propulsion system parameters.
It is designed to hold information related to the aircraft's name, description,
as well as different sets of parameters used for analysis and optimization.

Overloads Base.summary to print a summary of the `aircraft` model.

# Fields:
- `name::String` : Aircraft name (eg: "Boeing 777")      
- `description::String` : A brief description of the aircraft
- `options::TASOPT.Options` : Configuration options for the aircraft
- `parg::AbstractArray{Float64}` : Geometry parameters                   
- `parm::AbstractArray{Float64}` : Mission parameters                    
- `para::AbstractArray{Float64}` : Aero parameters                       
- `pare::AbstractArray{Float64}` : Engine parameters 
- `fuse_tank::fuselage_tank`: fuselage fuel tank object
- `fuselage::Fuselage`: fuselage fuel tank object
- `wing::Wing`: wing object
- `htail::Tail`: horizontal tail object
- `vtail::Tail`: vertical tail object
- `engine::Engine`: engine object
- `landing_gear::LandingGear`: landing gear object

- `fuselage::Fuselage` : Fuselage layout, data, and parameters
- `fuse_tank::fuselage_tank` : Fuselage tank data and parameters (when applicable)
- `wing::Wing` : Wing data and parameters
- `htail::Tail` : Horizontal tail data and parameters
- `vtail::Tail` : Vertical tail data and parameters
- `engine::Engine` : Engine models, data, and parameters

The indices for accessing specific data in the `par` arrays are defined in `/src/data_structs/index.inc`. 
Refer to the sample input file (`/example/defaults/default_input.toml` and `read_input.jl`) for usage.
Refer to the docs for a summary of the main `struct`s.
"""
@kwdef mutable struct aircraft
    name::String = "Untitled Aircraft"
    description::String = "Indescribable"
    options::TASOPT.Options

    parg::Vector{Float64}
    parm::Array{Float64, 2}
    para::Array{Float64, 3}
    pare::Array{Float64, 3}
    
    is_sized::Vector{Bool} = [false]

    fuselage::Fuselage = Fuselage()
    fuse_tank::fuselage_tank = fuselage_tank()
    wing::Wing = Wing()
    htail::Tail = Tail()
    vtail::Tail = Tail()
    engine::Engine = Engine()
    landing_gear::LandingGear = LandingGear()
    #TODO: update DOCSTRING for ANY NEW fields/sub-structures
end


function Base.getproperty(ac::aircraft, sym::Symbol)
    if sym === :parad #Design para
        return view(getfield(ac, :para), :, : , 1) 
    elseif sym === :pared #Design pare
        return view(getfield(ac, :pare), :, : , 1) 
    elseif sym === :parmd #Design parm
        return view(getfield(ac, :parm), :, 1) 
    else
        return getfield(ac, sym)
    end
end  # function getproperty

function Base.summary(ac::aircraft)
    println("\n----- TASOPT model summary -----")
    println(ac.name)
    println(ac.description)
    geometry(ac)
    weight_buildup(ac)
    aero(ac)
end
function Base.show(io::IO, ac::aircraft)
    print(io, 
    """Name: $(ac.name);
    Wpay = $(round(ac.parm[imWpay]/1e3, sigdigits = 3)) kN
    Des. Range  = $(round(ac.parm[imRange]/1e3, sigdigits = 3)) km
    Cruise Mach = $(round(ac.para[iaMach, ipcruise1, 1], sigdigits=3))""")
end