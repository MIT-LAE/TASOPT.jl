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
- `pari::AbstractVector{Int64}` : integer flag parameters               
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
- `sized::AbstractVector{1,Bool}`: flag if aircraft is sized (default is `[false]`)

For devs: the indices for accessing specific data are defined in `/src/misc/index.inc`. Refer to the sample input file (`/src/IO/default_input.toml` and `read_input.jl`) for usage.
"""
Base.@kwdef mutable struct aircraft #inner constructor
    name::String = "Untitled Aircraft"
    description::String = "Indescribable"
    aircraft_type::Symbol = :unknown
    pari::AbstractVector{Int64}
    parg::AbstractVector{Float64}
    parm::AbstractArray{Float64}
    para::AbstractArray{Float64}
    pare::AbstractArray{Float64}

    fuse_tank::fuselage_tank = fuselage_tank()
    fuselage::Fuselage = Fuselage()
    wing::Wing = Wing()
    htail::Tail = Tail()
    vtail::Tail = Tail()
    engine::Engine = Engine()
    landing_gear::LandingGear = LandingGear()

    sized::AbstractVector{Bool} = [false]

    #TODO: update DOCSTRING for ANY NEW fields/sub-structures
end

# #TODO: sort out a robust meta-structure such that new individual constructors aren't required
# #outer constructor for if `sized` and fuse_tank not given
# function aircraft(name::String, description::String, pari::AbstractVector{Int64}, parg::AbstractVector{Float64},
#         parm::AbstractArray{Float64}, para::AbstractArray{Float64}, pare::AbstractArray{Float64}) 
#         return aircraft(name, description, pari, parg, parm, para, pare, [false])
# end
# #constructor for if fuse_tank not given
function aircraft(name::String, description::String, aircraft_type::Symbol, pari::AbstractVector{Int64}, parg::AbstractVector{Float64},
        parm::AbstractArray{Float64}, para::AbstractArray{Float64}, pare::AbstractArray{Float64}, 
        sized::AbstractVector{Bool}) 
        return aircraft(name, description, pari, parg, parm, para, pare, fuselage_tank(), Fuselage(), Wing(), Tail(), Tail(), LandingGear(), sized)
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