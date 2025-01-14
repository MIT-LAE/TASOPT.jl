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
- `sized::AbstractVector{1,Bool}`: flag if aircraft is sized (default is `[false]`)

For devs: the indices for accessing specific data are defined in `/src/misc/index.inc`. Refer to the sample input file (`/src/IO/default_input.toml` and `read_input.jl`) for usage.
"""
Base.@kwdef mutable struct aircraft #inner constructor
    name::String = "Untitled Aircraft"
    description::String = "Indescribable"
    pari::AbstractVector{Int64}
    parg::AbstractVector{Float64}
    parm::AbstractArray{Float64}
    para::AbstractArray{Float64}
    pare::AbstractArray{Float64}
    
    sized::AbstractVector{Bool} = [false]
    fuse_tank::fuselage_tank = fuselage_tank()
    fuselage::Fuselage = Fuselage()
    wing::Wing = Wing()
    htail::Tail = Tail()
    vtail::Tail = Tail()
    engine::Engine = Engine()
end

# #TODO: sort out a robust meta-structure such that new individual constructors aren't required
# #outer constructor for if `sized` and fuse_tank not given
# function aircraft(name::String, description::String, pari::AbstractVector{Int64}, parg::AbstractVector{Float64},
#         parm::AbstractArray{Float64}, para::AbstractArray{Float64}, pare::AbstractArray{Float64}) 
#         return aircraft(name, description, pari, parg, parm, para, pare, [false])
# end
# #constructor for if fuse_tank not given
function aircraft(name::String, description::String, pari::AbstractVector{Int64}, parg::AbstractVector{Float64},
        parm::AbstractArray{Float64}, para::AbstractArray{Float64}, pare::AbstractArray{Float64}, 
        sized::AbstractVector{Bool}) 
        return aircraft(name, description, pari, parg, parm, para, pare, sized, fuselage_tank(), Fuselage(), Wing(), Tail(), Tail())
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

"""
    unpack_ac(ac, imission; ip = 0)
Helper function to unpack all aircraft parameters.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft` : aircraft object to unpack   
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index (optional)
    **Outputs:**
    - `pari::AbstractVector{Int64}` : integer flag parameters               
    - `parg::AbstractArray{Float64}` : Geometry parameters 
    - `parm::AbstractArray{Float64}` : Mission parameters                    
    - `para::AbstractArray{Float64}` : Aero parameters                       
    - `pare::AbstractArray{Float64}` : Engine parameters      
    - `fuse::Fuselage` : fuselage parameters             
    - `fuse_tank::fuselage_tank` : fuel tank in fuselage parameters
"""
function unpack_ac(ac, imission; ip = 0)
    pari = ac.pari
    parg = ac.parg
    parm = view(ac.parm, :, imission) 
    fuse_tank = ac.fuse_tank
    fuse = ac.fuselage 
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    engine = ac.engine

    if ip == 0 #If no point is given
        para = view(ac.para, :, :, imission)
        pare = view(ac.pare, :, :, imission)
    else #ip is given
        para = view(ac.para, :, ip, imission)
        pare = view(ac.pare, :, ip, imission)
    end

    return pari, parg, parm, para, pare, fuse, fuse_tank, wing, htail, vtail, engine
end

"""
    unpack_ac_components(ac)
Helper function to unpack aircraft physical components.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft` : aircraft object to unpack   
    **Outputs:**
    - `pari::AbstractVector{Int64}` : integer flag parameters               
    - `parg::AbstractArray{Float64}` : Geometry parameters      
    - `fuse::Fuselage` : fuselage parameters             
    - `fuse_tank::fuselage_tank` : fuel tank in fuselage parameters
    - `landing_gear::LandingGear`: landing gear parameters
"""
function unpack_ac_components(ac)
    pari = ac.pari
    parg = ac.parg
    fuse_tank = ac.fuse_tank
    fuse = ac.fuselage 
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail

    return pari, parg, fuse, fuse_tank, wing, htail, vtail
end