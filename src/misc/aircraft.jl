mutable struct fuselage_tank
    placement::String
    size_insulation::Bool
    Wfuelintank::Float64
    Rfuse::Float64
    dRfuse::Float64
    wfb::Float64
    nfweb::Float64
    clearance_fuse::Float64

    t_insul::Array{Float64}
    material_insul::Array{String}
    iinsuldes::Array{Int64}

    UTSinner::Float64
    rhoinner::Float64
    UTSouter::Float64
    rhoouter::Float64
    Eouter::Float64
    poissouter::Float64
    ARtank::Float64
    theta_inner::Float64
    theta_outer::Vector{Float64}
    Ninterm::Float64
    
    ptank::Float64
    rhofuel::Float64
    Tfuel::Float64
    boiloff_rate::Float64

    ftankadd::Float64
    ew::Float64
    ullage_frac::Float64
    qfac::Float64
    fuselage_tank() = new() 
end

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

For devs: the indices for accessing specific data are defined in `/src/misc/index.inc`. Refer to the sample input file (`/src/IO/default_input.toml` and `read_input.jl`) for usage.
"""
struct aircraft
    name::String
    description::String
    pari::AbstractVector{Int64}
    parg::AbstractVector{Float64}
    parm::AbstractArray{Float64}
    para::AbstractArray{Float64}
    pare::AbstractArray{Float64}
    fuse_tank::fuselage_tank
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
    Des. Range  = $(round(ac.parm[imRange], sigdigits = 3)) km
    Cruise Mach = $(round(ac.para[iaMach, ipcruise1, 1], sigdigits=3))""")
end


