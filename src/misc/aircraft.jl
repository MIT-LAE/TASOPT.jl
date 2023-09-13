"""
This is a structure representing an aircraft with various parameters.
It is designed to hold information related to the aircraft's name, description,
as well as different sets of parameters used for analysis and optimization.

Fields:
- `name::String` : Aircraft name (eg: "Boeing 777")      
- `description::String` : A brief description of the aircraft
- `pari::AbstractVector{Int64}` : integer flag parameters               
- `parg::AbstractArray{Float64}` : Geometry parameters                   
- `parm::AbstractArray{Float64}` : Mission parameters                    
- `para::AbstractArray{Float64}` : Aero parameters                       
- `pare::AbstractArray{Float64}` : Engine parameters                     
"""
struct aircraft
    name::String
    description::String
    pari::AbstractVector{Int64}
    parg::AbstractVector{Float64}
    parm::AbstractArray{Float64}
    para::AbstractArray{Float64}
    pare::AbstractArray{Float64}
end

function Base.summary(ac::aircraft)
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

function Base.getproperty(ac::aircraft, sym::Symbol)
    if sym == :description
        println(getfield(ac, sym))
    else
        getfield(ac, sym)
    end
end

