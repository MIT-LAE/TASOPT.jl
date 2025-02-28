export unpack_ac, unpack_ac_components
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
    - `wing::Wing` : wing object
    - `htail::Tail` : horizontal stabilizer object
    - `vtail::Tail` : vertical stabilizer object
    - `engine::Engine`: engine object
"""
function unpack_ac(ac, imission::Int64; ip::Int64 = 0)
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
    - `wing::Wing` : wing object
    - `htail::Tail` : horizontal stabilizer object
    - `vtail::Tail` : vertical stabilizer object
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
