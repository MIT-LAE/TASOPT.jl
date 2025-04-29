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
    - `parg::AbstractArray{Float64}` : Geometry parameters 
    - `parm::AbstractArray{Float64}` : Mission parameters                    
    - `para::AbstractArray{Float64}` : Aero parameters                       
    - `pare::AbstractArray{Float64}` : Engine parameters      
    - `options::Options` : aircraft configuration options
    - `fuse::Fuselage` : fuselage parameters             
    - `fuse_tank::fuselage_tank` : fuel tank in fuselage parameters
    - `wing::Wing` : wing object
    - `htail::Tail` : horizontal stabilizer object
    - `vtail::Tail` : vertical stabilizer object
    - `engine::Engine`: engine object
    - `landing_gear::LandingGear`: landing gear object
"""
function unpack_ac(ac, imission::Int64; ip::Int64 = 0)
    parg = ac.parg
    parm = view(ac.parm, :, imission) 
    options = ac.options

    fuse_tank = ac.fuse_tank
    fuse = ac.fuselage 
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    engine = ac.engine
    landing_gear = ac.landing_gear

    if ip == 0 #If no point is given
        para = view(ac.para, :, :, imission)
        pare = view(ac.pare, :, :, imission)
    else #ip is given
        para = view(ac.para, :, ip, imission)
        pare = view(ac.pare, :, ip, imission)
    end

    return parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear
end

"""
    unpack_ac_components(ac)
    
Helper function to unpack aircraft physical components.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft` : aircraft object to unpack   
    **Outputs:**
    - `parg::AbstractArray{Float64}` : Geometry parameters      
    - `options::Options` : aircraft configuration options
    - `fuse::Fuselage` : fuselage parameters             
    - `fuse_tank::fuselage_tank` : fuel tank in fuselage parameters
    - `wing::Wing` : wing object
    - `htail::Tail` : horizontal stabilizer object
    - `vtail::Tail` : vertical stabilizer object
    - `landing_gear::LandingGear`: landing gear object
"""
function unpack_ac_components(ac)
    parg = ac.parg
    options = ac.options
    fuse_tank = ac.fuse_tank
    fuse = ac.fuselage 
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    landing_gear = ac.landing_gear

    return parg, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear
end
