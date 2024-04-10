"""
    find_mdot_time(t, parg, para, pare)

This function outputs the fuel mass flow rate to the engines as a function of time for a TASOPT aircraft 
model.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: time in mission (s)
    - `parg::Vector{Float64}`: vector with aircraft geometric parameters
    - `para::Array{Float64}`: array with aircraft aerodynamic parameters
    - `pare::Array{Float64}`: array with aircraft engine parameters
    
    **Outputs:**
    - `t::mdot`: fuel mass flow rate out of the tank (kg/s)
"""
function find_mdot_time(t, parg, para, pare)
    mdots = pare[ieff, :,1] .* pare[iemcore, :,1] .* parg[igneng]

    times = para[iatime,:,1]

    mdot = 0.0
    for ip = 1:(length(times)-1) 
        if (t >= times[ip]) && (t< times[ip+1]) #If the point is the correct one
            t0 = times[ip] #Time at start of segment
            tf = times[ip+1] #Time at end of segment	

            if mdots[ip] > 0
                k = log(mdots[ip+1]/mdots[ip])/( t0 - tf)
                mdot = mdots[ip] * exp( - k * (t - times[ip])) #Assume fuel burn varies exponentially in mission
            else
                mdot = 0.0
            end
        end
    end
    return mdot
end

"""
    calc_Q_points(fuse_tank, pari, parg, para)

This function calculates the heat transfer rate into the tank at the design mission points.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuse_tank::FuselageTank`: struct with aircraft cryogenic tank parameters
    - `pari::Vector{Int64}`: vector with aircraft Boolean and integer parameters
    - `parg::Vector{Float64}`: vector with aircraft geometric parameters
    - `pare::Array{Float64}`: array with aircraft engine parameters
    
    **Outputs:**
    - `Qs::Vector{Float64}`: vector with heat transfer rate at mission points (W)
"""
function calc_Q_points(fuse_tank, pari, parg, para)
    #Extract tank parameters
    if fuse_tank.placement == "rear"
        xftank = parg[igxftankaft]
    else
        xftank = parg[igxftank]
    end
    ifuel = pari[iifuel]

    Npoints = size(para)[2]
    Qs = zeros(Npoints)
    for ip = 1:Npoints
        Mair = para[iaMach, ip, 1]
        z = para[iaalt, ip, 1]
        t = para[iatime, ip, 1]

        #Calculate heat rate at this point
        Qs[ip], _, _ = structures.tankWthermal(fuse_tank, z, Mair, xftank, t, ifuel)
       
    end
    return Qs
end

"""
    find_Q_time_interp(t, para, Qs)

This function estimates the heat transfer rate into the tank in a TASOPT model for a given time. It uses precomputed rates
at each mission point for speed.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: time in mission (s)
    - `para::Array{Float64}`: array with aircraft aerodynamic parameters
    - `Qs::Vector{Float64}`: vector with heat transfer rate at mission points (W)
    
    **Outputs:**
    - `Q::Float64`: heat transfer rate (W)
"""
function find_Q_time_interp(t, para, Qs)
    times = para[iatime,:,1]

    Q = 0.0
    for ip = 1:(length(times)-1)
        if (t >= times[ip]) && (t< times[ip+1]) #If the point is the correct one
            t0 = times[ip]
            tf = times[ip+1]
            Q0 = Qs[ip]
            Qf = Qs[ip+1]

            Q = Q0 + (Qf - Q0)/(tf-t0) * (t - t0) #Interpolate linearly between points
        end
    end
    return Q
end

"""
    find_Q_time(t, fuse_tank, pari, parg, para)

This function calculates the heat transfer rate into the tank in a TASOPT model for a given time.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: time in mission (s)
    - `fuse_tank::FuselageTank`: struct with aircraft cryogenic tank parameters
    - `pari::Vector{Int64}`: vector with aircraft Boolean and integer parameters
    - `parg::Vector{Float64}`: vector with aircraft geometric parameters
    - `para::Array{Float64}`: array with aircraft aerodynamic parameters
    
    **Outputs:**
    - `Q::Float64`: heat transfer rate (W)
"""
function find_Q_time(t, fuse_tank, pari, parg, para)
    #Extract tank parameters
    if ac.fuse_tank.placement == "rear"
        xftank = parg[igxftankaft]
    else
        xftank = parg[igxftank]
    end
    ifuel = pari[iifuel]
    times = para[iatime,:,1]

    Q = 0.0
    for ip = 1:(length(times)-1)
        if (t >= times[ip]) && (t< times[ip+1]) #If the point is the correct one
            t0 = times[ip]
            tf = times[ip+1]
            M0 = para[iaMach, ip, 1]
            Mf = para[iaMach, ip+1, 1]
            z0 = para[iaalt, ip, 1]
            zf = para[iaalt, ip+1, 1]

            #Interpolate Mach number and altitude linearly
            Mair = M0 + (Mf - M0)/(tf-t0) * (t - t0) 
            z = z0 + (zf - z0)/(tf-t0) * (t - t0)

            #Calculate heat rate at this point
            Q, _, _ = structures.tankWthermal(fuse_tank, z, Mair, xftank, t, ifuel)
        end
    end
    return Q
end

function analyze_TASOPT_tank(ac_orig, t_hold_orig::Float64 = 0.0, t_hold_dest::Float64 = 0.0, α::Float64 = 1.0, N::Int64 = 1000)
    ac = deepcopy(ac_orig) #Deepcopy original struct to avoid modifying it

    #Modify aircraft with holding times
    para_alt = zeros(size(ac.para)[1], size(ac.para)[2] + 3, size(ac.para)[3])
    ac.para[iatime, :, 1] .= ac.para[iatime, :, 1] .- minimum(ac.para[iatime, :, 1])
    para_alt[:, 3:(iptotal + 2), :] .= ac.para[:,:,:]
    para_alt[iatime, 2:(iptotal + 2), :] .= para_alt[iatime, 2:(iptotal + 2), :] .+ t_hold_orig #Apply hold at origin
    para_alt[iatime, 1, 1] = 0.0
    Np = size(para_alt)[2]
    para_alt[iatime, Np-1, 1] = maximum(para_alt[iatime, :, 1])
    para_alt[iatime, Np, 1] = para_alt[iatime, Np-1, 1] + t_hold_dest #Apply hold at destination

    pare_alt = zeros(size(ac.pare)[1], size(ac.pare)[2] + 3, size(ac.pare)[3])
    pare_alt[:, 3:(iptotal + 2), :] .= ac.pare[:,:,:]
    
    #Precompute heat transfer rate at each mission point for speed
    Qs_points = calc_Q_points(ac.fuse_tank, ac.pari, ac.parg, para_alt)

    #Define functions for heat and fuel burn profiles through mission 
    Q_calc(t) = find_Q_time_interp(t, para_alt, Qs_points)
    W_calc(t) = 0.0
    mdot_calc(t) = find_mdot_time(t, ac.parg, para_alt, pare_alt)

    #Store profiles in input struct
    u = pressure_inputs(Q_calc, W_calc, mdot_calc)

    # Original tank state in TASOPT
    if ac.pari[iifuel] == 11
        species = "CH4"
    elseif ac.pari[iifuel] == 40
        species = "H2"
    end

    ullage_frac = ac.fuse_tank.ullage_frac
    Wfuel = ac.fuse_tank.Wfuelintank
    ρfuel = ac.fuse_tank.rhofuel
    ρfgas = ac.fuse_tank.rhofuelgas
    ρfmix = (1 - ullage_frac) * ρfuel + ullage_frac * ρfgas #Density of saturated mixture in tank
    V = Wfuel / (gee * ρfmix) #Total tank volume taken by saturated mixture

    pmax = ac.fuse_tank.pvent
    p0 = ac.fuse_tank.pinitial
    βmaxp = (1.0 - ullage_frac) #beta at maximum allowable pressure
    β0 = convert_β_same_ρ(species, p0, pmax, βmaxp)

    mixture_init = SaturatedMixture(species, p0, β0) #Use the initial pressure for the intial saturated mixture

    M0 = V * mixture_init.ρ #Initial fluid tank in mass

    xout = 0.0 #Only liquid is being drawn from the tank
    xvent = 1.0 #Only vent gas
    #Store tank parameters in struct 
    params = pressure_params(mixture_init, V, pmax, xout, xvent, α)

    #Integrate profiles across mission
    y0 = [p0, β0, M0, 0.0, 0.0, 0.0] #Initial states
    ts = LinRange(0,maximum(para_alt[iatime,:,1]) - 1.0, N) #Vector with mission times, subtract one second in the end

    dy_dx(t, y, u, p) = TankDerivatives(t, y, u, p) #State derivatives
    y, _ = RK4(dy_dx, ts, y0, u, params) #Find the state evolution in time by integration

    #Produce outputs
    ps = y[1, :] #Pressure evolution
    βs = y[2, :] #Tank fill evolution
    Ms = y[3, :] #Total fluid mass in tank evolution
    Mburns = y[4, :] #Cumulative mass that has been burnt evolution
    Mvents = y[5, :] #Cumulative mass that has been vented
    Mboils = y[6, :] #Cumulative mass that has been boiled off evolution
    mdot_boils = calculate_boiloff_rate(ts, Mboils) #Calculate boiloff rate
    
    mdots = zeros(N)
    Qs = zeros(N)
    for (i,t) in enumerate(ts)
        mdots[i] = mdot_calc(t) #Fuel burn mass flow rate
        Qs[i] = Q_calc(t) #Heat transfer rate into tank
    end

    return ts, ps, βs, Ms, Mburns, Mboils, mdot_boils, Mvents, mdots, Qs
end