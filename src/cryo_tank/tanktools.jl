"""
    find_mdot_time(t, pari, parg, para, pare)

This function outputs the fuel mass flow rate to the engines as a function of time for a TASOPT aircraft 
model.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: time in mission (s)
    - `parg::Vector{Float64}`: vector with aircraft integer parameters
    - `parg::Vector{Float64}`: vector with aircraft geometric parameters
    - `para::Array{Float64}`: array with aircraft aerodynamic parameters
    - `pare::Array{Float64}`: array with aircraft engine parameters
    
    **Outputs:**
    - `t::mdot`: fuel mass flow rate out of the tank (kg/s)
"""
function find_mdot_time(t::Float64, pari::Vector{Int64}, parg::Vector{Float64}, para::Array{Float64}, pare::Array{Float64})
    nftanks = pari[iinftanks]

    #Mass flow rate out of tank is total mass flow to engines divided by number of tanks
    mdots = pare[ieff, :] .* pare[iemcore, :] .* parg[igneng] / nftanks

    times = para[iatime,:]

    # Handle cases where t is exactly one of the sample points
    if t in times
        return mdots[findfirst(isequal(t), times)]
    end

    #Otherwise interpolate exponentially
    for i in 1:(length(times)-1)
        if times[i] <= t < times[i+1]
            t0, tf = times[i], times[i+1]
            m0, mf = mdots[i], mdots[i+1]

            if m0 > 0
                k = log(mf / m0) / (tf - t0)
                return m0 * exp(k * (t - t0))
            else
                return 0.0
            end
        end
    end

end

"""
    calc_Q_points(fuse, fuse_tank, pari, parg, para)

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
function calc_Q_points(fuse::Fuselage, fuse_tank::fuselage_tank, pari::Vector{Int64}, parg::Vector{Float64}, para::Array{Float64})
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
        Mair = para[iaMach, ip]
        z = para[iaalt, ip]

        #Calculate heat rate at this point
        Qs[ip] = tankWthermal(fuse, fuse_tank, z, Mair, xftank, ifuel)
       
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
function find_Q_time_interp(t::Float64, para::Array{Float64}, Qs::Vector{Float64})
    times = para[iatime,:]

    # Handle cases where t is exactly one of the sample points
    if t in times
        return Qs[findfirst(isequal(t), times)]
    end

    # Find the interval where t belongs
    for i in 1:(length(times)-1)
        if times[i] <= t < times[i+1]
            t0, tf = times[i], times[i+1]
            Q0, Qf = Qs[i], Qs[i+1]

            # Linear interpolation formula
            return Q0 + (Qf - Q0) / (tf - t0) * (t - t0)
        end
    end
end

"""
    find_Q_time(t, fuse_tank, pari, parg, para)

This function calculates the heat transfer rate into the tank in a TASOPT model for a given time.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: time in mission (s)
    - `fuse::Fuselage`: fuselage object.
    - `fuse_tank::fuselage_tank`: fuselage tank object.
    - `pari::Vector{Int64}`: vector with aircraft Boolean and integer parameters
    - `parg::Vector{Float64}`: vector with aircraft geometric parameters
    - `para::Array{Float64}`: array with aircraft aerodynamic parameters
    
    **Outputs:**
    - `Q::Float64`: heat transfer rate (W)
"""
function find_Q_time(t::Float64, fuse::Fuselage, fuse_tank::fuselage_tank, pari::Vector{Int64}, parg::Vector{Float64}, para::Array{Float64})
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
        if t ≈ times[ip]
            M0 = para[iaMach, ip, 1]
            z0 = para[iaalt, ip, 1]

            Q = tankWthermal(fuse, fuse_tank, z0, M0, xftank, ifuel)
        elseif (t >= times[ip]) && (t< times[ip+1]) #If the point is the correct one
            t0 = times[ip]
            tf = times[ip+1]
            M0 = para[iaMach, ip]
            Mf = para[iaMach, ip+1]
            z0 = para[iaalt, ip]
            zf = para[iaalt, ip+1]

            #Interpolate Mach number and altitude linearly
            Mair = M0 + (Mf - M0)/(tf-t0) * (t - t0) 
            z = z0 + (zf - z0)/(tf-t0) * (t - t0)

            #Calculate heat rate at this point
            Q = tankWthermal(fuse, fuse_tank, z, Mair, xftank, ifuel)
        end
    end
    return Q
end

"""
    analyze_TASOPT_tank(ac_orig, t_hold_orig::Float64 = 0.0, t_hold_dest::Float64 = 0.0, N::Int64 = 50)

This function analyses the evolution in time of a cryogenic tank inside a TASOPT aircraft model.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac_orig::aicraft`: TASOPT aircraft model
    - `t_hold_orig::Float64`: hold at origin (s)
    - `t_hold_dest::Float64`: hold at destination (s)
    - `im::Int64`: mission index
    
    **Outputs:**
    - `ts::Vector{Float64}`: vector with times (s)
    - `ps::Vector{Float64}`: vector with pressure evolution in time (Pa)
    - `βs::Vector{Float64}`: vector with tank fill fraction evolution in time
    - `Ms::Vector{Float64}`: vector with tank fuel mass evolution in time (kg)
    - `Mburns::Vector{Float64}`: vector with cumulative mass burnt in engine (kg)
    - `Mboils::Vector{Float64}`: vector with cumulative mass that has boiled off (kg)
    - `mdot_boils::Vector{Float64}`: vector with evolution of boiloff mass flow rate in time (kg/s)
    - `Mvents::Vector{Float64}`: vector with cumulative mass that has been vented (kg)
    - `mdots::Vector{Float64}`: vector with fuel mass flow rate to engines (kg/s)
    - `Qs::Vector{Float64}`: vector with heat rate to tank (W)
    
"""
function analyze_TASOPT_tank(ac_orig::aircraft, t_hold_orig::Float64 = 0.0, t_hold_dest::Float64 = 0.0, im::Int64 = 1)
    ac = deepcopy(ac_orig) #Deepcopy original struct to avoid modifying it

    #Modify aircraft with holding times
    para_alt = zeros(size(ac.para)[1], size(ac.para)[2] + 3)
    ac.para[iatime, :, im] .= ac.para[iatime, :, im] .- minimum(ac.para[iatime, :, im])
    para_alt[:, 3:(iptotal + 1)] .= ac.para[:, 1:(size(ac.para)[2] - 1), im] #Do not copy iptest
    para_alt[iatime, 2:(iptotal + 1)] .= para_alt[iatime, 2:(iptotal + 1)] .+ t_hold_orig #Apply hold at origin
    para_alt[iatime, 1] = 0.0
    Np = size(para_alt)[2]
    para_alt[iatime, Np-1] = maximum(para_alt[iatime, :])
    para_alt[iatime, Np] = para_alt[iatime, Np-1] + t_hold_dest #Apply hold at destination

    pare_alt = zeros(size(ac.pare)[1], size(ac.pare)[2] + 3)
    pare_alt[:, 3:(iptotal + 2)] .= ac.pare[:,:, im]
    
    #Precompute heat transfer rate at each mission point for speed
    Qs_points = calc_Q_points(ac.fuselage, ac.fuse_tank, ac.pari, ac.parg, para_alt)

    #Define functions for heat and fuel burn profiles through mission 
    Q_calc(t) = find_Q_time_interp(t, para_alt, Qs_points)
    W_calc(t) = 0.0
    mdot_calc(t) = find_mdot_time(t, ac.pari, ac.parg, para_alt, pare_alt)

    #Store profiles in input struct
    u = tank_inputs(Q_calc, W_calc, mdot_calc)

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
    α = ac.fuse_tank.pfac
    #Store tank parameters in struct 
    params = tank_params(mixture_init, V, pmax, xout, xvent, α)

    ODEparams = (u, params)

    #Integrate profiles across mission
    y0 = [p0, β0, M0, 0.0, 0.0, 0.0] #Initial states

    tspan = (0.0, para_alt[iatime,end]) #start and end times

    dy_dt(y, p, t) = TankDerivatives(t, y, p[1], p[2]) #State derivatives
    #ODE problem; specify changes in mission segments for better speed
    prob = ODEProblem(dy_dt, y0, tspan, ODEparams, tstops = para_alt[iatime,:], reltol = 1e-8)
    sol = solve(prob, Tsit5()) #Solve ODE problem using recommended solver

    #Extract solution vector
    ts = sol.t 
    y = sol.u
    y = hcat(y...) #Concatenate into array

    #Produce outputs
    ps = y[1, :] #Pressure evolution
    βs = y[2, :] #Tank fill evolution
    Ms = y[3, :] #Total fluid mass in tank evolution
    Mburns = y[4, :] #Cumulative mass that has been burnt evolution
    Mvents = y[5, :] #Cumulative mass that has been vented
    Mboils = y[6, :] #Cumulative mass that has been boiled off evolution
    mdot_boils = calculate_boiloff_rate(ts, Mboils) #Calculate boiloff rate

    mdots = zeros(length(ts))
    Qs = zeros(length(ts))
    for (i,t) in enumerate(ts)
        mdots[i] = mdot_calc(t) #Fuel burn mass flow rate
        Qs[i] = Q_calc(t) #Heat transfer rate into tank
    end

    return ts, ps, βs, Ms, Mburns, Mboils, mdot_boils, Mvents, mdots, Qs
end