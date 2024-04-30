"""
        tanksize!(fuse_tank, z, TSL, Mair, xftank,
        time_flight, ifuel)

`tanksize` sizes a cryogenic fuel tank for a cryogenic-fuel aircraft

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `fuse_tank::Struct`: structure of type `fuselage_tank` with parameters.
        - `z::Float64`: flight altitude (m)
        - `TSL::Float64`: sea-level temperature (K)
        - `Mair::Float64`: external air Mach number
        - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
        - `ifuel::Int64`: fuel index.
        
        **Outputs:**
        - `mdot_boiloff::Float64`: Mass boiled off during the mission flight (kg).
        - `Vfuel::Float64`: Volume of fuel (m^3).
        - `Rtank::Float64`: Radius of the tank (m).
        - `Winsul_sum::Float64`: Sum of the insulation weight (N).
        - `l_tank::Float64`: Length of the tank (m).
        - `Wtank::Float64`: Weight of the tank structure (N).

See [here](@ref fueltanks).
"""
function tanksize!(fuse_tank, z::Float64, TSL::Float64, Mair::Float64, xftank::Float64,
                      time_flight::Float64, ifuel::Int64)

        #Unpack variables in fuse_tank
        t_cond = fuse_tank.t_insul
        Wfuel = fuse_tank.Wfuelintank
        Tfuel = fuse_tank.Tfuel
        flag_size = fuse_tank.size_insulation #Boolean for whether to size for a boiloff rate

        if flag_size #If insulation is sized for a given boiloff rate
                boiloff_percent = fuse_tank.boiloff_rate
                iinsuldes = fuse_tank.iinsuldes

                #Thermal calculations
                _, _, Taw = freestream_heat_coeff(z, TSL, Mair, xftank) #Find adiabatic wall temperature 

                ΔT = Taw - Tfuel

                #Create inline function with residuals as a function of x
                residual(x) = res_MLI_thick(x, fuse_tank, z, TSL, Mair, xftank, ifuel) #Residual in boiloff rate as a function of Δt
                #Assemble guess for non linear solver
                #x[1] = Δt; x[2] = T_tank; x[3:(end-1)]: T at edge of insulation layer; x[end] = T at fuselage wall
                guess = zeros(length(t_cond) + 2) 
                guess[1] = 0.0
                guess[2] = Tfuel + 1.0

                n_layers = length(t_cond)
                for i = 1:n_layers
                        guess[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ sum(t_cond)
                end
                #Solve non-linear problem with NLsolve.jl
                sol = nlsolve(residual, guess, ftol = 1e-7) 
                Δt = sol.zero[1] 

                for ind in iinsuldes #For every segment whose thickness can be changed
                        t_cond[ind] = t_cond[ind] + Δt #This will modify fuse_tank.t_insul
                end

                mdot_boiloff = boiloff_percent *  Wfuel / (gee * 100) /3600
                m_boiloff = boiloff_percent *  Wfuel / (gee * 100)*time_flight/3600 #initial value of boil-off mass
                
        else #If insulation does not need to be sized
                m_boiloff, mdot_boiloff = tankWthermal(fuse_tank, z, TSL, Mair, xftank, time_flight, ifuel)
        end
        
        thickness_insul = sum(t_cond)
        
        #Evaluate tank weight
        Winner_tot, lcyl1, tskin, Rinnertank, Vfuel, Winnertank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Wstiff, Winsul,
        Sinternal, Shead_insul, l_inner = size_inner_tank(fuse_tank, fuse_tank.t_insul)

        if ("vacuum" in fuse_tank.material_insul) || ("Vacuum" in fuse_tank.material_insul) #If tank is double-walled
                Routertank = fuse_tank.Rfuse - fuse_tank.clearance_fuse
                lcyl2 = lcyl1 * Routertank / Rinnertank #Scale outer vessel length for geometric similarity
                
                Ninterm = optimize_outer_tank(fuse_tank, Winner_tot, lcyl2) #Find optimal number of intermediate stiffeners
                
                fuse_tank.Ninterm = Ninterm #Store in fuse_tank to use as guess in next wsize iteration

                Wtank2, Wcyl2, Whead2, Wstiff2, Souter, Shead2, Scyl2, 
                t_cyl2, t_head2, l_outer = size_outer_tank(fuse_tank, Winner_tot, lcyl2, Ninterm)

                Wtank_total = Winner_tot + Wtank2
                Wtank = Winnertank + Wtank2
                l_tank = l_outer #If there is an outer vessel, the total length is the length of this tank
                Rtank = Routertank

        else
                Wtank_total = Winner_tot
                Wtank = Winnertank
                l_tank = l_inner #Tank length when there is only an inner vessel
                Rtank = Rinnertank
        end
        
        return mdot_boiloff, Vfuel, Rtank, Winsul_sum, l_tank, Wtank
end

"""
        res_MLI_thick(x, fuse_tank, z, TSL, Mair, xftank, ifuel)

This function evaluates the residual vector for a given state containing change in wall thickness, heat transfer rate and 
insulation interface temperatures.

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `x::Float64`: vector with states
        - `fuse_tank::Struct`: structure of type `fuselage_tank` with parameters.
        - `z::Float64`: flight altitude (m)
        - `TSL::Float64`: sea-level temperature (K)
        - `Mair::Float64`: external air Mach number
        - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
        - `ifuel::Int64`: fuel index.

        **Outputs:**
        - `res::Vector{Float64}`: residuals vector.
"""
function res_MLI_thick(x::Vector{Float64}, fuse_tank, z::Float64, TSL::Float64, Mair::Float64, xftank::Float64, ifuel::Int64)

        #Extract parameters from fuse_tank
        boiloff_percent = fuse_tank.boiloff_rate
        t_cond = fuse_tank.t_insul
        qfac = fuse_tank.qfac
        iinsuldes = fuse_tank.iinsuldes
        Tfuel = fuse_tank.Tfuel
        Wfuel = fuse_tank.Wfuelintank
        h_v = fuse_tank.hvap #heat of vaporization

        # Extract states
        Δt = x[1]
        x_thermal = x[2:end]

        #Prepare to call residual_Q for thermal-related residuals
        t_all = deepcopy(t_cond) #copy to avoid modifying input

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_all[ind] = t_all[ind] + Δt
        end

        Wtank_total, l_cyl, tskin, r_tank, Vfuel, Wtank, Wfuel_tot,
        Winsul_sum, t_head, Whead, Wcyl, Wstiff, Winsul, Sinternal, Shead, l_tank = size_inner_tank(fuse_tank, t_all)

        mdot_boiloff = boiloff_percent *  Wfuel / (gee * 100) / 3600  
        # Boil-off rate equals the heat rate divided by heat of vaporization
        Q_net = mdot_boiloff * h_v  # Heat rate from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
        Q = Q_net / qfac

        #Assemble struct with parameters for residual_Q
        p = thermal_params()
        p.Q = Q #Store heat rate as it is known
        p.l_cyl = l_cyl
        p.l_tank = l_tank
        p.r_tank = r_tank
        p.Shead = Shead
        p.t_cond = t_all
        p.material = fuse_tank.material_insul
        p.Tfuel = Tfuel
        p.z = z
        p.TSL = TSL
        p.Mair = Mair
        p.xftank = xftank
        p.ifuel = ifuel

        res = residuals_Q(x_thermal, p, "Q_known") #Find thermal-related residuals
        return res
end 