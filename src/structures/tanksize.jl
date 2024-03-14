"""
        tanksize(gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, z, Mair, xftank,
        t_cond, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel, qfac)

`tanksize` sizes a cryogenic fuel tank for a cryogenic-fuel aircraft

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `gee::Float64`: Gravitational acceleration (m/s^2).
        - `rhoFuel::Float64`: Density of fuel (kg/m^3).
        - `deltap::Float64`: Allowed pressure difference in vessel (Pa).
        - `Rfuse::Float64`: Fuselage radius (m).
        - `dRfuse::Float64`: Accounts for flatness at the bottom of the fuselage (m).
        - `hconvgas::Float64`: Convective coefficient of insulating purged gas (e.g., N2) (W/m2*K).
        - `Tfuel::Float64`: Fuel temperature (K).
        - `z::Float64`: flight altitude (m)
        - `Mair::Float64`: external air Mach number
        - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
        - `t_cond::Vector{Float64}`: Thickness array t (m) for each MLI layer.
        - `time_flight::Float64`: total flight time (s)
        - `fstring::Float64`: mass factor to account for stiffening material.
        - `ffadd::Float64`: Additional mass factor for the tank.
        - `wfb::Float64`: parameter for multi-bubble configuration.
        - `nfweb::Float64`: Number of bubbles.
        - `sigskin::Float64`: Material property.
        - `material_insul::Vector{String}`: material name for each MLI layer.
        - `rhoskin::Float64`: Material property.
        - `Wfuel::Float64`: Weight of fuel (N).
        - `threshold_percent::Float64`: Max allowed percentage of fuel that is allowed to boil off (%/hour).
        - `clearance_fuse::Float64`: Clearance for the fuselage (m).
        - `AR::Float64`: Aspect ratio.
        - `iinsuldes::Vector{Int64}`: indices for insulation layers to be sized.
        - `ifuel::Int64`: fuel index.
        - `qfac::Float64`: Factor to multiply heat tranfer rate by to account for heat leakae through structure, piping, etc

        
        **Outputs:**
        - `Wtank_total::Float64`: Total weight of the tank including fuel (N).
        - `thickness_insul::Float64`: Total thickness of the insulation (m).
        - `lshell::Float64`: Length of the tank (m).
        - `mdot_boiloff::Float64`: Mass boiled off during the mission flight (kg).
        - `Vfuel::Float64`: Volume of fuel (m^3).
        - `Wfuel_tot::Float64`: Weight of fuel (N).
        - `m_boiloff::Float64`: Mass boiled off (kg).
        - `tskin::Float64`: Thickness of the tank's skin (m).
        - `t_head::Float64`: Thickness of the tank's head (m).
        - `Rtank::Float64`: Radius of the tank (m).
        - `Whead::Float64`: Weight of the tank's head (N).
        - `Wcyl::Float64`: Weight of the tank's cylinder (N).
        - `Winsul_sum::Float64`: Sum of the insulation weight (N).
        - `Winsul::Float64`: Weight of insulation (N).
        - `l_tank::Float64`: Length of the tank (m).
        - `Wtank::Float64`: Weight of the tank structure (N).

See [here](@ref fueltanks).
"""
function tanksize(gee, rhoFuel, deltap,
                      Rfuse, dRfuse, hconvgas, Tfuel, z, Mair, xftank,
                      t_cond, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
                      iinsuldes, ifuel, qfac)

        Wfuel_init = Wfuel
        m_boiloff = threshold_percent *  Wfuel / (gee * 100)*time_flight/3600 #initial value of boil-off mass

        #Create inline function with residuals as a function of x
        residual(x) = res_MLI_thick(x, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, z, Mair, xftank,
        t_cond, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel, qfac) #Residual in boiloff rate as a function of Δt

        _, _, Taw = freestream_heat_coeff(z, Mair, xftank, 200) #Find adiabatic wall temperature with dummy wall temperature

        ΔT = Taw - Tfuel

        #Assemble guess for non linear solver
        #x[1] = Δt; x[2] = T_tank; x[3:(end-1)]: T at edge of insulation layer; x[end] = T at fuselage wall
        guess = zeros(length(t_cond) + 2) 
        guess[1] = 0.0
        guess[2] = Tfuel + 1

        n_layers = length(t_cond)
        for i = 1:n_layers
                guess[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ sum(t_cond)
        end
        #Solve non-linear problem with NLsolve.jl
        sol = nlsolve(residual, guess, ftol = 1e-7) 
        Δt = sol.zero[1] 

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_cond[ind] = t_cond[ind] + Δt  
        end
        thickness_insul = sum(t_cond)

        #Evaluate tank weight
        mdot_boiloff = threshold_percent *  Wfuel / (gee * 100) /3600
        Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul,
        Shead_insul, l_tank = tankWmech(gee, rhoFuel,
                                fstring, ffadd, deltap,
                                Rfuse, dRfuse, wfb, nfweb,
                                sigskin, material_insul, rhoskin,
                                Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

        return Wtank_total, thickness_insul, lshell, mdot_boiloff, 
        Vfuel, Wfuel_tot, m_boiloff, tskin, t_head, Rtank, Whead,
        Wcyl, Winsul_sum, Winsul, l_tank, Wtank 
end

"""
        res_MLI_thick(x, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, z, Mair, xftank,
        t_cond, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel, qfac)

This function evaluates the residual vector for a given state containing change in wall thickness, heat transfer rate and 
insulation interface temperatures.

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `x::Float64`: vector with states
        - `gee::Float64`: Gravitational acceleration (m/s^2).
        - `rhoFuel::Float64`: Density of fuel (kg/m^3).
        - `deltap::Float64`: Allowed pressure difference in vessel (Pa).
        - `Rfuse::Float64`: Fuselage radius (m).
        - `dRfuse::Float64`: Accounts for flatness at the bottom of the fuselage (m).
        - `hconvgas::Float64`: Convective coefficient of insulating purged gas (e.g., N2) (W/m2*K).
        - `Tfuel::Float64`: Fuel temperature (K).
        - `z::Float64`: flight altitude (m)
        - `Mair::Float64`: external air Mach number
        - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
        - `t_cond::Vector{Float64}`: Thickness array t (m) for each MLI layer.
        - `time_flight::Float64`: total flight time (s)
        - `fstring::Float64`: mass factor to account for stiffening material.
        - `ffadd::Float64`: Additional mass factor for the tank.
        - `wfb::Float64`: parameter for multi-bubble configuration.
        - `nfweb::Float64`: Number of bubbles.
        - `sigskin::Float64`: Material property.
        - `material_insul::Vector{String,1}`: material name for each MLI layer.
        - `rhoskin::Float64`: Material property.
        - `Wfuel::Float64`: Weight of fuel (N).
        - `threshold_percent::Float64`: Max allowed percentage of fuel that is allowed to boil off (%/hour).
        - `clearance_fuse::Float64`: Clearance for the fuselage (m).
        - `AR::Float64`: Aspect ratio.
        - `iinsuldes::Vector{Int64}`: indices for insulation layers to be sized.
        - `ifuel::Int64`: fuel index.
        - `qfac::Float64`: Factor to multiply heat tranfer rate by to account for heat leakae through structure, piping, etc

        **Outputs:**
        - `res::Vector{Float64}`: residuals vector.
"""
function res_MLI_thick(x, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, z, Mair, xftank,
        t_cond, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel, qfac)

        # Extract states
        Δt = x[1]
        x_thermal = x[2:end]

        #Prepare to call residual_Q for thermal-related residuals
        t_all = deepcopy(t_cond) #copy to avoid modifying input

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_all[ind] = t_all[ind] + Δt
        end

        m_boiloff = threshold_percent *  Wfuel / (gee * 100) *time_flight/3600 #value of boil-off mass

        Wtank_total, l_cyl, tskin, r_tank, Vfuel, Wtank, Wfuel_tot,
        Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead, l_tank = tankWmech(gee, rhoFuel,
                        fstring, ffadd, deltap,
                        Rfuse, dRfuse, wfb, nfweb,
                        sigskin, material_insul, rhoskin,
                        Wfuel, m_boiloff, t_all, clearance_fuse, AR)

        _, h_v = tank_heat_coeffs(Tfuel, ifuel, Tfuel, l_tank) #Liquid heat of vaporizatio

        mdot_boiloff = threshold_percent *  Wfuel / (gee * 100) / 3600  
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
        p.hconvgas = hconvgas
        p.t_cond = t_all
        p.material = material_insul
        p.Tfuel = Tfuel
        p.z = z
        p.Mair = Mair
        p.xftank = xftank
        p.ifuel = ifuel

        res = residuals_Q(x_thermal, p, "Q_known") #Find thermal-related residuals
        return res
end 