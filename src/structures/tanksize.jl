"""
        tanksize(gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel)

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
        - `Tair::Float64`: Ambient temperature (K).
        - `t_cond::Array{Float64}`: Thickness array t (m) for each MLI layer.
        - `hconvair::Float64`: Convective coefficient of ambient air (W/m2*K).
        - `time_flight::Float64`: total flight time (s)
        - `fstring::Float64`: mass factor to account for stiffening material.
        - `ffadd::Float64`: Additional mass factor for the tank.
        - `wfb::Float64`: parameter for multi-bubble configuration.
        - `nfweb::Float64`: Number of bubbles.
        - `sigskin::Float64`: Material property.
        - `material_insul::Array{String,1}`: material name for each MLI layer.
        - `rhoskin::Float64`: Material property.
        - `Wfuel::Float64`: Weight of fuel (N).
        - `threshold_percent::Float64`: Max allowed percentage of fuel that is allowed to boil off (%/hour).
        - `clearance_fuse::Float64`: Clearance for the fuselage (m).
        - `AR::Float64`: Aspect ratio.
        - `iinsuldes::Array{Int64}`: indices for insulation layers to be sized.
        - `ifuel::Int64`: fuel index.

        
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
                      Rfuse, dRfuse, hconvgas, Tfuel, Tair,
                      t_cond, hconvair, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
                      iinsuldes, ifuel)

        Wfuel_init = Wfuel
        m_boiloff = threshold_percent *  Wfuel / (gee * 100) #initial value of boil-off mass

        residual(x) = res_MLI_thick(x, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel) #Residual in boiloff rate as a function of Δt

        ΔT = Tair - Tfuel
        guess = zeros(length(t_cond) + 3) 
        guess[1] = 0.0
        guess[2] = 100.0
        guess[3] = Tfuel + 1
        for i = 1:length(t_cond)
                guess[i + 3] = Tfuel + ΔT * sum(t_cond[1:i])/ sum(t_cond)
        end

        sol = nlsolve(residual, guess, ftol = 1e-7)
        Δt = sol.zero[1] #Solve for change in layer thickness with NLsolve.jl

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_cond[ind] = t_cond[ind] + Δt  
        end
        thickness_insul = sum(t_cond)

        #Evaluate tank weight
        mdot_boiloff = threshold_percent
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
        res_MLI_thick(Δt, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel)

`tanksize` sizes a fuel tank for a cryogenic-fuel aircraft

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `Δt::Float64`: thickness change to each insulation layer in iinsuldes (m).
        - `gee::Float64`: Gravitational acceleration (m/s^2).
        - `rhoFuel::Float64`: Density of fuel (kg/m^3).
        - `deltap::Float64`: Allowed pressure difference in vessel (Pa).
        - `Rfuse::Float64`: Fuselage radius (m).
        - `dRfuse::Float64`: Accounts for flatness at the bottom of the fuselage (m).
        - `hconvgas::Float64`: Convective coefficient of insulating purged gas (e.g., N2) (W/m2*K).
        - `Tfuel::Float64`: Fuel temperature (K).
        - `Tair::Float64`: Ambient temperature (K).
        - `t_cond::Array{Float64}`: Thickness array t (m) for each MLI layer.
        - `hconvair::Float64`: Convective coefficient of ambient air (W/m2*K).
        - `time_flight::Float64`: total flight time (s)
        - `fstring::Float64`: mass factor to account for stiffening material.
        - `ffadd::Float64`: Additional mass factor for the tank.
        - `wfb::Float64`: parameter for multi-bubble configuration.
        - `nfweb::Float64`: Number of bubbles.
        - `sigskin::Float64`: Material property.
        - `material_insul::Array{String,1}`: material name for each MLI layer.
        - `rhoskin::Float64`: Material property.
        - `Wfuel::Float64`: Weight of fuel (N).
        - `threshold_percent::Float64`: Max allowed percentage of fuel that is allowed to boil off (%/hour).
        - `clearance_fuse::Float64`: Clearance for the fuselage (m).
        - `AR::Float64`: Aspect ratio.
        - `iinsuldes::Array{Int64}`: indices for insulation layers to be sized.
        - `ifuel::Int64`: fuel index.

        **Outputs:**
        - `res::Float64`: difference between desired boiloff rate and current boiloff rate (%/hour).
"""
function res_MLI_thick(x, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel)

        # Extract states
        Δt = x[1]
        Q = x[2]
        x_thermal = x[2:end]

        #Prepare to call residual_Q for thermal-related residuals
        t_all = deepcopy(t_cond) #copy to avoid modifying input

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_all[ind] = t_all[ind] + Δt  
        end

        m_boiloff = threshold_percent *  Wfuel / (gee * 100) #value of boil-off mass

        Wtank_total, l_cyl, tskin, r_tank, Vfuel, Wtank, Wfuel_tot,
        Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead, l_tank = tankWmech(gee, rhoFuel,
                        fstring, ffadd, deltap,
                        Rfuse, dRfuse, wfb, nfweb,
                        sigskin, material_insul, rhoskin,
                        Wfuel, m_boiloff, t_all, clearance_fuse, AR)

        #Assemble struct with parameters for residual_Q
        p = thermal_params()
        p.l_cyl = l_cyl
        p.l_tank = l_tank
        p.r_tank = r_tank
        p.Shead = Shead
        p.hconvgas = hconvgas
        p.hconvair = hconvair
        p.t_cond = t_all
        p.material = material_insul
        p.Tfuel = Tfuel
        p.Tair = Tair
        p.ifuel = ifuel

        F_thermal = residuals_Q(x_thermal, p) #Find thermal-related residuals
                        
        _, h_v = tank_heat_coeffs(Tfuel, ifuel, Tfuel, l_tank) #Liquid heat of vaporization
      
        qfac = 1.3
        Q_net = qfac * Q    # Heat rate from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
        mdot_boiloff = Q_net / h_v  # Boil-off rate equals the heat rate divided by heat of vaporization
        m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation

        #Create array with all residuals
        res = zeros(length(x))
        res[1] = mdot_boiloff * gee * 3600 / Wfuel * 100 - threshold_percent #difference between current boiloff and desired
        res[2:end] .= F_thermal 
        return res
end 
