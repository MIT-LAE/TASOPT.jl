"""
        tanksize(gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, k, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
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
        - `k::Array{Float64}`: Conductivity array k (W/(m*K)) for each MLI layer.
        - `hconvair::Float64`: Convective coefficient of ambient air (W/m2*K).
        - `time_flight::Float64`: total flight time (s)
        - `fstring::Float64`: mass factor to account for stiffening material.
        - `ffadd::Float64`: Additional mass factor for the tank.
        - `wfb::Float64`: parameter for multi-bubble configuration.
        - `nfweb::Float64`: Number of bubbles.
        - `sigskin::Float64`: Material property.
        - `rho_insul::Array{Float64}`: Array of insulation layer densities (kg/m3).
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
                      t_cond, k, hconvair, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
                      iinsuldes, ifuel)

        Wfuel_init = Wfuel
        m_boiloff = threshold_percent *  Wfuel / (gee * 100) #initial value of boil-off mass

        residual(Δt) = res_MLI_thick(Δt, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, k, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel) #Residual as a function of Δt

        Δt = Roots.find_zero(residual, t_cond[1]*1e-2) #Find root with Roots.jl

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_cond[ind] = t_cond[ind] + Δt  
        end
        thickness_insul = sum(t_cond)

        mdot_boiloff = threshold_percent
        Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul,
        Shead_insul, l_tank = tankWmech(gee, rhoFuel,
                                fstring, ffadd, deltap,
                                Rfuse, dRfuse, wfb, nfweb,
                                sigskin, rho_insul, rhoskin,
                                Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

        return Wtank_total, thickness_insul, lshell, mdot_boiloff, 
        Vfuel, Wfuel_tot, m_boiloff, tskin, t_head, Rtank, Whead,
        Wcyl, Winsul_sum, Winsul, l_tank, Wtank 
end

"""
        res_MLI_thick(Δt, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, k, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel)

`tanksize` sizes a cryogenic fuel tank for a cryogenic-fuel aircraft

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
        - `k::Array{Float64}`: Conductivity array k (W/(m*K)) for each MLI layer.
        - `hconvair::Float64`: Convective coefficient of ambient air (W/m2*K).
        - `time_flight::Float64`: total flight time (s)
        - `fstring::Float64`: mass factor to account for stiffening material.
        - `ffadd::Float64`: Additional mass factor for the tank.
        - `wfb::Float64`: parameter for multi-bubble configuration.
        - `nfweb::Float64`: Number of bubbles.
        - `sigskin::Float64`: Material property.
        - `rho_insul::Array{Float64}`: Array of insulation layer densities (kg/m3).
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
function res_MLI_thick(Δt, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, k, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
        iinsuldes, ifuel)

        t_all = deepcopy(t_cond) #copy to avoid modifying input

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_all[ind] = t_all[ind] + Δt  
        end

        m_boiloff = threshold_percent *  Wfuel / (gee * 100) #initial value of boil-off mass

        Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot,
        Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead_insul, l_tank = tankWmech(gee, rhoFuel,
                        fstring, ffadd, deltap,
                        Rfuse, dRfuse, wfb, nfweb,
                        sigskin, rho_insul, rhoskin,
                        Wfuel, m_boiloff, t_all, clearance_fuse, AR)

        m_boiloff, mdot_boiloff = tankWthermal(lshell, l_tank, Rtank, Shead_insul,
        hconvgas,  hconvair, 
        t_all, k,
        Tfuel, Tair, 
        time_flight, ifuel)

        res = mdot_boiloff * gee * 3600 / Wfuel * 100 - threshold_percent #difference between current boiloff and desired
        return res
end 
