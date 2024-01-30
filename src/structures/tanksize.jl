"""
        tanksize(gee, rhoFuel, deltap,
                Rfuse, dRfuse, hconvgas, h_LH2, Tfuel, Tair,
                h_v, t_cond, k, hconvair, time_flight, fstring,ffadd,
                wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR)

`tanksize` sizes a cryogenic fuel tank for a LH-fueled aircraft

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `gee::Float64`: Gravitational acceleration (m/s^2).
        - `rhoFuel::Float64`: Density of fuel (kg/m^3).
        - `deltap::Float64`: Allowed pressure difference in vessel (Pa).
        - `Rfuse::Float64`: Fuselage radius (m).
        - `rho_insul::Array{Float64}`: Array of insulation layer densities (kg/m3).
        - `dRfuse::Float64`: Accounts for flatness at the bottom of the fuselage (m).
        - `xshell1::Float64`: Start x-coordinate of the tank (m).
        - `xshell2::Float64`: End x-coordinate of the tank (m).
        - `Thermal::Array{Float64}`: Conductivity array k (W/(m*K)) for each MLI layer.
        - `hconvgas::Float64`: Convective coefficient of insulating purged gas (e.g., N2) (W/m2*K).
        - `hconvair::Float64`: Convective coefficient of ambient air (W/m2*K).
        - `Thickness::Array{Float64}`: Thickness array t (m) for each MLI layer.
        - `h_LH2::Float64`: LH2 convective coefficient (W/m2*K).
        - `Tfuel::Float64`: Fuel temperature (K).
        - `Tair::Float64`: Ambient temperature (K).
        - `r_tank::Float64`: Tank outer radius (m).
        - `h_v::Float64`: Heat of vaporization of liquid hydrogen (J/kg).
        - `r_gas::Float64`: Inner radius of the gas-purged chamber (m).
        - `threshold_percent::Float64`: Max allowed percentage of fuel that is allowed to boil off.
        - `mode::Char`: '1' means optimize m_boiloff, '2' means find m_boiloff based on given thickness.

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
                      wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, ifuel)

        Wfuel_init = Wfuel
        m_boiloff = threshold_percent *  Wfuel / (gee * 100) #initial value of boil-off mass
        thickness_insul = sum(t_cond)

        residual(Δt) = res_MLI_thick(Δt, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, k, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, ifuel)

        Δt = Roots.find_zero(residual, 0.0) #Find root with Roots.jl

        t_cond[1] = t_cond[1] + Δt  #increase foam insulation thickness of layers 1 and 3
        t_cond[3] = t_cond[3] + Δt

        mdot_boiloff = threshold_percent
        Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead_insul, l_tank = tankWmech(gee, rhoFuel,
                                fstring, ffadd, deltap,
                                Rfuse, dRfuse, wfb, nfweb,
                                sigskin, rho_insul, rhoskin,
                                Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

        return Wtank_total, thickness_insul, lshell, mdot_boiloff, Vfuel, Wfuel_tot, m_boiloff, tskin, t_head, Rtank, Whead, Wcyl, Winsul_sum, Winsul, l_tank, Wtank #boiloff rate output
end

function res_MLI_thick(Δt, gee, rhoFuel, deltap,
        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
        t_cond, k, hconvair, time_flight, fstring,ffadd,
        wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, ifuel)

        t_all = deepcopy(t_cond)

        t_all[1] = t_all[1] + Δt  #increase foam insulation thickness of layers 1 and 3
        t_all[3] = t_all[3] + Δt

        m_boiloff = threshold_percent *  Wfuel / (gee * 100) #initial value of boil-off mass

        Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot,
        Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead_insul, l_tank = tankWmech(gee, rhoFuel,
                        fstring, ffadd, deltap,
                        Rfuse, dRfuse, wfb, nfweb,
                        sigskin, rho_insul, rhoskin,
                        Wfuel, m_boiloff, t_all, clearance_fuse, AR)

        m_boiloff, mdot_boiloff = tankWthermal(lshell, Rtank, Shead_insul,
        hconvgas,  hconvair, 
        t_all, k,
        Tfuel, Tair, 
        time_flight, ifuel)

        res = mdot_boiloff * gee * 3600 / Wfuel * 100 - threshold_percent #difference between current boiloff and desired
        return res
end 
