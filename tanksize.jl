"""
tanksize sizes the fuel tank of aircraft

Inputs:
NOTE: Every parameter is in SI units
-gee is gravitational acceleration (m/s^2)
-rhoFuel is density of fuel (kg/m^3)
-deltap is allowed pressure difference in vessel (Pa)
-Rfuse is fuselage radius (m)
-rho_insul (kg/m3) is an array of insulation layer densities
-dRfuse accounts for flatness at bottom of fuselage (m)
-xshell1 and xshell2 are start and end x-coordinates of tank (m)
-Thermal conductivity array k (W/(m*K)) comprising of k value for each MLI layer
-hconvgas (W/m2*K)  is convective coefficient of insulating purged gas (e.g. N2)
-hconvair (W/m2*K) is convective coefficient of ambient air
-Thickness array t (m) corresponds to thickness value of each layer in MLI
-h_LH2 ((W/m2*K) corresponds to LH2 convective coefficient
-Tfuel (K) is fuel temperature
-Tair (K) is ambient temperature
-r_tank (m) is tank outer radius
-h_v (J/kg) is heat of vaporization of liquid hydrogen (from Hydrogen tank design paper)
-r_gas is inner radius of gas-purged chamber (m)
-threshold_percent is the max allowed percentage of fuel that is allowed to boil off
-mode: '1' means optimize m_boiloff, '2' means find m_boiloff based on given thickness

Outputs:
- Total thickness of the insulation (m)
- Total weight of tank Wtank including fuel (N)
- Length of tank lshell (m)
- Volume of fuel Vfuel (m^3)
- Wfuel weight of fuel (N)
- m_boiloff mass (kg)
"""

function tanksize(gee, rhoFuel, deltap,
                      Rfuse, dRfuse, hconvgas, h_LH2, Tfuel, Tair,
                      h_v, t_cond, k, hconvair, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, mode, clearance_fuse, AR)

       #include("tankWmech.jl")
       #include("tankWthermal.jl")
       Wfuel_init = Wfuel
       m_boiloff = 0 #initial value of boil-off mass
       thickness_insul = sum(t_cond)

       ##To find tank weight, tank length and skin thickness of tank wall (non-insulator part)
       Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul = tankWmech(gee, rhoFuel,
                             fstring, ffadd, deltap,
                             Rfuse, dRfuse, wfb, nfweb,
                             sigskin, rho_insul, rhoskin,
                             Wfuel, m_boiloff, thickness_insul, t_cond, clearance_fuse, AR)


       m_boiloff, mdot_boiloff = tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, Rtank,
                             h_v, t_cond, k, hconvair, time_flight)


       #Wfuel = (m_boiloff * gee) + Wfuel

       Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul = tankWmech(gee, rhoFuel,
                             fstring, ffadd, deltap,
                             Rfuse, dRfuse, wfb, nfweb,
                             sigskin, rho_insul, rhoskin,
                             Wfuel, m_boiloff, thickness_insul, t_cond, clearance_fuse, AR)

        if mode == 1
                for n = 1:500 #optimize boil off mass according to threshold
                        m_boiloff, mdot_boiloff = tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, Rtank,
                                              h_v, t_cond, k, hconvair, time_flight)


                        Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul = tankWmech(gee, rhoFuel,
                                              fstring, ffadd, deltap,
                                              Rfuse, dRfuse, wfb, nfweb,
                                              sigskin, rho_insul, rhoskin,
                                              Wfuel, m_boiloff, thickness_insul, t_cond, clearance_fuse, AR)

                        Wfuel = Wfuel_init

                        if((m_boiloff / (time_flight/3600)) > (threshold_percent *  Wfuel / (gee * 100))) || break
                        end
                        t_cond[1] = t_cond[1] + 0.01 * t_cond[1]  #increase foam insulation thickness and try again
                        t_cond[3] = t_cond[3] + 0.01 * t_cond[3]
                end
        end
        #Wfuel = (m_boiloff * gee) + Wfuel
        #Wtank = Wtank #+ m_boiloff * gee + Wfuel #weight of tank including fuel

return Wtank_total, thickness_insul, lshell, mdot_boiloff, Vfuel, Wfuel_tot, m_boiloff, tskin, t_head, Rtank, Whead, Wcyl, Winsul_sum, Winsul #boiloff rate output
end
