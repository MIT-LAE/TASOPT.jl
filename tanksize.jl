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
"""

function tanksize(gee, rhoFuel, deltap,
                      Rfuse, dRfuse, hconvgas, h_LH2, Tfuel, Tair,
                      h_v, t_cond, r_gas, k, hconvair, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent, mode)

       include("tankWmech.jl")
       include("tankWthermal.jl")

       m_boiloff = 0 #initial value of boil-off mass
       thickness_insul = sum(t_cond)

       ##To find tank weight, tank length and skin thickness of tank wall (non-insulator part)
       Wtank, lshell, tskin, r_tank, Vfuel = tankWmech(gee, rhoFuel,
                             fstring, ffadd, deltap,
                             Rfuse, dRfuse, wfb, nfweb,
                             sigskin, rho_insul, rhoskin,
                             Wfuel, m_boiloff, thickness_insul, t_cond)


       m_boiloff, mdot_boiloff = tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                             h_v, t_cond, k, hconvair, time_flight)


       Wfuel = (m_boiloff * gee) + Wfuel

       Wtank, lshell, tskin, r_tank, Vfuel = tankWmech(gee, rhoFuel,
                             fstring, ffadd, deltap,
                             Rfuse, dRfuse, wfb, nfweb,
                             sigskin, rho_insul, rhoskin,
                             Wfuel, m_boiloff, thickness_insul, t_cond)

        if mode == 1
                for n=1:500 #optimize boil off mass according to threshold
                        m_boiloff, mdot_boiloff = tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                                              h_v, t_cond, k, hconvair, time_flight)

                        Wtank, lshell, tskin, r_tank, Vfuel = tankWmech(gee, rhoFuel,
                                              fstring, ffadd, deltap,
                                              Rfuse, dRfuse, wfb, nfweb,
                                              sigskin, rho_insul, rhoskin,
                                              Wfuel, m_boiloff, thickness_insul, t_cond)

                        if(m_boiloff > (threshold_percent *  Wfuel / (gee * 100))) || break
                        end
                        t = t + 0.01 * t  #increase insulation thickness and try again
                end
        end

        Wtank = Wtank + m_boiloff * gee + Wfuel #weight of tank including fuel

return Wtank, thickness_insul, lshell, mdot_boiloff, Vfuel #boiloff rate output
end
