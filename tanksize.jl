"""
tanksize sizes the fuel tank of aircraft

Inputs:
NOTE: Every parameter is in SI units
-gee is gravitational acceleration (m/s^2)
-rhoFuel is density of fuel (kg/m^3)
-deltap is allowed pressure difference in vessel (Pa)
-Rfuse is fuselage radius (m)
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
-h_e (J/kg) is heat of vaporization of liquid hydrogen (from Hydrogen tank design paper)
-r_gas is inner radius of gas-purged chamber (m)
-threshold_percent is the max allowed percentage of fuel that is allowed to boil off
-mode: '1' means optimize m_boiloff, '2' means find m_boiloff based on given thickness

Outputs:
- Total thickness of the insulation (m)
- Total weight of tank Wtank including fuel (N)
"""

function tanksize(gee, rhoFuel, deltap,
                      Rfuse, dRfuse,
                      xshell1, xshell2, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                      h_e, t_cond, r_gas, k, hconvair, time_flight, fstring,fframe,ffadd,
                      wfb, nfweb, sigskin, Wppinsul, rhoskin, Wfuel, threshold_percent, mode)

       include("tankWmech.jl")
       include("tankWthermal.jl")

       m_boiloff = 0 #initial value of boil-off mass
       thickness_insul = sum(t_cond)

       ##First guess of tank weight, tank length and skin thickness of tank wall (non-insulator part)
       result = tankWmech(gee, rhoFuel,
                             fstring, fframe, ffadd, deltap,
                             Rfuse, dRfuse, wfb, nfweb,
                             sigskin, Wppinsul, rhoskin,
                             Wfuel, m_boiloff, thickness_insul)

       Wtank = result[1]
       lshell = result[2]
       tskin = result[3]
       r_tank = result[4]

       therm = tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                             h_e, t_cond, k, hconvair, time_flight)
       m_boiloff = therm[1]
       mdot_boiloff = therm[2]

       Wfuel = (m_boiloff * gee) + Wfuel

       mech = tankWmech(gee, rhoFuel,
                             fstring, fframe, ffadd, deltap,
                             Rfuse, dRfuse, wfb, nfweb,
                             sigskin, Wppinsul, rhoskin,
                             Wfuel, m_boiloff, thickness_insul)

        Wtank = mech[1]
        Wfuel = mech[2]

        if mode == 1
                for n=1:500 #optimize boil off mass according to threshold
                        m_boiloff = tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                                              h_e, t_cond, k, hconvair, time_flight)
                        if(m_boiloff > (threshold_percent *  Wfuel / (gee * 100))) || break
                        end



                        result = tankWmech(gee, rhoFuel,
                                              fstring, fframe, ffadd, deltap,
                                              Rfuse, dRfuse, wfb, nfweb,
                                              sigskin, Wppinsul, rhoskin,
                                              Wfuel, m_boiloff, thickness_insul)
                        Wtank = mech[1]
                        Wfuel = mech[2]
                        t = t + 0.01 * t  #increase insulation thickness and try again
                end
        end


        Wtank = Wtank + m_boiloff * gee + Wfuel #weight of tank including fuel

return Wtank, thickness_insul, lshell, mdot_boiloff #boiloff rate output
end
