"""
tankWthermal calculates boil-off rate of LH2

Inputs:

NOTE: Everything is in SI units. MLI: Multi-layer insulation
-Thermal conductivity array k (W/(m*K)) comprising of k value for each MLI layer
-hconvgas (W/m2*K)  is convective coefficient of insulating purged gas (e.g. N2)
-hconvair (W/m2*K) is convective coefficient of ambient air
-Thickness array t_cond (m) corresponds to thickness value of each layer in MLI
-h_LH2 ((W/m2*K) corresponds to LH2 convective coefficient
-Tfuel (K) is fuel temperature
-Tair (K) is ambient temperature
-r_tank (m) is tank outer radius
-h_v (J/kg) is heat of vaporization of liquid hydrogen (from Hydrogen tank design paper)
-lshell is the length of tank (m)
-Time of flight in given segment under analysis time_flight (s)


Outputs:
- m_boiloff (kg) is the boiloff LH2 mass
- mdot_boiloff (kg/s) is the boiloff rate of LH2
"""
function tankWthermal(lshell, hconvgas, h_LH2, Tfuel, Tair, r_tank,
                      h_v, t_cond, k, hconvair, time_flight, Shead_insul)

      N = length(t_cond) #Number of layers in insulation
      thickness = sum(t_cond) #total thickness of insulation

#--- Heat flux and resistances
      ΔT = Tair - Tfuel  # Overall temperature drop between ambient and LH2
      
      qfac = 1.3         # Account for heat leak from pipes and valves
      
      # Move constants to constant file
      σ = 5.67e-8
      ε = 0.95    # white aircraft (Verstraete)

      hradair = σ * ε * ((Tair^2) + (Tfuel^2)) * (Tair + Tfuel)
      h_air = hconvair + hradair
      Rair_conv_rad = 1 / (h_air * (2π*r_tank*lshell + 2*Shead_insul[end]))  #thermal resistance of ambient air
      #r_inner = r_tank - thickness  #inner radius of tank
      r_inner = r_tank

      #Not needed for MLI. May add later for purged He etc. Rgas = 1 / (hconvgas * 2 * pi * r_inner * lshell)  #thermal resistance of purged gas

      R_LH2 = 1 / (h_LH2 * (2*π*(r_inner - thickness) * lshell) + 2*Shead_insul[1]) #thermal resistance of LH2

      R_mli     = zeros(N)  #size of MLI resistance array (Based on number of layers)
      R_mli_lat = zeros(N)
      R_mli_cyl = zeros(N)

      for n in 1:N
            R_mli_cyl[n] = log((r_inner  + t_cond[n])/ (r_inner)) / (2 *π*lshell * k[n]) #Resistance of each MLI layer
            R_mli_lat[n] = t_cond[n] / (2*k[n] * Shead_insul[n])
            R_mli[n] = (R_mli_lat[n] * R_mli_cyl[n]/(R_mli_lat[n] + R_mli_cyl[n])) # || parallel addition
            r_inner = r_inner + t_cond[n]  #Inner layer w.r.t the nth MLI layer being evaluated

      end

      R_mli_tot = sum(R_mli)  #Total thermal resistance of MLI

      Req = R_mli_tot + Rair_conv_rad + R_LH2  # Total equivalent resistance of thermal circuit

      q = qfac * ΔT / Req     # Heat flux from ambient to LH2 30% extra as in eq 3.20 by Verstraete
      mdot_boiloff = q / h_v  # Boil-off rate equals the heat flux divided by heat of vaporization
      m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation

return  m_boiloff, mdot_boiloff
end
