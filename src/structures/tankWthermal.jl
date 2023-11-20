"""
      tankWthermal(l_cyl, r_tank, Shead,
            hconvgas, h_LH2, hconvair, 
            t_cond, k, Tfuel, Tair, 
            h_v:, time_flight)

`tankWthermal` calculates the boil-off rate of LH2 for a given insulation thickness.

This subroutine does **not** size the thermal insulation layers
but rather calculates the boil-off rate of fLH2, 
for a given insulation thickness
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `l_cyl::Float64`: Length of the tank (m).
      - `r_tank::Float64`: Tank outer radius (m).
      - `Shead::Array{Float64,1}`: Array of surface areas of each layer of the end/head of the tank [m²].
      - `hconvgas::Float64`: Convective coefficient of insulating purged gas (W/m²*K).
      - `h_LH2::Float64`: LH2 convective coefficient (W/m²*K).
      - `hconvair::Float64`: Convective coefficient of ambient air (W/m²*K).
      - `t_cond::Array{Float64,1}`: Array of thickness of each layer in MLI (m).
      - `k::Array{Float64,1}`: Thermal conductivity array (W/(m*K)) comprising k values for each MLI layer.
      - `Tfuel::Float64`: Fuel temperature (K).
      - `Tair::Float64`: Ambient temperature (K).
      - `h_v::Float64`: Heat of vaporization of liquid hydrogen (J/kg).
      - `time_flight::Float64`: Time of flight (s).

      **Outputs:**
      - `m_boiloff::Float64`: Boil-off LH2 mass (kg).
      - `mdot_boiloff::Float64`: Boil-off rate of LH2 (kg/s).

See [here](@ref fueltanks).
"""
function tankWthermal(l_cyl::Float64  , r_tank::Float64, Shead::Array{Float64,1},
                      hconvgas::Float64, h_LH2::Float64,  hconvair::Float64, 
                      t_cond::Array{Float64,1}, k::Array{Float64,1},
                      Tfuel::Float64 , Tair::Float64, 
                      h_v::Float64, time_flight::Float64)

      N = length(t_cond)       # Number of layers in insulation
      thickness = sum(t_cond)  # total thickness of insulation

#--- Heat flux and resistances
      ΔT = Tair - Tfuel  # Overall temperature drop between ambient and LH2
      
      qfac = 1.3         # Account for heat leak from pipes and valves
      
      # [TODO] Move constants to constant file
      σ = 5.67e-8
      ε = 0.95    # white aircraft (Verstraete)

      hradair = σ * ε * ((Tair^2) + (Tfuel^2)) * (Tair + Tfuel)
      h_air = hconvair + hradair # Combines radiative and convective heat transfer at outer end
      Rair_conv_rad = 1 / (h_air * (2π*r_tank*l_cyl + 2*Shead[end]))  # thermal resistance of ambient air (incl. conv and rad)

      r_inner = r_tank #- thickness

      # Not needed for MLI. May add later for purged He etc. Rgas = 1 / (hconvgas * 2 * pi * r_inner * l_cyl)  #thermal resistance of purged gas

      R_LH2 = 1 / (h_LH2 * (2*π*(r_inner - thickness) * l_cyl) + 2*Shead[1]) #thermal resistance of LH2

      R_mli      = zeros(Float64, N)  #size of MLI resistance array (Based on number of layers)
      R_mli_ends = zeros(Float64, N)
      R_mli_cyl  = zeros(Float64, N)

      for i in 1:N
            R_mli_cyl[i]  = log((r_inner  + t_cond[i])/ (r_inner)) / (2π*l_cyl * k[i]) #Resistance of each MLI layer
            R_mli_ends[i] = t_cond[i] / (k[i] * 2*Shead[i])
            # Parallel addition of resistance
            R_mli[i]  = (R_mli_ends[i] * R_mli_cyl[i]/(R_mli_ends[i] + R_mli_cyl[i])) 
            
            # Update r_inner
            r_inner   = r_inner + t_cond[i]  
      end

      R_mli_tot = sum(R_mli)  #Total thermal resistance of MLI

      Req = R_mli_tot + Rair_conv_rad + R_LH2  # Total equivalent resistance of thermal circuit

      q = qfac * ΔT / Req     # Heat flux from ambient to LH2, including extra heat leak from valves etc as in eq 3.20 by Verstraete
      mdot_boiloff = q / h_v  # Boil-off rate equals the heat flux divided by heat of vaporization
      m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation

      return  m_boiloff, mdot_boiloff
end
