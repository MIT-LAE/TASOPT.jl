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
      - `h_liq::Float64`: LH2 convective coefficient (W/m²*K).
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
                      hconvgas::Float64, h_liq::Float64,  hconvair::Float64, 
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
      Rair_conv_rad = 1 / (h_air * (2π * r_tank * l_cyl + 2*Shead[end]))  # thermal resistance of ambient air (incl. conv and rad)

      r_inner = r_tank #- thickness

      # Not needed for MLI. May add later for purged He etc. Rgas = 1 / (hconvgas * 2 * pi * r_inner * l_cyl)  #thermal resistance of purged gas

      R_liq = 1 / (h_liq * (2π * (r_inner - thickness) * l_cyl) + 2*Shead[1]) #thermal resistance of LH2

      R_mli      = zeros(Float64, N)  #size of MLI resistance array (Based on number of layers)
      R_mli_ends = zeros(Float64, N)
      R_mli_cyl  = zeros(Float64, N)

      for i in 1:N #TODO: check this carefully
            R_mli_cyl[i]  = log((r_inner  + t_cond[i])/ (r_inner)) / (2π*l_cyl * k[i]) #Resistance of each MLI layer
            R_mli_ends[i] = t_cond[i] / (k[i] * 2*Shead[i])
            # Parallel addition of resistance
            R_mli[i]  = (R_mli_ends[i] * R_mli_cyl[i]/(R_mli_ends[i] + R_mli_cyl[i])) 
            
            # Update r_inner
            r_inner   = r_inner + t_cond[i]  
      end

      R_mli_tot = sum(R_mli)  #Total thermal resistance of MLI

      Req = R_mli_tot + Rair_conv_rad + R_liq  # Total equivalent resistance of thermal circuit

      q = qfac * ΔT / Req     # Heat flux from ambient to LH2, including extra heat leak from valves etc as in eq 3.20 by Verstraete
      mdot_boiloff = q / h_v  # Boil-off rate equals the heat flux divided by heat of vaporization
      m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation

      return  m_boiloff, mdot_boiloff
end

"""
      tank_heat_coeffs(z, M, ifuel, Tfuel, ltank, xftank)

This function calculates the latent heat of vaporization of a liquid fuel and its liquid-side heat transfer coefficient. 
The heat leaks through some material towards an outside gas, which is assumed to be a freestream in forced convection. The
function also returns the gas-side heat transfer coefficient and freestream temperature.

This subroutine does **not** size the thermal insulation layers
but rather calculates the boil-off rate of fLH2, 
for a given insulation thickness
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `z::Float64`: flight altitude (m).
      - `M::Float64`: freestream Mach number.
      - `ifuel::Int64`: fuel index.
      - `h_liq::Float64`: LH2 convective coefficient (W/m²*K).
      - `Tfuel::Float64`: temperature of fuel in fuel tank (K).
      - `ltank::Float64`: fuel tank length (m).
      - `Tair::Float64`: Ambient temperature (K).
      - `xftank::Float64`: longitudinal position of the fuel tank CG (m).

      **Outputs:**
      - `h_liq::Float64`: liquid-side heat tansfer coefficient (W/m^2/K).
      - `h_convair::Float64`: air-side heat transfer coefficient (W/m^2/K).
      - `h_v::Float64`: liquid's enthalpy of vaporization (J/kg).
      - `Tair::Float64`: air-side temperature (K).
"""
function tank_heat_coeffs(z, M, ifuel, Tfuel, ltank, xftank)
      #Get fuel properties
      if ifuel == 11 #CH4
            h_v =  510e3 #J/kg, latent heat of vaporozation
            Pr_l = 2.0 #Prandtl number, from NIST at 2atm and 120 K
            β = 3.5e-3 #K^(-1), thermal expansion coefficient https://aiche.onlinelibrary.wiley.com/doi/10.1002/aic.14254
            ν_l = 2.4e-7 #m^2/s, kinematic viscosity from NIST at 2atm and 120K 
            k = 0.17185	#W/m/K, thermal conductivity from NIST at 2atm and 120K 
        
      elseif ifuel == 40 #LH2
            h_v =  447e3
            Pr_l = 1.3 #from NIST at 2atm and 20K
            β = 15e-3 #https://nvlpubs.nist.gov/nistpubs/jres/089/jresv89n4p317_A1b.pdf
            ν_l = 2.0e-7 #m^2/s, from NIST at 2atm and 20K 
            k = 0.10381
      end

      #Use ISA function to calculate freestream conditions
      Tair, p, ρ, a, μ = atmos(z / 1e3)

      Ra_l = gee * β * (Tair - Tfuel) * ltank^3 * Pr_l / ν_l^2 #Tank-length-based Rayleigh number

      Nu_l = 0.0605 * Ra_l^(1/3) #Length-based Nusselt number for liquid side, from https://doi.org/10.2514/6.1986-1253

      h_liq = Nu_l * k / ltank #heat transfer coefficient for liquid

      #Find h for air
      cp_air = 1004 #J/(kg K) specific heat at constant pressure for ambient air
      Pr_t = 0.85 #Prandlt number of turbulent air
      u = M * a #freestrean velocity

      Re_xftank = ρ * u * xftank / μ
      Cf_xftank = 0.0576 / Re_xftank^0.2 #Turbulent flat plate skin friction coefficient       
      #Calculate Stanton number using Reynolds analogy
      St_air = Cf_xftank / (2 * Pr_t^(2/3)) #Chilton-Colburn analogy
      hconvair = St_air * ρ *u* cp_air #In W/(m^2 K)

      return h_liq, hconvair, h_v, Tair
end
