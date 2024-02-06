"""
      tankWthermal(l_cyl::Float64, l_tank::Float64, r_tank::Float64, Shead::Array{Float64,1},
      hconvgas::Float64,  hconvair::Float64, 
      t_cond::Array{Float64,1}, k::Array{Float64,1},
      Tfuel::Float64 , Tair::Float64, 
      time_flight::Float64, ifuel::Int64)

`tankWthermal` calculates the boil-off rate of a cryogenic liquid for a given insulation thickness.

This subroutine does **not** size the thermal insulation layers
but rather calculates the boil-off rate of the fuel, 
for a given insulation thickness
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `l_cyl::Float64`: Length of cylindrical portion of the tank (m).
      - `l_tank::Float64`: Tank total length (m).
      - `r_tank::Float64`: Tank outer radius (m).
      - `Shead::Array{Float64,1}`: Array of surface areas of each layer of the end/head of the tank [m²].
      - `hconvgas::Float64`: Convective coefficient of insulating purged gas (W/m²*K).
      - `hconvair::Float64`: Convective coefficient of ambient air (W/m²*K).
      - `t_cond::Array{Float64,1}`: Array of thickness of each layer in MLI (m).
      - `k::Array{Float64,1}`: Thermal conductivity array (W/(m*K)) comprising k values for each MLI layer.
      - `Tfuel::Float64`: Fuel temperature (K).
      - `Tair::Float64`: Ambient temperature (K).
      - `time_flight::Float64`: Time of flight (s).
      - `ifuel::Int64`: fuel index.

      **Outputs:**
      - `m_boiloff::Float64`: Boil-off LH2 mass (kg).
      - `mdot_boiloff::Float64`: Boil-off rate of LH2 (kg/s).

See [here](@ref fueltanks).
"""
function tankWthermal(l_cyl::Float64, l_tank::Float64, r_tank::Float64, Shead::Array{Float64,1},
                      hconvgas::Float64,  hconvair::Float64, 
                      t_cond::Array{Float64,1}, k::Array{Float64,1},
                      Tfuel::Float64 , Tair::Float64, 
                      time_flight::Float64, ifuel::Int64)

      p = thermal_params()
      p.l_cyl = l_cyl
      p.l_tank = l_tank
      p.r_tank = r_tank
      p.Shead = Shead
      p.hconvgas = hconvgas
      p.hconvair = hconvair
      p.t_cond = t_cond
      p.material = "rohacell"
      p.Tfuel = Tfuel
      p.Tair = Tair
      p.ifuel = ifuel
      
      thickness = sum(t_cond)  # total thickness of insulation
      ΔT = Tair - Tfuel
      qfac = 1.3         # Account for heat leak from pipes and valves
      
      fun(x) = residuals_Q(x, p)
      
      guess = zeros(length(t_cond) + 2)
      guess[1] = 100
      guess[2] = Tfuel + 1
      
      for i = 1:length(t_cond)
            guess[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ thickness
      end
      sol = nlsolve(fun, guess)
      
      _, h_v = tank_heat_coeffs(Tfuel, ifuel, Tfuel, l_tank) #Liquid side h and heat of vaporization
      
      Q = qfac * sol.zero[1]    # Heat flux from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
      mdot_boiloff = Q / h_v  # Boil-off rate equals the heat flux divided by heat of vaporization
      m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation
      
      return  m_boiloff, mdot_boiloff
end


"""
      res_q_tank(T_w, ΔT, S_int, R_eq_ext, ifuel, Tfuel, ltank)

This function calculates the difference between the liquid-side heat transfer and the overall heat transfer for
a given wall temperature. This residual should be 0 at the correct wall temperature. 
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `T_w::Float64`: wall temperature (K).
      - `ΔT::Float64`: difference between liquid in tank and outside air.
      - `S_int::Float64`: liquid-side surface area (m^2).
      - `R_eq_ext::Float64`: total thermal resistance of wall and outside air (m^2 K / W).
      - `ifuel::Int64`: fuel index.
      - `Tfuel::Float64`: temperature of fuel in fuel tank (K).
      - `ltank::Float64`: fuel tank length (m).

      **Outputs:**
      - `res::Float64`: residual (W/m^2).
"""
function residuals_Q(x, p)

      Q = x[1]
      T_w = x[2]
      T_mli = x[3:end]
  
      #Unpack parameters
      l_cyl = p.l_cyl
      l_tank = p.l_tank
      r_tank = p.r_tank
      Shead = p.Shead
      hconvgas = p.hconvgas
      hconvair = p.hconvair
      t_cond = p.t_cond
      material = p.material
      Tfuel = p.Tfuel
      Tair = p.Tair
      ifuel = p.ifuel
  
      r_inner = r_tank #- thickness
      ΔT = Tair - Tfuel
      thickness = sum(t_cond)  # total thickness of insulation
  
      # Radiation
      σ = 5.6704e-8 #W/(m^2 K^4), Stefan-Boltzmann constant
      ε = 0.95    # white aircraft (Verstraete)
  
      hradair = σ * ε * ((Tair^2) + (Tfuel^2)) * (Tair + Tfuel) #Radiative heat transfer coefficient; Eq. (2.28) in https://ahtt.mit.edu/
      h_air = hconvair + hradair # Combines radiative and convective heat transfer at outer end
      Rair_conv_rad = 1 / (h_air * (2π * r_tank * l_cyl + 2*Shead[end]))  # thermal resistance of ambient air (incl. conv and rad)
  
      S_int = (2π * (r_inner - thickness) * l_cyl) + 2*Shead[1] #liquid side surface area
      h_liq, _ = tank_heat_coeffs(T_w, ifuel, Tfuel, l_tank) #Find liquid-side heat transfer coefficient
      R_liq = 1 / (h_liq * S_int) #Liquid-side thermal resistance
  
      N = length(t_cond)       # Number of layers in insulation
      R_mli      = zeros(Float64, N)  #size of MLI resistance array (Based on number of layers)
      R_mli_ends = zeros(Float64, N)
      R_mli_cyl  = zeros(Float64, N)
  
      T_prev = T_w
      for i in 1:N
          k = insulation_conductivity_calc((T_mli[i] + T_prev)/2, material)
          R_mli_cyl[i] = log((r_inner  + t_cond[i])/ (r_inner)) / (2π*l_cyl * k) #Resistance of each MLI layer; from integration of Fourier's law in cylindrical coordinates
          R_mli_ends[i] = t_cond[i] / (k * 2*Shead[i])
          # Parallel addition of resistance
          R_mli[i] = (R_mli_ends[i] * R_mli_cyl[i]/(R_mli_ends[i] + R_mli_cyl[i])) 
          
          # Update r_inner
          r_inner = r_inner + t_cond[i]  
          T_prev = T_mli[i]
      end
  
      R_mli_tot = sum(R_mli)  #Total thermal resistance of MLI
      Req = R_mli_tot + R_liq + Rair_conv_rad # Total equivalent resistance of thermal circuit
  
      F = zeros(length(x))
      F[1] = Q - ΔT / Req
  
      T_calc = Tfuel + R_liq * Q
      F[2] = T_w - T_calc
  
      for i = 1:length(T_mli)
          T_calc = T_calc + R_mli[i] * Q
          F[i + 2] = T_mli[i] - T_calc
          
      end
      return F
end  

function insulation_conductivity_calc(T, material)
      if material == "rohacell"
          k = 0.00235 + 8.824e-5 * T  # W m/(K) Linear fit to Fig. 4.78 in Brewer (1991)
      end
      return k
end
  
mutable struct thermal_params
      l_cyl::Float64
      l_tank::Float64
      r_tank::Float64
      Shead::Array{Float64}
      hconvgas::Float64
      hconvair::Float64
      t_cond::Array{Float64} 
      material::String
      Tfuel::Float64
      Tair::Float64
      ifuel::Int64
      thermal_params() = new() 
end

"""
      tank_heat_coeffs(T_w, ifuel, Tfuel, ltank)

This function calculates the latent heat of vaporization of a liquid fuel and its liquid-side heat transfer coefficient. 
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `T_w::Float64`: wall temperature (K).
      - `ifuel::Int64`: fuel index.
      - `Tfuel::Float64`: temperature of fuel in fuel tank (K).
      - `ltank::Float64`: fuel tank length (m).

      **Outputs:**
      - `h_liq::Float64`: liquid-side heat tansfer coefficient (W/m^2/K).
      - `h_v::Float64`: liquid's enthalpy of vaporization (J/kg).
"""
function tank_heat_coeffs(T_w, ifuel, Tfuel, ltank)
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

      Ra_l = gee * β * abs(T_w - Tfuel) * ltank^3 * Pr_l / ν_l^2 #Tank-length-based Rayleigh number

      Nu_l = 0.0605 * Ra_l^(1/3) #Length-based Nusselt number for liquid side, from https://doi.org/10.2514/6.1986-1253

      h_liq = Nu_l * k / ltank #heat transfer coefficient for liquid

      return h_liq, h_v
end

"""
      freestream_heat_coeff(z, M, xftank)

This function calculates the air-side heat transfer coefficient, which is assumed to be that of a freestream 
in forced convection. The freestream temperature is also returned.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `z::Float64`: flight altitude (m).
      - `M::Float64`: freestream Mach number.
      - `xftank::Float64`: longitudinal position of the fuel tank CG (m).

      **Outputs:**
      - `h_convair::Float64`: air-side heat transfer coefficient (W/m^2/K).
      - `Tair::Float64`: air-side temperature (K).
"""
function freestream_heat_coeff(z, M, xftank)
      #Use ISA function to calculate freestream conditions
      Tair, p, ρ, a, μ = atmos(z / 1e3)

      #Find h for air
      cp_air = 1004 #J/(kg K) specific heat at constant pressure for ambient air
      Pr_t = 0.85 #Prandlt number of turbulent air
      u = M * a #freestrean velocity

      Re_xftank = ρ * u * xftank / μ
      Cf_xftank = 0.455 / (log10(Re_xftank)^2.58 * (1 + 0.144 * M^2)^0.65) #Turbulent flat plate skin friction coefficient (12.27 in Raymer)
      #Calculate Stanton number using Reynolds analogy
      St_air = Cf_xftank / (2 * Pr_t^(2/3)) #Chilton-Colburn analogy
      hconvair = St_air * ρ *u* cp_air #In W/(m^2 K)

      return hconvair, Tair
end
