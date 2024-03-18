"""
      tankWthermal(l_cyl::Float64, l_tank::Float64, r_tank::Float64, Shead::Array{Float64,1}, material_insul::Array{String,1},
      hconvgas::Float64,
      t_cond::Array{Float64,1},
      Tfuel::Float64, z::Float64, Mair::Float64, xftank::Float64,
      time_flight::Float64, ifuel::Int64, qfac::Float64)

`tankWthermal` calculates the boil-off rate of a cryogenic liquid for a given insulation thickness.

This function does **not** size the thermal insulation layers
but rather calculates the boil-off rate of the fuel, 
for a given insulation thickness
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `l_cyl::Float64`: Length of cylindrical portion of the tank (m).
      - `l_tank::Float64`: Tank total length (m).
      - `r_tank::Float64`: Tank outer radius (m).
      - `Shead::Array{Float64,1}`: Array of surface areas of each layer of the end/head of the tank [m²].
      - `material_insul::Array{String,1}`: material name for each MLI layer.
      - `hconvgas::Float64`: Convective coefficient of insulating purged gas (W/m²*K).
      - `t_cond::Array{Float64,1}`: Array of thickness of each layer in MLI (m).
      - `Tfuel::Float64`: Fuel temperature (K).
      - `z::Float64`: flight altitude (m)
      - `Mair::Float64`: external air Mach number
      - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
      - `time_flight::Float64`: Time of flight (s).
      - `ifuel::Int64`: fuel index.
      - `qfac::Float64`: Factor to multiply heat tranfer rate by to account for heat leakae through structure, piping, etc

      **Outputs:**
      - `m_boiloff::Float64`: Boil-off LH2 mass (kg).
      - `mdot_boiloff::Float64`: Boil-off rate of LH2 (kg/s).

See [here](@ref fueltanks).
"""
function tankWthermal(l_cyl::Float64, l_tank::Float64, r_tank::Float64, Shead::Array{Float64,1}, material_insul::Array{String,1},
                      hconvgas::Float64,
                      t_cond::Array{Float64,1},
                      Tfuel::Float64, z::Float64, Mair::Float64, xftank::Float64,
                      time_flight::Float64, ifuel::Int64, qfac::Float64)

      p = thermal_params()
      p.l_cyl = l_cyl
      p.l_tank = l_tank
      p.r_tank = r_tank
      p.Shead = Shead
      p.hconvgas = hconvgas
      p.t_cond = t_cond
      p.material = material_insul
      p.Tfuel = Tfuel
      p.z = z
      p.Mair = Mair
      p.xftank = xftank
      p.ifuel = ifuel

      _, _, Taw = freestream_heat_coeff(z, Mair, xftank, 200) #Find adiabatic wall temperature with dummy wall temperature
      
      thickness = sum(t_cond)  # total thickness of insulation
      ΔT = Taw - Tfuel
      
      fun(x) = residuals_Q(x, p, "Q_unknown") #Create function handle to be zeroed
      
      #Initial guess for function
      guess = zeros(length(t_cond) + 2) 
      guess[1] = 1000
      guess[2] = Tfuel + 1
      
      for i = 1:length(t_cond)
            guess[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ thickness
      end
      sol = nlsolve(fun, guess, xtol = 1e-7, ftol = 1e-6) #Solve non-linear problem with NLsolve.jl
      
      _, h_v = tank_heat_coeffs(Tfuel, ifuel, Tfuel, l_tank) #Liquid heat of vaporization
      
      Q = qfac * sol.zero[1]    # Heat rate from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
      mdot_boiloff = Q / h_v  # Boil-off rate equals the heat rate divided by heat of vaporization
      m_boiloff = mdot_boiloff * time_flight # Boil-off mass calculation
      
      return  m_boiloff, mdot_boiloff
end

"""
      residuals_Q(x, p)

This function calculates the residual for a non-linear solver. The states are the heat transfer rate (optional), 
the tank wall temperature, and the temperatures at the interfaces of MLI insulation layers.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `x::Vector{Float64}`: vector with unknowns.
      - `p::Struct`: structure of type `thermal_params`.
      - `mode::String`: mode to find residual, options are "Q_known" and "Q_unknown"

      **Outputs:**
      - `F::Vector{Float64}`: vector with residuals.
"""
function residuals_Q(x, p, mode)
      #Unpack states
      if mode == "Q_known" #If the heat is known because the max boiloff is an input
            Q = p.Q

            T_w = x[1]
            T_mli = x[2:end]
            Tfuse = x[end] #fuselage wall temperature
            F = zeros(length(x)+1) #initialize residual vector

      elseif mode == "Q_unknown" #heat is not known and we are solving for it too
            Q = x[1] #Heat rate is first input
            
            T_w = x[2]
            T_mli = x[3:end]
            Tfuse = x[end]
            F = zeros(length(x))
      end
      
      #Unpack parameters
      l_cyl = p.l_cyl
      l_tank = p.l_tank
      r_tank = p.r_tank
      Shead = p.Shead
      hconvgas = p.hconvgas
      t_cond = p.t_cond
      material = p.material
      Tfuel = p.Tfuel
      z = p.z
      Mair = p.Mair
      xftank = p.xftank
      ifuel = p.ifuel    
      
      #Calculate heat transfer coefficient, freestream temperature and adiabatic wall temperature
      hconvair, _, Taw = freestream_heat_coeff(z, Mair, xftank, Tfuse)
  
      r_inner = r_tank #- thickness
      ΔT = Taw - Tfuel #Heat transfer is driven by difference between external adiabatic wall temperature and fuel temperature
      thickness = sum(t_cond)  # total thickness of insulation
  
      # Radiation
      ε = 0.95    # white aircraft (Verstraete)
  
      hradair = σ_SB * ε * ((Taw^2) + (Tfuse^2)) * (Taw + Tfuse) #Radiative heat transfer coefficient; Eq. (2.28) in https://ahtt.mit.edu/
      h_air = hconvair + hradair # Combines radiative and convective heat transfer at outer end
      Rair_conv_rad = 1 / (h_air * (2π * (r_tank + thickness) * l_cyl + 2*Shead[end]))  # thermal resistance of ambient air (incl. conv and rad)
  
      S_int = (2π * (r_inner) * l_cyl) + 2*Shead[1] #liquid side surface area
      h_liq, _ = tank_heat_coeffs(T_w, ifuel, Tfuel, l_tank) #Find liquid-side heat transfer coefficient
      R_liq = 1 / (h_liq * S_int) #Liquid-side thermal resistance

      N = length(t_cond)       # Number of layers in insulation
      R_mli      = zeros(Float64, N)  #size of MLI resistance array (Based on number of layers)
      R_mli_ends = zeros(Float64, N)
      R_mli_cyl  = zeros(Float64, N)
  
      #Find resistance of each MLI layer
      T_prev = T_w
      for i in 1:N
            if lowercase(material[i]) == "vacuum"
                  S_inner = 2π * l_cyl * r_inner
                  S_outer = 2π * l_cyl * (r_inner + t_cond[i])
                  R_mli[i] = vacuum_resistance(T_prev, T_mli[i], S_inner, S_outer)
            else #If insulation layer is not a vacuum
                  k = insulation_conductivity_calc((T_mli[i] + T_prev)/2, material[i])
                  R_mli_cyl[i] = log((r_inner  + t_cond[i])/ (r_inner)) / (2π*l_cyl * k) #Resistance of each MLI layer; from integration of Fourier's law in cylindrical coordinates
                  
                  Area_coeff = Shead[i] / r_inner^2 #Proportionality parameter in ellipsoidal area, approximately a constant
                  R_mli_ends[i] = t_cond[i] / (k * (Shead[i+1] + Shead[i] - Area_coeff * t_cond[i]^2)) #From closed-form solution for hemispherical caps
                  # Parallel addition of resistance
                  R_mli[i] = (R_mli_ends[i] * R_mli_cyl[i]/(R_mli_ends[i] + R_mli_cyl[i])) 
                  
                  # Update r_inner
                  r_inner = r_inner + t_cond[i]  
                  T_prev = T_mli[i]
            end
      end
  
      R_mli_tot = sum(R_mli)  #Total thermal resistance of MLI
      Req = R_mli_tot + R_liq + Rair_conv_rad # Total equivalent resistance of thermal circuit

      #Assemble array with residuals
      F[1] = Q - ΔT / Req #Heat transfer rate residual
  
      T_calc = Tfuel + R_liq * Q #Wall temperature residual
      F[2] = T_w - T_calc
      for i = 1:length(T_mli)
          T_calc = T_calc + R_mli[i] * Q 
          F[i + 2] = T_mli[i] - T_calc #Residual at the edge of each MLI layer
      end

      return F
end  

"""
      insulation_conductivity_calc(T, material)

This function calculates the thermal conductivity of different insulation materials as a function of temperature.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `T::Float64`: temperature (K).
      - `material::String`: material name.

      **Outputs:**
      - `k::Float64`: thermal conductivity (W/(m K)).
"""
function insulation_conductivity_calc(T, material)
      if material == "rohacell41S"
            #Note: Brewer (1991) only had one point for rohacell thermal conductivity. They assumed the same curve as PVC
            k = 0.001579 + 1.283e-4 * T - 3.353e-7*T^2 + 8.487e-10 * T^3 # W/(m K), polynomial fit to Fig. 4.78 in Brewer (1991) between 20 and 320 K
      elseif material == "polyurethane27" #polyurethane with density of 27 kg/m^3
            k = 2.114e-13 * T^5 - 1.639e-10 *T^4 + 4.438e-8 * T^3 - 5.222e-6*T^2 + 3.740e-4*T - 2.192e-3
            # W/(m K), polynomial fit to Fig. 4.78 in Brewer (1991) between 20 and 320 K
      elseif material == "polyurethane32" #polyurethane with density of 32 kg/m^3
            k = 2.179E-13 * T^5 - 1.683E-10* T^4 + 4.542E-08* T^3 - 5.341E-06* T^2 + 3.816E-04* T - 2.367E-03
            # W/(m K), polynomial fit to Fig. 4.78 in Brewer (1991) between 20 and 320 K
      elseif material == "polyurethane35" #polyurethane with density of 35 kg/m^3
            k = 2.104E-13* T^5 - 1.695E-10* T^4 + 4.746E-08* T^3 - 5.662E-06* T^2 + 3.970E-04* T - 2.575E-03
            # W/(m K), polynomial fit to Fig. 4.78 in Brewer (1991) between 20 and 320 K
      end
      return k
end

"""
      thermal_params

This structure stores the material and thermal properties of a cryogenic tank insulation layer.
      
!! details "💾 Data fields"
    **Inputs:**
    - `Q::Float64`: heat rate (W)
    - `l_cyl::Float64`: length of cylindrical portion of tank (m)
    - `l_tank::Float64`: full tank length (m)
    - `r_tank::Float64`: tank radius (m)
    - `Shead::Array{Float64}`: surface area of elliptical caps at different cross-sections (m^2)
    - `hconvgas::Float64`: convective heat transfer coefficient across purged gas layer (W / (m^2 K))
    - `hconvair::Float64 `: convective heat transfer coefficient on fuselage (W / (m^2 K))
    - `material::Array{String} `: array with material names for different insulation layers
    - `Tfuel::Float64`: fuel temperature (K)
    - `z::Float64`: flight altitude (m)
    - `Mair::Float64`: external air Mach number
    - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
    - `ifuel::Int64`: fuel species index
"""
mutable struct thermal_params
      Q::Float64
      l_cyl::Float64
      l_tank::Float64
      r_tank::Float64
      Shead::Array{Float64}
      hconvgas::Float64
      t_cond::Array{Float64} 
      material::Array{String}
      Tfuel::Float64
      z::Float64
      Mair::Float64
      xftank::Float64
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
in forced convection at a given altitude. The freestream temperature is also returned. Heat transfer is modeled via the
Meador-Smart reference temperature model with the Chilton-Colburn analogy, described on p. 1056 in Anderson, Fundamentals
of Aerodynamics.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `z::Float64`: flight altitude (m).
      - `M::Float64`: freestream Mach number.
      - `xftank::Float64`: longitudinal position of the fuel tank CG (m).
      - `Tw::Float64`: wall temperature (K).

      **Outputs:**
      - `h_convair::Float64`: air-side heat transfer coefficient (W/m^2/K).
      - `Tair::Float64`: air-side temperature (K).
      - `Taw::Float64`: adiabatic wall temperature (K).
"""
function freestream_heat_coeff(z, M, xftank, Tw)
      #Use ISA function to calculate freestream conditions
      Tair, p, _, a, _ = atmos(z / 1e3)
      u = M * a #freestrean velocity

      #Assumed parameters for air
      Pr = 0.71
      γ = 1.4
      cp = 1005 #J/(kg K)
      R = 287 #J/(kg K)
      
      r = Pr^(1/3) #recovery factor for turbulent air
      Taw = Tair * (1 + r*M^2*(γ - 1)/2)  #K, adiabatic wall temperature
      #Find h for air
      # This code uses the reference temperature model and the Chilton-Colburn analogy
      T_s = Tair * (0.5 * (1 + Tw/Tair) + 0.16 * r * (γ - 1)/2 * M^2) #Reference temperature
      
      #Find viscosity from Sutherland's law
      μ0 = 1.716e-5
      S_μ = 111
      T0 = 273
      μ_s = μ0 * (T_s / T0)^(3 / 2) * ( (T0 + S_μ) / (T_s + S_μ) )

      ρ_s = p / (R * T_s) #density at reference temperature

      Re_xftank = ρ_s * u * xftank / μ_s
      cf_xftank = 0.02296 / (Re_xftank)^0.139 #From Meador-Smart method
      #Calculate Stanton number using Reynolds analogy
      St_air = cf_xftank / (2 * Pr^(2/3)) #Chilton-Colburn analogy
      hconvair = St_air * ρ_s *u* cp #In W/(m^2 K)

      return hconvair, Tair, Taw
end

function vacuum_resistance(Tcold, Thot, S_inner, S_outer)
      #Assumed tank and gas properties
      a_outer = 0.9 
      a_inner = 1.0
      ε = 0.04    # highly polished aluminum
      p_vacuum = 1e-3 #Assumed vacuum pressure #TODO maybe make this an input?

      Rgas = 287.05  # specific gas constant
      gamma = 1.4
      
      #Emissivity factor
      Fe = 1 / (1/ε + S_inner/S_outer * (1/ε - 1))

      #Find heat transfer coeff. for radiation and corresponding resistance
      hrad = σ_SB * Fe * ((Tcold^2) + (Thot^2)) * (Tcold + Thot) #Radiative heat transfer coefficient; Eq. (2.28) in https://ahtt.mit.edu/
      R_rad = 1/(hrad * S_inner)  # radiative resistance

      #Find resistance due to convection by residual air
      Fa = 1 / (1/a_inner + S_inner/S_outer * (1/a_outer - 1)) #accomodation factor
      G = (gamma + 1)/(gamma - 1) * sqrt(Rgas / (8*pi * Thot)) * Fa
      R_conv = 1 / (G * p_vacuum * S_inner)  # convective resistance due to imperfect vacuum
      
      #Parallel addition of resistance
      R_eq = R_conv * R_rad / (R_conv + R_rad) 
      return R_eq
end