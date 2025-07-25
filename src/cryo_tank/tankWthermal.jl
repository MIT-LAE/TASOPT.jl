"""
      thermal_params

This structure stores the material and thermal properties of a cryogenic tank insulation layer.
"""
@kwdef mutable struct thermal_params{T<:AbstractCrossSection}
      """Heat rate (W)"""
      Q::Float64 = 0.0
      """Length of cylindrical portion of tank (m)"""
      l_cyl::Float64 = 0.0
      """Full tank length (m)"""
      l_tank::Float64 = 0.0
      """Tank radius (m)"""
      r_tank::Float64 = 0.0
      """Surface area of elliptical caps at different cross-sections (m^2)"""
      Shead::Vector{Float64} = Vector{Float64}()
      """Thickness of each insulation layer (m)"""
      t_cond::Vector{Float64} = Vector{Float64}()
      """Array with materials for different insulation layers"""
      material::Vector{ThermalInsulator} = Vector{ThermalInsulator}()
      """Fuel temperature (K)"""
      Tfuel::Float64 = 0.0
      """Flight altitude (m)"""
      z::Float64 = 0.0
      """Sea-level temperature (K)"""
      TSL::Float64 = 0.0
      """External air Mach number"""
      Mair::Float64 = 0.0
      """Longitudinal coordinate of fuel tank centroid from nose (m)"""
      xftank::Float64 = 0.0
      """Fuel species index"""
      ifuel::Int64 = 0
      """Fuselage cross section"""
      fuse_cs::T = T()
end

"""
      tankWthermal(fuse::Fuselage, fuse_tank::fuselage_tank, z::Float64, TSL::Float64, Mair::Float64, xftank::Float64, ifuel::Int64)

`tankWthermal` calculates the heat rate to a cryogenic tank for a given insulation thickness.

This function does **not** size the thermal insulation layers
but rather calculates the heat rate to the tank, 
for a given insulation thickness
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `fuse::Fuselage`: fuselage object.
      - `fuse_tank::fuselage_tank`: fuselage tank object.
      - `z::Float64`: flight altitude (m)
      - `TSL::Float64`: sea-level temperature (K)
      - `Mair::Float64`: external air Mach number
      - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
      - `ifuel::Int64`: fuel index.
      
      **Outputs:**
      - `Q::Float64`: Heat transfer rate into the tank (W).

See [here](@ref fueltanks).
"""
function tankWthermal(fuse::Fuselage, fuse_tank::fuselage_tank, z::Float64, TSL::Float64, Mair::Float64, xftank::Float64, ifuel::Int64)
      qfac = fuse_tank.qfac
      t_cond = fuse_tank.t_insul
      Tfuel = fuse_tank.Tfuel

      #Create struct with thermal parameters
      p = thermal_params{typeof(fuse.layout.cross_section)}()
      p.l_cyl = fuse_tank.l_cyl_inner
      p.l_tank = fuse_tank.l_inner
      p.r_tank = fuse_tank.Rinnertank
      p.Shead = fuse_tank.Shead_insul #Surface area of elliptical caps at different cross-sections (m^2)
      p.t_cond = t_cond
      p.material = fuse_tank.material_insul
      p.Tfuel = fuse_tank.Tfuel
      p.z = z
      p.TSL = TSL
      p.Mair = Mair
      p.xftank = xftank
      p.ifuel = ifuel
      p.fuse_cs = fuse.layout.cross_section

      _, _, Taw = freestream_heat_coeff(z, TSL, Mair, xftank) #Find adiabatic wall temperature
      
      thickness = sum(t_cond)  # total thickness of insulation
      ΔT = Taw - Tfuel
      
      #Function to zero in solver
      function residual(x) 
            try
                  return residuals_Q(x, p, "Q_unknown") #This should be zeroed
            catch #Return some high residual if it fails
                  return ones(length(x))*1e3
            end
      end

      #Initial guess for function
      guess = zeros(length(t_cond) + 2) 

      Rguess = 0.01
      guess[1] = ΔT / Rguess
      guess[2] = Tfuel + 1.0
      
      for i = 1:length(t_cond)
            guess[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ thickness
      end
      guess[end] = guess[end] - 1.0 #fuselage wall temperature
      sol = nlsolve(residual, guess, xtol = 1e-7, ftol = 1e-6) #Solve non-linear problem with NLsolve.jl
      
      Q = qfac * sol.zero[1]    # Heat rate from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
      
      return  Q
end

"""
      residuals_Q(x, p)

This function calculates the residual for a non-linear solver. The states are the heat transfer rate (optional), 
the tank wall temperature, and the temperatures at the interfaces of each insulation layers.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `x::Vector{Float64}`: vector with unknowns.
      - `p::Struct`: structure of type `thermal_params`.
      - `mode::String`: mode to find residual, options are "Q_known" and "Q_unknown"

      **Outputs:**
      - `F::Vector{Float64}`: vector with residuals.
"""
function residuals_Q(x::Vector{Float64}, p::thermal_params, mode::String)
      #Unpack states
      if mode == "Q_known" #If the heat is known because the max boiloff is an input
            Q = p.Q

            T_w = x[1]
            T_ins = x[2:end]
            Tfuse = x[end] #fuselage wall temperature
            F = zeros(Float64, length(x)+1) #initialize residual vector

      elseif mode == "Q_unknown" #heat is not known and we are solving for it too
            Q = x[1] #Heat rate is first input
            
            T_w = x[2]
            T_ins = x[3:end]
            Tfuse = x[end]
            F = zeros(Float64, length(x))
      end
      
      #Unpack parameters
      l_cyl = p.l_cyl
      l_tank = p.l_tank
      r_tank = p.r_tank
      Shead = p.Shead
      t_cond = p.t_cond
      material = p.material
      Tfuel = p.Tfuel
      z = p.z
      TSL = p.TSL
      Mair = p.Mair
      xftank = p.xftank
      fuse_cs = p.fuse_cs
      Rfuse = fuse_cs.radius
      ifuel = p.ifuel    
      
      #Calculate heat transfer coefficient, freestream temperature and adiabatic wall temperature
      h_air, _, Taw = freestream_heat_coeff(z, TSL, Mair, xftank, Tfuse, Rfuse)

      r_inner = r_tank #- thickness
      ΔT = Taw - Tfuel #Heat transfer is driven by difference between external adiabatic wall temperature and fuel temperature
      thickness = sum(t_cond)  # total thickness of insulation

      #Geometry
      perim_inner, _ = scaled_cross_section(fuse_cs, r_inner) #Tank perimeter and cross-sectional area

      S_cyl = perim_inner*l_cyl #Surface area of cylindrical portion
      S_int = S_cyl + 2*Shead[1] #liquid side surface area
      perim_R = perim_inner/r_inner #ratio of geometric shape perimeter to radius; if a circle this is 2π
  
      perim_fuse = fuse_cs.perimeter #Fuselage perimeter

      Rair = 1 / (h_air * (perim_fuse * l_tank ))  # thermal resistance of ambient air

      #Liquid-side resistance
      h_liq = tank_heat_coeff(T_w, ifuel, Tfuel, l_tank) #Find liquid-side heat transfer coefficient
      R_liq = 1 / (h_liq * S_int) #Liquid-side thermal resistance

      N = length(t_cond)       # Number of layers in insulation
      R_ins      = zeros(Float64, N)  #size of resistance vector (Based on number of layers)
      R_ins_ends = zeros(Float64, N)
      R_ins_cyl  = zeros(Float64, N)
  
      #Find resistance of each insulation layer
      T_prev = T_w
      for i in 1:N
            if lowercase(material[i].name) == "vacuum"
                  S_inner = perim_R * l_cyl * r_inner + 2*Shead[i]
                  S_outer = perim_R * l_cyl * (r_inner + t_cond[i]) + 2*Shead[i+1]
                  R_ins[i] = vacuum_resistance(T_prev, T_ins[i], S_inner, S_outer)

            else #If insulation layer is not a vacuum
                  k = thermal_conductivity(material[i], (T_ins[i] + T_prev)/2)
                  R_ins_cyl[i] = log((r_inner  + t_cond[i])/ (r_inner)) / (perim_R*l_cyl * k) #Resistance of each layer; from integration of Fourier's law in cylindrical coordinates
                  
                  Area_coeff = Shead[i] / r_inner^2 #Proportionality parameter in ellipsoidal area, approximately a constant
                  R_ins_ends[i] = t_cond[i] / (k * (Shead[i+1] + Shead[i] - Area_coeff * t_cond[i]^2)) #From closed-form solution for hemispherical caps
                  # Parallel addition of resistance
                  R_ins[i] = (R_ins_ends[i] * R_ins_cyl[i]/(R_ins_ends[i] + R_ins_cyl[i])) 
                  
                  # Update r_inner
                  r_inner = r_inner + t_cond[i]  
                  T_prev = T_ins[i]
            end
      end
  
      R_ins_tot = sum(R_ins)  #Total thermal resistance of insulation
      Req = R_ins_tot + R_liq + Rair # Total equivalent resistance of thermal circuit

      #Assemble array with residuals
      F[1] = Q - ΔT / Req #Heat transfer rate residual
  
      T_calc = Tfuel + R_liq * Q #Wall temperature residual
      F[2] = T_w - T_calc
      for i = 1:length(T_ins)
          T_calc = T_calc + R_ins[i] * Q 
          F[i + 2] = T_ins[i] - T_calc #Residual at the edge of each layer
      end

      return F
end  

"""
      tank_heat_coeff(T_w, ifuel, Tfuel, ltank)

This function calculates the liquid-side heat transfer coefficient of a cryogenic fuel in a tank. 
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `T_w::Float64`: wall temperature (K).
      - `ifuel::Int64`: fuel index.
      - `Tfuel::Float64`: temperature of fuel in fuel tank (K).
      - `ltank::Float64`: fuel tank length (m).

      **Outputs:**
      - `h_liq::Float64`: liquid-side heat tansfer coefficient (W/m^2/K).
"""
function tank_heat_coeff(T_w::Float64, ifuel::Int64, Tfuel::Float64, ltank::Float64)
      #Get fuel properties
      if ifuel == 11 #CH4
            Pr_l = 2.0 #Prandtl number, from NIST at 2atm and 120 K
            β = 3.5e-3 #K^(-1), thermal expansion coefficient https://aiche.onlinelibrary.wiley.com/doi/10.1002/aic.14254
            ν_l = 2.4e-7 #m^2/s, kinematic viscosity from NIST at 2atm and 120K 
            k = 0.17185	#W/m/K, thermal conductivity from NIST at 2atm and 120K 
        
      elseif ifuel == 40 #LH2
            Pr_l = 1.3 #from NIST at 2atm and 20K
            β = 15e-3 #https://nvlpubs.nist.gov/nistpubs/jres/089/jresv89n4p317_A1b.pdf
            ν_l = 2.0e-7 #m^2/s, from NIST at 2atm and 20K 
            k = 0.10381
      end

      Ra_l = gee * β * abs(T_w - Tfuel) * ltank^3 * Pr_l / ν_l^2 #Tank-length-based Rayleigh number

      Nu_l = 0.0605 * Ra_l^(1/3) #Length-based Nusselt number for liquid side, from https://doi.org/10.2514/6.1986-1253

      h_liq = Nu_l * k / ltank #heat transfer coefficient for liquid

      return h_liq
end

"""
      freestream_heat_coeff(z, TSL, M, xftank, Tw)

This function calculates the air-side heat transfer coefficient, which is assumed to be that of a freestream 
at a given altitude. Depending on the Mach number, either a natural or forced convection model is used.
The freestream temperature is also returned. In the forced convection case, heat transfer is modeled via the Meador-Smart 
reference temperature model with the Chilton-Colburn analogy, described on p. 1056 in Anderson, Fundamentals
of Aerodynamics. For the free convection case, the heat transfer is modeled as that on a horizontal cylinder, described on 
p. 334 in Holman.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `z::Float64`: flight altitude (m).
      - `TSL::Float64`: sea-level temperature (K)
      - `M::Float64`: freestream Mach number.
      - `xftank::Float64`: longitudinal position of the fuel tank CG (m).
      - `Tw::Float64`: wall temperature (K).

      **Outputs:**
      - `h_convair::Float64`: air-side heat transfer coefficient (W/m^2/K).
      - `Tair::Float64`: air-side temperature (K).
      - `Taw::Float64`: adiabatic wall temperature (K).
"""
function freestream_heat_coeff(z::Float64, TSL::Float64, M::Float64, xftank::Float64, Tw::Float64 = Tref, Rfuse::Float64 = 1.0)
      #Use ISA function to calculate freestream conditions
      Tair, p, _, a, _ = atmos(z / 1e3, TSL - Tref)
      u = M * a #freestrean velocity

      #Parameters for air
      R, Pr, γ, cp, _, _ = gasPr("air_simple", Tair) #This saves some computational time by using a constant cp for air
      
      r = Pr^(1/3) #recovery factor for turbulent air
      Taw = Tair * (1 + r*M^2*(γ - 1)/2)  #K, adiabatic wall temperature
      #Find h for air
      # This code uses the reference temperature model and the Chilton-Colburn analogy
      T_s = Tair * (0.5 * (1 + Tw/Tair) + 0.16 * r * (γ - 1)/2 * M^2) #Reference temperature
      
      #Find properties at reference temperature
      _, Pr_s, _, cp, μ_s, k_s = gasPr("air_simple", T_s)
      #This uses a "simple air" model with constant R and cp for speed

      ρ_s = p / (R * T_s) #density at reference temperature

      if M ≈ 0 #Natural convection
            L = 2 * Rfuse #reference length
            β = 1 / Tair #thermal expansion coefficient
            ν = μ_s / ρ_s

            Gr = gee * β * abs(Tair - Tw) * L^3 / ν^2 #Grasshoff number
            Ra = Gr * Pr_s #Rayleigh number

            # Parameters from Table 7-1 in Holton (2010)
            C = 0.13
            m = 0.333
            Nu = C * Ra^m #Nusselt number for 1e9 < Ra < 12
            hconvair = Nu * k_s / L

      else #Forced convection
            Re_xftank = ρ_s * u * xftank / μ_s
            cf_xftank = 0.02296 / (Re_xftank)^0.139 #From Meador-Smart method
            #Calculate Stanton number using Reynolds analogy
            St_air = cf_xftank / (2 * Pr_s^(2/3)) #Chilton-Colburn analogy
            hconvair = St_air * ρ_s *u* cp #In W/(m^2 K)
      end

      return hconvair, Tair, Taw
end

"""
      vacuum_resistance(Tcold, Thot, S_inner, S_outer)

This function calculates the thermal resistance of a vacuum layer using the methods of Barron (1985).
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `Tcold::Float64`: cold wall temperature (K).
      - `Thot::Float64`: hot wall temperature (K).
      - `S_inner::Float64`: area of tank inner surface (m^2).
      - `S_outer::Float64`: area of tank outer surface (m^2).

      **Outputs:**
      - `R_eq::Float64`: equivalent thermal resistance of vacuum gap (K/W).
"""
function vacuum_resistance(Tcold::Float64, Thot::Float64, S_inner::Float64, S_outer::Float64)
      #Assumed tank and gas properties
      a_outer = 0.9 
      a_inner = 1.0
      ε = 0.04    # highly polished aluminum
      p_vacuum = 1e-2 #Assumed vacuum pressure, approximately 1e-4 Torr as in Brewer (1991) #TODO maybe make this an input?

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