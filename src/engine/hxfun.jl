# hxfun.jl
# These functions can be used to size and model a heat exchanger with involute staggered tubes in a crossflow
# The design method is based on the effectiveness-NTU method, described in many sources such as 
# https://www.mathworks.com/help/hydro/ref/entuheattransfer.html

const alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127] #Air composition

"""
    HX_gas

Structure containing the fluid properties of the process (`_p`) and coolant (`_c`) streams.
"""
@kwdef mutable struct HX_gas
      "Process fluid name"
      fluid_p :: String = ""
      "Coolant fluid name"
      fluid_c :: String = "" 
      "Process gas composition"
      alpha_p :: Vector{Float64} = []
      "Coolant gas index, if coolant is a gas"
      igas_c :: Float64 = 0.0
      "Mass flow rate of process gas (kg/s)"
      mdot_p :: Float64 = 0.0
      "Mass flow rate of coolant gas (kg/s)"
      mdot_c :: Float64 = 0.0
      "Process gas inlet temperature (K)"
      Tp_in :: Float64 = 0.0
      "Coolant gas inlet temperature (K)"
      Tc_in :: Float64 = 0.0
      "Process gas inlet pressure (Pa)"
      pp_in :: Float64 = 0.0
      "Coolant gas inlet pressure (Pa)"
      pc_in :: Float64 = 0.0
      "Process gas inlet Mach number"
      Mp_in  :: Float64 = 0.0
      "Coolant gas inlet Mach number"
      Mc_in :: Float64 = 0.0
      "Process gas outlet temperature (K)"
      Tp_out :: Float64 = 0.0
      "Coolant gas outlet temperature (K)"
      Tc_out :: Float64 = 0.0
      "Wall temperature (K)"
      Tw :: Float64 = 0.0
      "Enthalpy change across HX (J/kg)"
      Δh_p :: Float64 = 0.0
      "Enthalpy change across coolant (J/kg)"
      Δh_c :: Float64 = 0.0
      "Pressure drop of process gas across heat exchanger (Pa)"
      Δp_p :: Float64 = 0.0
      "Pressure drop of coolant gas across tubes (Pa)"
      Δp_c :: Float64 = 0.0
      "Power loss due to pressure drop in process stream (W)"
      Pl_p :: Float64 = 0.0
      "Power loss due to pressure drop in coolant stream (W)"
      Pl_c :: Float64 = 0.0
      "Desired heat exchanger effectiveness"
      ε :: Float64 = 0.0 
      "Temperature of recirculating flow at HX inlet (K)"
      recircT :: Float64 = 0.0 
      "Recirculating flow mass flow rate (kg/s)"
      mdot_r :: Float64 = 0.0 
      "Latent heat capacity in freestream coolant liquid (J/kg)"
      h_lat :: Float64 = 0.0 
      "Power required to pump recirculating flow (W)"
      P_recirc :: Float64 = 0.0
end
# Overload Base.getproperty for convenience
function Base.getproperty(HXgas::HX_gas, sym::Symbol)
      if (sym === :Q) 
            return abs(getfield(HXgas, :mdot_c) * getfield(HXgas, :Δh_c))
      else
         return getfield(HXgas, sym)
      end
end

"""
    HX_tubular

Structure containing the heat exchanger geometric and material properties.
"""
@kwdef mutable struct HX_tubular
      "Flag for concentric geometry (true: concentric ; false: rectangular)"
      is_concentric :: Bool = false
      "Flag for recirculation (true: recirculation ; false: no recirculation)"
      has_recirculation :: Bool = false
      "Flag for whether HX contains shaft (true: shaft ; false: no shaft)"
      has_shaft :: Bool = false
      "Number of tubes per row"
      N_t :: Float64 = 0.0
      "Number of different coolant stages with different coolant flows"
      n_stages :: Float64 = 0.0
      "Number of coolant passes"
      n_passes:: Float64 = 0.0
      "Process side freestream cross-sectional area (m^2)"
      A_cs:: Float64  = 0.0
      "Length of tubes (m)"
      l :: Float64 = 0.0
      "Cooling tube wall thickness (m)"
      t :: Float64 = 0.0
      "Tube outer diameter (m)"
      tD_o :: Float64 = 0.0
      "Circumferential pitch between tubes at the root over tube outer diameter"
      xt_D :: Float64 = 0.0
      "Longitudinal pitch between rows over tube outer diameter"
      xl_D :: Float64 = 0.0
      "Process-side fouling factor (m^2 K/W)"
      Rfp :: Float64 = 0.0
      "Coolant-side fouling factor (m^2 K/W)"
      Rfc :: Float64 = 0.0
      "Inner diameter of core (m)"
      D_i :: Float64 = 0.0
      "Material"
      material :: StructuralAlloy = StructuralAlloy("Al-2219-T87")
      "Design pressure difference between tube and outside (Pa)"
      Δpdes::Float64 = 0.0
      "Maximum allowable HEX length (m)"
      maxL::Float64 = 0.0
end

# Overload Base.getproperty for convenience
function Base.getproperty(HXgeom::HX_tubular, sym::Symbol)
      if (sym === :D_o) && getfield(HXgeom, :is_concentric)
            A_cs = getfield(HXgeom, :A_cs)
            D_i = getfield(HXgeom, :D_i)
            D_o = sqrt(4 * A_cs / pi +D_i^2)
         return D_o
      elseif sym === :N_tubes_tot
            return getfield(HXgeom, :n_stages) * getfield(HXgeom, :n_passes) * getfield(HXgeom, :N_t)
      elseif sym === :tD_i
            return getfield(HXgeom, :tD_o) - 2* getfield(HXgeom, :t) 
      elseif sym === :L
            return getfield(HXgeom, :n_stages) * getfield(HXgeom, :n_passes) * getfield(HXgeom, :xl_D) * getfield(HXgeom, :tD_o)
      else
         return getfield(HXgeom, sym)
      end
end

"""
    HeatExchanger

Structure containing all the heat exchanger geometry and operational information.
"""
@kwdef mutable struct HeatExchanger
      "HX type"
      type :: String = ""
      "Geometry object"
      HXgeom :: HX_tubular = HX_tubular()
      "Array of HX_gas objects for each mission point"
      HXgas_mission :: Array{HX_gas} = HX_gas[]
      "Order"
      order::Int64 = 0
      "Design-point effectiveness"
      design_effectiveness::Float64 = 0.0 #Design effectiveness, used for sizing
      "Design point inlet Mach number"
      design_Mach::Float64 = 0.0 #Design Mach number, used for sizing
      "Maximum allowable length"
      maximum_length::Float64 = 0.0 #Maximum allowable HX length (m)
      "Added mass fraction"
      added_mass_fraction::Float64 = 0.0 #Added mass fraction due to hubs and structures
      "Flag for recirculation"
      has_recirculation::Bool = false #Flag for recirculation
      "Recirculating fluid temperature"
      recirculation_temperature::Float64 = 0.0 #temperature of recirculating flow at HX inlet (K)
      "Minimum wall temperature in the sizing mission"
      min_wall_temperature::Float64 = 0.0 #Minimum wall temperature (K)
end

# Outer constructor with custom size for HXgas_mission
function make_HeatExchanger(nmis::Int; kwargs...)
    obj = HeatExchanger(; kwargs...)  # initialize with other kwargs
    obj.HXgas_mission = Array{HX_gas}(undef, iptotal, nmis)
    return obj
end
  
"""
    hxsize!(HXgas, HXgeom)

Sizes a crossflow heat exchanger and calculates the pressure drop. Uses the ε-NTU method to size the heat exchanger
from a prescribed ε. For representative fouling factors see Standards of the Tubular Exchanger Manufacturers Association
or `https://powderprocess.net/Tools_html/Data_Diagrams/Heat_Exchanger_Fouling_Factor.html`

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    
    **Outputs:**
    No direct outputs. Input structures are modified with outlet gas properties and HX design geometry.
"""
function hxsize!(HXgas::HX_gas, HXgeom::HX_tubular)
      #---------------------------------
      # Extract inputs
      #---------------------------------
      #Gas parameters
      fluid_p = HXgas.fluid_p #Process-side fluid name
      fluid_c = HXgas.fluid_c #coolant fluid name
      alpha_p = HXgas.alpha_p #Process-side gas composition
      igas_c = HXgas.igas_c #Coolant-side gas index
      mdot_p = HXgas.mdot_p #mass flow rate of process gas (kg/s)
      mdot_c_inf = HXgas.mdot_c #mass flow rate of coolant gas (kg/s)
      ε = HXgas.ε 
      Tp_in = HXgas.Tp_in
      Tc_inf = HXgas.Tc_in 
      Mp_in  = HXgas.Mp_in
      Mc_in = HXgas.Mc_in
      pp_in = HXgas.pp_in
      pc_in = HXgas.pc_in

      #Flags 
      is_concentric = HXgeom.is_concentric
      has_recirculation = HXgeom.has_recirculation

      #HX geometry
      t = HXgeom.t #Initial wall thickness, may be overwritten
      n_stages = HXgeom.n_stages
      xt_D = HXgeom.xt_D
      xl_D = HXgeom.xl_D
      Rfp = HXgeom.Rfp
      Rfc = HXgeom.Rfc 
      l = HXgeom.l

      tol = 1e-10 #convergence tolerance

      if is_concentric #If geometry is concentric
            D_i = HXgeom.D_i
      end

      if has_recirculation
            recircT = HXgas.recircT
            h_lat = HXgas.h_lat

      end

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------
      _, _, hp_in, _, cp_p_in, Rp = gassum(alpha_p, length(alpha_p), Tp_in)

      if occursin("liquid", fluid_c)
            cp_c_inf = liquid_properties(fluid_c, Tc_inf).cp
            hc_inf = cp_c_inf * Tc_inf #TODO: replace with enthalpy calculation for non-constant cp
      else
            _, _, hc_inf, _, cp_c_inf, Rc = gasfun(igas_c, Tc_inf)
      end
      #---------------------------------
      # Thermal calculations
      #---------------------------------

      C_p = mdot_p * cp_p_in #Process-side heat capacity rate
      mdot_c = mdot_c_inf
      modot_c_prev = mdot_c

      if has_recirculation #If there is recirculation
            _, _, hc_in, _, cp_c_in, _ = gasfun(igas_c, recircT)
            Tc_in = recircT
            N_iters_rec = 15 #Iterations in recirculation loop
            mdot_r = 0 #Initialize

            for j =1:N_iters_rec
                  #Calculate assuming that C_c = C_min
                  mdot_r = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * cp_c_in * (Tp_in -  recircT) )
            
                  C_check = (mdot_c_inf + mdot_r) * cp_c_in #Coolant heat capacity rate, to check validity of above assumption

                  C_max = max(C_check, C_p)

                  if C_check == C_max #If the calculation above is incorrect because C_c = C_max
                        A = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * C_p * (Tp_in -  recircT) )
                        
                        if A > 1
                              error("Insufficient heat capacity in process stream")
                              return
                        end
                        mdot_r = A * mdot_c_inf / (1 - A)
                  end

                  mdot_c = mdot_c_inf + mdot_r #update coolant mass flow rate

                  if (abs(modot_c_prev - mdot_c)/mdot_c < tol)
                        break #Break for loop if convergence has been reached
                  end 
                  modot_c_prev = mdot_c #otherwise store current value for comparison
            end
      
      else #No recirculation
            Tc_in = Tc_inf
            hc_in = hc_inf
            cp_c_in = cp_c_inf
      end
      
      C_c = mdot_c * cp_c_in #Coolant heat capacity rate

      C_min = min(C_c, C_p)
      C_max = max(C_c, C_p)
      C_r = C_min / C_max

      if C_c == C_min
            ε_max = 1 / C_r * (1 - exp(-C_r)) #At ε = ε_max, NTU tends to infinity
            if ε > ε_max
                  ε = 0.99 * ε_max #Limit effectiveness to 99% of maximum possible one
                                    #effectiveness is limited to avoid an error in the NTU formula
                                    #User is warned if this happens at the end of the weight loop
            end
            NTU = -log(1 + log(1 - C_r * ε) / C_r) # For cross-flow with C_max mixed and C_min unmixed
      else
            ε_max = 1 - exp(-1 / C_r)#At ε = ε_max, NTU tends to infinity
            if ε > ε_max
                  ε = 0.99 * ε_max #Limit effectiveness to 99% of maximum possible one
            end
            NTU = -1 / C_r * log(1 + C_r * log(1 - ε) ) # For cross-flow with C_max unmixed and C_min mixed
      end

      # Calculate total heat transfer and exit temperatures
      Qmax = C_min * (Tp_in - Tc_in) #Maximum heat transfer rate
      Q = ε * Qmax #Actual heat transfer rate

      Tp_out_guess = Tp_in - Q / C_p
      Tc_out_guess = Tc_in + Q / C_c

      hp_out = hp_in - Q / mdot_p
      hc_out = hc_in + Q / mdot_c

      Δh_p = hp_out - hp_in
      Δh_c = hc_out - hc_in

      Tp_out = gas_tset(alpha_p, length(alpha_p), hp_out, Tp_out_guess)

      if occursin("liquid", fluid_c)
            Tc_out = Tc_out_guess
      else
            Tc_out = gas_tset_single(igas_c, hc_out, Tc_out_guess)
      end

      #---------------------------------
      # Fluid calculations
      #---------------------------------
      γ_p_in = cp_p_in / (cp_p_in - Rp)
      ρ_p_in = pp_in / (Rp * Tp_in)
      Vp_in = Mp_in * sqrt(γ_p_in * Rp * Tp_in)

      if occursin("liquid", fluid_c)
            lc_in = liquid_properties(fluid_c, Tc_inf)
            ρ_c_in = lc_in.ρ
            Vc_in  = Mc_in * lc_in.a
      else
            γ_c_in = cp_c_in / (cp_c_in - Rc)
            ρ_c_in = pc_in / (Rc * Tc_in)
            Vc_in = Mc_in * sqrt(γ_c_in * Rc * Tc_in)
      end

      #---------------------------------
      # Mean fluid properties
      #---------------------------------
      # Evaluate gas properties at a mean temperature and use these to find heat transfer coefficients 
      # See Kays. Compact Heat Exchangers (1984), p. 106

      Tp_m = (Tp_out + Tp_in) / 2 #Mean temperature of process stream
      Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of coolant stream

      _, Pr_p_m, _, _, μ_p_m, k_p_m = gasPr(fluid_p, Tp_m)
      ρ_p_m = pp_in / (Rp * Tp_m) #Assume pressure is constant TODO: consider adding Rayleigh flow model

      if occursin("liquid", fluid_c) #if coolant is liquid
            lc_m = liquid_properties(fluid_c, Tc_inf)
            ρ_c_m  = lc_m.ρ
            μ_c_m  = lc_m.μ
            k_c_m  = lc_m.k
            Pr_c_m = lc_m.Pr
            Vc_m = ρ_c_in * Vc_in / ρ_c_m
      else
            _, Pr_c_m, _, cp_c_m, μ_c_m, k_c_m = gasPr(fluid_c, Tc_m)

            ρ_c_m = pc_in / (Rc * Tc_m)
            Vc_m = ρ_c_in * Vc_in / ρ_c_m #conservation of mass
      end
      ν_c_m = μ_c_m / ρ_c_m

      #---------------------------------
      # Geometry calculations
      #---------------------------------
      A_cs = mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream

      if is_concentric #If channel is concentric, e.g., engine core
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter
            b = pi * D_i #Inner circumference
      else #If channel is rectangular
            b = A_cs / l #Side length
      end

      A_cc = mdot_c / (ρ_c_in * Vc_in) #total coolant cross-sectional area

      K = pi * b * n_stages / (4 * xt_D * A_cc) #Constant for tubesize!
      tubesize!(K, HXgeom) #size for design pressure difference

      #Extract outputs from tubesize!
      tD_o = HXgeom.tD_o
      t = HXgeom.t
      kw = HXgeom.material.k

      tD_i = tD_o - 2 * t #tube inner diameter

      N_t = 4 * A_cc / (pi * tD_i^2 * n_stages) #number of different coolant tubes per row

      xtm_D = A_cs / (N_t * tD_o * l) #Mean tangential pitch to diameter ratio

      A_D = N_t * l * tD_o #Total cross-sectional area of tubes in row, for drag calculation
      A_min = A_cs - A_D
      G = mdot_p / A_min #mass flow rate per unit area at minimum area

      #---------------------------------
      # Iterative loop
      #---------------------------------
      N_iter = 15 #Expect fast convergence

      n_passes = 4.0 #Initialize number of passes
      n_passes_prev = n_passes
      Ah = 0.0
      Cf = 0.0
      Tw_c = Tw_p = (Tp_out + Tp_in + Tc_out + Tc_in) / 4 #guess wall temperature
      for i = 1:N_iter
            N_L = n_passes * n_stages #total number of rows

            #Calculate heat transfer coefficient for coolant
            Re_D_c = Vc_m * tD_i / ν_c_m #Reynolds number based on pipe diameter
            jc, Cf = jcalc_pipe(Re_D_c) #Colburn j-factor and skin-friction coefficient
            Nu_cm =  Re_D_c * jc * Pr_c_m ^ (1/3) #Nusselt number in mean flow
            if ~occursin("liquid", fluid_c) #if fluid is a gas
                  Nu_c = Nu_cm * (Tw_c/Tc_m)^(-0.5) #Eq.(4.1) in Kays and London (1998)
            else
                  Nu_c = Nu_cm
            end
            h_c = Nu_c * k_c_m / tD_i

            # Calculate heat transfer coefficient for process side
            Re_D_p = G * tD_o / μ_p_m #Reynolds number based on minimum free flow and tube outer diameter
            Nu_pm = Nu_calc_staggered_cyl(Re_D_p, Pr_p_m, N_L, xtm_D, xl_D) #Nusselt number based on mean flow
            Nu_p = Nu_pm * (Tw_p/Tp_m)^0.0 #Eq.(4.1) in Kays and London (1998)
            h_p = Nu_p * k_p_m / tD_o

            #Overall thermal resistance times area (m^2 K / W)
            RA = 1 / ( h_c * (tD_i/tD_o) ) + Rfp + tD_o / tD_i * Rfc + t / kw + 1 / h_p 

            # Size heat exchanger
            Ah = NTU * C_min * RA   #Find required process-side cooling area from NTU
            n_passes = Ah / (N_t * n_stages * pi * tD_o * l)

            #Wall temperature (on process side)
            Tw_p = Tc_m + ((Tp_m - Tc_m)/RA) * (1 / ( h_c * (tD_i/tD_o) ) + tD_o / tD_i * Rfc + Rfp + t / kw)

            #Wall temperature (on coolant side)
            Tw_c = Tc_m + ((Tp_m - Tc_m)/RA) * (1 / ( h_c * (tD_i/tD_o) ))

            if (abs((n_passes_prev - n_passes)/n_passes) < tol)
                  break #Break for loop if convergence has been reached
            end 
            n_passes_prev = n_passes #otherwise store current value for comparison
      end

      #---------------------------------
      # Compute pressure drops
      #---------------------------------
      _, _, _, _, μ_p_w, _ = gasPr(fluid_p, Tw_p)
      μ_μw = μ_p_m / μ_p_w #Ratio of free flow viscosity to wall viscosity

      N_tubes_tot = N_t * n_passes * n_stages #total number of tubes across all rows
      N_L = n_passes * n_stages #total number of rows
      L = N_L * xl_D * tD_o #total axial length

      #Volumetric hydraulic diameter; Dv = 4 * (Net free volume) / (Friction surface)
      NFV = A_cs * L - N_tubes_tot * pi * tD_o^2 * l / 4 #Net free volume
      Dv = 4 * NFV / Ah

      Re_Dv = Dv * G / μ_p_m #Reynolds number based on hydraulic diameter and minimum free area

      if Re_Dv < 0 #If whole section is blocked
            Δp_p = Inf
      else
            Δp_p = Δp_calc_staggered_cyl(Re_Dv, G, L, ρ_p_m, Dv, tD_o, xtm_D, xl_D, μ_μw) #Calculate using the method of Gunter and Shaw (1945)
      end

      Pl_p = Δp_p * mdot_p / ρ_p_m #Power loss due to pressure drop in process stream

      # Compute coolant pressure drop
      τw = ρ_c_m * Vc_m^2 / 2 * Cf
      A_s_c = pi * tD_i * l * n_passes #Surface area on one coolant stream
      A_cs_c = pi * tD_i^2 / 4 #cross-sectional area of one coolant stream
      Δp_c = τw * A_s_c / A_cs_c
      Δp_c = Δp_c * (Tw_c/Tc_m)^(-0.1) #Eq.(4.2) in Kays and London (1998)

      Pl_c = Δp_c * mdot_c / ρ_c_m #Power loss due to pressure drop in coolant stream

      #---------------------------------
      # Output structs
      #---------------------------------
      #These lines rewrite the input structs

      #Gas parameters
      HXgas.Tp_out = Tp_out
      HXgas.Tc_out = Tc_out
      HXgas.Δh_p = Δh_p
      HXgas.Δh_c = Δh_c
      HXgas.Δp_p = Δp_p
      HXgas.Δp_c = Δp_c
      HXgas.Pl_p = Pl_p
      HXgas.Pl_c = Pl_c
      HXgas.ε = ε

      #Geometry parameters
      HXgeom.N_t = N_t
      HXgeom.n_passes = n_passes
      HXgeom.A_cs = A_cs

      if has_recirculation
            HXgas.mdot_r = mdot_r
      end

end #hxsize!

"""
    hxoper!(HXgas, HXgeom)

Evaluates crossflow heat exchanger performance for off-design operation. Uses the ε-NTU 
method to calculate effectiveness from prescribed geometry.       

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    
    **Outputs:**
    No direct outputs. Input structures are modified with outlet gas properties.
"""
function hxoper!(HXgas::HX_gas, HXgeom::HX_tubular)
      #---------------------------------
      # Extract inputs
      #---------------------------------
      #Gas parameters
      fluid_p = HXgas.fluid_p
      fluid_c = HXgas.fluid_c
      alpha_p = HXgas.alpha_p
      igas_c = HXgas.igas_c 
      mdot_p = HXgas.mdot_p
      mdot_c_inf = HXgas.mdot_c 
      Tp_in = HXgas.Tp_in
      Tc_inf = HXgas.Tc_in 
      pp_in = HXgas.pp_in
      pc_in = HXgas.pc_in

      #Flags 
      is_concentric = HXgeom.is_concentric
      has_recirculation = HXgeom.has_recirculation

      #HX geometry
      t = HXgeom.t
      N_t = HXgeom.N_t
      n_stages = HXgeom.n_stages
      n_passes = HXgeom.n_passes
      xl_D = HXgeom.xl_D
      tD_o = HXgeom.tD_o 
      kw = HXgeom.material.k
      Rfp = HXgeom.Rfp
      Rfc = HXgeom.Rfc 
      l = HXgeom.l
      A_cs = HXgeom.A_cs

      tol = 1e-10

      if has_recirculation
            recircT = HXgas.recircT
            h_lat = HXgas.h_lat
      end

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------

      _, _, hp_in, _, cp_p_in, Rp = gassum(alpha_p, length(alpha_p), Tp_in)

      if occursin("liquid", fluid_c)
            cp_c_inf = liquid_properties(fluid_c, Tc_inf).cp
            hc_inf = cp_c_inf * Tc_inf #TODO: replace with enthalpy calculation for non-constant cp
      else
            _, _, hc_inf, _, cp_c_inf, Rc = gasfun(igas_c, Tc_inf)
      end

      #---------------------------------
      # Fluid calculations
      #---------------------------------

      if has_recirculation
            _, _, hc_in, _, cp_c_in, _ = gasfun(igas_c, recircT)
            Tc_in = recircT
      
      else #No recirculation
            mdot_c = mdot_c_inf
            Tc_in = Tc_inf
            hc_in = hc_inf
            cp_c_in = cp_c_inf
      end

      ρ_p_in = pp_in / (Rp * Tp_in)

      if occursin("liquid", fluid_c)
            ρ_c_in = liquid_properties(fluid_c, Tc_inf).ρ
      else
            ρ_c_in = pc_in / (Rc * Tc_in)
      end
      
      # Calculate process-side velocity from geometry
      Vp_in = mdot_p / (ρ_p_in * A_cs) #process freestream velocity
      N_hyd_ways = N_t * n_stages #number of different hydraulic pathways

      #---------------------------------
      # Geometry calculations
      #---------------------------------
      tD_i = tD_o - 2 * t #tube inner diameter
      A_cs_tube = pi * tD_i^2 / 4 #Tube cross-sectional area

      A_D = N_t * l * tD_o #Total cross-sectional area of tubes in row, for drag calculation
      A_min = A_cs - A_D
      G = ρ_p_in * Vp_in * A_cs / A_min #mass flow rate per unit area at minimum area      

      N_tubes_tot = HXgeom.N_tubes_tot #total number of tubes across all rows
      N_L = n_passes * n_stages #total number of rows
      L =  HXgeom.L #total axial length

      Ah = N_tubes_tot * pi * tD_o * l #Total external surface area of cooling tubes

      xtm_D = A_cs / (N_t * tD_o * l) #Mean tangential pitch to diameter ratio
      
      #---------------------------------
      # Thermal calculations
      #---------------------------------

      # Compute C_p and guess C_min
      C_p = mdot_p * cp_p_in #process-side heat capacity rate
      C_c = mdot_c_inf * cp_c_in #Coolant heat capacity rate, guess

      C_min = min(C_c, C_p)

      Qmax = C_min * (Tp_in - Tc_in) #Maximum possible heat transfer rate, guess

      #---------------------------------
      # Iterative loop
      #---------------------------------
      ε = 0.95 #guess for effectiveness

      Qg = ε * Qmax #Actual heat transfer rate, guess
      Qprev = Qg

      # Guess outlet and wall temperatures
      Tp_out = Tp_in - Qg / C_p 
      Tc_out = Tc_in + Qg / C_c
      Tw_c = Tw_p = (Tp_out + Tp_in + Tc_out + Tc_in) / 4

      N_iter = 20 #Rapid convergence expected

      ρ_p_m = 0.0 #Initiliaze because of Julia scope
      μ_p_m = 0.0
      Δh_p = 0.0
      Δh_c = 0.0
      Vc_m = 0.0
      ρ_c_m = 0.0
      Cf = 0.0
      Tc_m = 0.0

      mdot_r = 0.0 #Initialize
      for i = 1 : N_iter
            
            if has_recirculation
                  #Calculate assuming that C_c = C_min
                  #Note that this ignores the enthalpy increase in the recirculation stream when it is recompressed
                  mdot_r = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * cp_c_in * (Tp_in -  recircT) )
                  
                  C_check = (mdot_c_inf + mdot_r) * cp_c_in #Coolant heat capacity rate, to check validity of above assumption

                  C_max = max(C_check, C_p)

                  if C_check == C_max #If the calculation above is incorrect because C_c = C_max
                        A = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * C_p * (Tp_in -  recircT) )
                        if A > 1
                              error("Insufficient heat capacity in process stream")
                        end
                        mdot_r = A * mdot_c_inf / (1 - A)
                  end
                  mdot_c = mdot_c_inf + mdot_r
            end

            #Heat capacity rates
            C_c = mdot_c * cp_c_in #Coolant heat capacity rate

            C_min = min(C_c, C_p)
            C_max = max(C_c, C_p)
            C_r = C_min / C_max

            Qmax = C_min * (Tp_in - Tc_in) #Maximum heat transfer rate
            
            #Properties at mean temperature
            Tp_m = (Tp_out + Tp_in) / 2 #Mean temperature of process stream
            Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of coolant stream

            _, Pr_p_m, _, _, μ_p_m, k_p_m = gasPr(fluid_p, Tp_m)
            ρ_p_m = pp_in / (Rp * Tp_m) #Assume pressure is constant TODO: consider adding Rayleigh flow model
            
            if occursin("liquid", fluid_c)
                  lc_m = liquid_properties(fluid_c, Tc_inf)
                  ρ_c_m  = lc_m.ρ
                  μ_c_m  = lc_m.μ
                  k_c_m  = lc_m.k
                  Pr_c_m = lc_m.Pr
      
            else
                  _, Pr_c_m, _, cp_c_m, μ_c_m, k_c_m = gasPr(fluid_c, Tc_m)
      
                  ρ_c_m = pc_in / (Rc * Tc_m)
            end
            ν_c_m = μ_c_m / ρ_c_m

            Vc_m =  mdot_c / (N_hyd_ways * ρ_c_m * A_cs_tube) #coolant velocity at mean temperature

            # Calculate thermal resistance
            # Calculate heat transfer coefficient for process side
            Re_D_p = G * tD_o / μ_p_m #Reynolds number based on minimum free flow and tube outer diameter
            Nu_pm = Nu_calc_staggered_cyl(Re_D_p, Pr_p_m, N_L, xtm_D, xl_D) #Nusselt number
            Nu_p = Nu_pm * (Tw_p/Tp_m)^0.0 #Eq.(4.1) in Kays and London (1998)
            h_p = Nu_p * k_p_m / tD_o

            #Calculate heat transfer coefficient for coolant
            Re_D_c = Vc_m * tD_i / ν_c_m #Reynolds number based on pipe diameter
            jc, Cf = jcalc_pipe(Re_D_c) #Colburn j-factor
            Nu_cm =  Re_D_c * jc * Pr_c_m ^ (1/3) #Nusselt number in mean flow
            if ~occursin("liquid", fluid_c) #if fluid is a gas
                  Nu_c = Nu_cm * (Tw_c/Tc_m)^(-0.5) #Eq.(4.1) in Kays and London (1998)
            else
                  Nu_c = Nu_cm
            end
            h_c = Nu_c * k_c_m / tD_i

            #Overall thermal resistance times area (m^2 K / W)
            RA = 1 / ( h_c * (tD_i/tD_o) ) + Rfp + tD_o / tD_i * Rfc + t / kw + 1 / h_p 

            # Find NTU and use it to calculate effectiveness
            NTU = Ah / (C_min * RA)

            if C_c == C_min
                  ε = 1 / C_r *(1 - exp(-C_r * (1 - exp(-NTU)) ) ) #effectiveness for cross flow with C_max mixed and C_min unmixed
            else
                  ε = 1 - exp(-1 / C_r * (1 - exp(-C_r * NTU) ) ) #effectiveness for cross flow with C_max unmixed and C_min mixed
            end

            # Calculate total heat transfer and exit temperatures
            Q = ε * Qmax #Actual heat transfer rate

            #Wall temperature (on process side)
            Tw_p = Tc_m + ((Tp_m - Tc_m)/RA) * (1 / ( h_c * (tD_i/tD_o) ) + tD_o / tD_i * Rfc + Rfp + t / kw)

            #Wall temperature (on coolant side)
            Tw_c = Tc_m + ((Tp_m - Tc_m)/RA) * (1 / ( h_c * (tD_i/tD_o) ))

            #Outlet properties
            Tp_out_guess = Tp_in - Q / C_p 
            Tc_out_guess = Tc_in + Q / C_c

            hp_out = hp_in - Q / mdot_p
            hc_out = hc_in + Q / mdot_c

            Δh_p = hp_out - hp_in
            Δh_c = hc_out - hc_in

            Tp_out = gas_tset(alpha_p, length(alpha_p), hp_out, Tp_out_guess)

            if occursin("liquid", fluid_c) #if coolant is liquid
                  Tc_out = Tc_out_guess
            else
                  Tc_out = gas_tset_single(igas_c, hc_out, Tc_out_guess)
            end

            if (abs((Q-Qprev)/Q) < tol)
                  break #break loop if convergence has been reached
            end
            Qprev = Q #else update previous heat
      end
      Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of coolant stream
      Tp_m = (Tp_out + Tp_in) / 2 #Mean temperature of process stream
      #---------------------------------
      # Compute pressure drop
      #---------------------------------
      _, _, _, _, μ_p_w, _ = gasPr(fluid_p, Tw_p)
      _, _, _, _, μ_p_m, _ = gasPr(fluid_p, Tp_m)
      μ_μw = μ_p_m / μ_p_w #Ratio of free flow viscosity to wall viscosity

      #Volumetric hydraulic diameter; Dv = 4 * (Net free volume) / (Friction surface)
      NFV = A_cs * L - N_tubes_tot * pi * tD_o^2 * l / 4 #Net free volume
      Dv = 4 * NFV / Ah

      Re_Dv = Dv * G / μ_p_m #Reynolds number based on hydraulic diameter and minimum free area

      if Re_Dv < 0 #If whole section is blocked
            Δp_p = Inf
      else
            Δp_p = Δp_calc_staggered_cyl(Re_Dv, G, L, ρ_p_m, Dv, tD_o, xtm_D, xl_D, μ_μw) #Calculate using the method of Gunter and Shaw (1945)
      end

      # Compute coolant pressure drop
      τw = ρ_c_m * Vc_m^2 / 2 * Cf
      A_s_c = pi * tD_i * l * n_passes #Surface area on one coolant stream
      A_cs_c = pi * tD_i^2 / 4 #cross-sectional area of one coolant stream
      Δp_c = τw * A_s_c / A_cs_c
      Δp_c = Δp_c * (Tw_c/Tc_m)^(-0.1) #Eq.(4.2) in Kays and London (1998)

      #---------------------------------
      # Output structs
      #---------------------------------
      #Gas parameters
      HXgas.Tp_out = Tp_out
      HXgas.Tc_out = Tc_out
      HXgas.Tw = Tw_p #Store only process side wall temperature
      HXgas.Δh_p = Δh_p
      HXgas.Δh_c = Δh_c
      HXgas.Δp_p = Δp_p
      HXgas.Δp_c = Δp_c
      HXgas.ε = ε

      if has_recirculation
            HXgas.mdot_r = mdot_r
      end

end #hxoper!

"""
      radiator_coolant_mass(HXgas::HX_gas, Q::Float64)

Calculates the coolant mass flow rate that a heat exchanger needs to reject a desired amount of heat. It assumes
that the process side has maximum heat capacity rate.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::HX_gas`: structure with the gas properties
    - `Q::Float64`: required heat transfer rate (W)
    
    **Outputs:**
    - `mdot_c::Float64`: coolant mass flow rate (kg/s)
"""
function radiator_coolant_mass(HXgas::HX_gas, Q::Float64)

      #Fluid parameters
      Tp_in = HXgas.Tp_in
      Tc_in = HXgas.Tc_in
      ε = HXgas.ε

      #specific heats
      if occursin("liquid", HXgas.fluid_c)
            cp_c = liquid_properties(HXgas.fluid_c, HXgas.Tc_in).cp
      else
            _, _, _, _, cp_c, _ = gasfun(HXgas.igas_c, HXgas.Tc_in)
      end
      _, _, _, _, cp_p, _ = gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)
      
      #Calculate minimum heat capacity rate from desired effectiveness and temperature difference
      C_min = abs(Q / (ε * (Tp_in - Tc_in)))

      #Design for C_min being C_c
      mdot_c = C_min / cp_c

      C_r = mdot_c * cp_c / (HXgas.mdot_p * cp_p)
      if C_r > 1 #If the coolant has the maximum heat capacity rate
            error("In radiator sizing, process stream has minimum heat capacity rate")
      end
            
      return mdot_c
end

"""
      RadiatorOffDesignCalc!(HXgas::HX_gas, HXgeom::HX_tubular, Q::Float64)

Evaluates the off-design performance of a heat exchanger for a given process-side mass flow rate and required heat transfer rate.
The function assumes that the minimum heat capacity rate is in the coolant stream, and calculates the coolant mass flow rate required to 
meet the heat transfer requirement.    

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::HX_gas`: structure with the gas properties
    - `HXgeom::HX_tubular`: structure with the HX geometric properties
    - `Q::Float64`: required heat transfer rate (W)
    
    **Outputs:**
    No direct outputs. Input structures are modified with outlet gas properties.
"""
function RadiatorOffDesignCalc!(HXgas::HX_gas, HXgeom::HX_tubular, Q::Float64)

      #TODO: consider case with recirculation
      if occursin("liquid", HXgas.fluid_c)
            cp_c = liquid_properties(HXgas.fluid_c, HXgas.Tc_in).cp
      else
            _, _, _, _, cp_c, _ = gasfun(HXgas.igas_c, HXgas.Tc_in)
      end
      _, _, _, _, cp_p, _ = gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)

      HXod_res(C_r) = HXheating_residual!(HXgas, HXgeom, Q, C_r) #Should be 0 at correct C_r

      Crg = 0.5 #guess for heat capacity rate
      C_r = find_zero(HXod_res, Crg) #Find root with Roots.jl

      if abs(HXod_res(C_r)) > 1e-4 #If the non-linear solver did not find a suitable solution
            error("Failed to find coolant mass flow rate in radiator off-design model")
      end
      HXgas.mdot_c = C_r * HXgas.mdot_p * cp_p / cp_c

      hxoper!(HXgas, HXgeom)

end
  
  """
      HXheating_residual!(HXgas::HX_gas, HXgeom::HX_tubular, Q::Float64, C_r::Float64)

Calculates the difference between the heat transfer rate of a heat exchanger and its required heat capacity rate for a given
ratio of heat capacity rates.  

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    - `Q::Float64`: required heat transfer rate (W)
    - `C_r::Float64`: ratio of heat capacity rate in process-side to heat capacity rate in coolant-side
    
    **Outputs:**
    - `res::Float64`: relative difference between desired heat rate and actual heat rate
"""
function HXheating_residual!(HXgas::HX_gas, HXgeom::HX_tubular, Q::Float64, C_r::Float64)

      if occursin("liquid", HXgas.fluid_c)
            cp_c = liquid_properties(HXgas.fluid_c, HXgas.Tc_in).cp
      else
            _, _, _, _, cp_c, _ = gasfun(HXgas.igas_c, HXgas.Tc_in)
      end
      _, _, _, _, cp_p, _ = gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)

      HXgas.mdot_c = C_r * (HXgas.mdot_p * cp_p) / cp_c

      try
            hxoper!(HXgas, HXgeom)
      catch
            return 1.0 #return something other than 0 if it fails
      end

      Q_HX = HXgas.Δh_p * HXgas.mdot_p

      res = 1 - Q_HX / Q
      return res
end

"""
    hxoptim!(HXgas, HXgeom, initial_x)

Optimizes heat exchanger design parameters for a given set of inputs. Uses the NLopt.jl package. The optimization
variables are `Mc_in`, `n_stages`, `xt_D` and `l`. The length of `initial_x` is the flag to determine how many parameters 
to optimize: if it is 4, all parameters are optimized; if it is 3, the tube length `l` is assumed to be an input and is not 
optimized.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    - `initial_x::Vector{Float64}`: vector with the initial guess for the optimization
    
    **Outputs:**
    No direct outputs. Input structures are modified with HX design geometry.
"""
function hxoptim!(HXgas::HX_gas, HXgeom::HX_tubular, initial_x::Vector{Float64})
      #Parameters to optimize: x[1]: 100 * Mc_in; x[2]: n_stages; x[3]: xt_D; x[4]: l
      #Set function to minimize
      obj(x, grad) =  hxobjf(x, HXgas, HXgeom) #Minimize objective function

      #Calculate minimum tube length
      alpha_p = HXgas.alpha_p #process gas composition
      mdot_p = HXgas.mdot_p #mass flow rate of process gas (kg/s)
      Tp_in = HXgas.Tp_in
      Mp_in  = HXgas.Mp_in
      pp_in  = HXgas.pp_in

      #Flags 
      is_concentric = HXgeom.is_concentric

      _, _, _, _, cp_p_in, Rp = gassum(alpha_p, length(alpha_p), Tp_in)
      γ_p_in = cp_p_in / (cp_p_in - Rp)
      ρ_p_in = pp_in / (Rp * Tp_in)
      Vp_in = Mp_in * sqrt(γ_p_in * Rp * Tp_in)

      A_cs = mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream
      
      if is_concentric #Flow is concentric
            D_i = HXgeom.D_i
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

            lmin = (D_o - D_i) / 2 #minimum tube length
            lmax = Inf
            
      else #square cross-section
            AR_max = 10 #Maximum aspect ratio
            #lmax = sqrt(AR_max * A_cs)
            lmax = Inf
            lmin = sqrt(A_cs / AR_max)

      end

      #Set bounds
      lower = [1e-9, 1.0, 1.0, lmin]
      upper = [30.0, 20.0, 6.0, lmax]
      initial_dx = [0.1, -0.1, -0.1, 0.1]

      #Use NLopt.jl to minimize function 
      opt = NLopt.Opt(:LN_COBYLA, length(initial_x)) #COBYLA can handle constraints natively

      # Set Optimization parameters
      opt.lower_bounds = lower
      opt.upper_bounds = upper
      opt.ftol_rel = 1e-9
      opt.initial_step = initial_dx
      opt.maxeval = 1000  # Set the maximum number of function evaluations

      opt.min_objective = obj
      
      tol = 1e-4
      inequality_constraint!(opt, (x, grad) -> MinPassesCstr(HXgeom), tol)
      inequality_constraint!(opt, (x, grad) -> MaxPassesCstr(HXgeom), tol)
      inequality_constraint!(opt, (x, grad) -> MinTubesCstr(HXgeom), tol)
      inequality_constraint!(opt, (x, grad) -> MaxTubesCstr(HXgeom), tol)
      inequality_constraint!(opt, (x, grad) -> MaxProcessΔPCstr(HXgas), tol)
      inequality_constraint!(opt, (x, grad) -> MaxCoolantΔPCstr(HXgas), tol)
      inequality_constraint!(opt, (x, grad) -> MaxLengthCstr(HXgeom), tol)

      (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
      #println(opt.numevals)
      #xopt_round = round.(xopt) #elements 2 and 3 must be integers
      
      #Modify structs with output
      HXgas.Mc_in = xopt[1] / 100 #x[1] has 100*Mc_in

      HXgeom.n_stages = xopt[2]
      HXgeom.xt_D = xopt[3]
      HXgeom.l = xopt[4]

      #Return optimum parameters by modifying input structs

end #hxoptim!

 
"""
      hxobjf(x, HXgas, HXgeom)

Objective function for HX optimization in hxoptim!(). It returns the sum of the power dissipated due to pressure
drops in the process and coolant streams.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x::Vector{Float64}`: state vector with [`100 * Mc_in`, `l`, `n_stages`, `xt_D`]
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
 
    **Outputs:**
    - `Iobj::Float64`: objective function (W)
"""
function hxobjf(x::Vector{Float64}, HXgas::HX_gas, HXgeom::HX_tubular)

      #Apply states
      HXgas.Mc_in = x[1] / 100

      HXgeom.n_stages = x[2]
      HXgeom.xt_D = x[3]
      HXgeom.l = x[4]

      #Size HX
      Iobj = Inf #Start with very high value of objective function
      try 
            hxsize!(HXgas, HXgeom)

            #Extract outputs
            Pl_p = HXgas.Pl_p
            Pl_c = HXgas.Pl_c

            Iobj = (Pl_p + Pl_c) #Initialize objective function

      catch #Do nothing if it errors
      end

      return Iobj
end #hxobjf

#Dictionaries with indices for HXs in engine
PreCDict = Dict{String, Int64}(
      "iTp_in" => ieTt19,
      "ipp_in" => iept19,
      "ipc_in" => iept3,
      "iDh"    => iePreCDeltah,
      "iDp"    => iePreCDeltap
)

InterCDict = Dict{String, Int64}(
      "iTp_in" => ieTt25,
      "ipp_in" => iept25,
      "ipc_in" => iept3,
      "iDh"    => ieInterCDeltah,
      "iDp"    => ieInterCDeltap
)

RegenDict = Dict{String, Int64}(
      "iTp_in" => ieTt49,
      "ipp_in" => iept49,
      "ipc_in" => iept3,
      "iDh"    => ieRegenDeltah,
      "iDp"    => ieRegenDeltap
)

TurbCDict = Dict{String, Int64}(
      "iTp_in" => ieTt3,
      "ipp_in" => iept3,
      "ipc_in" => iept3,
      "iDh"    => ieTurbCDeltah,
      "iDp"    => ieTurbCDeltap
)

rad_dict = Dict(
      "iTp_in" => ieTt21,
      "ipp_in" => iept21,
      "iTc_in" => ieRadiatorCoolantT,
      "ipc_in" => ieRadiatorCoolantP,
      "iDh"    => ieRadiatorDeltah,
      "iDp"    => ieRadiatorDeltap,
      "imp_in" => iemfan,
      "iQheat" => ieRadiatorHeat
      )

HXsDict = Dict{String, Dict}(
      "PreC" => PreCDict,
      "InterC" => InterCDict,
      "Regen" => RegenDict,
      "TurbC" =>  TurbCDict,
      "Radiator" => rad_dict
)

"""
      hxdesign!(pare, igas, ipdes, HXs_prev)

Heat exchanger design and operation function. It calls hxoptim!() to optimize the heat exchanger at the design point and 
then evaluates performance for all missions and points with hxoper!().      

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pare::Array{Float64 , 3}`: array with engine parameters
    - `igas::Float64`: coolant gas index
    - `ipdes::Float64`: index for design mission segment
    - `HXs_prev::Vector{Any}`: vector with heat exchanger data from the previous wsize iteration; elements are `HeatExchanger` structures
    - `rlx::Float64`: relaxation factor for pare update
    **Outputs:**
    - `HeatExchangers::Vector{Any}`: vector with heat exchanger data; elements are `HeatExchanger` structures
    - Also modifies `pare` with the fuel temperature and the HX enthalpy and pressure changes
"""
function hxdesign!(ac, ipdes, imission; rlx = 1.0)
      #Unpack aircraft
      pare = view(ac.pare,:, :, imission)
      pare_sl = view(pare,:, ipdes)
      igas = ac.options.ifuel
      HXs = ac.engine.heat_exchangers
      
      #Initialize Heat Exchanger vector
      Mc_opts = []

      for (i,HX) in enumerate(HXs) #For every desired type of heat exchanger (skipped if empty)
            type = HX.type
            #---------------------------------
            # Design exchangers
            #---------------------------------
            HXgeom, HXgas = PrepareHXobjects(HXs, i, ipdes, imission, igas, pare_sl, type, "sizing")
            HXgeom.Δpdes = max(maximum(ac.pare[iept3,:,:]), maximum(ac.pare[ieRadiatorCoolantP,:,:])) #size wall thickness for maximum HPC or coolant pressure
            HXgeom.maxL = HX.maximum_length #HXgeom.maxL is redundant with HX.maximum_length. TODO: make more elegant

            # Guess starting point for optimization
            #First calculate minimum tube length
            lmin, linit = calculate_min_tube_length(HXgeom, HXgas) #Minimum tube lenght and initial guess

            #Now set starting point
            if ~isassigned(HX.HXgas_mission, 1, 1) #If there is no previous heat exchanger design point
                  #Calculate initial length

                  initial_x = [3, 4, 4, linit] #Initial guess
            else 
                  #x[1]: 100 * Mc_in; x[2]: n_stages; x[3]: xt_D; x[4]: l;
                  initial_x = [100 * HX.HXgas_mission[ipdes, 1].Mc_in, 
                  HX.HXgeom.n_stages, HX.HXgeom.xt_D, max(HX.HXgeom.l, lmin)] #guess is previous iteration design point
            end

            hxoptim!(HXgas, HXgeom, initial_x) #Optimize heat exchanger geometry
            hxsize!(HXgas, HXgeom) #Evaluate all geometry properties at design point

            push!(Mc_opts, HXgas.Mc_in)
            
            HX.HXgeom = HXgeom #Store design point geometry
            HX.HXgas_mission[ipdes,imission] = HXgas #Store design point fluid state
      end
      #---------------------------------
      # Analyze off-design performance
      #---------------------------------
     
      HXOffDesign!(HXs, pare, igas, imission, rlx = rlx)

      for i in 1:length(HXs)
            HXs[i].HXgas_mission[ipdes,imission].Mc_in = Mc_opts[i] #Store optimum Mc_in
      end

      return HXs
end #hxdesign!

"""
      PrepareHXobjects(HeatExchangers, idx, ip, imission, igas, pare_sl, type, mode = "off_design")

This function prepares the gas and geometry heat exchanger objects to be used in design or 
off design analysis.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HeatExchangers::Vector{HeatExchanger}`: vector with heat exchanger data
    - `idx::Int64`: index for the heat exchanger number
    - `ip::Int64`: mission point index
    - `igas::Int64`: gas index
    - `pare_sl::Vector{Float64}`: sliced engine array with engine parameters
    - `type::String`: heat exchanger type
    - `mode::String`: design or off design indicator

    **Outputs:**
    - `HXgeom::HX_tubular`: structure with the HX geometric properties
    - `HXgas::HX_gas`: structure with the gas properties
"""  
function PrepareHXobjects(HeatExchangers, idx, ip, imission, igas, pare_sl, type, mode = "off_design")
      #Initiliaze design geometry and gas property as empty structs
      HXgas = HX_gas()
      HXgeom = HX_tubular()

      #Extract some inputs
      D_i = pare_sl[ieDi]
      Tc_ft = pare_sl[ieTft]
      has_recirculation = HeatExchangers[1].has_recirculation
      recircT = HeatExchangers[1].recirculation_temperature
      h_lat = pare_sl[iehvap]

      if igas == 11 #TODO: add more options
            coolant_name = "ch4"
      elseif igas == 40
            coolant_name = "h2"
      else
            coolant_name = ""
      end

      # Heat exchanger materials and wall properties
      HXgeom.xl_D = 1
      HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W

      mcore = pare_sl[iemcore]
      mofft = pare_sl[iemofft]
      fc = pare_sl[iefc] #Extract cooling gas factor
      ff = pare_sl[ieff]

      if mcore > 0
            fo = mofft / mcore
      else #Avoid divided by 0
            fo = 0
      end

      HXgas.fluid_p = "air"
      HXgas.alpha_p = alpha #Use alpha by default, except for Regen

      #Find indices in pare corresponding to HX variables
      iTp_in = HXsDict[type]["iTp_in"]
      ipp_in = HXsDict[type]["ipp_in"]
      ipc_in = HXsDict[type]["ipc_in"]
 
      #Store inputs
      if mode == "sizing"
            HXgas.ε = HeatExchangers[idx].design_effectiveness
            HXgas.Mp_in = HeatExchangers[idx].design_Mach
      end

      HXgas.Tp_in = pare_sl[iTp_in]
      HXgas.pp_in = pare_sl[ipp_in]
      HXgas.pc_in = pare_sl[ipc_in]
      HXgas.fluid_c = coolant_name
      HXgas.igas_c = igas
      HXgas.mdot_c = mcore * ff #Fuel fraction times core mass flow rate

      if type == "PreC" #Compressor Precooler
            HXgeom.is_concentric = true #Concentric
            HXgeom.D_i = D_i 
            HXgeom.Rfp = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 
            HXgeom.material = StructuralAlloy("Al-2219-T87")
            HXgeom.has_shaft = true #HX contains a shaft

            HXgas.mdot_p = mcore   #Core mass flow 

      elseif type == "InterC" #Compressor Intercooler
            HXgeom.is_concentric = true #Concentric
            HXgeom.D_i = D_i
            HXgeom.Rfp = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 
            HXgeom.material = StructuralAlloy("Al-2219-T87")
            HXgeom.has_shaft = true

            HXgas.mdot_p = mcore * (1 - fo) #Core mass flow minus offtake

      elseif type == "Regen" #Regenerative cooling
            HXgeom.is_concentric = true 
            HXgeom.D_i = D_i
            HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W 
            HXgeom.material = StructuralAlloy("SS-304") #use stainless steel for regenerative cooler as temp is above melting for Al

            HXgas.mdot_p = mcore * (1 - fo) #Core mass flow minus offtake
            
            HXgas.alpha_p = lambdap_calc(pare_sl, alpha, igas) #Calculate postcombustion and mixing composition

      elseif type =="TurbC" #Cooling of turbine cooling flow
            HXgeom.is_concentric = false 
            HXgeom.Rfp = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 
            HXgeom.material = StructuralAlloy("Al-2219-T87")

            HXgas.mdot_p = mcore * fc #Only cooling mass flow rate

      elseif type == "Radiator" #Radiator to reject residual heat
            HXgeom.is_concentric = true 
            HXgeom.D_i = D_i
            HXgeom.Rfp = 1.761e-4 #Compressed air fouling resistance, m^2*K/W 
            HXgeom.Rfc = 9E-05 #Distilled water fouling resistance, m^2*K/W
            HXgeom.material = StructuralAlloy("Al-2219-T87")

            iTc_in = HXsDict[type]["iTc_in"]
            imp_in =  HXsDict[type]["imp_in"]

            #Calculate coolant mass rate to release required heat
            HXgas.Tc_in = pare_sl[iTc_in]
            HXgas.mdot_p = pare_sl[imp_in]

            findLiquidCoolant!(HXgas) #Find coolant based on temperature
            if mode == "sizing"
                  iQheat = HXsDict[type]["iQheat"]
                  Q = pare_sl[iQheat]
                  HXgas.mdot_c = radiator_coolant_mass(HXgas, Q) 
            end
      end

      if type != "Radiator" #If HX is in the engine TODO make more robust to future HEX options
            if idx == 1 #At first heat exchanger
                  HXgas.Tc_in = Tc_ft #Coolant temperature is the tank temperature
                  if has_recirculation #There can only be recirculation in the first heat exchanger
                        HXgeom.has_recirculation = true 
                        HXgas.recircT = recircT
                        HXgas.h_lat = h_lat
                  else
                        HXgeom.has_recirculation = false 
                  end
                        
            else # For subsequent exchangers
                  HXprev = HeatExchangers[idx - 1] #Get previous heat exchanger struct
                  all_gas_prev = HXprev.HXgas_mission
                  HXgas.Tc_in = all_gas_prev[ip, imission].Tc_out #The inlet temperature is the outlet of previous HX at design point
                  
            end
      end

      return HXgeom, HXgas
end

"""
     HXOffDesign!(HeatExchangers, pare, igas, imission; rlx = 1.0)
This function runs the heat exchangers through an aircraft mission and calculates performance at every
mission point.      

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HeatExchangers::Vector{HeatExchanger}`: vector with heat exchanger data
    - `pare::Array{Float64 , 3}`: array with engine parameters
    - `igas::Int64`: gas index
    - `imission::Int64`: mission index
    - `rlx::Float64`: relaxation factor for pare update
    **Outputs:**
    Modifies `pare` with the fuel temperature and the HX enthalpy and pressure changes and
    `HeatExchangers` with the gas properties at every mission point.
"""
function HXOffDesign!(HeatExchangers, pare, igas, imission; rlx = 1.0)
      if length(HeatExchangers) == 0 #Skip if no HXs
            return
      end
      has_recirculation = HeatExchangers[1].has_recirculation #Check if there is recirculation in the first HX
      #Operate off-design for engine-integrated HEXs
      for (i,HX) in enumerate(HeatExchangers)
            HXgas_mis = Vector{Any}(undef, size(pare)[2]) #Vector to store gas properties across missions and segments

            type = HX.type

            Dh_i = HXsDict[type]["iDh"]
            Dp_i = HXsDict[type]["iDp"]
      
            for ip = 1:size(pare)[2] #For every point

                  pare_sl = pare[:, ip] #Slice pare with the parameters for the current point

                  _, HXgasp = PrepareHXobjects(HeatExchangers, i, ip, imission, igas, pare_sl, type)

                  if HXgasp.mdot_p == 0 #If the mass flow rate in this mission is 0, nothing happens
                        HXgasp.Tp_out = HXgasp.Tp_in
                        HXgasp.Tc_out = HXgasp.Tc_in
                        HXgasp.Δh_p = 0
                        HXgasp.Δp_p = 0
                        HXgasp.ε = 0
                  else #Otherwise, call HX off-design routine
                        if type == "Radiator" #If HEX is a radiator
                              iQheat = HXsDict[type]["iQheat"]
                              Q = pare_sl[iQheat]
                              
                              if Q > 0 #Radiator with non-zero heat
                                    RadiatorOffDesignCalc!(HXgasp, HX.HXgeom, Q)

                              else #Radiator with zero heat
                                    HXgasp.Tp_out = HXgasp.Tp_in
                                    HXgasp.Tc_out = HXgasp.Tc_in
                                    HXgasp.Δh_p = 0
                                    HXgasp.Δp_p = 0
                                    HXgasp.ε = 0
                              end
                        else #Not a radiator
                              hxoper!(HXgasp, HX.HXgeom)
                        end
                  end

                  if i == 1 && has_recirculation #If there is recirculation in the HX
                        #Store power drawn by recirculation to use it in the engine
                        P_recirc = find_recirculation_power(HXgasp)
                        pare[ieHXrecircP, ip] =  P_recirc
                        HXgasp.P_recirc = P_recirc
                  end

                  HXgas_mis[ip] = HXgasp

                  #Store output in pare
                  pare[Dh_i, ip] =  (1 - rlx) * pare[Dh_i, ip] + rlx * HXgasp.Δh_p
                  pare[Dp_i, ip] = (1 - rlx) * pare[Dp_i, ip] + rlx * HXgasp.Δp_p
                  
            end
            HeatExchangers[i].HXgas_mission[:,imission] = HXgas_mis
      end

      #---------------------------------
      # Update fuel temperature and heat of vaporization
      #---------------------------------

      if HeatExchangers[end].type != "Radiator" #If HX is in the engine TODO make more robust to future HEX options
            for ip = 1:size(pare)[2] #For every mission point
                  if length(HeatExchangers) > 0
                        lastHX = HeatExchangers[end]
                        HXgas = lastHX.HXgas_mission[ip, imission]
                        Tf = HXgas.Tc_out

                        pare[ieTfuel, ip] = (1 - rlx) * pare[ieTfuel, ip] + rlx * Tf
                  end
            end

            if (has_recirculation)  #Currently, non-zero heat of vaporization is only accounted for if there is recirculation
                  pare[iehvapcombustor, :, :] .= 0.0 #Fuel is vaporized in HX
            end

            findMinWallTemperature!(HeatExchangers) #Store minimum wall temperature at each mission point to check for freezing
      end
end

"""
      findLiquidCoolant!(HXgas)

This function calculates the coolant liquid based on its temperature.      

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::HX_gas`: heat exchanger gas object

    **Outputs:**
    Modifies `HXgas` with the coolant name
"""
function findLiquidCoolant!(HXgas)
      if HXgas.Tc_in > 373.15 #TODO add more coolant options
            coolant_name = "liquid ethylene glycol"
      else
            coolant_name = "liquid water"
      end
      HXgas.fluid_c = coolant_name
end

"""
      resetHXs(pare)

This function sets the fuel temperature to the temperature in the tank and sets the enthalpy and pressure changes
across the HXs to zero.      

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pare::Array{Float64 , 3}`: array with engine parameters

    **Outputs:**
    Modifies `pare` with the fuel temperature and the HX enthalpy and pressure changes
"""
function resetHXs(pare)
      #Reset fuel temperature
      pare[ieTfuel, :] = pare[ieTft, :] #Fuel tank temperature

      #Reset enthalpy differences and pressure differences in engine
      pare[iePreCDeltah, :] .= 0.0
      pare[iePreCDeltap, :] .= 0.0
      pare[ieInterCDeltah, :] .= 0.0
      pare[ieInterCDeltap, :] .= 0.0
      pare[ieTurbCDeltah, :] .= 0.0
      pare[ieTurbCDeltap, :] .= 0.0
      pare[ieRegenDeltah, :] .= 0.0
      pare[ieRegenDeltap, :] .= 0.0
      pare[ieRadiatorDeltah, :] .= 0.0
      pare[ieRadiatorDeltap, :] .= 0.0
      pare[ieHXrecircP, :] .= 0.0

      #Reset heat of vaporization in combustor
      pare[iehvapcombustor, :, :] = pare[iehvap, :, :]

end

"""
      check_HX_overwriting(HXs)

This function checks if a heat exchanger design effectiveness has been overwritten because it was too high. 
The heat transfer area tends to infinity as the effectiveness tends to the maximum possible one,
so a limit is set to 99% of the maximum possible one.     

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXs::Vector{HeatExchanger}`: vector with heat exchanger data

    **Outputs:**
    Produces a warning if the effectiveness has been overwritten.
"""
function check_HX_overwriting(HXs)
      flag = false
      for HX in HXs 
            if abs(HX.design_effectiveness - HX.HXgas_mission[ipcruise1, 1].ε) > 1e-5
                  flag = true
                  break
            end
      end
      if flag
            @warn "Heat-exchanger design effectiveness limited to 99% of maximum possible one"
      end
end

"""
      calculate_min_tube_length(HXgeom::HX_tubular, HXgas::HX_gas)

This function calculates the minimum tube length in a tubular HX pass. It also returns an initial guess 
for its optimization.  

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgeom::HX_tubular`: structure with the HX geometric properties
    - `HXgas::HX_gas`: structure with the gas properties
    
    **Outputs:**
    - `lmin::Float64`: minimum tube length (m)
    - `linit::Float64`: initial guess for tube length (m)
"""
function calculate_min_tube_length(HXgeom::HX_tubular, HXgas::HX_gas)
      _, _, _, _, cp_p_in, Rp = gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)
      γ_p_in = cp_p_in / (cp_p_in - Rp)
      ρ_p_in = HXgas.pp_in / (Rp * HXgas.Tp_in)
      Vp_in = HXgas.Mp_in * sqrt(γ_p_in * Rp * HXgas.Tp_in)

      A_cs = HXgas.mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream
      
      if HXgeom.is_concentric #Flow is concentric
            D_i = HXgeom.D_i
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

            lmin = (D_o - D_i) / 2 #minimum tube length
            linit = 1.1 * lmin
      else #square cross-section
            AR_min = 0.1 #Minimum aspect ratio
            lmin = sqrt(AR_min * A_cs)
            linit = sqrt(A_cs)
      end
      return lmin, linit
end

"""
    jcalc_pipe(Re_D)

Calculates the Colburn j-factor and skin-friction coefficient for flow inside a circular pipe, 
assuming flow is fully developed and turbulent. Uses the 1913 Blasius correlation.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Re_D::Float64`: Reynolds number based on pipe diameter

    **Outputs:**
    - `j::Float64`: Colburn j-factor for heat calculations
    - `Cf::Float64`: skin-friction coefficient

"""
function jcalc_pipe(Re_D::Float64)
      #turbulent flow
      Cf = 0.0791 * Re_D^(-0.25) #Blasius solution for smooth pipes 
      
      j = Cf / 2 #Colburn j-factor by Reynolds analogy
      return j, Cf
end #jcalc_pipe

"""
    Nu_calc_staggered_cyl(Re_D, Pr, N_L, xt_D, xl_D)

Calculates the Nusselt number for cross flow on a set of staggered circular cylinders. Based on
the model in A. Žkauskas. Heat Transfer from Tubes in Crossflow. Advances in Heat Transfer v.18 (1987).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Re_D::Float64`: Reynolds number based on cylinder diameter with velocity at minimum free flow area
    - `Pr::Float64`: Prandtl number
    - `N_L::Float64`: number of cylinder rows
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
    - `xt_D::Float64`: circumferential pitch between tubes over tube outer diameter
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
 
    **Outputs:**
    - `Nu::Float64`: Nusselt number based on cylinder diameter 
"""
function Nu_calc_staggered_cyl(Re_D::Float64, Pr::Float64, N_L::Float64, xt_D::Float64, xl_D::Float64)

      #First calculate C_2
      if (Re_D > 1000)
            C2 = 1 - exp(-N_L^(1 / sqrt(3)))
      else
            C2 = 1 - exp(-sqrt(3 * N_L^(1 / sqrt(2))))
      end

      if (Re_D < 40)
            C1 = 1.04
            m = 0.4
            n = 0.36
      elseif (Re_D >= 40) & (Re_D < 1000)
            C1 = 0.71
            m = 0.5
            n = 0.36
      elseif (Re_D >= 1000) & (Re_D <= 2e5) & (xt_D / xl_D < 2)
            C1 = 0.35 * (xt_D / xl_D) ^ 0.2
            m = 0.6
            n = 0.36
      elseif  (Re_D >= 1000) & (Re_D <= 2e5) & (xt_D / xl_D >= 2)
            C1 = 0.4
            m = 0.6
            n = 0.36
      else
            C1 = 0.031 * (xt_D / xl_D) ^ 0.2
            m = 0.8
            n = 0.4
      end
      Nu = C1 * C2 * Re_D^m * Pr^n

      return Nu
end #Nu_calc_staggered_cyl

"""
    Δp_calc_staggered_cyl(Re, G, L, ρ, Dv, tD_o, xt_D, xl_D, μ_μw)

Calculates the pressure drop across a set of staggered cylinders in cross flow. Uses the method of Gunter and Shaw. 
A General Correlation of Friction Factors for Various Types of Surfaces in Crossflow. Journal of Fluids Engineering, 1945.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Re::Float64`: Reynolds number based on hydraulic diameter and minimum free flow area: Re = Dv G /μ
    - `G::Float64`: mass flow rate divided by minimum free flow area. `G = mdot / (A_min)`, `A_min` is the minimum free-flow area (kg/s/m^2)
    - `L::Float64`: axial channel length (m)
    - `ρ::Float64`: density (kg/m^3)
    - `Dv::Float64`: volumetric hydraulic diameter. Dv = 4 * (Net free volume) / (Friction surface)
    - `tD_o::Float64`: cylinder outer diameter (m)
    - `xt_D::Float64`: circumferential pitch between tubes over tube outer diameter
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
    - `μ_μw::Float64`: ratio of free flow viscosity to wall viscosity
 
    **Outputs:**
    - `Δp::Float64`: pressure drop across staggered cylinders (Pa)
"""
function Δp_calc_staggered_cyl(Re::Float64, G::Float64, L::Float64, ρ::Float64, Dv::Float64, tD_o::Float64, xt_D::Float64, xl_D::Float64, μ_μw::Float64)

      #Compute friction factor, f_2 = f/2
      if (Re <= 200)
            f_2 = 90 / Re 
      else
            f_2 = 0.96 * Re ^ (-0.145)
      end

      #Calculate pressure drop
      Δp = G^2 * L / (Dv * ρ) * f_2 * (Dv / (xt_D * tD_o) )^0.4 * (xl_D / xt_D)^0.6 * μ_μw^(-0.14)
      
      return Δp
end #Δp_calc_staggered_cyl

"""
      tubesize!(K, HXgeom)

Calculates the tube diameter and thickness from flow and hoop stress balance.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `K::Float64`: constant in equation for `tD_o`; `K = pi * b * n_stages / (4 * xt_D * A_cc)`
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric and material properties
    
    **Outputs:**
    Modifies `HXgeom`.
"""
function tubesize!(K, HXgeom)
      Δp = HXgeom.Δpdes #design pressure difference
      safety_factor = 2
      σy = HXgeom.material.YTS

      tmin = 3e-4 #m, from Brewer 1991. Corresponds to 30 BWG.
      C = safety_factor * Δp / (2 * σy) #t = C * tD_o, from hoop stress balance
      
      thoop = C * K / (K - 2 * K * C)^2

      t = max(thoop, tmin) #Set a minimum thickness
      
      tD_o = (4 * K * t + sqrt(8 * K * t + 1) + 1) / (2 * K) #Compute tube outer diameter

      #Store outputs
      HXgeom.t = t
      HXgeom.tD_o = tD_o

end #tubesize!

"""
      hxweight(gee, HXgeom, fouter)

Calculates the weight of a heat exchanger with involute tubes.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gee::Float64`: gravitational acceleration (m/s^2)
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric and material properties
    - `fouter::Float64`: ratio of HX external mass to tube mass
 
    **Outputs:**
    - `W_hx::Float64`: weight of heat exchanger (N)
"""
function hxweight(gee, HXgeom, fouter)
      #Extract inputs
      tD_o = HXgeom.tD_o 
      ρ = HXgeom.material.ρ
      l = HXgeom.l
      t = HXgeom.t
      N_tubes_tot = HXgeom.N_tubes_tot #total number of tubes across all rows
      tD_i = HXgeom.tD_i

      V_t = N_tubes_tot * pi * (tD_o^2 - tD_i^2) / 4 * l #total tube volume
      m_t = ρ * V_t #total tube mass

      W_hx = gee * m_t * (1 + fouter)

      if HXgeom.has_shaft #If the HEX is in the core with a shaft through it, add additional mass for the shaft
            shaft_material = StructuralAlloy("AISI-4340") #Assume shaft is made of 4340 steel
            W_shaft = gee * shaft_material.ρ * HXgeom.L * HXgeom.D_i^2 * pi / 4 #Weight of the extra shaft length because of the HEX
            W_hx = W_hx + W_shaft
      end
      return W_hx
end #hxweight

"""
      findMinWallTemperature!(HXs)

Calculates and stores the minimum tube wall temperature at each mission point.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXs_prev::Vector{HeatExchanger}`: vector with heat exchanger data

    **Outputs:**
    Modifies `pare` with the fuel temperature and the HX enthalpy and pressure changes
"""
function findMinWallTemperature!(HXs)
      minT = Inf #Start with infinite temperature
      for ip in 1:iptotal
            for HX in HXs
                  if HX.HXgas_mission[ip, 1].mdot_c > 0  #TODO extend to off-design missions
                        minT = min(minT,  HX.HXgas_mission[ip, 1].Tw)
                  end
            end
      end
      if minT == Inf
            minT = 0.0 #Replace Inf with 0 when there is no mass flow rate
      end
      HXs[1].min_wall_temperature = minT
end

"""
      find_recirculation_power(HXgas, η_i = 1.0)

Calculates and stores the power needed to drive recirculation in a HX.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::HX_gas`: structure with HX fluid data
    - `η_i::Float64`: isentropic efficiency of the recirculation compressor (default 1.0, no losses)

    **Outputs:**
   - `P::Float64`: power needed to drive recirculation (W)
"""
function find_recirculation_power(HXgas, η_i = 1.0)
      
      mdot_r = HXgas.mdot_r
      if mdot_r ≈ 0.0 #If there is no recirculation, return 0
            return 0.0
      end
      Tout = HXgas.Tc_out
      igas_c = HXgas.igas_c
      Δp_c = HXgas.Δp_c
      pc_in = HXgas.pc_in

      pc_out = pc_in - Δp_c #Outlet pressure of coolant in HX

      #Find specific heat capacity at outlet
      _, _, _, _, cp_out, R = gasfun(igas_c, Tout)
      γ = cp_out / (cp_out - R)

      #Isentropic compression power
      P = 1/η_i * mdot_r * cp_out * Tout * ((pc_in/pc_out)^(1 - 1 / γ) - 1)
      return P
end

"""
      lambdap_calc(pare_sl, alpha_in, ifuel)

Calculates the mass fractions of the gases in post-combustion air.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pare_sl::Vector{Float64 }`: array slice with engine parameters
    - `alpha_in::Vector{Float64}`: vector with process-side gas composition before combustion and mixing
    - `ifuel::Float64`: fuel gas index
    
    **Outputs:**
    - `lambda_p::Vector{Float64}`: vector with process-side gas composition after combustion and mixing
"""
function lambdap_calc(pare_sl, alpha_in, ifuel)
      #Extract inputs
      Tt3 = pare_sl[ieTt3]
      Ttf = pare_sl[ieTfuel]
      Tt4 = pare_sl[ieTt4]
      hvap = pare_sl[iehvapcombustor]

      etab = pare_sl[ieetab]
      beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
      alpha = deepcopy(alpha_in)
      alpha = push!(alpha, 0.0)

      n = 6
      nair = n - 1
      gamma = gasfuel(ifuel, n)

      # Create buffer
      buf = Zygote.Buffer(gamma, length(gamma))
      for i = 1:length(gamma)
            buf[i] = gamma[i]
      end

      # Zygote can handle this
      for i = 1:nair
            buf[i] = etab * buf[i]
      end
      buf[n] = 1.0 - etab

      gamma = copy(buf)
      #

      ffb, lambda = gas_burn(alpha, beta, gamma, n, ifuel, Tt3, Ttf, Tt4, hvap)

      lambdap = zeros(nair)

      mcore = pare_sl[iemcore]
      mofft = pare_sl[iemofft]
      fc = pare_sl[iefc] #Extract cooling gas factor

      if mcore > 0
            fo = mofft / mcore
      else #Avoid divided by 0
            fo = 0
      end

      ff = pare_sl[ieff]

      #----- IGV exit mixing
      frac4 = (1.0 - fo - fc + ff) / (1.0 - fo + ff)
      fracm = fc / (1.0 - fo + ff)

      #----- mixed constituent fraction vector from mass equation
      # for i = 1:nair
      #       lambdap[i] = frac4 * lambda[i] + fracm * alpha[i]
      # end

      buf = Zygote.Buffer(lambdap, length(lambdap))
      for i = 1:nair
            buf[i] = frac4 * lambda[i] + fracm * alpha[i]
      end

      lambdap = copy(buf)

      return lambdap
end #lambdap_calc

"""
Thermodynamic transport properties of a liquid coolant.

Returned by `liquid_properties`.
"""
struct LiquidCoolantProps
      ρ::Float64   # density [kg/m³]
      cp::Float64  # specific heat capacity [J/(kg·K)]
      μ::Float64   # dynamic viscosity [Pa·s]
      k::Float64   # thermal conductivity [W/(m·K)]
      Pr::Float64  # Prandtl number [-]
      a::Float64   # speed of sound [m/s]
end

function liquid_properties(fluid, T)
      #TODO: account for temperature dependence of properties
      if (fluid == "liquid water") #properties at 100 degrees Celsius and 3 bar pressure from NIST
            ρ = 958.44
            cp = 4215.2
            μ = 0.00028161
            k = 0.67733
            a = 1500 #m/s, speed of sound in water

      elseif (fluid == "liquid ethylene glycol") #data for 180-200 degrees Celsius
            ρ = 1060.1
            cp = 2746.5
            μ = 0.0024
            k =  0.2362
            a = 1500 #m/s, speed of sound in water (only used as scaling for coolant speed)

      end
      Pr = cp * μ / k

      return LiquidCoolantProps(ρ, cp, μ, k, Pr, a)
end #liquid_properties

#Constraint functions for HX optimization
#TODO consider making the minima and maxima inputs rather than hardcoded params
function MinPassesCstr(HXgeom)
      n_passes = HXgeom.n_passes
      min_passes = 1.0
      return min_passes/n_passes - 1.0
end

function MaxPassesCstr(HXgeom)
      n_passes = HXgeom.n_passes
      max_passes = 20.0
      return n_passes/max_passes - 1.0
end

function MinTubesCstr(HXgeom)
      N_t = HXgeom.N_t
      min_tubes = 1.0
      return min_tubes/N_t - 1.0
end

function MaxTubesCstr(HXgeom)
      N_t = HXgeom.N_t
      max_tubes = 200.0
      return N_t/max_tubes - 1.0
end

function MaxProcessΔPCstr(HXgas)
      pp_in = HXgas.pp_in
      p_thres = 0.5
      Δp_p = HXgas.Δp_p
      Δp_max = p_thres * pp_in
      return Δp_p/Δp_max - 1.0
end

function MaxCoolantΔPCstr(HXgas)
      pc_in = HXgas.pc_in
      p_thres = 0.5 
      Δp_c = HXgas.Δp_c
      Δp_max = p_thres * pc_in
      return Δp_c/Δp_max - 1.0
end

function MaxLengthCstr(HXgeom)
      L = HXgeom.L
      Lmax = HXgeom.maxL
      return L/Lmax - 1.0
end