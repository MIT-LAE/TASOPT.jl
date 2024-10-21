# hxfun.jl
# These functions can be used to size and model a heat exchanger with involute staggered tubes in a crossflow
# The design method is based on the effectiveness-NTU method, described in many sources such as 
# https://www.mathworks.com/help/hydro/ref/entuheattransfer.html
# Nicolas Gomez Vega, Oct 2023

"""
    HX_gas

Structure containing the gas properties of the process and coolant streams.

!!! details "💾 Data fields"
    **Inputs:**
    - `fluid_p::String`: process fluid name
    - `fluid_c::String`: coolant fluid name
    - `alpha_p::Vector{Float64}`: process gas composition
    - `igas_c::Float64`: coolant gas index, if coolant is a gas
    - `mdot_p::Float64`: mass flow rate of process gas (kg/s)
    - `mdot_c::Float64`: mass flow rate of coolant gas (kg/s)
    - `Tp_in::Float64`: process gas inlet temperature (K)
    - `Tc_in::Float64`: coolant gas inlet temperature (K)
    - `pp_in::Float64`: process gas inlet pressure (Pa)
    - `pc_in::Float64`: coolant gas inlet pressure (Pa)
    - `Mp_in::Float64`: process gas inlet Mach number
    - `Mc_in::Float64`: coolant gas inlet Mach number
    - `Tp_out::Float64`: process gas outlet temperature
    - `Tc_out::Float64`: coolant gas outlet temperature
    - `Δh_p::Float64`: enthalpy change across HX (J/kg)
    - `Δp_p::Float64`: pressure drop of process gas across heat exchanger (Pa)
    - `Δp_c::Float64`: pressure drop of coolant gas across tubes (Pa)
    - `Pl_p::Float64`: power loss due to pressure drop in process stream (W)
    - `Pl_c::Float64`: power loss due to pressure drop in coolant stream (W)
    - `ε::Float64`: desired heat exchanger effectiveness
    - `recircT::Float64`: temperature of recirculating flow at HX inlet (K)
    - `mdot_r::Float64`: recirculating flow mass flow rate (kg/s)
    - `h_lat::Float64`: latent heat capacity in freestream coolant liquid (J/kg)
"""
@kwdef mutable struct HX_gas
      fluid_p :: String = ""
      fluid_c :: String = "" 
      alpha_p :: Vector{Float64} = []
      igas_c :: Float64 = 0.0
      mdot_p :: Float64 = 0.0
      mdot_c :: Float64 = 0.0
      Tp_in :: Float64 = 0.0
      Tc_in :: Float64 = 0.0
      pp_in :: Float64 = 0.0
      pc_in :: Float64 = 0.0
      Mp_in  :: Float64 = 0.0
      Mc_in :: Float64 = 0.0
      Tp_out :: Float64 = 0.0
      Tc_out :: Float64 = 0.0
      Δh_p :: Float64 = 0.0
      Δh_c :: Float64 = 0.0
      Δp_p :: Float64 = 0.0
      Δp_c :: Float64 = 0.0
      Pl_p :: Float64 = 0.0
      Pl_c :: Float64 = 0.0
      ε :: Float64 = 0.0 
      recircT :: Float64 = 0.0 
      mdot_r :: Float64 = 0.0 
      h_lat :: Float64 = 0.0 
end

"""
    HX_tubular

Structure containing the heat exchanger geometric and material properties.

!!! details "💾 Data fields"
    **Inputs:**
    - `fconc::Bool`: flag for concentric geometry (true: concentric ; false: rectangular)
    - `frecirc::Bool`: flag for recirculation (true: recirculation ; false: no recirculation)
    - `fshaf::Bool`: flag for whether HX contains shaf(true: shaft ; false: no shaft)
    - `N_t::Float64`: number of tubes per row
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
    - `n_passes::Float64`: number of coolant passes
    - `A_cs::Float64`: process side freestream cross-sectional area (m^2)
    - `l::Float64`: length of tubes (m)
    - `t::Float64`: cooling tube wall thickness (m)
    - `tD_o::Float64`: tube outer diameter (m)
    - `xt_D::Float64`: circumferential pitch between tubes at the root over tube outer diameter 
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
    - `Rfp::Float64`: process-side fouling factor (m^2 K/W)
    - `Rfc::Float64`: coolant-side fouling factor (m^2 K/W)
    - `D_i::Float64`: inner diameter of core (m)
    - `Δpdes::Float64`: design pressure difference between tube and outside (Pa)
"""
@kwdef mutable struct HX_tubular
      fconc :: Bool = false
      frecirc :: Bool = false
      fshaft :: Bool = false
      N_t :: Float64 = 0.0
      n_stages :: Float64 = 0.0
      n_passes:: Float64 = 0.0
      A_cs:: Float64  = 0.0
      l :: Float64 = 0.0
      t :: Float64 = 0.0
      tD_o :: Float64 = 0.0
      xt_D :: Float64 = 0.0
      xl_D :: Float64 = 0.0
      Rfp :: Float64 = 0.0
      Rfc :: Float64 = 0.0
      D_i :: Float64 = 0.0
      material :: StructuralAlloy = StructuralAlloy("Al-2219-T87")
      Δpdes::Float64 = 0.0
end

# Overload Base.getproperty for convenience
function Base.getproperty(HXgeom::HX_tubular, sym::Symbol)
      if (sym === :D_o) && getfield(HXgeom, :fconc)
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
    HX_struct

Structure containing all the heat exchanger geometry and operational information.

!!! details "💾 Data fields"
    **Inputs:**
    - `type::String`: type of heat exchanger ("PreC": precooler; "InterC": intercooler; "Regen": regenerative; "TurbC": turbine cooling)
    - `HXgeom::HX_tubular`: structure containing the HX geometric information
    - `HXgas_mission::Array{Any}`: array containing the gas properties, of type `HX_gas` for each mission and segment
"""
mutable struct HX_struct
      type :: String 
      HXgeom :: HX_tubular  
      HXgas_mission :: Array{Any}
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
      fconc = HXgeom.fconc
      frecirc = HXgeom.frecirc

      #HX geometry
      t = HXgeom.t #Initial wall thickness, may be overwritten
      n_stages = HXgeom.n_stages
      xt_D = HXgeom.xt_D
      xl_D = HXgeom.xl_D
      Rfp = HXgeom.Rfp
      Rfc = HXgeom.Rfc 
      l = HXgeom.l

      tol = 1e-10 #convergence tolerance

      if fconc #If geometry is concentric
            D_i = HXgeom.D_i
      end

      if frecirc
            recircT = HXgas.recircT
            h_lat = HXgas.h_lat

      end

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------
      _, _, hp_in, _, cp_p_in, Rp = gassum(alpha_p, length(alpha_p), Tp_in)

      if occursin("liquid", fluid_c)
            _, cp_c_inf, _, _, _, _  = liquid_properties(fluid_c, Tc_inf)
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

      if frecirc #If there is recirculation
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
                  error("Effectiveness exceeds maximum possible one")
            end
            NTU = -log(1 + log(1 - C_r * ε) / C_r) # For cross-flow with C_max mixed and C_min unmixed
      else
            ε_max = 1 - exp(-1 / C_r)#At ε = ε_max, NTU tends to infinity
            if ε > ε_max
                  error("Effectiveness exceeds maximum possible one")
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
            ρ_c_in, _, _, _, _, ac_in = liquid_properties(fluid_c, Tc_inf)
            Vc_in = Mc_in * ac_in
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
            ρ_c_m, cp_c_m, μ_c_m, k_c_m, Pr_c_m, ac_m  = liquid_properties(fluid_c, Tc_inf)
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

      if fconc #If channel is concentric, e.g., engine core
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
      Tw = (Tp_out + Tp_in + Tc_out + Tc_in) / 4 #guess wall temperature
      for i = 1:N_iter
            N_L = n_passes * n_stages #total number of rows

            #Calculate heat transfer coefficient for coolant
            Re_D_c = Vc_m * tD_i / ν_c_m #Reynolds number based on pipe diameter
            jc, Cf = jcalc_pipe(Re_D_c) #Colburn j-factor and skin-friction coefficient
            Nu_cm =  Re_D_c * jc * Pr_c_m ^ (1/3) #Nussel number in mean flow
            if ~occursin("liquid", fluid_c) #if fluid is a gas
                  Nu_c = Nu_cm * (Tw/Tc_m)^(-0.5) #Eq.(4.1) in Kays and London (1998)
            end
            h_c = Nu_c * k_c_m / tD_i

            # Calculate heat transfer coefficient for process side
            Re_D_p = G * tD_o / μ_p_m #Reynolds number based on minimum free flow and tube outer diameter
            Nu_pm = Nu_calc_staggered_cyl(Re_D_p, Pr_p_m, N_L, xtm_D, xl_D) #Nusselt number based on mean flow
            Nu_p = Nu_pm * (Tw/Tc_m)^0.0 #Eq.(4.1) in Kays and London (1998)
            h_p = Nu_p * k_p_m / tD_o

            #Overall thermal resistance times area (m^2 K / W)
            RA = 1 / ( h_c * (tD_i/tD_o) ) + Rfp + tD_o / tD_i * Rfc + t / kw + 1 / h_p 

            # Size heat exchanger
            Ah = NTU * C_min * RA   #Find required process-side cooling area from NTU
            n_passes = Ah / (N_t * n_stages * pi * tD_o * l)

            #Wall temperature (neglect change across wall)
            Tw = Tc_m + ((Tp_m - Tc_m)/RA) * (1 / ( h_c * (tD_i/tD_o) ) + tD_o / tD_i * Rfc)

            if (abs(n_passes_prev - n_passes)/n_passes < tol)
                  break #Break for loop if convergence has been reached
            end 
            n_passes_prev = n_passes #otherwise store current value for comparison
      end

      #---------------------------------
      # Compute pressure drops
      #---------------------------------
      _, _, _, _, μ_p_w, _ = gasPr(fluid_p, Tw)
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
      Δp_c = Δp_c * (Tw/Tc_m)^(-0.1) #Eq.(4.2) in Kays and London (1998)

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

      if frecirc
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
      fconc = HXgeom.fconc
      frecirc = HXgeom.frecirc

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

      if frecirc
            recircT = HXgas.recircT
            h_lat = HXgas.h_lat
      end

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------

      _, _, hp_in, _, cp_p_in, Rp = gassum(alpha_p, length(alpha_p), Tp_in)

      if occursin("liquid", fluid_c)
            _, cp_c_inf, _, _, _, _  = liquid_properties(fluid_c, Tc_inf)
            hc_inf = cp_c_inf * Tc_inf #TODO: replace with enthalpy calculation for non-constant cp
      else
            _, _, hc_inf, _, cp_c_inf, Rc = gasfun(igas_c, Tc_inf)
      end

      #---------------------------------
      # Fluid calculations
      #---------------------------------

      if frecirc
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
            ρ_c_in, _, _, _, _, _ = liquid_properties(fluid_c, Tc_inf)
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
      Tw = (Tp_out + Tp_in + Tc_out + Tc_in) / 4

      N_iter = 20 #Rapid convergence expected

      ρ_p_m = 0.0 #Initiliaze because of Julia scope
      μ_p_m = 0.0
      Δh_p = 0.0
      Δh_c = 0.0

      mdot_r = 0.0 #Initialize
      for i = 1 : N_iter
            
            if frecirc
                  #Calculate assuming that C_c = C_min
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
                  ρ_c_m, cp_c_m, μ_c_m, k_c_m, Pr_c_m, _  = liquid_properties(fluid_c, Tc_inf)
      
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
            Nu_p = Nu_pm * (Tw/Tc_m)^0.0 #Eq.(4.1) in Kays and London (1998)
            h_p = Nu_p * k_p_m / tD_o

            #Calculate heat transfer coefficient for coolant
            Re_D_c = Vc_m * tD_i / ν_c_m #Reynolds number based on pipe diameter
            jc, Cf = jcalc_pipe(Re_D_c) #Colburn j-factor
            Nu_cm =  Re_D_c * jc * Pr_c_m ^ (1/3) #Nussel number in mean flow
            if ~occursin("liquid", fluid_c) #if fluid is a gas
                  Nu_c = Nu_cm * (Tw/Tc_m)^(-0.5) #Eq.(4.1) in Kays and London (1998)
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

            #Wall temperature (neglect change across wall)
            Tw = Tc_m + ((Tp_m - Tc_m)/RA) * (1 / ( h_c * (tD_i/tD_o) ) + tD_o / tD_i * Rfc)

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

            if (abs(Q-Qprev)/Q < tol)
                  break #break loop if convergence has been reached
            end
            Qprev = Q #else update previous heat
      end
      Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of coolant stream
      Tp_m = (Tp_out + Tp_in) / 2 #Mean temperature of process stream
      #---------------------------------
      # Compute pressure drop
      #---------------------------------
      _, _, _, _, μ_p_w, _ = gasPr(fluid_p, Tw)
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

      #---------------------------------
      # Output structs
      #---------------------------------
      #Gas parameters
      HXgas.Tp_out = Tp_out
      HXgas.Tc_out = Tc_out
      HXgas.Δh_p = Δh_p
      HXgas.Δh_c = Δh_c
      HXgas.Δp_p = Δp_p
      HXgas.ε = ε

      if frecirc
            HXgas.mdot_r = mdot_r
      end

end #hxoper!

"""
      radiator_design!(HXgas, HXgeom, Q)

    Evaluates the off-design performance of a heat exchanger for a given process-side mass flow rate and required heat transfer rate.
    The function assumes that the minimum heat capacity rate is in the coolant stream, and calculates the coolant mass flow rate required to 
    meet the heat transfer requirement. 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    - `Q::Float64`: required heat transfer rate (W)
    
    **Outputs:**
    No direct outputs. Input structures are modified with outlet gas and HX properties.

"""
function radiator_design!(HXgas::HX_gas, HXgeom::HX_tubular, Q::Float64)

      #Fluid parameters
      Tp_in = HXgas.Tp_in
      Tc_in = HXgas.Tc_in
      pp_in = HXgas.pp_in
      alpha_p = HXgas.alpha_p

      #specific heats
      if occursin("liquid", HXgas.fluid_c)
            _, cp_c, _, _, _, _ = liquid_properties(HXgas.fluid_c, HXgas.Tc_in)
      else
            _, _, _, _, cp_c, _ = gasfun(HXgas.igas_c, HXgas.Tc_in)
      end
      _, _, _, _, cp_p, Rp = gassum(alpha_p, length(alpha_p), Tp_in)

      #TODO: this calculation does not work when there is recirculation or if the process side has C_min

      #Calculate minimum heat capacity rate from desired effectiveness and temperature difference
      C_min = abs(Q / (ε * (Tp_in - Tc_in)))

      #Design for C_min being C_c
      mdot_c = C_min / cp_c

      #Find tube length, assuming square cross section
      ρ_p = pp_in / (Rp * Tp_in)
      γ = cp_p / (cp_p - Rp)
      A_cs = mdot_p / (ρ_p * Mp_in * sqrt(γ * Rp * Tp_in))
      l = sqrt(A_cs)

      HXgas.mdot_c = mdot_c

      HXgeom.l = l
      HXgeom.xl_D = 1

      initial_x = [0.1, 6, 4] #do not optimize tube length for this radiator

      hxoptim!(HXgas, HXgeom, initial_x)
      hxsize!(HXgas, HXgeom)

end

"""
      HXoffDesignCalc!(HXgas, HXgeom, Q)

Evaluates the off-design performance of a heat exchanger for a given process-side mass flow rate and required heat transfer rate.
The function assumes that the minimum heat capacity rate is in the coolant stream, and calculates the coolant mass flow rate required to 
meet the heat transfer requirement.    

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    - `Q::Float64`: required heat transfer rate (W)
    
    **Outputs:**
    No direct outputs. Input structures are modified with outlet gas properties.
"""
function HXoffDesignCalc!(HXgas::HX_gas, HXgeom::HX_tubular, Q::Float64)

      #TODO: consider case with recirculation
      if occursin("liquid", HXgas.fluid_c)
            _, cp_c, _, _, _, _ = liquid_properties(HXgas.fluid_c, HXgas.Tc_in)
      else
            _, _, _, _, cp_c, _ = gasfun(HXgas.igas_c, HXgas.Tc_in)
      end
      _, _, _, _, cp_p, _ = gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)

      HXod_res(C_r) = HXheating_residual!(HXgas, HXgeom, Q, C_r) #Should be 0 at correct C_r

      Crg = 1.5 #guess for heat capacity rate
      C_r = find_zero(HXod_res, Crg) #Find root with Roots.jl

      HXgas.mdot_c = 1 / C_r * HXgas.mdot_p * cp_p / cp_c

      hxoper!(HXgas, HXgeom)

end
  
  """
      HXheating_residual!(HXgas, HXgeom, Q, C_r)

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
            _, cp_c, _, _, _, _ = liquid_properties(HXgas.fluid_c, HXgas.Tc_in)
      else
            _, _, _, _, cp_c, _ = gasfun(HXgas.igas_c, HXgas.Tc_in)
      end
      _, _, _, _, cp_p, _ = gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)

      HXgas.mdot_c = 1 / C_r * (HXgas.mdot_p * cp_p) / cp_c

      hxoper!(HXgas, HXgeom)

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
      #Parameters to optimize: x[1]: 100 * Mc_in; x[2]: n_stages; x[3]: xt_D; x[4]: l (optional)
      #Set function to minimize
      obj(x, grad) =  hxobjf(x, HXgas, HXgeom) #Minimize objective function

      #Calculate minimum tube length
      alpha_p = HXgas.alpha_p #process gas composition
      mdot_p = HXgas.mdot_p #mass flow rate of process gas (kg/s)
      Tp_in = HXgas.Tp_in
      Mp_in  = HXgas.Mp_in
      pp_in  = HXgas.pp_in

      #Flags 
      fconc = HXgeom.fconc

      _, _, _, _, cp_p_in, Rp = gassum(alpha_p, length(alpha_p), Tp_in)
      γ_p_in = cp_p_in / (cp_p_in - Rp)
      ρ_p_in = pp_in / (Rp * Tp_in)
      Vp_in = Mp_in * sqrt(γ_p_in * Rp * Tp_in)

      A_cs = mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream
      
      if fconc #Flow is concentric
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
      if length(initial_x) == 4
            lower = [0.0, 1.0, 1.0, lmin]
            upper = [30.0, 20.0, 6.0, lmax]
            initial_dx = [0.1, -0.1, -0.1, 0.1]
      else #Only 3 optimization variables
            lower = [0.0, 1.0, 1.0]
            upper = [30.0, 20.0, 6.0]
            initial_dx = [0.1, -0.1, -0.1]
      end
      
      #Use NLopt.jl to minimize function 
      opt = Opt(:LN_NELDERMEAD, length(initial_x))
      opt.lower_bounds = lower
      opt.upper_bounds = upper
      opt.ftol_rel = 1e-9
      opt.initial_step = initial_dx
      opt.maxeval = 500  # Set the maximum number of function evaluations

      opt.min_objective = obj
      
      (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
      #println(opt.numevals)
      #xopt_round = round.(xopt) #elements 2 and 3 must be integers

      #Modify structs with output
      HXgas.Mc_in = xopt[1] / 100 #x[1] has 100*Mc_in

      HXgeom.n_stages = xopt[2]
      HXgeom.xt_D = xopt[3]

      if length(initial_x) == 4 #only add length if it is being optimized
            HXgeom.l = xopt[4]
      end

      #Return optimum parameters by modifying input structs

end #hxoptim!

"""
      hxobjf(x, HXgas, HXgeom)

Objective function for HX optimization in hxoptim!(). It returns the sum of the power dissipated due to pressure
drops in the process and coolant streams, with penalty factors to enforce constraints.

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

      if length(x) == 4 #only add length if it is being optimized
            HXgeom.l = x[4]
      end

      #Size HX
      Iobj = 1e9 #Start with very high value of objective function
      try 
            hxsize!(HXgas, HXgeom)

            #Extract outputs
            Pl_p = HXgas.Pl_p
            Pl_c = HXgas.Pl_c

            n_passes = HXgeom.n_passes
            N_t = HXgeom.N_t
            Δp_p = HXgas.Δp_p
            Δp_c = HXgas.Δp_c
            L = HXgeom.L #HX length

            #Inlet pressures (pressure drops should not exceed these)
            pp_in = HXgas.pp_in
            pc_in = HXgas.pc_in
            p_thres = 0.5 #start applying penalty function is pressure drops exceed this fraction of the inlet pressure

            #Constraints are applied to these variables
            vars = [n_passes, N_t, Δp_p, Δp_c, L]

            lower = [1.0, 1.0, 1.0, 1.0, 1e-3] #desired lower limits
            upper = [20.0, 200.0, p_thres * pp_in, p_thres * pc_in, 0.25]
            #TODO: the HEX length (L) constraint should be an input (maybe a global optimization variable?)

            Iobj = (Pl_p + Pl_c) #Initialize objective function

            # Apply penalty function so that outputs are within desired range
            for (i,var) in enumerate(vars)
                  vmin = lower[i]
                  vmax = upper[i]

                  pmin = ( vmin / min(vmin, var) )^2 #Penalty function
                  pmax = ( max(vmax, var) / vmax )^2
                  p = max(pmin, pmax)

                  Iobj = Iobj * p
            end
      catch #Do nothing if it errors; this returns the high Iobj
      end

      return Iobj
end #hxobjf

"""
      hxdesign!(pare, pari, ipdes, HXs_prev)

Heat exchanger design and operation function. It calls hxoptim!() to optimize the heat exchanger at the design point and 
then evaluates performance for all missions and points with hxoper!().      

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pare::Array{Float64 , 3}`: array with engine parameters
    - `pari::Vector{Int}`: vector with integer parameters
    - `ipdes::Float64`: index for design mission segment
    - `HXs_prev::Vector{Any}`: vector with heat exchanger data from the previous wsize iteration; elements are `HX_struct` structures
    - `rlx::Float64`: relaxation factor for pare update
    **Outputs:**
    - `HeatExchangers::Vector{Any}`: vector with heat exchanger data; elements are `HX_struct` structures
    - Also modifies `pare` with the fuel temperature and the HX enthalpy and pressure changes
"""
function hxdesign!(pare, pari, ipdes, HXs_prev; rlx = 1.0)
      
      #---------------------------------
      # Extract inputs
      #---------------------------------
      pare_sl = pare[:, ipdes] #Slice pare at design point

      D_i = pare_sl[ieDi]
      Tc_ft = pare_sl[ieTft]
      frecirc = Bool(pare_sl[iefrecirc])
      recircT = pare_sl[ierecircT]
      h_lat = pare_sl[iehvap]
      igas = pari[iifuel]
      PreCorder = pare_sl[iePreCorder]
      PreCepsilon = pare_sl[iePreCepsilon]
      PreCMp = pare_sl[iePreCMp]
      InterCorder = pare_sl[ieInterCorder]
      InterCepsilon = pare_sl[ieInterCepsilon]
      InterCMp = pare_sl[ieInterCMp]
      Regenorder = pare_sl[ieRegenorder]
      Regenepsilon = pare_sl[ieRegenepsilon]
      RegenMp = pare_sl[ieRegenMp]
      TurbCorder = pare_sl[ieTurbCorder]
      TurbCepsilon = pare_sl[ieTurbCepsilon]
      TurbCMp = pare_sl[ieTurbCMp]

      if igas == 11 #TODO: add more options
            coolant_name = "ch4"
      elseif igas == 40
            coolant_name = "h2"
      end

      # Sort heat exchangers
      all_types = ["PreC", "InterC", "Regen", "TurbC"]
      all_orders = [PreCorder, InterCorder, Regenorder, TurbCorder]
      all_Mp = [PreCMp, InterCMp, RegenMp, TurbCMp]
      all_eps = [PreCepsilon, InterCepsilon, Regenepsilon, TurbCepsilon]

      HXtypes = []
      Mp_in = []
      ε_des = []
      sort_i = sortperm(all_orders) #Sort according to order

      for ind in sort_i
            if (all_eps[ind] > 0) && (all_eps[ind] <= 1) #If effectiveness is between 0 and 1
                  push!(HXtypes, all_types[ind])
                  push!(Mp_in, all_Mp[ind])
                  push!(ε_des, all_eps[ind])
            end
      end

      alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127] #Air composition
      #Initialize Heat Exchanger vector
      HeatExchangers = []

      for (i,type) in enumerate(HXtypes) #For every desired type of heat exchanger (skipped if empty)
            #---------------------------------
            # Design exchangers
            #---------------------------------
            pare_sl = pare[:, ipdes] #Slice pare at design point
            #Initiliaze design geometry and gas property as empty structs
            HXgas = HX_gas()
            HXgeom = HX_tubular()

            HXgas.ε = ε_des[i]

            # Heat exchanger materials and wall properties
            HXgeom.xl_D = 1
            HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W
            HXgeom.Δpdes = maximum(pare[iept3,:]) #size wall thickness for maximum HPC pressure

            mcore = pare_sl[iemcore]
            mofft = pare_sl[iemofft]
            fc = pare_sl[iefc] #Extract cooling gas factor
            fo = mofft / mcore
            ff = pare_sl[ieff]

            HXgas.fluid_p = "air"
            HXgas.alpha_p = alpha #Use alpha by default, except for Regen

            if type == "PreC" #Compressor Precooler
                  HXgeom.fconc = true #Concentric
                  HXgeom.D_i = D_i 
                  HXgeom.Rfp = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 
                  HXgeom.material = StructuralAlloy("Al-2219-T87")
                  HXgeom.fshaft = true #HX contains a shaft

                  HXgas.mdot_p = mcore   #Core mass flow 
                  iTp_in = ieTt19
                  ipp_in = iept19
                  ipc_in = iept3

                  Dh_i = iePreCDeltah
                  Dp_i = iePreCDeltap

            elseif type == "InterC" #Compressor Intercooler
                  HXgeom.fconc = true #Concentric
                  HXgeom.D_i = D_i
                  HXgeom.Rfp = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 
                  HXgeom.material = StructuralAlloy("Al-2219-T87")
                  HXgeom.fshaft = true

                  HXgas.mdot_p = mcore * (1 - fo) #Core mass flow minus offtake
                  iTp_in = ieTt25
                  ipp_in = iept25
                  ipc_in = iept3

                  Dh_i = ieInterCDeltah
                  Dp_i = ieInterCDeltap

            elseif type == "Regen" #Regenerative cooling
                  HXgeom.fconc = true 
                  HXgeom.D_i = D_i
                  HXgeom.Rfp = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W 
                  HXgeom.material = StructuralAlloy("SS-304") #use stainless steel for regenerative cooler as temp is above melting for Al

                  HXgas.mdot_p = mcore * (1 - fo) #Core mass flow minus offtake
                  iTp_in = ieTt49
                  ipp_in = iept49
                  ipc_in = iept3

                  Dh_i = ieRegenDeltah
                  Dp_i = ieRegenDeltap

                  HXgas.alpha_p = lambdap_calc(pare, alpha, igas, ipdes) #Calculate postcombustion and mixing composition

            elseif type =="TurbC" #Cooling of turbine cooling flow
                  HXgeom.fconc = false 
                  HXgeom.Rfp = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 
                  HXgeom.material = StructuralAlloy("Al-2219-T87")

                  HXgas.mdot_p = mcore * fc #Only cooling mass flow rate
                  iTp_in = ieTt3
                  ipp_in = iept3
                  ipc_in = iept3

                  Dh_i = ieTurbCDeltah
                  Dp_i = ieTurbCDeltap
            end

            HXgas.Tp_in = pare_sl[iTp_in]
            HXgas.pp_in = pare_sl[ipp_in]
            HXgas.Mp_in = Mp_in[i]

            HXgas.fluid_c = coolant_name
            HXgas.igas_c = igas
            HXgas.mdot_c = mcore * ff #Fuel fraction times core mass flow rate
            
            HXgas.pc_in = pare_sl[ipc_in]

            if i == 1 #At first heat exchanger
                  HXgas.Tc_in = Tc_ft #Coolant temperature is the tank temperature
                  if frecirc #There can only be recirculation in the first heat exchanger
                        HXgeom.frecirc = true 
                        HXgas.recircT = recircT
                        HXgas.h_lat = h_lat
                  else
                        HXgeom.frecirc = false 
                  end
                        
            else # For subsequent exchangers
                  HXprev = HeatExchangers[i - 1] #Get previous heat exchanger struct
                  all_gas_prev = HXprev.HXgas_mission
                  HXgas.Tc_in = all_gas_prev[ipdes].Tc_out #The inlet temperature is the outlet of previous HX at design point
                  
            end

            # Guess starting point for optimization
            #First calculate minimum tube length
            _, _, _, _, cp_p_in, Rp = gassum(HXgas.alpha_p, length(HXgas.alpha_p), HXgas.Tp_in)
            γ_p_in = cp_p_in / (cp_p_in - Rp)
            ρ_p_in = HXgas.pp_in / (Rp * HXgas.Tp_in)
            Vp_in = HXgas.Mp_in * sqrt(γ_p_in * Rp * HXgas.Tp_in)

            A_cs = HXgas.mdot_p / (ρ_p_in * Vp_in) #Cross-sectional area of freestream
            
            if HXgeom.fconc #Flow is concentric
                  D_i = HXgeom.D_i
                  D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

                  lmin = (D_o - D_i) / 2 #minimum tube length
                  linit = 1.1 * lmin
            else #square cross-section
                  AR_min = 0.1 #Minimum aspect ratio
                  lmin = sqrt(AR_min * A_cs)
                  linit = sqrt(A_cs)
            end

            #Now set starting point
            if isempty(HXs_prev) #If there is no previous heat exchanger design point
                  #Calculate initial length

                  initial_x = [3, 4, 4, linit] #Initial guess
            else 
                  #x[1]: 100 * Mc_in; x[2]: n_stages; x[3]: xt_D; x[4]: l;
                  initial_x = [100 * HXs_prev[i].HXgas_mission[ipdes].Mc_in, 
                  HXs_prev[i].HXgeom.n_stages, HXs_prev[i].HXgeom.xt_D, max(HXs_prev[i].HXgeom.l, lmin)] #guess is previous iteration design point
            end

            hxoptim!(HXgas, HXgeom, initial_x) #Optimize heat exchanger geometry
            hxsize!(HXgas, HXgeom) #Evaluate all geometry properties at design point

            #---------------------------------
            # Analyze off-design performance
            #---------------------------------
            HXgas_mis = Vector{Any}(undef, size(pare)[2]) #Vector to store gas properties across missions and segments

            for ip = 1:size(pare)[2] #For every point

                  pare_sl = pare[:, ip] #Slice pare with the parameters for the current point
                  HXgasp = deepcopy(HXgas) #Initialize gas property struct for this segment

                  mcore = pare_sl[iemcore]
                  mofft = pare_sl[iemofft]
                  fc = pare_sl[iefc] #Extract cooling gas factor

                  if mcore > 0
                        fo = mofft / mcore
                  else #Avoid divided by 0
                        fo = 0
                  end

                  ff = pare_sl[ieff]

                  if (type == "PreC") 
                        HXgasp.mdot_p = mcore 

                  elseif (type == "InterC") 
                        HXgasp.mdot_p = mcore * (1 - fo)

                  elseif (type == "Regen")
                        HXgasp.mdot_p = mcore * (1 - fo)
                        HXgasp.alpha_p = lambdap_calc(pare, alpha, igas, ip) #Calculate postcombustion and mixing composition
                        
                  elseif type =="TurbC"
                        HXgasp.mdot_p = mcore * fc

                  end        

                  HXgasp.Tp_in = pare_sl[iTp_in] #The indices come from the design process above as the HX is the same
                  HXgasp.pp_in = pare_sl[ipp_in]

                  HXgasp.mdot_c = mcore * ff #Fuel fraction times core mass flow rate
                  HXgasp.pc_in = pare_sl[ipc_in]

                  if i == 1 #If this is the first heat exchanger
                        HXgas.Tc_in = Tc_ft #Temperature is the tank temperature
                        if frecirc #If there is recirculation, it can only happen at fist HX
                              HXgasp.recircT = recircT
                              HXgasp.h_lat = h_lat
                        end
                  else
                        HXprev = HeatExchangers[i - 1]
                        all_gas_prev = HXprev.HXgas_mission
                        HXgasp.Tc_in = all_gas_prev[ip].Tc_out
                  end

                  if HXgasp.mdot_p == 0 #If the mass flow rate in this mission is 0, nothing happens
                        HXgasp.Tp_out = HXgasp.Tp_in
                        HXgasp.Tc_out = HXgasp.Tc_in
                        HXgasp.Δh_p = 0
                        HXgasp.Δp_p = 0
                        HXgasp.ε = 0
                  else #Otherwise, call HX off-design routine
                        hxoper!(HXgasp, HXgeom)

                  end

                  HXgas_mis[ip] = HXgasp

                  #Store output in pare
                  pare[Dh_i, ip] =  (1 - rlx) * pare[Dh_i, ip] + rlx * HXgasp.Δh_p
                  pare[Dp_i, ip] = (1 - rlx) * pare[Dp_i, ip] + rlx * HXgasp.Δp_p
                  
            end
            push!(HeatExchangers, HX_struct(type, HXgeom, HXgas_mis)) #Store HX struct in overall array
      end

      #---------------------------------
      # Update fuel temperature and heat of vaporization
      #---------------------------------

      for ip = 1:size(pare)[2] #For every mission point
            if length(HeatExchangers) > 0
                  lastHX = HeatExchangers[end]
                  HXgas = lastHX.HXgas_mission[ip]
                  Tf = HXgas.Tc_out

                  pare[ieTfuel, ip] = (1 - rlx) * pare[ieTfuel, ip] + rlx * Tf
            end
      end

      if (frecirc) && (length(HeatExchangers) > 0) #Currently, non-zero heat of vaporization is only accounted for if there is recirculation
            pare[iehvapcombustor, :, :] .= 0.0 #Fuel is vaporized in HX
      end
     
      return HeatExchangers
end #hxdesign!

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

      #Reset heat of vaporization in combustor
      pare[iehvapcombustor, :, :] = pare[iehvap, :, :]

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

      if HXgeom.fshaft #If the HEX is in the core with a shaft through it, add additional mass for the shaft
            shaft_material = StructuralAlloy("AISI-4340") #Assume shaft is made of 4340 steel
            W_shaft = gee * shaft_material.ρ * HXgeom.L * HXgeom.D_i^2 * pi / 4 #Weight of the extra shaft length because of the HEX
            W_hx = W_hx + W_shaft
      end
      return W_hx
end #hxweight


"""
      lambdap_calc(pare, alpha_in, ifuel, ip)

Calculates the mass fractions of the gases in post-combustion air.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pare::Array{Float64 , 3}`: array with engine parameters
    - `alpha_in::Vector{Float64}`: vector with process-side gas composition before combustion and mixing
    - `ifuel::Float64`: fuel gas index
    - `ip::Float64`: index for mission point
    
    **Outputs:**
    - `lambda_p::Vector{Float64}`: vector with process-side gas composition after combustion and mixing
"""
function lambdap_calc(pare, alpha_in, ifuel, ip)
      #Extract inputs
      pare_sl = pare[:, ip]

      Tt3 = pare_sl[ieTt3]
      Ttf = pare_sl[ieTfuel]
      Tt4 = pare_sl[ieTt4]
      hvap = pare_sl[iehvapcombustor]

      etab = pare[ieetab]
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

      return ρ, cp, μ, k, Pr, a 
end #liquid_properties