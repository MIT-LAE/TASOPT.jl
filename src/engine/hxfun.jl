# hxfun.jl
# These functions can be used to size and model a heat exchanger with involute staggered tubes in a crossflow
# The design method is based on the effectiveness-NTU method, described in many sources such as 
# https://www.mathworks.com/help/hydro/ref/entuheattransfer.html
# Nicolas Gomez Vega, Oct 2023
using NLopt


"""
    HX_gas

Structure containing the gas properties of the hot and cold streams.

!!! details "💾 Data fields"
    **Inputs:**
    - `gas_h::String`: hot gas name
    - `gas_c::String`: coolant gas name
    - `alpha_h::Vector{Float64}`: hot gas composition
    - `igas_c::Float64`: cold gas index
    - `mdot_h::Float64`: mass flow rate of hot gas (kg/s)
    - `mdot_c::Float64`: mass flow rate of cold gas (kg/s)
    - `Th_in::Float64`: hot gas inlet temperature (K)
    - `Tc_in::Float64`: cold gas inlet temperature (K)
    - `ph_in::Float64`: hot gas inlet pressure (Pa)
    - `pc_in::Float64`: cold gas inlet pressure (Pa)
    - `Mh_in::Float64`: hot gas inlet Mach number
    - `Mc_in::Float64`: cold gas inlet Mach number
    - `Th_out::Float64`: hot gas outlet temperature
    - `Tc_out::Float64`: cold gas outlet temperature
    - `Δh_h::Float64`: enthalpy change across HX (J/kg)
    - `Δp_h::Float64`: pressure drop of hot gas across heat exchanger (Pa)
    - `Δp_c::Float64`: pressure drop of cold gas across tubes (Pa)
    - `Pl_h::Float64`: power loss due to pressure drop in hot stream (W)
    - `Pl_c::Float64`: power loss due to pressure drop in cold stream (W)
    - `ε::Float64`: desired heat exchanger effectiveness
    - `recircT::Float64`: temperature of recirculating flow at HX inlet (K)
    - `mdot_r::Float64`: recirculating flow mass flow rate (kg/s)
    - `h_lat::Float64`: latent heat capacity in freestream cold liquid (J/kg)
"""
mutable struct HX_gas
      gas_h :: String 
      gas_c :: String 
      alpha_h :: Vector{Float64} 
      igas_c :: Float64 
      mdot_h :: Float64
      mdot_c :: Float64
      Th_in :: Float64
      Tc_in :: Float64
      ph_in :: Float64
      pc_in :: Float64
      Mh_in  :: Float64
      Mc_in :: Float64
      Th_out :: Float64
      Tc_out :: Float64 
      Δh_h :: Float64
      Δp_h :: Float64 
      Δp_c :: Float64
      Pl_h :: Float64
      Pl_c :: Float64
      ε :: Float64 
      recircT :: Float64 
      mdot_r :: Float64 
      h_lat :: Float64 
end

"""
    HX_tubular

Structure containing the heat exchanger geometric and material properties.

!!! details "💾 Data fields"
    **Inputs:**
    - `fconc::Bool`: flag for concentric geometry (1: concentric ; 0: rectangular)
    - `frecirc::Bool`: flag for recirculation (1: recirculation ; 0: no recirculation)
    - `N_t::Float64`: number of tubes per row
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
    - `n_passes::Float64`: number of coolant passes
    - `A_cs::Float64`: hot side freestream cross-sectional area (m^2)
    - `l::Float64`: length of tubes (m)
    - `t::Float64`: cooling tube wall thickness (m)
    - `tD_o::Float64`: tube outer diameter (m)
    - `xt_D::Float64`: circumferential pitch between tubes at the root over tube outer diameter 
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
    - `Rfh::Float64`: hot-side fouling factor (m^2 K/W)
    - `Rfc::Float64`: cold-side fouling factor (m^2 K/W)
    - `kw::Float64`: thermal conductivity of wall material (W/m/K)
    - `ρw::Float64`: mean density of HE (kg/m^3)
    - `D_i::Float64`: inner diameter of core (m)
"""
mutable struct HX_tubular
      fconc :: Bool
      frecirc :: Bool 
      N_t :: Float64 
      n_stages :: Float64 
      n_passes:: Float64 
      A_cs:: Float64 
      l :: Float64
      t :: Float64
      tD_o :: Float64
      xt_D :: Float64
      xl_D :: Float64
      Rfh :: Float64
      Rfc :: Float64
      kw :: Float64
      ρw :: Float64
      D_i :: Float64
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
from a prescribed ε.  For representative fouling factors see Standards of the Tubular Exchanger Manufacturers Association
or `https://powderprocess.net/Tools_html/Data_Diagrams/Heat_Exchanger_Fouling_Factor.html`

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    
    **Outputs:**
    No direct outputs. Input structures are modified with outlet gas properties and HX design geometry.
"""
function hxsize!(HXgas, HXgeom)
      #---------------------------------
      # Extract inputs
      #---------------------------------
      #Gas parameters
      gas_h = HXgas.gas_h #Hot gas name
      gas_c = HXgas.gas_c #Cold gas name
      alpha_h = HXgas.alpha_h #Hot gas composition
      igas_c = HXgas.igas_c #Cold gas index
      mdot_h = HXgas.mdot_h #mass flow rate of hot gas (kg/s)
      mdot_c_inf = HXgas.mdot_c #mass flow rate of cold gas (kg/s)
      ε = HXgas.ε 
      Th_in = HXgas.Th_in
      Tc_inf = HXgas.Tc_in 
      Mh_in  = HXgas.Mh_in
      Mc_in = HXgas.Mc_in
      ph_in = HXgas.ph_in
      pc_in = HXgas.pc_in

      #Flags 
      fconc = HXgeom.fconc
      frecirc = HXgeom.frecirc

      #HX geometry
      t = HXgeom.t
      n_stages = HXgeom.n_stages
      xt_D = HXgeom.xt_D
      xl_D = HXgeom.xl_D
      kw = HXgeom.kw
      Rfh = HXgeom.Rfh
      Rfc = HXgeom.Rfc 
      l = HXgeom.l

      if fconc #If geometry is concentric
            D_i = HXgeom.D_i
      end

      if frecirc
            recircT = HXgas.recircT
            h_lat = HXgas.h_lat

            if isnan(h_lat) #If the latent heat was not specified, assume it is 0
                  h_lat = 0
            end
      end

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------
      _, _, hh_in, _, cp_h_in, Rh = gassum(alpha_h, length(alpha_h), Th_in)
      _, _, hc_inf, _, cp_c_inf, Rc = gasfun(igas_c, Tc_inf)

      #---------------------------------
      # Thermal calculations
      #---------------------------------

      C_h = mdot_h * cp_h_in #Hot-side heat capacity rate
      mdot_c = mdot_c_inf

      if frecirc #If there is recirculation
            _, _, hc_in, _, cp_c_in, _ = gasfun(igas_c, recircT)
            Tc_in = recircT
            N_iters_rec = 15 #Iterations in recirculation loop
            mdot_r = 0 #Initialize

            for j =1:N_iters_rec
                  #Calculate assuming that C_c = C_min
                  mdot_r = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * cp_c_in * (Th_in -  recircT) )
            
                  C_check = (mdot_c_inf + mdot_r) * cp_c_in #Coolant heat capacity rate, to check validity of above assumption

                  C_max = max(C_check, C_h)

                  if C_check == C_max #If the calculation above is incorrect because C_c = C_max
                        A = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * C_h * (Th_in -  recircT) )
                        
                        if A > 1
                              println("Insufficient heat capacity in hot stream")
                              return
                        end
                        mdot_r = A * mdot_c_inf / (1 - A)
                  end

                  mdot_c = mdot_c_inf + mdot_r
                  C_c = mdot_c * cp_c_in #Coolant heat capacity rate
                  C_min = min(C_c, C_h)
                  C_r = C_min / C_max

                  if C_h == C_max
                        ε_max = 1 / C_r * (1 - exp(-C_r)) #At ε = ε_max, NTU tends to infinity
                        ε = min(ε, 0.95 * ε_max) #Limit effectiveness to 95% of maximum possible
                  else
                        ε_max = 1 - exp(-1 / C_r) #At ε = ε_max, NTU tends to infinity
                        ε = min(ε, 0.95 * ε_max) #Limit effectiveness to 95% of maximum possible
                  end
            end
      
      else #No recirculation
            Tc_in = Tc_inf
            hc_in = hc_inf
            cp_c_in = cp_c_inf
      end
      
      C_c = mdot_c * cp_c_in #Coolant heat capacity rate

      C_min = min(C_c, C_h)
      C_max = max(C_c, C_h)
      C_r = C_min / C_max

      if C_c == C_min
            ε_max = 1 / C_r * (1 - exp(-C_r)) #At ε = ε_max, NTU tends to infinity
            ε = min(ε, 0.95 * ε_max) #Limit effectiveness to 95% of maximum possible
            NTU = -log(1 + log(1 - C_r * ε) / C_r) # For cross-flow with C_max mixed and C_min unmixed
      else
            ε_max = 1 - exp(-1 / C_r)#At ε = ε_max, NTU tends to infinity
            ε = min(ε, 0.95 * ε_max) #Limit effectiveness to 95% of maximum possible
            NTU = -1 / C_r * log(1 + C_r * log(1 - ε) ) # For cross-flow with C_max unmixed and C_min mixed
      end

      # Calculate total heat transfer and exit temperatures
      Qmax = C_min * (Th_in - Tc_in) #Maximum heat transfer rate
      Q = ε * Qmax #Actual heat transfer rate

      Th_out_guess = Th_in - Q / C_h 
      Tc_out_guess = Tc_in + Q / C_c

      hh_out = hh_in - Q / mdot_h
      hc_out = hc_in + Q / mdot_c

      Δh_h = hh_out - hh_in

      Th_out = gas_tset(alpha_h, length(alpha_h), hh_out, Th_out_guess)
      Tc_out = gas_tset_single(igas_c, hc_out, Tc_out_guess)

      #---------------------------------
      # Fluid calculations
      #---------------------------------
      γ_h_in = cp_h_in / (cp_h_in - Rh)
      γ_c_in = cp_c_in / (cp_c_in - Rc)

      ρ_h_in = ph_in / (Rh * Th_in)
      ρ_c_in = pc_in / (Rc * Tc_in)

      Vh_in = Mh_in * sqrt(γ_h_in * Rh * Th_in)
      Vc_in = Mc_in * sqrt(γ_c_in * Rc * Tc_in)

      #---------------------------------
      # Mean fluid properties
      #---------------------------------
      # Evaluate gas properties at a mean temperature and use these to find heat transfer coefficients 
      # See Kays. Compact Heat Exchangers (1984), p. 106

      Th_m = (Th_out + Th_in) / 2 #Mean temperature of hot stream
      Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of cold stream

      _, Pr_h_m, _, _, μ_h_m, k_h_m = gasPr(gas_h, Th_m)
      _, Pr_c_m, _, cp_c_m, μ_c_m, _ = gasPr(gas_c, Tc_m)

      ρ_h_m = ph_in / (Rh * Th_m) #Assume pressure is constant TODO: consider adding Rayleigh flow model
      ρ_c_m = pc_in / (Rc * Tc_m)

      ν_c_m = μ_c_m / ρ_c_m

      Vc_m = ρ_c_in * Vc_in / ρ_c_m #conservation of mass

      #---------------------------------
      # Geometry calculations
      #---------------------------------
      A_cs = mdot_h / (ρ_h_in * Vh_in) #Cross-sectional area of freestream

      if fconc #If channel is concentric, e.g., engine core
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter
            b = pi * D_i #Inner circumference
      else #If channel is rectangular
            b = A_cs / l #Side length
      end

      A_cc = mdot_c / (ρ_c_in * Vc_in) #total coolant cross-sectional area

      K = pi * b * n_stages / (4 * xt_D * A_cc) #Constant for next equation
      tD_o = (4 * K * t + sqrt(8 * K *t + 1) + 1) / (2 * K) #Compute tube outer diameter
      tD_i = tD_o - 2 * t #tube inner diameter

      N_t = 4 * A_cc / (pi * tD_i^2 * n_stages) #number of different coolant tubes per row

      xtm_D = A_cs / (N_t * tD_o * l) #Mean tangential pitch to diameter ratio

      A_D = N_t * l * tD_o #Total cross-sectional area of tubes in row, for drag calculation
      A_min = A_cs - A_D
      G = mdot_h / A_min #mass flow rate per unit area at minimum area

      #---------------------------------
      # Calculate thermal resistance 
      #---------------------------------

      #Calculate heat transfer coefficient for coolant
      Re_D_c = Vc_m * tD_i / ν_c_m #Reynolds number based on pipe diameter
      jc, Cf = jcalc_pipe(Re_D_c) #Colburn j-factor and skin-friction coefficient
      h_c = ρ_c_m * Vc_m * cp_c_m * jc / Pr_c_m^(2/3)

      #---------------------------------
      # Iterative loop
      #---------------------------------
      N_iter = 15 #Expect fast convergence

      n_passes = 4 #Initialize number of passes
      Ah = 0
      for i = 1:N_iter
            N_L = n_passes * n_stages #total number of rows

            # Calculate heat transfer coefficient for air
            Re_D_h = G * tD_o / μ_h_m #Reynolds number based on minimum free flow and tube outer diameter
            Nu_h = Nu_calc_staggered_cyl(Re_D_h, Pr_h_m, N_L, xtm_D, xl_D) #Nusselt number
            h_h = Nu_h * k_h_m / tD_o

            #Overall thermal resistance times area (m^2 K / W)
            RA = 1 / ( h_c * (tD_i/tD_o) ) + Rfh + tD_o / tD_i * Rfc + t / kw + 1 / h_h 

            # Size heat exchanger
            Ah = NTU * C_min * RA   #Find required hot-side cooling area from NTU
            n_passes = Ah / (N_t * n_stages * pi * tD_o * l)
      end

      #---------------------------------
      # Compute pressure drops
      #---------------------------------
      N_tubes_tot = N_t * n_passes * n_stages #total number of tubes across all rows
      N_L = n_passes * n_stages #total number of rows
      L = N_L * xl_D * tD_o #total axial length

      #Volumetric hydraulic diameter; Dv = 4 * (Net free volume) / (Friction surface)
      NFV = A_cs * L - N_tubes_tot * pi * tD_o^2 * l / 4 #Net free volume
      Dv = 4 * NFV / Ah

      Re_Dv = Dv * G / μ_h_m #Reynolds number based on hydraulic diameter and minimum free area

      if Re_Dv < 0 #If whole section is blocked
            Δp_h = Inf
      else
            Δp_h = Δp_calc_staggered_cyl(Re_Dv, G, L, ρ_h_m, Dv, tD_o, xtm_D, xl_D) #Calculate using the method of Gunter and Shaw (1945)
      end

      Pl_h = Δp_h * mdot_h / ρ_h_m #Power loss due to pressure drop in hot stream

      # Compute coolant pressure drop
      τw = ρ_c_m * Vc_m^2 / 2 * Cf
      A_s_c = pi * tD_i * l * n_passes #Surface area on one coolant stream
      A_cs_c = pi * tD_i^2 / 4 #cross-sectional area of one coolant stream
      Δp_c = τw * A_s_c / A_cs_c

      Pl_c = Δp_c * mdot_c / ρ_c_m #Power loss due to pressure drop in cold stream

      #---------------------------------
      # Output structs
      #---------------------------------
      #These lines rewrite the input structs

      #Gas parameters
      HXgas.Th_out = Th_out
      HXgas.Tc_out = Tc_out
      HXgas.Δh_h = Δh_h
      HXgas.Δp_h = Δp_h
      HXgas.Δp_c = Δp_c
      HXgas.Pl_h = Pl_h
      HXgas.Pl_c = Pl_c
      HXgas.ε = ε

      #Gas parameters
      HXgeom.N_t = N_t
      HXgeom.n_passes = n_passes
      HXgeom.tD_o = tD_o
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
function hxoper!(HXgas, HXgeom)
      #---------------------------------
      # Extract inputs
      #---------------------------------
      #Gas parameters
      gas_h = HXgas.gas_h
      gas_c = HXgas.gas_c
      alpha_h = HXgas.alpha_h
      igas_c = HXgas.igas_c 
      mdot_h = HXgas.mdot_h
      mdot_c_inf = HXgas.mdot_c 
      ε = HXgas.ε
      Th_in = HXgas.Th_in
      Tc_inf = HXgas.Tc_in 
      ph_in = HXgas.ph_in
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
      kw = HXgeom.kw
      Rfh = HXgeom.Rfh
      Rfc = HXgeom.Rfc 
      l = HXgeom.l
      A_cs = HXgeom.A_cs

      if frecirc
            recircT = HXgas.recircT
            h_lat = HXgas.h_lat

            if isnan(h_lat) #If the latent heat was not specified, assume it is 0
                  h_lat = 0
            end
      end

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------

      _, _, hh_in, _, cp_h_in, Rh = gassum(alpha_h, length(alpha_h), Th_in)
      _, _, hc_inf, _, cp_c_inf, Rc = gasfun(igas_c, Tc_inf)

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

      ρ_h_in = ph_in / (Rh * Th_in)
      ρ_c_in = pc_in / (Rc * Tc_in)

      # Calculate hot-side velocity from geometry
      Vh_in = mdot_h / (ρ_h_in * A_cs) #hot freestream velocity
      N_hyd_ways = N_t * n_stages #number of different hydraulic pathways

      #---------------------------------
      # Geometry calculations
      #---------------------------------
      tD_i = tD_o - 2 * t #tube inner diameter
      A_cs_tube = pi * tD_i^2 / 4 #Tube cross-sectional area

      A_D = N_t * l * tD_o #Total cross-sectional area of tubes in row, for drag calculation
      A_min = A_cs - A_D
      G = ρ_h_in * Vh_in * A_cs / A_min #mass flow rate per unit area at minimum area

      N_tubes_tot = N_t * n_passes * n_stages #total number of tubes across all rows
      N_L = n_passes * n_stages #total number of rows
      L = N_L * xl_D * tD_o #total axial length

      Ah = N_tubes_tot * pi * tD_o * l #Total external surface area of cooling tubes

      xtm_D = A_cs / (N_t * tD_o * l) #Mean tangential pitch to diameter ratio
      
      #---------------------------------
      # Thermal calculations
      #---------------------------------

      # Compute C_h and guess C_min
      C_h = mdot_h * cp_h_in #Hot-side heat capacity rate
      C_c = mdot_c_inf * cp_c_in #Coolant heat capacity rate, guess

      C_min = min(C_c, C_h)

      Qmax = C_min * (Th_in - Tc_in) #Maximum possible heat transfer rate, guess

      #---------------------------------
      # Iterative loop
      #---------------------------------
      ε = 0.95 #guess for effectiveness

      Qg = ε * Qmax #Actual heat transfer rate, guess

      # Guess outlet temperatures
      Th_out = Th_in - Qg / C_h 
      Tc_out = Tc_in + Qg / C_c

      N_iter = 20 #Rapid convergence expected

      ρ_h_m = 0 #Initiliaze because of annoying Julia scope
      μ_h_m = 0
      Δh_h = 0

      for i = 1 : N_iter
            
            if frecirc
                  #Calculate assuming that C_c = C_min
                  mdot_r = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * cp_c_in * (Th_in -  recircT) )
                  
                  C_check = (mdot_c_inf + mdot_r) * cp_c_in #Coolant heat capacity rate, to check validity of above assumption

                  C_max = max(C_check, C_h)

                  if C_check == C_max #If the calculation above is incorrect because C_c = C_max
                        A = mdot_c_inf * (hc_in - hc_inf + h_lat) / (ε * C_h * (Th_in -  recircT) )
                        if A > 1
                              println("Insufficient heat capacity in hot stream")
                              return
                        end
                        mdot_r = A * mdot_c_inf / (1 - A)
                  end
                  mdot_c = mdot_c_inf + mdot_r
            end

            #Heat capacity rates
            C_c = mdot_c * cp_c_in #Coolant heat capacity rate

            C_min = min(C_c, C_h)
            C_max = max(C_c, C_h)
            C_r = C_min / C_max

            Qmax = C_min * (Th_in - Tc_in) #Maximum heat transfer rate
            
            #Properties at mean temperature
            Th_m = (Th_out + Th_in) / 2 #Mean temperature of hot stream
            Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of cold stream

            _, Pr_h_m, _, _, μ_h_m, k_h_m = gasPr(gas_h, Th_m)
            _, Pr_c_m, _, cp_c_m, μ_c_m, _ = gasPr(gas_c, Tc_m)

            ρ_h_m = ph_in / (Rh * Th_m) #Assume pressure is constant TODO: consider adding Rayleigh flow model
            ρ_c_m = pc_in / (Rc * Tc_m)

            ν_c_m = μ_c_m / ρ_c_m

            Vc_m =  mdot_c / (N_hyd_ways * ρ_c_m * A_cs_tube) #cold velocity at mean temperature

            # Calculate thermal resistance
            # Calculate heat transfer coefficient for air
            Re_D_h = G * tD_o / μ_h_m #Reynolds number based on minimum free flow and tube outer diameter
            Nu_h = Nu_calc_staggered_cyl(Re_D_h, Pr_h_m, N_L, xtm_D, xl_D) #Nusselt number
            h_h = Nu_h * k_h_m / tD_o

            #Calculate heat transfer coefficient for coolant
            Re_D_c = Vc_m * tD_i / ν_c_m #Reynolds number based on pipe diameter
            jc, Cf = jcalc_pipe(Re_D_c) #Colburn j-factor
            h_c = ρ_c_m * Vc_m * cp_c_m * jc / Pr_c_m^(2/3)

            #Overall thermal resistance times area (m^2 K / W)
            RA = 1 / ( h_c * (tD_i/tD_o) ) + Rfh + tD_o / tD_i * Rfc + t / kw + 1 / h_h 

            # Find NTU and use it to calculate effectiveness
            NTU = Ah / (C_min * RA)

            if C_c == C_min
                  ε = 1 / C_r *(1 - exp(-C_r * (1 - exp(-NTU)) ) ) #effectiveness for cross flow with C_max mixed and C_min unmixed
            else
                  ε = 1 - exp(-1 / C_r * (1 - exp(-C_r * NTU) ) ) #effectiveness for cross flow with C_max unmixed and C_min mixed
            end

            # Calculate total heat transfer and exit temperatures
            Q = ε * Qmax #Actual heat transfer rate

            Th_out_guess = Th_in - Q / C_h 
            Tc_out_guess = Tc_in + Q / C_c

            hh_out = hh_in - Q / mdot_h
            hc_out = hc_in + Q / mdot_c

            Δh_h = hh_out - hh_in

            Th_out = gas_tset(alpha_h, length(alpha_h), hh_out, Th_out_guess)
            Tc_out = gas_tset_single(igas_c, hc_out, Tc_out_guess)

      end

      #---------------------------------
      # Compute pressure drop
      #---------------------------------

      #Volumetric hydraulic diameter; Dv = 4 * (Net free volume) / (Friction surface)
      NFV = A_cs * L - N_tubes_tot * pi * tD_o^2 * l / 4 #Net free volume
      Dv = 4 * NFV / Ah

      Re_Dv = Dv * G / μ_h_m #Reynolds number based on hydraulic diameter and minimum free area

      if Re_Dv < 0 #If whole section is blocked
            Δp_h = Inf
      else
            Δp_h = Δp_calc_staggered_cyl(Re_Dv, G, L, ρ_h_m, Dv, tD_o, xtm_D, xl_D) #Calculate using the method of Gunter and Shaw (1945)
      end

      #---------------------------------
      # Output structs
      #---------------------------------
      #Gas parameters
      HXgas.Th_out = Th_out
      HXgas.Tc_out = Tc_out
      HXgas.Δh_h = Δh_h
      HXgas.Δp_h = Δp_h
      HXgas.ε = ε

end #hxoper!

"""
    hxoptim!(HXgas, HXgeom, initial_x)

Optimizes heat exchanger design parameters for a given set of inputs. Uses the NLopt.jl package. The optimization
variables are `Mc_in`, `l`, `n_stages` and `xt_D`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
    - `initial_x::Vector{Float64}`: vector with the initial guess for the optimization
    
    **Outputs:**
    No direct outputs. Input structures are modified with HX design geometry.
"""
function hxoptim!(HXgas, HXgeom, initial_x)
      #Parameters to optimize: x[1]: 100 * Mc_in; x[2]: l; x[3]: n_stages; x[4]: xt_D
      #Set function to minimize
      obj(x, grad) =  hxobjf(x, HXgas, HXgeom) #Minimize objective function

      #Calculate minimum tube length
      alpha_h = HXgas.alpha_h #Hot gas composition
      mdot_h = HXgas.mdot_h #mass flow rate of hot gas (kg/s)
      Th_in = HXgas.Th_in
      Mh_in  = HXgas.Mh_in
      ph_in  = HXgas.ph_in

      #Flags 
      fconc = HXgeom.fconc

      _, _, _, _, cp_h_in, Rh = gassum(alpha_h, length(alpha_h), Th_in)
      γ_h_in = cp_h_in / (cp_h_in - Rh)
      ρ_h_in = ph_in / (Rh * Th_in)
      Vh_in = Mh_in * sqrt(γ_h_in * Rh * Th_in)

      A_cs = mdot_h / (ρ_h_in * Vh_in) #Cross-sectional area of freestream
      
      if fconc #Flow is concentric
            D_i = HXgeom.D_i
            D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

            lmin = (D_o - D_i) / 2 #minimum tube length
            
      else #square cross-section
            AR_min = 0.1 #Minimu aspect ratio
            lmin = sqrt(AR_min * A_cs)

      end

      #Set bounds
      lower = [0, lmin, 1, 1]
      upper = [30, Inf, 10, 6]
      
      #Use NLopt.jl to minimize function 
      opt = Opt(:LN_NELDERMEAD, length(initial_x))
      opt.lower_bounds = lower
      opt.upper_bounds = upper
      opt.ftol_rel = 1e-10
      
      opt.min_objective = obj
      
      (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
      
      #xopt_round = round.(xopt) #elements 2 and 3 must be integers

      #Modify structs with output
      HXgas.Mc_in = xopt[1] / 100 #x[1] has 100*Mc_in

      HXgeom.l = xopt[2]
      HXgeom.n_stages = xopt[3]
      HXgeom.xt_D = xopt[4]

      #Return optimum parameters by modifying input structs

end #hxoptim!

"""
      hxobjf(x, HXgas, HXgeom)

Objective function for HX optimization in hxoptim!(). It returns the sum of the power dissipated due to pressure
drops in the hot and cold streams, with penalty factors to enforce constraints.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x::Vector{Float64}`: state vector with [`100 * Mc_in`, `l`, `n_stages`, `xt_D`]
    - `HXgas::Struct`: structure of type HX_gas with the gas properties
    - `HXgeom::Struct`: structure of type HX_tubular with the HX geometric properties
 
    **Outputs:**
    - `I::Float64`: objective function (W)
"""
function hxobjf(x, HXgas, HXgeom)

      # Create local copy of structs
      HXg = deepcopy(HXgas)
      HXgeo = deepcopy(HXgeom)

      #Apply states
      HXg.Mc_in = x[1] / 100

      HXgeo.l = x[2]
      HXgeo.n_stages = x[3]
      HXgeo.xt_D = x[4]
      
      #Size HX
      hxsize!(HXg, HXgeo)

      #Extract outputs
      Pl_h = HXg.Pl_h
      Pl_c = HXg.Pl_c

      n_passes = HXgeo.n_passes
      N_t = HXgeo.N_t
      fconc = HXgeo.fconc

      vars = [n_passes, N_t]
      lower = [1, 1] #desired lower limits

      if fconc
            upper = [10, 200] #desired upper limits for concentric case
      else
            upper = [20, 200]  #allow more passes in rectangular case
      end
      I = (Pl_h + Pl_c) #Initialize objective function

      # Apply penalty function so that outputs are within desired range
      for (i,var) in enumerate(vars)
            vmin = lower[i]
            vmax = upper[i]

            pmin = ( vmin / min(vmin, var) )^2 #Penalty function
            pmax = ( max(vmax, var) / vmax )^2
            p = max(pmin, pmax)

            I = I * p
      end

      return I
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
 
    **Outputs:**
    - `HeatExchangers::Vector{Any}`: vector with heat exchanger data; elements are `HX_struct` structures
    - Also modifies `pare` with the fuel temperature and the HX enthalpy and pressure changes
"""
function hxdesign!(pare, pari, ipdes, HXs_prev)
      
      #---------------------------------
      # Extract inputs
      #---------------------------------
      pare_sl = pare[:, ipdes] #Slice pare at design point

      D_i = pare_sl[ieDi]
      Tc_ft = pare_sl[ieTft]
      frecirc = Bool(pare_sl[iefrecirc])
      recircT = pare_sl[ierecircT]
      h_lat = pare_sl[iehlat]
      igas = pari[iifuel]
      PreCorder = pare_sl[iePreCorder]
      PreCepsilon = pare_sl[iePreCepsilon]
      PreCMh = pare_sl[iePreCMh]
      InterCorder = pare_sl[ieInterCorder]
      InterCepsilon = pare_sl[ieInterCepsilon]
      InterCMh = pare_sl[ieInterCMh]
      Regenorder = pare_sl[ieRegenorder]
      Regenepsilon = pare_sl[ieRegenepsilon]
      RegenMh = pare_sl[ieRegenMh]
      TurbCorder = pare_sl[ieTurbCorder]
      TurbCepsilon = pare_sl[ieTurbCepsilon]
      TurbCMh = pare_sl[ieTurbCMh]

      if igas == 11 #TODO: add more options
            coolant_name = "ch4"
      elseif igas == 40
            coolant_name = "h2"
      end

      # Sort heat exchangers
      all_types = ["PreC", "InterC", "Regen", "TurbC"]
      all_orders = [PreCorder, InterCorder, Regenorder, TurbCorder]
      all_Mh = [PreCMh, InterCMh, RegenMh, TurbCMh]
      all_eps = [PreCepsilon, InterCepsilon, Regenepsilon, TurbCepsilon]

      HXtypes = []
      Mh_in = []
      ε_des = []
      sort_i = sortperm(all_orders) #Sort according to order

      for ind in sort_i
            if (all_eps[ind] > 0) && (all_eps[ind] <= 1) #If effectiveness is between 0 and 1
                  push!(HXtypes, all_types[ind])
                  push!(Mh_in, all_Mh[ind])
                  push!(ε_des, all_eps[ind])
            end
      end

      alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127] #Air composition

      #Initialize Heat Exchanger vector
      HeatExchangers = []

      #Initiliaze structures with NaNs
      HXgas_NaN = HX_gas("0","0", [NaN], NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
      HXgeom_NaN = HX_tubular(0, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)

      for (i,type) in enumerate(HXtypes) #For every desired type of heat exchanger (skipped if empty)
            
            #---------------------------------
            # Design exchangers
            #---------------------------------
            pare_sl = pare[:, ipdes] #Slice pare at design point
            #Initiliaze design geometry and gas property as NaNs
            HXgas = deepcopy(HXgas_NaN)
            HXgeom = deepcopy(HXgeom_NaN)

            HXgas.ε = ε_des[i]

            # Heat exchanger materials and wall properties
            HXgeom.t = 0.03e-2 #m, wall thicknesss
            HXgeom.xl_D = 1
            HXgeom.kw = 45 #thermal conductivity of steel, W/m/K
            HXgeom.ρw = 7930
            HXgeom.Rfc = 8.815E-05 #Hydrogen gas fouling resistance, m^2*K/W

            mcore = pare_sl[iemcore]
            mofft = pare_sl[iemofft]
            fc = pare_sl[iefc] #Extract cooling gas factor
            fo = mofft / mcore
            ff = pare_sl[ieff]

            HXgas.gas_h = "air"
            HXgas.alpha_h = alpha #Use alpha by default, except for Regen

            if type == "PreC" #Compressor Precooler
                  HXgeom.fconc = 1 #Concentric
                  HXgeom.D_i = D_i 
                  HXgeom.Rfh = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 

                  HXgas.mdot_h = mcore   #Core mass flow 
                  iTh_in = ieTt19
                  iph_in = iept19
                  ipc_in = iept3

                  Dh_i = iePreCDeltah
                  Dp_i = iePreCDeltap

            elseif type == "InterC" #Compressor Intercooler
                  HXgeom.fconc = 1 #Concentric
                  HXgeom.D_i = D_i
                  HXgeom.Rfh = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 

                  HXgas.mdot_h = mcore * (1 - fo) #Core mass flow minus offtake
                  iTh_in = ieTt25
                  iph_in = iept25
                  ipc_in = iept3

                  Dh_i = ieInterCDeltah
                  Dp_i = ieInterCDeltap

            elseif type == "Regen" #Regenerative cooling
                  HXgeom.fconc = 1 
                  HXgeom.D_i = D_i
                  HXgeom.Rfh = 0.01*0.1761 #Engine exhaust air fouling resistance, m^2*K/W 

                  HXgas.mdot_h = mcore * (1 - fo) #Core mass flow minus offtake
                  iTh_in = ieTt49
                  iph_in = iept49
                  ipc_in = iept3

                  Dh_i = ieRegenDeltah
                  Dp_i = ieRegenDeltap

                  HXgas.alpha_h = lambdap_calc(pare, alpha, igas, ipdes) #Calculate postcombustion and mixing composition

            elseif type =="TurbC" #Cooling of turbine cooling flow
                  HXgeom.fconc = 0 
                  HXgeom.Rfh = 0.001*0.1761 #Compressed air fouling resistance, m^2*K/W 

                  HXgas.mdot_h = mcore * fc #Only cooling mass flow rate
                  iTh_in = ieTt3
                  iph_in = iept3
                  ipc_in = iept3

                  Dh_i = ieTurbCDeltah
                  Dp_i = ieTurbCDeltap
            end

            HXgas.Th_in = pare_sl[iTh_in]
            HXgas.ph_in = pare_sl[iph_in]
            HXgas.Mh_in = Mh_in[i]

            HXgas.gas_c = coolant_name
            HXgas.igas_c = igas
            HXgas.mdot_c = mcore * ff #Fuel fraction times core mass flow rate
            
            HXgas.pc_in = pare_sl[ipc_in]

            if i == 1 #At first heat exchanger
                  HXgas.Tc_in = Tc_ft #Coolant temperature is the tank temperature
                  if frecirc #There can only be recirculation in the first heat exchanger
                        HXgeom.frecirc = 1 
                        HXgas.recircT = recircT
                        HXgas.h_lat = h_lat
                  end
                        
            else # For subsequent exchangers
                  HXprev = HeatExchangers[i - 1] #Get previous heat exchanger struct
                  all_gas_prev = HXprev.HXgas_mission
                  HXgas.Tc_in = all_gas_prev[ipdes].Tc_out #The inlet temperature is the outlet of previous HX at design point
                  
            end

            # Guess starting point for optimization
            #First calculate minimum tube length
            _, _, _, _, cp_h_in, Rh = gassum(HXgas.alpha_h, length(HXgas.alpha_h), HXgas.Th_in)
            γ_h_in = cp_h_in / (cp_h_in - Rh)
            ρ_h_in = HXgas.ph_in / (Rh * HXgas.Th_in)
            Vh_in = HXgas.Mh_in * sqrt(γ_h_in * Rh * HXgas.Th_in)

            A_cs = HXgas.mdot_h / (ρ_h_in * Vh_in) #Cross-sectional area of freestream
            
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

                  initial_x = [3, linit, 4, 4] #Initial guess
            else 
                  #x[1]: 100 * Mc_in; x[2]: l; x[3]: n_stages; x[4]: xt_D
                  initial_x = [100 * HXs_prev[i].HXgas_mission[ipdes].Mc_in, max(HXs_prev[i].HXgeom.l, lmin), 
                  HXs_prev[i].HXgeom.n_stages, HXs_prev[i].HXgeom.xt_D] #guess is previous iteration design point
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
                        HXgasp.mdot_h = mcore 

                  elseif (type == "InterC") 
                        HXgasp.mdot_h = mcore * (1 - fo)

                  elseif (type == "Regen")
                        HXgasp.mdot_h = mcore * (1 - fo)
                        HXgasp.alpha_h = lambdap_calc(pare, alpha, igas, ip) #Calculate postcombustion and mixing composition
                        
                  elseif type =="TurbC"
                        HXgasp.mdot_h = mcore * fc

                  end        

                  HXgasp.Th_in = pare_sl[iTh_in] #The indices come from the design process above as the HX is the same
                  HXgasp.ph_in = pare_sl[iph_in]

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

                  if HXgasp.mdot_h == 0 #If the mass flow rate in this mission is 0, nothing happens
                        HXgasp.Th_out = HXgasp.Th_in
                        HXgasp.Tc_out = HXgasp.Tc_in
                        HXgasp.Δh_h = 0
                        HXgasp.Δp_h = 0
                        HXgasp.ε = 0
                  else #Otherwise, call HX off-design routine
                        hxoper!(HXgasp, HXgeom)

                  end

                  HXgas_mis[ip] = HXgasp

                  #Store output in pare
                  pare[Dh_i, ip] = HXgasp.Δh_h
                  pare[Dp_i, ip] = HXgasp.Δp_h
                  
            end
            push!(HeatExchangers, HX_struct(type, HXgeom, HXgas_mis)) #Store HX struct in overall array
      end

      #---------------------------------
      # Update fuel temperature
      #---------------------------------

      for ip = 1:size(pare)[2] #For every mission point
            if length(HeatExchangers) > 0
                  lastHX = HeatExchangers[end]
                  HXgas = lastHX.HXgas_mission[ip]
                  Tf = HXgas.Tc_out
                  _, _, hf, _, _, _ = gasfun(igas, Tf)

                  pare[ieTfuel, ip] = Tf
                  pare[iehfuel, ip] = hf
            end
      end
     
      return HeatExchangers
end #hxdesign!

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
function jcalc_pipe(Re_D)
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
function Nu_calc_staggered_cyl(Re_D, Pr, N_L, xt_D, xl_D)

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
    Δp_calc_staggered_cyl(Re, G, L, ρ, Dv, tD_o, xt_D, xl_D)

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
 
    **Outputs:**
    - `Δp::Float64`: pressure drop across staggered cylinders (Pa)
"""
function Δp_calc_staggered_cyl(Re, G, L, ρ, Dv, tD_o, xt_D, xl_D)

      #Compute friction factor, f_2 = f/2
      if (Re <= 200)
            f_2 = 90 / Re 
      else
            f_2 = 0.96 * Re ^ (-0.145)
      end

      #Calculate pressure drop
      Δp = G^2 * L / (Dv * ρ) * f_2 * (Dv / (xt_D * tD_o) )^0.4 * (xl_D / xt_D)^0.6
      
      return Δp
end #Δp_calc_staggered_cyl

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
      N_t = HXgeom.N_t
      n_stages = HXgeom.n_stages
      n_passes = HXgeom.n_passes
      tD_o = HXgeom.tD_o 
      ρ = HXgeom.ρw
      l = HXgeom.l
      t = HXgeom.t

      N_tubes_tot = N_t * n_passes * n_stages #total number of tubes across all rows
      tD_i = tD_o - 2 * t
      V_t = N_tubes_tot * pi * (tD_o^2 - tD_i^2) / 4 * l #total tube volume
      m_t = ρ * V_t #total tube mass

      W_hx = gee * m_t * (1 + fouter)

      return W_hx
end #hxweight

"""
      lambdap_calc(pare, alpha_in, ifuel, ip, im)

Calculates the weight of a heat exchanger with involute tubes.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pare::Array{Float64 , 3}`: array with engine parameters
    - `alpha_in::Vector{Float64}`: vector with hot-side gas composition before combustion and mixing
    - `ifuel::Float64`: fuel gas index
    - `ip::Float64`: index for mission point
    
    **Outputs:**
    - `lambda_p::Vector{Float64}`: vector with hot-side gas composition after combustion and mixing
"""
function lambdap_calc(pare, alpha_in, ifuel, ip)
      #Extract inputs
      pare_sl = pare[:, ip]

      Tt3 = pare_sl[ieTt3]
      Ttf = pare_sl[ieTfuel]
      Tt4 = pare_sl[ieTt4]

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

      ffb, lambda = gas_burn(alpha, beta, gamma, n, ifuel, Tt3, Ttf, Tt4)

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