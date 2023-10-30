# hxfun.jl
# These functions can be used to size and model a heat exchanger with involute staggered tubes in a crossflow
# The design method is based on the effectiveness-NTU method, described in many sources such as 
# https://www.mathworks.com/help/hydro/ref/entuheattransfer.html
# Nicolas Gomez Vega, Oct 2023
using NLopt

"""
    hxsize(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, Mc_in, ph_in, pc_in, D_i, t, 
    n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)

Sizes a crossflow heat exchanger and calculates the pressure drop. Uses the ε-NTU method to size the heat exchanger
from a prescribed ε.  For representative fouling factors see Standards of the Tubular Exchanger Manufacturers Association
or https://powderprocess.net/Tools_html/Data_Diagrams/Heat_Exchanger_Fouling_Factor.html 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gas_h::Char`: hot gas name
    - `gas_c::Char`: coolant gas name
    - `ε::Float64`: desired heat exchanger effectiveness
    - `mdot_h::Float64`: mass flow rate of hot gas (kg/s)
    - `mdot_c::Float64`: mass flow rate of cold gas (kg/s)
    - `Th_in::Float64`: hot gas inlet temperature (K)
    - `Tc_in::Float64`: cold gas inlet temperature (K)
    - `Mh_in::Float64`: hot gas inlet Mach number
    - `Mc_in::Float64`: cold gas inlet Mach number
    - `ph_in::Float64`: hot gas inlet pressure (Pa)
    - `pc_in::Float64`: cold gas inlet pressure (Pa)
    - `D_i::Float64`: inner diameter of core (m)
    - `t::Float64`: cooling tube thickness (m)
    - `n_passes::Float64`: number of coolant passes
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
    - `xt_D::Float64`: circumferential pitch between tubes at the root over tube outer diameter 
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
    - `kw::Float64`: thermal conductivity of wall material (W/m/K)
    - `Rfh::Float64`: hot-side fouling factor (m^2 K/W)
    - `Rfc::Float64`: cold-side fouling factor (m^2 K/W)
 
    **Outputs:**
    - `Th_out::Float64`: hot gas outlet temperature
    - `Tc_out::Float64`: cold gas outlet temperature
    - `Δp_h::Float64`: pressure drop of hot gas across heat exchanger (Pa)
    - `tD_o::Float64`: tube outer diameter (m)
    - `N_t::Float64`: number of tubes per row
    - `l::Float64`: length of tubes (m)
    - `A_cs::Float64`: cross-sectional area of HX freestream (m^2)
    - `D_i::Float64`: outer diameter of core (m)
    - `Pl_h::Float64`: power loss due to pressure drop in hot stream (W)
    - `Pl_c::Float64`: power loss due to pressure drop in cold stream (W)
"""
function hxsize(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, Mc_in, ph_in, pc_in, D_i, t, 
      n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------

      Rh, _, γ_h_in, cp_h_in, _, _ = gasPr(gas_h, Th_in)
      Rc, _, γ_c_in, cp_c_in, _, _ = gasPr(gas_c, Tc_in)

      #---------------------------------
      # Fluid calculations
      #---------------------------------

      ρ_h_in = ph_in / (Rh * Th_in)
      ρ_c_in = pc_in / (Rc * Tc_in)

      Vh_in = Mh_in * sqrt(γ_h_in * Rh * Th_in)
      Vc_in = Mc_in * sqrt(γ_c_in * Rc * Tc_in)

      #---------------------------------
      # Thermal calculations
      #---------------------------------

      # Compute C_min, C_max and C_r
      C_c = mdot_c * cp_c_in #Coolant heat capacity rate
      C_h = mdot_h * cp_h_in #Air heat capacity rate

      C_min = min(C_c, C_h)
      C_max = max(C_c, C_h)
      C_r = C_min / C_max

      NTU = -log(1 + log(1 - C_r * ε) / C_r) # For cross-flow with C_max mixed and C_min unmixed

      # Calculate total heat transfer and exit temperatures
      Qmax = C_min * (Th_in - Tc_in) #Maximum possible heat transfer rate
      Q = ε * Qmax #Actual heat transfer rate

      Th_out = Th_in - Q / C_h #TODO: replace with gas_calc calls
      Tc_out = Tc_in + Q / C_c

      #---------------------------------
      # Geometry calculations
      #---------------------------------

      A_cs = mdot_h / (ρ_h_in * Vh_in) #Cross-sectional area of freestream
      D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

      A_cc = mdot_c / (ρ_c_in * Vc_in) #total coolant cross-sectional area

      K = pi^2 * D_i * n_stages / (4 * xt_D * A_cc) #Constant for next equation
      tD_o = (4 * K * t + sqrt(8 * K *t + 1) + 1) / (2 * K) #Compute tube outer diameter
      tD_i = tD_o - 2 * t #tube inner diameter

      N_t = 4 * A_cc / (pi * tD_i^2 * n_stages) #number of different coolant tubes per row
      N_tubes_tot = N_t * n_passes * n_stages #total number of tubes across all rows
      N_L = n_passes * n_stages #total number of rows
      L = N_L * xl_D * tD_o #total axial length

      #---------------------------------
      # Mean fluid properties
      #---------------------------------
      # Evaluate gas properties at a mean temperature and use these to find heat transfer coefficients 
      # See Kays. Compact Heat Exchangers (1984), p. 106

      Th_m = (Th_out + Th_in) / 2 #Mean temperature of hot stream
      Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of cold stream

      _, Pr_h_m, _, cp_h_m, μ_h_m, k_h_m = gasPr(gas_h, Th_m)
      _, Pr_c_m, _, cp_c_m, μ_c_m, k_c_m = gasPr(gas_c, Tc_m)

      ρ_h_m = ph_in / (Rh * Th_m) #Assume pressure is constant TODO: consider adding Rayleigh flow model
      ρ_c_m = pc_in / (Rc * Tc_m)

      ν_h_m = μ_h_m / ρ_h_m
      ν_c_m = μ_c_m / ρ_c_m

      Vh_m = ρ_h_in * Vh_in / ρ_h_m #conservation of mass
      Vc_m = ρ_c_in * Vc_in / ρ_c_m

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
      N_iter = 5 #Expect fast convergence

      lmin = (D_o - D_i) / 2 #minimum tube length, core width
      l = lmin #initial guess
      G = 0 
      Ah = 0
      xtm_D = 0
      for i = 1:N_iter
            #TODO: think of a better way to model this 
            lxt = max(l, lmin) #Do not account for the xtm change if the length is smaller than the core width
            xtm_D = A_cs / (N_t * tD_o * lxt) #Mean tangential pitch to diameter ratio
            
            A_D = N_t * l * tD_o #Total cross-sectional area of tubes in row, for drag calculation
            A_min = A_cs - A_D
            G = mdot_h / A_min #mass flow rate per unit area at minimum area

            # Calculate heat transfer coefficient for air
            Re_D_h = G * tD_o / μ_h_m #Reynolds number based on minimum free flow and tube outer diameter
            Nu_h = Nu_calc_staggered_cyl(Re_D_h, Pr_h_m, N_L, xtm_D, xl_D) #Nusselt number
            h_h = Nu_h * k_h_m / tD_o

            #Overall thermal resistance times area (m^2 K / W)
            RA = 1 / ( h_c * (tD_i/tD_o) ) + Rfh + tD_o / tD_i * Rfc + t / kw + 1 / h_h 
      
            # Size heat exchanger
            Ah = NTU * C_min * RA   #Find required hot-side cooling area from NTU
            l = Ah / (N_tubes_tot * pi * tD_o)
      end

      #---------------------------------
      # Compute pressure drops
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

      Pl_h = Δp_h * mdot_h / ρ_h_m #Power loss due to pressure drop in hot stream

      # Compute coolant pressure drop
      τw = ρ_c_m * Vc_m^2 / 2 * Cf
      A_s_c = pi * tD_i * l * n_passes #Surface area on one coolant stream
      A_cs_c = pi * tD_i^2 / 4 #cross-sectional area of one coolant stream
      Δp_c = τw * A_s_c / A_cs_c
      
      Pl_c = Δp_c * mdot_c / ρ_c_m #Power loss due to pressure drop in cold stream

      return Th_out, Tc_out, Δp_h, tD_o, N_t, l, A_cs, D_o, Pl_h, Pl_c
end #hxsize


"""
    hxoper(gas_h, gas_c, mdot_h, mdot_c, Th_in, Tc_in, ph_in, pc_in, A_cs, tD_o, t, 
    l, N_t, n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)

Evaluates crossflow heat exchanger performance for off-design operation. Uses the ε-NTU 
method to calculate effectiveness from prescribed geometry. 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gas_h::Char`: hot gas name
    - `gas_c::Char`: coolant gas name
    - `mdot_h::Float64`: mass flow rate of hot gas (kg/s)
    - `mdot_c::Float64`: mass flow rate of cold gas (kg/s)
    - `Th_in::Float64`: hot gas inlet temperature (K)
    - `Tc_in::Float64`: cold gas inlet temperature (K)
    - `ph_in::Float64`: hot gas inlet pressure (Pa)
    - `pc_in::Float64`: cold gas inlet pressure (Pa)
    - `A_cs::Float64`: cross-sectional area of HX freestream (m^2)
    - `tD_o::Float64`: tube outer diameter (m)
    - `t::Float64`: cooling tube thickness (m)
    - `l::Float64`: length of tubes (m)
    - `N_t::Float64`: number of tubes per row
    - `n_passes::Float64`: number of coolant passes
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
    - `kw::Float64`: thermal conductivity of wall material (W/m/K)
    - `Rfh::Float64`: hot-side fouling factor (m^2 K/W)
    - `Rfc::Float64`: cold-side fouling factor (m^2 K/W)
 
    **Outputs:**
    - `Th_out::Float64`: hot gas outlet temperature
    - `Tc_out::Float64`: cold gas outlet temperature
    - `Δp_h::Float64`: pressure drop of hot gas across heat exchanger (Pa)
"""
function hxoper(gas_h, gas_c, mdot_h, mdot_c, Th_in, Tc_in, ph_in, pc_in, A_cs, tD_o, t, 
      l, N_t, n_passes, n_stages, xl_D, kw, Rfh, Rfc)

      #---------------------------------
      # Inlet gas parameters
      #---------------------------------

      Rh, _, _, cp_h_in, _, _ = gasPr(gas_h, Th_in)
      Rc, _, _, cp_c_in, _, _ = gasPr(gas_c, Tc_in)

      #---------------------------------
      # Fluid calculations
      #---------------------------------

      ρ_h_in = ph_in / (Rh * Th_in)
      ρ_c_in = pc_in / (Rc * Tc_in)

      tD_i = tD_o - 2 * t #tube inner diameter
      A_cs_tube = pi *tD_i^2 / 4 #Tube cross-sectional area

      # Calculate velocities from geometry
      Vh_in = mdot_h / (ρ_h_in * A_cs) #hot freestream velocity
      N_hyd_ways = N_t * n_stages #number of different hydraulic pathways
      Vc_in =  mdot_c / (N_hyd_ways * ρ_c_in * A_cs_tube) #cold freestream velocity

      #---------------------------------
      # Geometry calculations
      #---------------------------------
      D_o = sqrt(4 * (A_cs + pi * D_i^2 / 4) / pi) #Core outer diameter

      A_D = N_t * l * tD_o #Total cross-sectional area of tubes in row, for drag calculation
      A_min = A_cs - A_D
      G = ρ_h_in * Vh_in * A_cs / A_min #mass flow rate per unit area at minimum area

      N_tubes_tot = N_t * n_passes * n_stages #total number of tubes across all rows
      N_L = n_passes * n_stages #total number of rows
      L = N_L * xl_D * tD_o #total axial length

      Ac = N_tubes_tot * pi * tD_i * l #Total internal surface area of cooling tubes
      Ah = N_tubes_tot * pi * tD_o * l #Total external surface area of cooling tubes

      lmin = (D_o - D_i) / 2 #minimum tube length
      lxt = max(l, lmin) 
      xtm_D = A_cs / (N_t * tD_o * lxt) 
      #---------------------------------
      # Thermal calculations
      #---------------------------------

      # Compute C_min, C_max and C_r
      C_c = mdot_c * cp_c_in #Coolant heat capacity rate
      C_h = mdot_h * cp_h_in #Air heat capacity rate

      C_min = min(C_c, C_h)
      C_max = max(C_c, C_h)
      C_r = C_min / C_max

      Qmax = C_min * (Th_in - Tc_in) #Maximum possible heat transfer rate

      #---------------------------------
      # Iterative loop
      #---------------------------------
      ε_guess = 0.5 #guess for effectiveness

      Qg = ε_guess * Qmax #Actual heat transfer rate

      # Guess outlet temperatures
      Th_out = Th_in - Qg / C_h 
      Tc_out = Tc_in + Qg / C_c

      N_iter = 5 #Rapid convergence expected

      ρ_h_m = 0 #Initiliaze because of annoying Julia scope
      μ_h_m = 0

      for i = 1 : N_iter

            #Properties at mean temperature
            Th_m = (Th_out + Th_in) / 2 #Mean temperature of hot stream
            Tc_m = (Tc_out + Tc_in) / 2 #Mean temperature of cold stream

            _, Pr_h_m, _, cp_h_m, μ_h_m, k_h_m = gasPr(gas_h, Th_m)
            _, Pr_c_m, _, cp_c_m, μ_c_m, k_c_m = gasPr(gas_c, Tc_m)

            ρ_h_m = ph_in / (Rh * Th_m) #Assume pressure is constant TODO: consider adding Rayleigh flow model
            ρ_c_m = pc_in / (Rc * Tc_m)

            ν_h_m = μ_h_m / ρ_h_m
            ν_c_m = μ_c_m / ρ_c_m

            Vh_m = ρ_h_in * Vh_in / ρ_h_m #conservation of mass
            Vc_m = ρ_c_in * Vc_in / ρ_c_m

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
            ε = 1 / C_r *(1 - exp(-C_r * (1 - exp(-NTU)))) #effectiveness for cross flow with C_max mixed and C_min unmixed

            # Calculate total heat transfer and exit temperatures
            
            Q = ε * Qmax #Actual heat transfer rate

            Th_out = Th_in - Q / C_h #TODO: replace with gas_calc calls
            Tc_out = Tc_in + Q / C_c

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

      return Th_out, Tc_out, Δp_h
end #hxoper

"""
    hxoptim(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, ph_in, pc_in, D_i, t, kw, Rfh, Rfc)

Optimizes heat exchanger design parameters for a given set of inputs. Uses the NLopt package.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gas_h::Char`: hot gas name
    - `gas_c::Char`: coolant gas name
    - `ε::Float64`: desired heat exchanger effectiveness
    - `mdot_h::Float64`: mass flow rate of hot gas (kg/s)
    - `mdot_c::Float64`: mass flow rate of cold gas (kg/s)
    - `Th_in::Float64`: hot gas inlet temperature (K)
    - `Tc_in::Float64`: cold gas inlet temperature (K)
    - `Mh_in::Float64`: hot gas inlet Mach number
    - `ph_in::Float64`: hot gas inlet pressure (Pa)
    - `pc_in::Float64`: cold gas inlet pressure (Pa)
    - `D_i::Float64`: inner diameter of core (m)
    - `t::Float64`: cooling tube thickness (m)
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
    - `kw::Float64`: thermal conductivity of wall material (W/m/K)
    - `Rfh::Float64`: hot-side fouling factor (m^2 K/W)
    - `Rfc::Float64`: cold-side fouling factor (m^2 K/W)
 
    **Outputs:**
    - `Mc_in::Float64`: cold gas inlet Mach number
    - `n_passes::Float64`: number of coolant passes
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
    - `xt_D::Float64`: circumferential pitch between tubes at root over tube outer diameter
"""
function hxoptim(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, ph_in, pc_in, D_i, t, xl_D, kw, Rfh, Rfc)

      #Parameters to optimize: x[1]: 100 * Mc_in; x[2]: n_passes; x[3]: n_stages; x[4]: xt_D
      #Set function to minimize
      obj(x, grad) =  hxobjf(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, x[1] / 100, ph_in, pc_in, D_i, t, 
      x[2], x[3], x[4], xl_D, kw, Rfh, Rfc) #Minimize objective function

      #Set bounds
      lower = [0, 0, 0, 1]
      upper = [30, 10, 10, 6]

      initial_x = [3, 8, 4, 4] #Initial guess
      
      #Use Optim.jl to minimize function 
      opt = Opt(:LN_NELDERMEAD, length(initial_x))
      opt.lower_bounds = lower
      opt.upper_bounds = upper
      opt.ftol_rel = 1e-4 #High tolerance to speed up convergence
      
      opt.min_objective = obj
      
      (minf,xopt,ret) = NLopt.optimize(opt, initial_x)
      
      #xopt_round = round.(xopt) #elements 2 and 3 must be integers

      Mc_in = xopt[1] / 100 #x[1] has 100*Mc_in

      n_passes = xopt[2]
      n_stages = xopt[3]
      xt_D = xopt[4]

      #Return optimum parameters
      return Mc_in, n_passes, n_stages, xt_D
end #hxoptim

"""
    hxobjf(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, Mc_in, ph_in, pc_in, D_i, t, 
    n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)

Objective function for HX optimization in hxoptim(). It returns the sum of the power dissipated due to pressure
drops in the hot and cold streams.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gas_h::Char`: hot gas name
    - `gas_c::Char`: coolant gas name
    - `ε::Float64`: desired heat exchanger effectiveness
    - `mdot_h::Float64`: mass flow rate of hot gas (kg/s)
    - `mdot_c::Float64`: mass flow rate of cold gas (kg/s)
    - `Th_in::Float64`: hot gas inlet temperature (K)
    - `Tc_in::Float64`: cold gas inlet temperature (K)
    - `Mh_in::Float64`: hot gas inlet Mach number
    - `Mc_in::Float64`: cold gas inlet Mach number
    - `ph_in::Float64`: hot gas inlet pressure (Pa)
    - `pc_in::Float64`: cold gas inlet pressure (Pa)
    - `D_i::Float64`: inner diameter of core (m)
    - `t::Float64`: cooling tube thickness (m)
    - `n_passes::Float64`: number of coolant passes
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
    - `xt_D::Float64`: circumferential pitch between tubes over tube outer diameter
    - `xl_D::Float64`: longitudinal pitch between rows over tube outer diameter
    - `kw::Float64`: thermal conductivity of wall material (W/m/K)
    - `Rfh::Float64`: hot-side fouling factor (m^2 K/W)
    - `Rfc::Float64`: cold-side fouling factor (m^2 K/W)
 
    **Outputs:**
    - `I::Float64`: objective function (W)
"""
function hxobjf(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, Mc_in, ph_in, pc_in, D_i, t, 
      n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)
      
      _, _, _, _, _, l, _, D_o, Pl_h, Pl_c = hxsize(gas_h, gas_c, ε, mdot_h, mdot_c, Th_in, Tc_in, Mh_in, Mc_in, 
      ph_in, pc_in, D_i, t, n_passes, n_stages, xt_D, xl_D, kw, Rfh, Rfc)

      lmin = (D_o - D_i) / 2 #minimum tube length

      #Introduce a penalty function to ensure that the tube length is not less 
      #than the distance between the inner and outer cylinders
      p = (lmin / min(lmin, l))^2 #Penalty function

      I = (Pl_h + Pl_c) * p #sum of power losses times penalty

      return I
end #hxobjf

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
    - `G::Float64`: mass flow rate divided by minimum free flow area. G = mdot / (A_min), A_min is the minimum free-flow area (kg/s/m^2)
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
    hxweight(gee, fouter, ρ_m, tD_o, t, l, N_t, n_passes, n_stages)

Calculates the weight of a heat exchanger with involute tubes.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gee::Float64`: gravitational acceleration (m/s^2)
    - `fouter::Float64`: ratio of mass of non-tube parts of HX to mass of the tubes
    - `ρ_m::Float64`: density of HX material (kg/m^3)
    - `tD_o::Float64`: cylinder outer diameter (m)
    - `t::Float64`: cooling tube thickness (m)
    - `l::Float64`: length of tubes (m)
    - `N_t::Float64`: number of tubes per row
    - `n_passes::Float64`: number of coolant passes
    - `n_stages::Float64`: number of different coolant stages with different coolant flows
 
    **Outputs:**
    - `W_hx::Float64`: weight of heat exchanger (N)
"""
function hxweight(gee, fouter, ρ_m, tD_o, t, l, N_t, n_passes, n_stages)
      # Calculates the weight of a heat exchanger with involute tubes

      # Input
      #-------
      # gee: gravitational acceleration (m/s^2)
      # fouter: ratio of mass of non-tube parts of HE to mass of the tubes
      # ρ_m: mean density of HE (kg/m^3)
      # tD_o: cooling tube outer diameter (m)
      # t: cooling tube thickness (m)
      # l: length of each tube (m)
      # N_t: number of cooling tubes per row
      # n_passes: number of coolant passes
      # n_stages: number of different coolant stages with different coolant flows

      # Output
      #-------
      # W_hx: weight of heat exchanger (N)

      N_tubes_tot = N_t * n_passes * n_stages #total number of tubes across all rows
      tD_i = tD_o - 2 * t
      V_t = N_tubes_tot * pi * (tD_o^2 - tD_i^2) / 4 * l #total tube volume
      m_t = ρ_m * V_t #total tube mass

      W_hx = gee * m_t * (1 + fouter)

      return W_hx
end #hxweight