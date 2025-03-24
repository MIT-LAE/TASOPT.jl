using DifferentialEquations

"""
    LT_PEMFC_voltage_simple(j, T, p_H2, p_air)

A simple 0D model of a low temperature PEM fuel cell, accounting for thermodynamic, activation, ohmic and concentration losses.      

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `j::Float64`: current density (A/m^2)
    - `T::Float64`: fuel cell temperature (K)
    - `p_H2::Float64`: inlet pressure of hydrogen at the anode (Pa)
    - `p_air::Float64`: inlet pressure of air at the cathode (Pa)
 
    **Outputs:**
    - `V::Float64`: voltage across cell (V)
"""
function LT_PEMFC_voltage_simple(j, T, p_H2, p_air)

    # Parameters
    R = 8.314 # J/(mol*K)
    F = 96485.3329 # C/mol
    T0 = 298.15  # Temperature in Kelvin
    p0 = 1.01325e5  # Pressure in Pascal
    n = 2 #number of electrons in half-reaction

    #Typical low-temperature PEMFC parameters (O'Hayre et al. 2016)
    j0_H2 = 1e3 #A/m^2
    α_H2 = 0.5
    j0_O2 = 1 #A/m^2
    α_O2 = 0.3
    ASR = 1e-6 #Ohm m^2
    j_leak = 1e2 #A/m^2
    j_L = 2e4 #A/m^2
    c = 0.1 #V

    # Calculations
    p_O2 = 0.2095 * p_air

    a_H2 = p_H2 / p0 #activity of hydrogen
    a_O2 = p_O2 / p0 #activity of oxygen

    T_vap = 383.15

    if T > T_vap #gaseous water product
        Δs0 = -44.34 # entropy change in chemical reaction, J/mol/K
        E0 = 1.185 # reversible voltage at STP, V
    else #liquid water product
        Δs0 = -163.23 # entropy change in chemical reaction, J/mol/K
        E0 = 1.229 # reversible voltage at STP, V
    end

    #Nernst equation for reversible voltage
    E_r = E0 + Δs0 / (n * F) * (T - T0) - R * T / (n * F) * log(1 / (a_H2 * a_O2^0.5))

    #Activation voltage from Tafel equation
    a_A = - R * T / (α_H2 * n * F) * log(j0_H2)
    b_A = R * T / (α_H2 * n * F)
    
    a_C = - R * T / (α_O2 * n * F) * log(j0_O2)
    b_C = R * T / (α_O2 * n * F)

    η_act = a_A + b_A * log(j + j_leak) + a_C + b_C * log(j + j_leak)

    #Ohmic losses
    η_ohm = j * ASR

    # Concentration losses
    η_conc = c * log(j_L / (j_L - (j + j_leak)))

    # Cell voltage
    V = E_r - η_act - η_ohm - η_conc
    return V
end

"""
    LT_PEMFC_voltage_OHayre(j, T, p_A, p_C, x_H2O_A, x_H2O_C, λ_O2)

A 1-D model of the voltage across a low-temperature fuel cell with a Nafion membrane, based on the model
in O'Hayre et al. (2016), which is a simplified version of that in Springer et al. (1991). 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `j::Float64`: current density (A/m^2)
    - `T::Float64`: fuel cell temperature (K)
    - `p_A::Float64`: inlet pressure of anode reactants (Pa)
    - `p_C::Float64`: inlet pressure of cathode reactants (Pa)
    - `x_H2O_A::Float64`: mole fraction of water vapour in the anode inlet
    - `x_H2O_A::Float64`: mole fraction of water vapour in the cathode inlet
    - `λ_O2::Float64`: stoichiometric ratio of oxygen in the cathode
 
    **Outputs:**
    - `V::Float64`: voltage across cell (V)
"""
function LT_PEMFC_voltage_OHayre(j, T, p_A, p_C, x_H2O_A, x_H2O_C, λ_O2)

    # Parameters from physics
    R = 8.314 # J/(mol*K), universal gas constant
    F = 96485.3329 # C/mol, Faraday contant
    n_O2 = 4 #number of electrons in reduction reaction
    n = 2 #number of electrons in oxidation reaction
    p0 = 101325 #Pa, reference pressure
    T0 = 298.15  #K, reference temperature
    n_drag = 2.5 #electro-osmotic drag coefficient in Nafion at saturation
    M_m = 1 #kg/mol, Nafion equaivalent weight, typically 1-1.1 kg/mol
    ρ_dry = 1970 #kg/m^3, dry density of Nafion
    x_ON = 0.21 #mole fraction of oxygen in dry air

    #Paremeters from model in O'Hayre (2016)
    λ_ref = 10 #reference water content for D_λ calculation
    α = 0.5 #symmetry parameter
    ε = 1 #porosity of electrodes
    τ = 1.5 #diffusion tortuosity
    j0 = 1 #A/m^2, reference current density at cathode
    t_M = 125e-6 #m, membrane thickness
    t_A = 350e-6 #m, anode thickness
    t_C = 350e-6 #m, cathode thickness

    #Calculate other mole fractions
    x_H2_in = 1 - x_H2O_A #Hydrogen mole fraction at anode inlet
    x_O2_in = x_ON * (1 - x_H2O_C) #Oxygen mole fraction at cathode inlet
    
    #Calculate other parameters
    D_H2H2O = binary_diffusion(T, p_A, ["H2", "H2O"]) #m^2/s, diffusion coefficient of H2 in H2O
    D_O2H2O = binary_diffusion(T, p_C, ["O2", "H2O"]) #m^2/s, diffusion coefficient of O2 in H2O

    Deff_H2H2O = porous_diffusion(D_H2H2O, ε, τ) #m^2/s, effective diffusion coefficient of H2 in H2O
    Deff_O2H2O = porous_diffusion(D_O2H2O, ε, τ) #m^2/s, effective diffusion coefficient of O2 in H2O

    D_λ = Nafion_diffusion(T, λ_ref) #m^2/s, water diffusivity in Nafion

    p_SAT = water_sat_pressure(T) #Pa, Water saturation pressure

    # Reversible voltage
    p_O2 = x_O2_in * p_C #Pa, partial pressure of oxygen in cathode
    p_H2 = x_H2_in * p_A  #Pa, partial pressure of hydrogen in anode

    a_H2 = p_H2 / p0 #activity of hydrogen
    a_O2 = p_O2 / p0 #activity of oxygen

    Δs0 = -163.23 #J/mol/K, entropy change in chemical reaction
    E0 = 1.229 #V,  reversible voltage at STP

    #Nernst equation for reversible voltage
    E_r = E0 + Δs0 / (n * F) * (T - T0) - R * T / (n * F) * log(1 / (a_H2 * a_O2^0.5))
 
    x_H2O_a = x_H2O_A
    x_H2O_d = x_H2O_C

    # Water balance calculation
    #λb = ξ1+ξ2*α_star = ξ5 * α_star + C
    #λc = ξ3+ξ4*α_star = ξ5 * α_star + ξ6 * C
    ξ1 = 14 * p_A / p_SAT * x_H2O_a
    ξ2 = - 14 * p_A / p_SAT * t_A * j * R * T / (2 * F * p_A * Deff_H2H2O)
    ξ3 = 12.6 + 1.4 * p_C / p_SAT * (x_H2O_d + t_C * j * R * T / (2 * F * p_C * Deff_O2H2O))
    ξ4 = 1.4 * p_C / p_SAT * t_C * j * R * T / (2 * F * p_C * Deff_O2H2O)
    ξ5 = 11 / n_drag
    ξ6 = exp(j * M_m * n_drag * t_M / (22 * F * ρ_dry * D_λ) )

    A = ξ2 - ξ5
    B = ξ4 - ξ5
    D = -ξ6
    α_star = 1 / (A * D + B) * (-ξ1 * D - ξ3) #ratio of water flux to proton flux
    C = 1 / (A * D + B) * (ξ1 * B - ξ3 * A) #integration constant in equation for λ

    # Ohmic losses 
    #σ = a + b * exp(c * z)
    #ASR = integral between 0 a t_M of dz/σ
    a = ( 0.005193 * (11 / n_drag * α_star) - 0.00326) * exp(1268 * (1/303 - 1/T) )
    b = ( 0.005193 * C ) * exp(1268 * (1/303 - 1/T) )
    c = j * M_m * n_drag / (22 * F * ρ_dry * D_λ) * 1e-2
    t_M_cm = t_M * 1e2 #cm, membrane thickness
    ASR = ((c * t_M_cm - log(a + b * exp(c * t_M_cm))) / (a * c) + log(a + b) / (a * c) ) * 1e-4 #Ohm*m^2
    η_ohm = j * ASR #V, ohmic losses

    # Cathode activation and concentration losses
    #Oxygen depletion
    x_O2_d = ((λ_O2 - 1) * (1 - x_H2O_C) * x_ON) / ( λ_O2 + (2 * α_star + 1) * (1 - x_H2O_C) * x_ON) #molar fraction of oxygen at outlet due to depletion
    
    x_O2_interf = x_O2_d - t_C * j * R * T / ( n_O2 * F * p_C * Deff_O2H2O ) #oxygen concentration at membrane interface

    if x_O2_interf <= 0 #If all oxygen has been depleted
        V = 0 #voltage is zero as there is not enough oxygen
    else
        η_C = R * T / (n_O2 * α * F) * log( j / ( j0 * p_C/p0 * x_O2_interf ) ) #cathode overvoltage

        #Overall voltage
        V =  max(E_r - η_ohm - η_C, 0) #limit voltage to 0
    end
    
    return V
end

"""
    PEMFC_params

Structure containing the LT-PEMFC problem parameters.

!!! details "💾 Data fields"
    **Inputs:**
    - `R::Float64`: universal gas constant (J/(mol*K))
    - `F::Float64`: Faraday contant (C/mol) 
    - `p0::Float64`: reference pressure (Pa)
    - `T0::Float64`: reference temperature (T)
    - `n_drag::Float64`: number of water molecules dragged per proton in fully saturated Nafion
    - `M_m::Float64`: Nafion equivalent weight (kg/mol)
    - `ρ_dry::Float64`: dry density of Nafion (kg/m^3)
    - `x_ON::Float64`: mole fraction of oxygen in dry air
    - `α::Float64`: symmetry parameter
    - `ε::Float64`: porosity of electrodes
    - `τ::Float64`: diffusion tortuosity
    - `t_M::Float64`: membrane thickness (m)
    - `t_A::Float64`: anode thickness (m)
    - `t_C::Float64`: cathode thickness (m)
    - `α_star::Float64`: ratio of water flux to proton flux
    - `Iflux::Float64`: proton flux through membrane (mol/m^2/s)
"""
mutable struct PEMFC_params
    R :: Float64 
    F :: Float64 
    p0 :: Float64
    T0 :: Float64
    n_drag :: Float64
    M_m :: Float64
    ρ_dry :: Float64
    x_ON :: Float64
    α  :: Float64
    ε :: Float64
    τ :: Float64
    α_star :: Float64
    Iflux :: Float64
    PEMFC_params() = new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

"""
    PEMFC_inputs

Structure containing the LT-PEMFC and HT-PEMFC model inputs.

!!! details "💾 Data fields"
    **Inputs:**
    - `j::Float64`: current density (A/m^2)
    - `T::Float64`: fuel cell temperature (K)
    - `p_A::Float64`: anode pressure (Pa) 
    - `p_C::Float64`: cathode pressure (Pa)
    - `x_H2O_A::Float64`: mole fraction of water in anode inlet 
    - `x_H2O_C::Float64`: mole fraction of water in cathode inlet 
    - `λ_H2::Float64`: stoichiometric ratio of hydrogen
    - `λ_O2::Float64`: stoichiometric ratio of oxygen
    - `t_M::Float64`: membrane thickness (m)
    - `t_A::Float64`: anode thickness (m)
    - `t_C::Float64`: cathode thickness (m)
    - `type::String`: type of fuel cell
"""
mutable struct PEMFC_inputs
    j :: Float64
    T :: Float64 
    p_A :: Float64 
    p_C :: Float64
    x_H2O_A :: Float64 
    x_H2O_C :: Float64
    λ_H2 :: Float64
    λ_O2 :: Float64
    t_M :: Float64
    t_A :: Float64
    t_C :: Float64
    type :: String
    PEMFC_inputs() = new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "")
end

"""
    LT_PEMFC_voltage(u, α_guess::Float64 = 0.25)

A 1-D model of the voltage across a low-temperature PEM fuel cell with a Nafion membrane, based on the model
in Springer et al. (1991), which captures the effect of reactant depletion, multispecies diffusion and water transport
in the membrane.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `u::Struct`: structure of type `PEMFC_inputs` with inputs
    - `α_guess::Float64`: guess for ratio of water flux to proton flux; default is 0.25
 
    **Outputs:**
    - `V::Float64`: voltage across cell (V)
    - `α_star::Float64`: ratio of water flux to proton flux
"""
function LT_PEMFC_voltage(u, α_guess::Float64 = 0.25)
    #---------------------------------
    # Parameters
    #---------------------------------
    # Parameters from physics
    R = 8.314 # J/(mol*K), universal gas constant
    F = 96485.3329 # C/mol, Faraday contant
    n = 2 #number of electrons in oxidation reaction
    n_C = 4 #number of electrons in reduction reaction
    p0 = 101325 #Pa, reference pressure
    T0 = 298.15  #K, reference temperature
    n_drag = 2.5 #electro-osmotic drag coefficient in Nafion at saturation
    M_m = 1.0 #kg/mol, Nafion equivalent weight, typically 1-1.1 kg/mol
    ρ_dry = 1970 #kg/m^3, dry density of Nafion
    x_ON = 0.21 #mole fraction of oxygen in dry air

    # Assumed fuel cell parameters
    α = 0.3 #symmetry parameter (typically 0.1-0.5)
    ε = 0.4 #porosity of electrodes
    τ = 1.5 #diffusion tortuosity
    Aeff_Ratio = 1000 #ratio of effective catalyst area to geometry area
    
    #---------------------------------
    # Extract inputs
    #---------------------------------
    j = u.j
    T = u.T
    p_A = u.p_A
    p_C = u.p_C
    x_H2O_A = u.x_H2O_A
    x_H2O_C = u.x_H2O_C
    λ_H2 = u.λ_H2
    λ_O2 = u.λ_O2
    t_M = u.t_M
    t_A = u.t_A
    t_C = u.t_C

    Iflux = j / (n * F)
    p_SAT = water_sat_pressure(T) #Pa, water saturation pressure

    # Store parameters in stuct
    p = PEMFC_params()

    p.R = R
    p.F = F
    p.p0 = p0
    p.T0 = T0
    p.n_drag = n_drag
    p.M_m = M_m
    p.ρ_dry = ρ_dry
    p.x_ON = x_ON
    p.α = α
    p.ε = ε
    p.τ = τ
    p.Iflux = Iflux

    #---------------------------------
    # Use Roots.jl solver to find α_star
    #---------------------------------
    f(x) = water_balance(x, u, p)
    α_star = find_zero(f, α_guess) #Find root with Roots.jl
    p.α_star = α_star

    #---------------------------------
    # Reversible voltage
    #---------------------------------
    #Calculate other mole fractions
    x_H2_in = 1 - x_H2O_A #Hydrogen mole fraction at anode inlet
    x_O2_in = x_ON * (1 - x_H2O_C) #Oxygen mole fraction at cathode inlet

    #Reversible voltage
    p_O2 = x_O2_in * p_C #Pa, partial pressure of oxygen in cathode
    p_H2 = x_H2_in * p_A  #Pa, partial pressure of hydrogen in anode

    a_H2 = p_H2 / p0 #activity of hydrogen
    a_O2 = p_O2 / p0 #activity of oxygen

    Δs0 = -163.23 #J/mol/K, entropy change in chemical reaction with water vapor as product
    E0 = 1.229 #V,  reversible voltage at STP with water vapor as product

    #Nernst equation for reversible voltage
    E_r = E0 + Δs0 / (n * F) * (T - T0) - R * T / (n * F) * log(1 / (a_H2 * a_O2^0.5))
    #E_r = 1.1 #In Springer

    #---------------------------------
    # Find concentrations at all stations
    #---------------------------------
    # Concentrations at cathode/fluid interface
    x_H2O_4 = (x_H2O_C * λ_O2 + 2 * (1 + α_star) * (1 - x_H2O_C) * x_ON) / (λ_O2 + (2 * α_star + 1) * (1 - x_H2O_C) * x_ON)
    x_O2_4 = ((λ_O2 - 1) * (1 - x_H2O_C) * x_ON) / ( λ_O2 + (2 * α_star + 1) * (1 - x_H2O_C) * x_ON) 

    # Concentration at anode/fluid interface
    x_H2O_1 = (λ_H2 * x_H2O_A - α_star * (1 - x_H2O_A)) / (x_H2O_A - α_star * (1 - x_H2O_A) + λ_H2 - 1)

    #Diffusion coefficients
    D_H2H2O = binary_diffusion(T, p_A, ["H2", "H2O"])
    D_O2N2 =  binary_diffusion(T, p_C, ["O2", "N2"])
    D_N2H2O =  binary_diffusion(T, p_C, ["N2", "H2O"])
    D_O2H2O =  binary_diffusion(T, p_C, ["O2", "H2O"])

    Deff_H2H2O = porous_diffusion(D_H2H2O, ε, τ)
    Deff_O2N2 =  porous_diffusion(D_O2N2, ε, τ)
    Deff_N2H2O =  porous_diffusion(D_N2H2O, ε, τ)
    Deff_O2H2O =  porous_diffusion(D_O2H2O, ε, τ)

    # Find water concentration at cathode/membrane interface from analytical integration of diffusion from cathode side
    # Problem is of the form dx/dz = M * x + B
    M = [(-1*((1 + α_star) / Deff_O2H2O - 1 / Deff_O2N2)) (-1*(0.5 / Deff_O2H2O - 1 / Deff_O2N2)) ;
        (-1*(1 + α_star) / Deff_N2H2O + (1 + α_star) / Deff_O2H2O) (-1*(1 + α_star) / Deff_N2H2O + 0.5 / Deff_O2H2O)] * (R * T * Iflux / p_C)

    B = [ - 1 / Deff_O2N2;
        (1 + α_star) / Deff_N2H2O]* (R * T * Iflux / p_C)

    x_C4 = [x_O2_4; x_H2O_4]

    x_C3 = solve_diffusion_ODE(M, B, x_C4, t_C) #solve diffusion problem analytically with eigendecomposition
    x_O2_3 = x_C3[1]
    x_H2O_3 = x_C3[2]

    xliq = x_H2O_3 - p_SAT / p_C #mole fraction of water in liquid phase

    # Find water concentration at cathode/membrane interface from integration from anode side
    x_H2O_2 = (x_H2O_1 - α_star / (1 + α_star)) * exp((1 + α_star) * R * T * Iflux * t_A / (p_A * Deff_H2H2O)) + α_star / (1 + α_star)

    aw2 = x_H2O_2 * p_A / p_SAT #water activity at 2
    λ2 = λ_calc(aw2) #water content at 2, for boundary condition for ODE

    #Integrate ODE for λ
    dy_dx(y, p, x) = [dλ_dz_membrane(y[1], p[1], p[2]), 1/conductivity_Nafion(T, y[1])] #Vector with dλ/dz and dASR/dz
    
    y0 = [λ2, 0.0]
    xspan = (0.0, t_M)
    ODEparams = (u, p)
    
    prob = ODEProblem(dy_dx, y0, xspan, ODEparams, reltol = 1e-6)
    sol = solve(prob, Tsit5()) #Solve ODE problem using recommended solver

    #Extract solution vector
    y = sol.u
    y = hcat(y...) #Concatenate into array

    ASR = y[2,end] #Extract final resistance
    η_ohm = j * ASR # Ohmic overvoltage

    #---------------------------------
    # Find cathode overvoltage and cell voltage
    #---------------------------------

    if (α_star == 0) ||  (x_H2O_1 > 1) ||  (x_H2O_2 > 1) ||  (x_O2_3 < 0) #Set voltage to 0 in unphysical scenarios
        V = 0
    else
        j0 = cathode_j0(T, p_C * x_O2_3 / (1 - xliq) , Aeff_Ratio) #A/m^2, exchange current density
        η_C = R * T / (n_C * α * F) * log( j / ( j0 * p_C/p0 * x_O2_3 / (1 - xliq) ) ) #V, cathode overvoltage
        #η_C = R * T / (0.5 * F) * log( j / ( j0 * p_C/p0 * x_O2_3 / (1 - xliq) ) ) #V, cathode overvoltage, Springer

        V = E_r - η_C - η_ohm
        V = max(V, 0) #limit voltage to 0
    end

    return V, α_star
end

"""
    HT_PEMFC_voltage(u)

A 1-D model of the voltage across a high-temperature PEM fuel cell with a PBI membrane. The code is based on the LT-PEMFC
model by Springer et al. (1991), modified to eliminate water transport across the membrane and with conductivity values for
PBI.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `u::Struct`: structure of type `PEMFC_inputs` with inputs
    
    **Outputs:**
    - `V::Float64`: voltage across cell (V)
"""
function HT_PEMFC_voltage(u)
    #---------------------------------
    # Parameters
    #---------------------------------
    # Parameters from physics
    R = 8.314 # J/(mol*K), universal gas constant
    F = 96485.3329 # C/mol, Faraday contant
    n = 2 #number of electrons in oxidation reaction
    n_C = 4 #number of electrons in reduction reaction
    p0 = 101325 #Pa, reference pressure
    T0 = 298.15  #K, reference temperature
    x_ON = 0.21 #mole fraction of oxygen in dry air

    # Assumed fuel cell parameters
    α = 0.3 #symmetry parameter (typically 0.1-0.5)
    ε = 0.4 #porosity of electrodes
    τ = 1.5 #diffusion tortuosity
    Aeff_Ratio = 1000 #ratio of effective catalyst area to geometry area
    DL = 6 # membrane doping level

    #---------------------------------
    # Extract inputs
    #---------------------------------
    j = u.j
    T = u.T
    p_A = u.p_A
    p_C = u.p_C
    x_H2O_A = u.x_H2O_A
    x_H2O_C = u.x_H2O_C
    λ_H2 = u.λ_H2
    λ_O2 = u.λ_O2
    t_M = u.t_M
    t_A = u.t_A
    t_C = u.t_C

    Iflux = j / (n * F)
    p_SAT = water_sat_pressure(T) #Pa, water saturation pressure

    #---------------------------------
    # Reversible voltage
    #---------------------------------
    #Calculate other mole fractions
    x_H2_in = 1 - x_H2O_A #Hydrogen mole fraction at anode inlet
    x_O2_in = x_ON * (1 - x_H2O_C) #Oxygen mole fraction at cathode inlet

    #Reversible voltage
    p_O2 = x_O2_in * p_C #Pa, partial pressure of oxygen in cathode
    p_H2 = x_H2_in * p_A  #Pa, partial pressure of hydrogen in anode

    a_H2 = p_H2 / p0 #activity of hydrogen
    a_O2 = p_O2 / p0 #activity of oxygen

    Δs0 = -44.34 #J/mol/K, entropy change in chemical reaction with water vapor as product
    E0 = 1.1847 #V,  reversible voltage at STP with water vapor as product

    #Nernst equation for reversible voltage
    E_r = E0 + Δs0 / (n * F) * (T - T0) - R * T / (n * F) * log(1 / (a_H2 * a_O2^0.5))
    #---------------------------------
    # Find concentrations at all stations
    #---------------------------------
    # Concentrations at cathode/fluid interface
    x_H2O_4 = (x_H2O_C * λ_O2 + 2 * (1 - x_H2O_C) * x_ON) / (λ_O2 + (1 - x_H2O_C) * x_ON)
    x_O2_4 = ((λ_O2 - 1) * (1 - x_H2O_C) * x_ON) / ( λ_O2 + (1 - x_H2O_C) * x_ON) 

    # Concentration at anode/fluid interface
    x_H2O_1 = (λ_H2 * x_H2O_A) / (x_H2O_A + λ_H2 - 1)

    #Diffusion coefficients
    D_H2H2O = binary_diffusion(T, p_A, ["H2", "H2O"])
    D_O2N2 =  binary_diffusion(T, p_C, ["O2", "N2"])
    D_N2H2O =  binary_diffusion(T, p_C, ["N2", "H2O"])
    D_O2H2O =  binary_diffusion(T, p_C, ["O2", "H2O"])

    Deff_H2H2O = porous_diffusion(D_H2H2O, ε, τ)
    Deff_O2N2 =  porous_diffusion(D_O2N2, ε, τ)
    Deff_N2H2O =  porous_diffusion(D_N2H2O, ε, τ)
    Deff_O2H2O =  porous_diffusion(D_O2H2O, ε, τ)

    # Find water concentration at cathode/membrane interface from analytical integration of diffusion from cathode side
    # Problem is of the form dx/dz = M * x + B
    M = [(-1*(1 / Deff_O2H2O - 1 / Deff_O2N2)) (-1*(0.5 / Deff_O2H2O - 1 / Deff_O2N2)) ;
        (-1 / Deff_N2H2O + 1 / Deff_O2H2O) (-1 / Deff_N2H2O + 0.5 / Deff_O2H2O)] * (R * T * Iflux / p_C)

    B = [ - 1 / Deff_O2N2
        1 / Deff_N2H2O]* (R * T * Iflux / p_C)

    x_C4 = [x_O2_4; x_H2O_4]

    x_C3 = solve_diffusion_ODE(M, B, x_C4, t_C) #solve diffusion problem analytically with eigendecomposition
    x_O2_3 = x_C3[1]
    x_H2O_3 = x_C3[2]

    # Find water concentration at cathode/membrane interface from integration from anode side
    x_H2O_2 = x_H2O_1 * exp( R * T * Iflux * t_A / (p_A * Deff_H2H2O)) 

    # Find conductivity from RH and doping level
    RH_2 = x_H2O_2 * p_A / p_SAT
    RH_3 = x_H2O_3 * p_C / p_SAT
    RH_avg = 0.5 * (RH_2 + RH_3) 

    σ = conductivity_PBI(T, DL, RH_avg)
    ## Ohmic loss
    ASR = t_M / σ #Ohm m^2, area-specific resistance

    η_ohm = j * ASR # Ohmic overvoltage

    #---------------------------------
    # Find cathode overvoltage and cell voltage
    #---------------------------------
    
    if (x_O2_3 < 0) #Set voltage to 0 in unphysical scenarios
        V = 0
    else
        j0 = cathode_j0(T, x_O2_3 * p_C, Aeff_Ratio) #A/m^2, exchange current density
        η_C = R * T / (n_C * α * F) * log( j / ( j0 * p_C/p0 * x_O2_3) ) #V, cathode overvoltage
        
        V = E_r - η_C - η_ohm
        V = max(V, 0) #limit voltage to 0
    end
    return V
end

"""
    water_sat_pressure(T)

Function to calculate the saturation partial pressure of water. It uses different models for temperatures above
or below 100 degrees Celsius. 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `T::Float64`: gas temperature (K)
 
    **Outputs:**
    - `p_SAT::Float64`: saturation pressure (Pa)
"""
function water_sat_pressure(T)
    
    t = T - 273.15 #C, temperature in Celsius

    if t < 100 #Huang, J. (2018)
        p_SAT = exp(34.494 - 4924.99 / (t + 237.1)) / (t + 105) ^ 1.57 #Pa, water saturation pressure
    else #Jiao, K., and Li, X. (2010)
        p_SAT = 0.68737 * T^3 - 732.39 * T^2 +263390 * T -31919000 #Pa
    end
    return p_SAT
end

"""
    binary_diffusion(T, p, sps)

This model estimates the binary diffusion coefficient of a gas species in a water vapor, to be used inside the electrodes. It uses
the method of Slattery and Bird (1958) for low pressures.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `T::Float64`: gas temperature (K)
    - `p::Float64`: gas pressure (Pa)
    - `sps::Vec{String}`: gas species vector ("H2", "H2", "N2" or "O2")
 
    **Outputs:**
    - `D::Float64`: diffusion coefficient (m^2/s)
"""
function binary_diffusion(T, p, sps)
    a = 2.745e-4 #Use this parameters if H2O is not involved
    b = 1.823

    p_atm = p / 101325 

    Tc = zeros(Float64, 2)
    pc = zeros(Float64, 2)
    M = zeros(Float64, 2)

    for (i, sp) in enumerate(sps)
        if sp == "H2"
            #hydrogen gas
            Tc[i] = 33.3
            pc[i] = 12.8
            M[i] = 2.016
            
        elseif sp == "O2"
            #oxygen gas
            Tc[i] = 154.4
            pc[i] = 49.7
            M[i] = 31.999

        elseif sp == "N2"
            #nitrogen gas
            Tc[i] = 126.2
            pc[i] = 33.5
            M[i] = 28.013

        elseif sp == "H2O"
            a = 3.64e-4 #Redefine parameters for H2O
            b = 2.334

            Tc[i] = 647.3
            pc[i] = 217.5
            M[i]  = 18.015
            
        end
    end

    Tc_i = Tc[1]
    pc_i = pc[1]
    M_i = M[1]

    Tc_j = Tc[2]
    pc_j = pc[2]
    M_j = M[2]

    D = a * (T / (sqrt(Tc_i * Tc_j)))^b * (pc_i * pc_j)^(1/3) * (Tc_i * Tc_j)^(5/12) * sqrt(1/M_i + 1/M_j) / p_atm
    D = D * 1e-4 #m^2/s
    return D
end

"""
    porous_diffusion(D, ε, τ)

This model estimates the effective diffusion coefficient of a gas in a porous environment, such as a PEM electrode.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `D::Float64`: diffusion coefficient outside porous material (m^2/s)
    - `ε::Float64`: porosity of material
    - `τ::String`: tortuosity of material
 
    **Outputs:**
    - `Deff::Float64`: effective diffusion coefficient (m^2/s)
"""
function porous_diffusion(D, ε, τ)
    
    Deff = ε^τ * D
    return Deff
end

"""
    Nafion_diffusion(T, λ)

This function estimates the diffusion coefficient of water in Nafion.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `T::Float64`: fuel cell temperature (K)
    - `λ::Float64`: water content; ratio of water molecules to SO3- sites
 
    **Outputs:**
    - `D_λ::Float64`: water diffusion coefficient in Nafion (m^2/s)
"""
function Nafion_diffusion(T, λ)
    
    if λ < 16.8
        D_λ = exp(2416 * (1/303 - 1/T)) * (2.563 - 0.33 * λ + 0.0264 * λ^2 - 0.000671 * λ^3) * 1e-10 
    else #Extrapolate with constant slope
        D0 = exp(2416 * (1/303 - 1/T)) * 1.2885009279999995e-10
        slope = exp(2416 * (1/303 - 1/T)) * (-1.1109120000000084e-12) * 1e-10 
        D_λ = D0 + slope * (λ - 16.8)
    end
    
    return D_λ
end
"""
    cathode_j0(T, p, Aeff_ratio)

This function calculates the exchange current density of a PEM with a platinum catalyst. 
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `T::Float64`: fuel cell temperature (K)
    - `p::Float64`: reactant partial pressure (Pa)
    - `Aeff_ratio::Float64`: ratio of catalyst surface area to geometric cross-sectional area
 
    **Outputs:**
    - `j0::Float64`: exchange current density (A/m^2)
"""
function cathode_j0(T, p, Aeff_ratio)
    #Parameters from physics
    R = 8.314 # J/(mol*K), universal gas constant
    E_c = 66e3 #J/mol for ORR on Pt
    T0 = 298.15 #K
    p0 = 101250 #Pa
    γ = 0.75 #typicall 0.5-1.0
    j0ref = 1e-5 #A/m^2, oxygen on Pt catalyst at T0 and p0

    j0 = j0ref * Aeff_ratio * (p / p0)^γ * exp(-E_c / (R * T) * (1 - T/T0))
    
    return j0
end

"""
    conductivity_Nafion(T, λ)

This function calculates the conductivity of Nafion as a function of water content.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `T::Float64`: fuel cell temperature (K)
    - `λ::Float64`: water content; ratio of water molecules to SO3- sites
 
    **Outputs:**
    - `σ::Float64`: conductivity (Ohm m)^-1
"""
function conductivity_Nafion(T, λ)
    σ_30 = (0.005139 * λ .- 0.00326) * 1e2
    σ = exp(1268 * (1/303 - 1/T)) * σ_30
    return σ
end

"""
    conductivity_PBI(T, DL, RH)

This function calculates the conductivity of a PBI membrane as a function of temperature, doping level and humidity.
Model from K. Jiao and X. Li (2009). A Three-Dimensional Non-isothermal Model of High Temperature Proton Exchange 
Membrane Fuel Cells with Phosphoric Acid Doped Polybenzimidazole Membranes.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `T::Float64`: fuel cell temperature (K)
    - `DL::Float64`: phosphoric acid doping level
    - `RH::Float64`: average relative humidity across membrane
 
    **Outputs:**
    - `σ::Float64`: conductivity (Ohm m)^-1
"""
function conductivity_PBI(T, DL, RH)
    R = 8.314 # J/(mol*K), universal gas constant

    E_a = -619.6 * DL + 21750 #activation energy
    a = 168 * DL^3 - 6324 * DL^2 + 65750 * DL + 8460
    if (T <= 413.15)
        b = 1 + (0.01704 * T - 4.767) * RH
    elseif (T <= 453.15)
        b = 1 + (0.1432 * T - 56.89) * RH
    else
        b = 1 + (0.7 * T -309.2) * RH
    end
    σ = a * b / T * exp(-E_a / (R * T)) # 1/(Ohm m)
    return σ
end

"""
    λ_calc(a)

This function calculates the water content at the Nafion/electrode interface based on water activity.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `a::Float64`: water activity
    
    **Outputs:**
    - `λ::Float64`: water content; ratio of water molecules to SO3- sites
"""
function λ_calc(a)
    if (a < 1)
        λ = 0.043 + 17.81 * a - 39.85 * a^2 + 36 * a^3
    elseif (a >= 1)
        λ = 14 + 1.4 * (a - 1)
    else #Extrapolated for robustness if a < 1
        λ = 14 * a
    end
    return λ
end

"""
    dλ_dz_membrane(λ, u, p)

This function evaluates the derivative in space of the water content in the membrane.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `λ::Float64`: water content; ratio of water molecules to SO3- sites
    - `u::Struct`: structure of type PEMFC_inputs with inputs 
    - `p::Struct`: structure of type PEMFC_params with parameters 
    
    **Outputs:**
    - `dλ_dz::Float64`: derivative of λ in space
"""
function dλ_dz_membrane(λ, u, p)
    #Extract inputs
    λ = λ[1]

    #Extract parameters and inputs
    Iflux = p.Iflux
    T = u.T

    α_star = p.α_star
    M_m = p.M_m
    ρ_dry = p.ρ_dry
    n_drag = p.n_drag

    #Diffusion coefficient
    D_λ = Nafion_diffusion(T, λ)

    #Oxygen concentration rate of change
    dλ_dz = (n_drag * λ / 11 - α_star) * Iflux * M_m / (ρ_dry * D_λ)

    return dλ_dz
end

"""
    solve_diffusion_ODE(M, B, x0, d)

This function uses eigendecomposition to solve a problem of the form dx/dz = M * x + B
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `M::Matrix{Float64}`: matrix with coefficients
    - `B::Vector{Float64}`: vector with right-hand parameters
    - `x0::Vector{Float64}`: vector with the boundary conditions
    - `d::Float64`: distance at which to evaluate x, z = d

    **Outputs:**
    - `x_end::Vector{Float64}`: vector with values of x at z = d
"""
function solve_diffusion_ODE(M, B, x0, d)
    #Find eigenvalues and eigenvectors
    U = eigvecs(M)
    Uinv = inv(U)
    eigsM = eigvals(M)

    Btrans = Uinv * B
    z_0 = Uinv * x0 #Boundary conditions

    #Integration constants
    C = z_0 + Btrans ./ eigsM

    #Solve problem analytically
    z_end = C .* exp.(eigsM * d) - Btrans ./ eigsM

    x_end = U * z_end
    return x_end
end

"""
    water_balance(α_star, u, p)

This function calculates the difference between the water content at 3 from integration from the anode or cathode sides.
This residual should be 0 if α_star is the correct one.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `α_star::Float64`: ratio of water flux to proton flux
    - `u::Struct`: structure of type `PEMFC_inputs` with inputs 
    - `p::Struct`: structure of type `PEMFC_params` with parameters 

    **Outputs:**
    - `x_end::Vector{Float64}`: vector with values of x at z = d
"""
function water_balance(α_star, u, p)
    p.α_star = α_star #store current α_star as parameter
    
    #Extract inputs
    j = u.j
    T = u.T
    p_A = u.p_A
    p_C = u.p_C
    x_H2O_A = u.x_H2O_A
    x_H2O_C = u.x_H2O_C
    λ_H2 = u.λ_H2
    λ_O2 = u.λ_O2
    t_M = u.t_M
    t_A = u.t_A
    t_C = u.t_C
    
    #Extract parameters
    R = p.R
    x_ON = p.x_ON 
    ε = p.ε
    τ = p.τ
    Iflux = p.Iflux

    p_SAT = water_sat_pressure(T)

    # Concentrations at cathode/fluid interface
    x_H2O_4 = (x_H2O_C * λ_O2 + 2 * (1 + α_star) * (1 - x_H2O_C) * x_ON) / (λ_O2 + (2 * α_star + 1) * (1 - x_H2O_C) * x_ON)
    x_O2_4 = ((λ_O2 - 1) * (1 - x_H2O_C) * x_ON) / ( λ_O2 + (2 * α_star + 1) * (1 - x_H2O_C) * x_ON) 

    # Concentrations at anode/fluid interface
    x_H2O_1 = (λ_H2 * x_H2O_A - α_star * (1 - x_H2O_A)) / (x_H2O_A - α_star * (1 - x_H2O_A) + λ_H2 - 1)

    #Diffusion coefficients
    D_H2H2O = binary_diffusion(T, p_A, ["H2", "H2O"])
    D_O2N2 =  binary_diffusion(T, p_C, ["O2", "N2"])
    D_N2H2O =  binary_diffusion(T, p_C, ["N2", "H2O"])
    D_O2H2O =  binary_diffusion(T, p_C, ["O2", "H2O"])

    Deff_H2H2O = porous_diffusion(D_H2H2O, ε, τ)
    Deff_O2N2 =  porous_diffusion(D_O2N2, ε, τ)
    Deff_N2H2O =  porous_diffusion(D_N2H2O, ε, τ)
    Deff_O2H2O =  porous_diffusion(D_O2H2O, ε, τ)

    # Find water content at cathode/membrane interface from analytical integration of diffusion from cathode side
    # Problem is of the form dx/dz = M * x + B
    M = [(-1*((1 + α_star) / Deff_O2H2O - 1 / Deff_O2N2)) (-1*(0.5 / Deff_O2H2O - 1 / Deff_O2N2)) ;
        (-1*(1 + α_star) / Deff_N2H2O + (1 + α_star) / Deff_O2H2O) (-1*(1 + α_star) / Deff_N2H2O + 0.5 / Deff_O2H2O)] * (R * T * Iflux / p_C)

    B = [ - 1 / Deff_O2N2
        (1 + α_star) / Deff_N2H2O]* (R * T * Iflux / p_C)

    x_C4 = [x_O2_4; x_H2O_4]

    x_C3 = solve_diffusion_ODE(M, B, x_C4, t_C) #solve diffusion problem analytically with eigendecomposition
    x_H2O_3 = x_C3[2]
    aw3 = x_H2O_3 * p_C / p_SAT
    λ3_prime = λ_calc(aw3) #water content at 3

    # Find water content at cathode/membrane interface from integration from anode side
    x_H2O_2 = (x_H2O_1 - α_star / (1 + α_star)) * exp((1 + α_star) * R * T * Iflux * t_A / (p_A * Deff_H2H2O)) + α_star / (1 + α_star)

    aw2 = x_H2O_2 * p_A / p_SAT
    λ2 = λ_calc(aw2)
   
    dy_dx(y, p, x) = [dλ_dz_membrane(y[1], p[1], p[2])]
    
    y0 = [λ2]
    xspan = (0.0, t_M)
    ODEparams = (u, p)
    
    prob = ODEProblem(dy_dx, y0, xspan, ODEparams, reltol = 1e-6)
    sol = solve(prob, Tsit5()) #Solve ODE problem using recommended solver

    #Extract solution vector
    y = sol.u
    y = hcat(y...) #Concatenate into array

    λ3 = y[1,end]

    res = λ3 - λ3_prime #difference between the two methods

    return res
end

"""
    PEMsize(P_des, V_des, u)

Designs the fuel cell stack for the design point conditions.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `P_des::Float64`: design stack output power, ideally maximum power in mission (W)
    - `V_des::Float64`: design stack voltage (V)
    - `u::Struct`: structure of type `PEMFC_inputs` with inputs 
 
    **Outputs:**
    - `n_cells::Float64`: number of cells in stack
    - `A_cell::Float64`: cell surface area (m^2)
    - `Q::Float64`: waste power produced by the fuel cell at design point (W)
"""
function PEMsize(P_des, V_des, u)
    #Extract inputs
    type = u.type
    j = u.j
    
    #Find heating voltage
    if type == "LT-PEMFC"
        V_heat = 1.482
    elseif type == "HT-PEMFC"
        V_heat = 1.254
    end
    #Calculate power density of a cell
    P2A = P2Acalc(u, j)
    V_cell = P2A / j

    #Size stack
    n_cells = V_des / V_cell #number of cells in series

    A_cell = P_des /(n_cells * P2A) #total area of a cell, including possible cells in parallel

    #Calculate heat to be dissipated
    Q = n_cells * A_cell * (V_heat - V_cell) * j

    return n_cells, A_cell, Q
end #PEMsize

"""
    PEMoper(P_stack, n_cells, A_cell, u)

Evaluates fuel cell stack performance in off-design conditions.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `P_stack::Float64`: stack output power (W)
    - `n_cells::Float64`: number of cells in stack
    - `A_cell::Float64`: cell surface area (m^2)
    - `u::Struct`: structure of type `PEMFC_inputs` with inputs 
 
    **Outputs:**
    - `V_stack::Float64`: stack voltage (V)
    - `Q::Float64`: waste power produced by the fuel cell (W)
"""
function PEMoper(P_stack, n_cells, A_cell, u)
    #Extract inputs
    type = u.type
    
    #Find heating voltage
    if type == "LT-PEMFC"
        V_heat = 1.482

    elseif type == "HT-PEMFC"
        V_heat = 1.254
    end

    P2A = P_stack / (n_cells * A_cell) #Power density of a cell

    f(x) = P2A - P2Acalc(u, x) #Residual function; it should be 0 if x = j
    j_g = 1 #A/m^2, very low current density as a guess
    j = find_zero(f, j_g) #Find root with Roots.jl. Careful! There are two roots, must use the smallest one

    #Find stack voltage
    V_cell = P2A / j
    V_stack = V_cell * n_cells

    #Calculate heat to be dissipated
    Q = n_cells * A_cell * (V_heat - V_cell) * j

    return V_stack, Q
end #PEMoper

"""
    P2Acalc(u, j)

Calculates the power density of a fuel cell.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `u::Struct`: structure of type `PEMFC_inputs` with inputs 
    - `j::Float64`: current density (A/m^2)
 
    **Outputs:**
    - `P2A::Float64`: power density (W/m^2)
"""
function P2Acalc(u, j)
    u.j = j

    if u.type == "LT-PEMFC"
        V_cell, _ = LT_PEMFC_voltage(u)

    elseif u.type == "HT-PEMFC"
        V_cell = HT_PEMFC_voltage(u)
    end
    
    P2A = j * V_cell
    return P2A
end

"""
    PEMstackweight(gee, u, A, n_cells, fouter)

Calculates the weight of a stack of PEM fuel cells.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gee::Float64`: gravitational acceleration (m/s^2)
    - `u::Struct`: structure of type `PEMFC_inputs` with inputs 
    - `n_cells::Float64`: number of cells in stack
    - `A_cell::Float64`: cell surface area (m^2)
    - `fouter::Float64`: ratio of stack structural mass (inc. bipolar plates) to membrane and electrode mass
 
    **Outputs:**
    - `W_stack::Float64`: weight of FC stack (N)
"""
function PEMstackweight(gee, u, n_cells, A_cell, fouter)
    #Extract inputs
    type = u.type
    t_M = u.t_M
    t_A = u.t_A
    t_C = u.t_C

    ρ_elect = 1.9e3 #kg/m^3, density of porous graphite plates

    if type == "LT-PEMFC"
        ρ_M = 1970 #kg/m^3, density of dry Nafion
    elseif type == "HT-PEMFC"
        ρ_M = 1300 #kg/m^3, density of PBI
    end

    W_stack = gee * n_cells * A_cell * ((t_A + t_C) * ρ_elect + t_M * ρ_M) * (1 + fouter)

    return W_stack
end #PEMstackweight
