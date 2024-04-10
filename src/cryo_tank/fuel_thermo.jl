"""
    gas_properties(species::String, p::Float64)

This function returns the thermodynamic properties of a saturated vapor.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: pressure (Pa)
    
    **Outputs:**
    - `Tsat::Float64`: temperature (K)
    - `ρ::Float64`: density (kg/m^3)
    - `ρ_p::Float64`: derivative of density with pressure (kg/m^3/Pa)
    - `h::Float64`: specific enthalpy (J/kg)
    - `u::Float64`: specific internal energy (J/kg)
    - `u_p::Float64`: derivative of internal energy with pressure (J/kg/Pa)
"""
function gas_properties(species::String, p::Float64)
    x = p / p_atm #pressure in atm

    #Fits to NIST data at p increments of 0.1 atm from 0.1 to 10 atm
    if uppercase(species) == "H2" || uppercase(species) == "LH2"
        #Saturation temperature (K)
        Tsat = -3.35772E-04*x^6 + 1.14637E-02*x^5 - 1.54871E-01*x^4 + 1.05803E+00*x^3 - 3.93770E+00*x^2 + 9.09205E+00*x + 1.42913E+01
        
        #Gas density (kg/m^3)
        ρ = 5.61446E-03*x^3 - 4.19141E-02*x^2 + 1.28761E+00*x + 6.91405E-02
        #Density-pressure derivative (kg/m^3/atm)
        ρ_p = 3*5.61446E-03*x^2 - 2*4.19141E-02*x + 1.28761E+00
        #Gas enthalpy (kJ/kg)
        h = -3.31730E-03*x^6 + 1.12617E-01*x^5 - 1.51499E+00*x^4 + 1.02766E+01*x^3 - 3.79360E+01*x^2 + 7.33546E+01*x + 4.04343E+02
        #Gas internal energy (kJ/kg)
        u = -2.01181E-03*x^6 + 6.82454E-02*x^5 - 9.18066E-01*x^4 + 6.22852E+00*x^3 - 2.30505E+01*x^2 + 4.44191E+01*x + 3.45858E+02
        #Internal energy-pressure derivative (kJ/kg/atm)
        u_p = -6*2.01181E-03*x^5 + 5*6.82454E-02*x^4 - 4*9.18066E-01*x^3 + 3*6.22852E+00*x^2 - 2*2.30505E+01*x + 4.44191E+01
    elseif uppercase(species) == "CH4"
        Tsat = -9.38765E-04*x^6 + 3.25637E-02*x^5 - 4.49376E-01*x^4 + 3.16285E+00*x^3 - 1.22945E+01*x^2 + 3.00925E+01*x + 9.09645E+01
        
        ρ = 1.66091E-03*x^3 - 2.82413E-02*x^2 + 1.69557E+00*x + 1.33041E-01
        ρ_p = 3*1.66091E-03*x^2 - 2*2.82413E-02*x + 1.69557E+00

        h = -1.88367E-03*x^6 + 6.52861E-02*x^5 - 8.99795E-01*x^4 + 6.31818E+00*x^3 - 2.44324E+01*x^2 + 5.58852E+01*x + 4.73575E+02
        u = -1.41615E-03*x^6 + 4.90837E-02*x^5 - 6.76548E-01*x^4 + 4.75172E+00*x^3 - 1.83908E+01*x^2 + 4.24831E+01*x + 4.26591E+02
        u_p = -6*1.41615E-03*x^5 + 5*4.90837E-02*x^4 - 4*6.76548E-01*x^3 + 3*4.75172E+00*x^2 - 2*1.83908E+01*x + 4.24831E+01
    end

    h = h * 1e3 #J/kg
    u = u * 1e3 #J/kg
    ρ_p = ρ_p / p_atm #kg/m^3/Pa
    u_p = u_p * 1e3 / p_atm #J/kg/Pa
    return Tsat, ρ, ρ_p, h, u, u_p
end

"""
    liquid_properties(species::String, p::Float64)

This function returns the thermodynamic properties of a saturated liquid.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `species::String`: Species name
    - `p::Float64`: pressure (Pa)
    
    **Outputs:**
    - `Tsat::Float64`: temperature (K)
    - `ρ::Float64`: density (kg/m^3)
    - `ρ_p::Float64`: derivative of density with pressure (kg/m^3/Pa)
    - `h::Float64`: specific enthalpy (J/kg)
    - `u::Float64`: specific internal energy (J/kg)
    - `u_p::Float64`: derivative of internal energy with pressure (J/kg/Pa)
"""
function liquid_properties(species::String, p::Float64)
    x = p / p_atm #pressure in atm

    #Fits to NIST data at p increments of 0.1 atm from 0.1 to 10 atm
    if uppercase(species) == "H2" || uppercase(species) == "LH2"
        #Saturation temperature (K)
        Tsat = -3.35772E-04*x^6 + 1.14637E-02*x^5 - 1.54871E-01*x^4 + 1.05803E+00*x^3 - 3.93770E+00*x^2 + 9.09205E+00*x + 1.42913E+01
        
        #Liquid density (kg/m^3)
        ρ = 2.58354E-04*x^6 - 8.90930E-03*x^5 + 1.21248E-01*x^4 - 8.37637E-01*x^3 + 3.14780E+00*x^2 - 8.42843E+00*x + 7.68629E+01
        #Density-pressure derivative (kg/m^3/atm)
        ρ_p = 6*2.58354E-04*x^5 - 5*8.90930E-03*x^4 + 4*1.21248E-01*x^3 - 3*8.37637E-01*x^2 + 2*3.14780E+00*x - 8.42843E+00
        #Liquid enthalpy (kJ/kg)
        h = -2.258517E-03*x^6 + 7.767483E-02*x^5 - 1.055739E+00*x^4 + 7.280899E+00*x^3 - 2.739412E+01*x^2 + 7.379258E+01*x - 5.277765E+01
        #Liquid internal energy (kJ/kg)
        u = -2.268218E-03*x^6 + 7.793319E-02*x^5 - 1.058739E+00*x^4 + 7.297785E+00*x^3 - 2.749429E+01*x^2 + 7.244553E+01*x - 5.277403E+01
        #Internal energy-pressure derivative (kJ/kg/atm)
        u_p = -6*2.268218E-03*x^5 + 5*7.793319E-02*x^4 - 4*1.058739E+00*x^3 + 3*7.297785E+00*x^2 - 2*2.749429E+01*x + 7.244553E+01
    elseif uppercase(species) == "CH4"
        Tsat = -9.38765E-04*x^6 + 3.25637E-02*x^5 - 4.49376E-01*x^4 + 3.16285E+00*x^3 - 1.22945E+01*x^2 + 3.00925E+01*x + 9.09645E+01
        
        ρ = 1.23474E-03*x^6 - 4.28482E-02*x^5 + 5.91776E-01*x^4 - 4.17205E+00*x^3 + 1.62724E+01*x^2 - 4.14826E+01*x + 4.51398E+02
        ρ_p = 6*1.23474E-03*x^5 - 5*4.28482E-02*x^4 + 4*5.91776E-01*x^3 - 3*4.17205E+00*x^2 + 2*1.62724E+01*x - 4.14826E+01

        h = -3.14958E-03*x^6 + 1.09246E-01*x^5 - 1.50774E+00*x^4 + 1.06164E+01*x^3 - 4.13089E+01*x^2 + 1.02732E+02*x - 7.11737E+01
        u = -3.14796E-03*x^6 + 1.09204E-01*x^5 - 1.50735E+00*x^4 + 1.06151E+01*x^3 - 4.13130E+01*x^2 + 1.02493E+02*x - 7.11704E+01
        u_p = -6*3.14796E-03*x^5 + 5*1.09204E-01*x^4 - 4*1.50735E+00*x^3 + 3*1.06151E+01*x^2 - 2*4.13130E+01*x + 1.02493E+02
    end

    h = h * 1e3 #J/kg
    u = u * 1e3 #J/kg
    ρ_p = ρ_p / p_atm #kg/m^3/Pa
    u_p = u_p * 1e3 / p_atm #J/kg/Pa

    return Tsat, ρ, ρ_p, h, u, u_p
end
