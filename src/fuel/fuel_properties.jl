"""
    cryo_fuel_properties(fuel::String, p::Float64)

Calculates the saturation temperature and density of a cryogenic fuel inside a tank, using fits to data from NIST.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuel::String`: name of the fuel. 
    - `p::String`: tank pressure (Pa).
    **Outputs:**
    - `Tsat::Float64`: saturation temperature (K) 
    - `ρl::Float64`: liquid density (kg/m^3).
    - `ρg::Float64`: gas density (kg/m^3).
    - `hv::Float64`: enthalpy of vaporization (J/kg).
"""
function cryo_fuel_properties(fuel::String, p::Float64)
    x = p / 101325 #pressure in atm

    #Fits to NIST data at p increments of 0.1 atm from 0.1 to 10 atm
    if uppercase(fuel) == "LH2"
        #Saturation temperature (K)
        Tsat = -3.35772E-04*x^6 + 1.14637E-02*x^5 - 1.54871E-01*x^4 + 1.05803E+00*x^3 - 3.93770E+00*x^2 + 9.09205E+00*x + 1.42913E+01
        
        #Liquid density (kg/m^3)
        ρl = 2.58354E-04*x^6 - 8.90930E-03*x^5 + 1.21248E-01*x^4 - 8.37637E-01*x^3 + 3.14780E+00*x^2 - 8.42843E+00*x + 7.68629E+01
        #Liquid enthalpy (kJ/kg)
        hl = -2.258517E-03*x^6 + 7.767483E-02*x^5 - 1.055739E+00*x^4 + 7.280899E+00*x^3 - 2.739412E+01*x^2 + 7.379258E+01*x - 5.277765E+01
    
        #Gas density (kg/m^3)
        ρg = 5.61446E-03*x^3 - 4.19141E-02*x^2 + 1.28761E+00*x + 6.91405E-02
        #Gas enthalpy (kJ/kg)
        hg = -3.31730E-03*x^6 + 1.12617E-01*x^5 - 1.51499E+00*x^4 + 1.02766E+01*x^3 - 3.79360E+01*x^2 + 7.33546E+01*x + 4.04343E+02
    elseif uppercase(fuel) == "CH4"
        Tsat = -9.38765E-04*x^6 + 3.25637E-02*x^5 - 4.49376E-01*x^4 + 3.16285E+00*x^3 - 1.22945E+01*x^2 + 3.00925E+01*x + 9.09645E+01
        
        ρl = 1.23474E-03*x^6 - 4.28482E-02*x^5 + 5.91776E-01*x^4 - 4.17205E+00*x^3 + 1.62724E+01*x^2 - 4.14826E+01*x + 4.51398E+02
        hl = -3.14958E-03*x^6 + 1.09246E-01*x^5 - 1.50774E+00*x^4 + 1.06164E+01*x^3 - 4.13089E+01*x^2 + 1.02732E+02*x - 7.11737E+01
    
        ρg = 1.66091E-03*x^3 - 2.82413E-02*x^2 + 1.69557E+00*x + 1.33041E-01
        hg = -1.88367E-03*x^6 + 6.52861E-02*x^5 - 8.99795E-01*x^4 + 6.31818E+00*x^3 - 2.44324E+01*x^2 + 5.58852E+01*x + 4.73575E+02
    end

    #hl and hg are in kJ/kg
    hv = (hg - hl) * 1e3 #J/kg, Enthalpy of vaporization is enthalpy difference at saturation

    return Tsat, ρl, ρg, hv
end
