"""
    cryo_fuel_properties(fuel::String, p::Float64)

Calculates the saturation temperature and density of a cryogenic fuel inside a tank, using fits to data from NIST.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fuel::String`: name of the fuel. 
    - `p::String`: tank pressure (Pa).
    **Outputs:**
    - `Tfuel::Float64`: fuel evaporation temperature (K) 
    - `ρfuel::String`: density of fuel (kg/m^3).
"""
function cryo_fuel_properties(fuel::String, p::Float64)
    x = p / 101325 #pressure in atm

    #Fits to NIST data at p increments of 0.1 atm from 0.1 to 10 atm
    if fuel == "LH2" 
        Tfuel = 3.35772E-04*x^6 + 1.14637E-02*x^5 - 1.54871E-01*x^4 + 1.05803E+00*x^3 - 3.93770E+00*x^2 + 9.09205E+00*x + 1.42913E+01
        ρfuel = 2.58354E-04*x^6 - 8.90930E-03*x^5 + 1.21248E-01*x^4 - 8.37637E-01*x^3 + 3.14780E+00*x^2 - 8.42843E+00*x + 7.68629E+01
    elseif fuel == "CH4"
        Tfuel = -9.38765E-04*x^6 + 3.25637E-02*x^5 - 4.49376E-01*x^4 + 3.16285E+00*x^3 - 1.22945E+01*x^2 + 3.00925E+01*x + 9.09645E+01
        ρfuel = 1.23474E-03*x^6 - 4.28482E-02*x^5 + 5.91776E-01*x^4 - 4.17205E+00*x^3 + 1.62724E+01*x^2 - 4.14826E+01*x + 4.51398E+02
    end
    return Tfuel, ρfuel
end
