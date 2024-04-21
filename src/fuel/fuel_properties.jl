export cryo_fuel_properties

"""
    cryo_fuel_properties(fuel::String, p::Float64)

Calculates the temperature and density of a cryogenic fuel inside a tank, using fits to data from NIST.
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

    #Fits to NIST data at p = [0.5,1,1.5,2,2.5,3,4,5,10] atm
    if fuel == "LH2" 
        Tfuel = 20.3839 * x^0.182022 
        ρfuel = -2.88304E-02*x^3 + 4.93319E-01*x^2 - 4.61334E+00*x + 7.51570E+01
    elseif fuel == "CH4"
        Tfuel = 111.69030 * x^0.12093
        ρfuel = -1.49198E-01*x^3 + 2.85956E+00*x^2 - 2.20464E+01*x + 4.42724E+02
    end
    return Tfuel, ρfuel
end
