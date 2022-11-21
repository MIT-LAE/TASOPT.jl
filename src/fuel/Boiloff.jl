"""
Boiloff calculates boil-off rate of LH2 given thermal resistances

Inputs:

NOTE: Everything is in SI units.
-Array R_thermal consisting of thermal resistances (K/W) of each layer
-Ambient temperature T_ambient (K)
-Fuel temperature T_fuel (K)

Outputs:
- m_boiloff (kg) is the boiloff LH2 mass
- mdot_boiloff (kg/s) is the boiloff rate of LH2
"""
function Boiloff(R_thermal, T_ambient, T_fuel)

      N = size(R_thermal) #Number of layers in insulation

#--- Heat flux and resistances

      deltaT = T_ambient - T_fuel  #Overall temperature drop between ambient and LH2
      R_insul = sum(R_thermal)  #Total equivalent thermal resistance of MLI
      q = deltaT / R_insul  #Heat flux from ambient to LH2
      mdot_boiloff = q / h_e  #Boil-off rate equals the heat flux divided by heat of vaporization
      m_boiloff = mdot_boiloff * time_flight #Boil-off mass calculation
  
return  mdot_boiloff
end
