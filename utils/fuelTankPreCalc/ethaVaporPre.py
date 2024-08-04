import numpy as np
#from DDBST, use Antoine Equation Parameter to calculate the vapor pressure at different temperature
##Constant for Ethanol
A = 8.20417
B = 1642.89
C = 230.3
##
T_test = 300 #218.85K temperature at 35000ft based on ISA
T_test = T_test - 273.15 #Celcius
P_varp = 10**(A-B/(C+T_test)) # mmHg
P_varp = P_varp*101325/760 #Pa
print("P_varp is : ",P_varp)