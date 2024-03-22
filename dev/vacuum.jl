Tcold = 78
Thot = 300
S_outer = pi*2.1^2
S_inner = pi*1.5^2

p_vacuum = 1e-3 #Assumed vacuum pressure

mu = 18.47*10^-6 #mu at 300K
Rgas = 287.05  # specific gas constant

gamma = 1.4
a_outer = 0.9 
a_inner = 1.0
ε = 0.04    # highly polished aluminum

#Emissivity factor
Fe = 1 / (1/0.05 + S_inner/S_outer * (1/0.1 - 1))

#Find heat transfer coeff for radiation and corresponding resistance
hrad = σ_SB * Fe * ((Tcold^2) + (Thot^2)) * (Tcold + Thot) #Radiative heat transfer coefficient; Eq. (2.28) in https://ahtt.mit.edu/
R_rad = 1/(hrad * S_inner)  # radiative resistance

#Find resistance due to convection by residual air
Fa = 1 / (1/a_inner + S_inner/S_outer * (1/a_outer - 1)) #accomodation factor
G = (gamma + 1)/(gamma - 1) * sqrt(Rgas / (8*pi * Thot)) * Fa
R_conv = 1 / (G * p_vacuum * S_inner)  # convective resistance due to imperfect vacuum
 
R_eq = R_conv * R_rad / (R_conv + R_rad) #parallel addition of resistance