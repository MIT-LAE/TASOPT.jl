#include("../src/structures/structures.jl")
include("../src/structures/tanksize.jl")
include("../src/structures/tankWmech.jl")
include("../src/structures/tankWthermal.jl")
using Roots

hconvgas = 0.0
Tfuel = 20.0
Tair = 288.0 #Heated cabin temp
h_LH2 = 210.0 #W/m^2/K, heat transfer coefficient of LH2 #TODO: replace by function
h_v = 447000.0 #TODO: replace by function
t_cond = [0.05, 1.524e-5, 0.05, 1.524e-5, 1.57e-2] #assumed from energies -- Total thickness is 11.6 cm ~ Brewer's Rigid closed cell foam tank type A pg194 
k = ones(length(t_cond)) .* 5.0e-3 #foam conductivities

H = 6
ifuel = 40

#Convective cooling
hconvair = 15.0 * (288-20) #In W/(m^2 K)

h_liq, hconvair, h_v, Tair = tank_heat_coeffs(11e3, 0.8, 11, 120, 6, 10)