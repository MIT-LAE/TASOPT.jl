include("../src/structures/structures.jl")
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

#Convective cooling
hconvair = 15.0 #In W/(m^2 K)

time_flight = 5 * 3600.0
sigskin = 172.4e6 #AL 2219 Brewer / energies stress for operating conditions (290e6 ultimate operation)
rho_insul = [35.24, 14764, 35.24, 14764, 83] #energies
rhoskintank = 2825.0 #Al 2219 / energies
max_boiloff = 0.1 #%/h, maximum percentage of full fuel tank boiling off per hour 
ARtank = 2.0 #TODO: why? doesn't this overconstrain?
clearance_fuse = 0.10 #TODO: why?
rhofuel = 70.0
ptank = 2.0 #atm #TODO: why? maybe write as input
ftankstiff = 0.1
ftankadd = 0.1

cargotank = false #TODO: figure out why this is here

if cargotank
    Wfmaintank = parg[igWfuel] * 2 / 3
    Wfcargotank = parg[igWfuel] * 1 / 3
else
    Wfmaintank = 48512 * 9.81 / 2 #Main tank carries 1/nftanks of the fuel
    Wfcargotank = 0.0
end

gee = 9.81

Wtank_total, thickness_insul, ltank, mdot_boiloff, Vfuel, Wfuel_tot,
m_boiloff, tskin, t_head, Rtank, Whead, Wcyl,
Winsul_sum, Winsul, l_tank, Wtank = tanksize(gee, rhofuel, ptank * 101325.0,
    3.0, 0.38, hconvgas, h_LH2, Tfuel, Tair,
    h_v, t_cond, k, hconvair, time_flight, ftankstiff, ftankadd,
    0.0, 1.0, sigskin, rho_insul, rhoskintank,
    Wfmaintank, max_boiloff, clearance_fuse, ARtank)