include("tanksize.jl")
include("tankWmech.jl")
include("tankWthermal.jl")
include("hydrogen.jl")


y = 0.1
ρmean = ρmix(y, 1.5)
rhoFuel = ρmean
println("ρ fuel = $rhoFuel")

# gee = 9.81

deltap = 200_000.0 #150000#150000#101325  #based on 0.5 bar cabin pressure, 2 or 1.5 bar LH2 pressure (Verstraete / Brewer)

Rfuse = 2.07 #radius of fuselage of A320

dRfuse = 0.0 #extended portion of fuselage at bottom (take as zero due to cylindrical tank)

hconvgas = 0.0 #no purged gas in foam insulated tank

h_LH2 = 550.0 #based on elsevier pub film boiling heat transfer properties of liq hydrogen 

Tfuel = 20.0 #20 K

Tair = 288.0 #cabin heated temperature

h_v = 447000.0 #Engineering toolbox, Hydrogen - Thermophysical properties

t_cond = [0.1, 1.524e-5, 0.1, 1.524e-5, 1.57e-2] #assumed from energies

k = [5e-3, 5e-3, 5e-3, 5e-3, 5e-3] #foam conductivities

#http://www.hysafe.org/download/1196/BRHS_Chap1_V1p2.pdf
hconvair = 15.0 #from sciencedirect.com https://www.sciencedirect.com/topics/engineering/convection-heat-transfer-coefficient

time_flight = 18000.0 #seconds, 5 hr flight

fstring = 0.34

ffadd = 0.2

wfb = 0.0

nfweb = 0.0

sigskin = 172.4e6 #AL 2219 Brewer / energies stress for operating conditions (290e6 ultimate operatoin)

rho_insul = [35.24, 14764, 35.24, 14764, 83] #energies

rhoskin = 2825.0 #Al 2219 / energies
AR = 2 #ellipsoidal head aspect ratio
#Wfuel = 8000
use_factor = 0.8
Wfuel = 14000 * 9.81 / use_factor #8000 * 9.81
# Wfuel = 9.3e4
#threshold_percent = 0.1:0.1:1

#r_tank = 1.722

mode = 1
Winsulation = zeros(5)
Wtot_tank = zeros(5)
gravimetric_eff = zeros(5)
boiloff = zeros(5)
Wmetal = zeros(5)
Wins_to_tot = zeros(5)
Wmet_to_tot = zeros(5)
boiloff_percentage = zeros(5)

clearance_fuse = 0.13 #m , Brewer suggests 5 in clearance
#n = 1
threshold_percent = 0.2:-0.05:0.1
for n in 1:1:length(threshold_percent)

#deltap = 50000:50000:250000
Wtank_total, thickness_insul, lshell, mdot_boiloff, Vfuel, Wfuel_tot,
 m_boiloff, tskin, t_head, Rtank, Whead, Wcyl,
  Winsul_sum, Winsul, l_tank, Wtank = tanksize(gee, rhoFuel, deltap,
                      Rfuse, dRfuse, hconvgas, h_LH2, Tfuel, Tair,
                      h_v, t_cond, k, hconvair, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, rho_insul, rhoskin, Wfuel, threshold_percent[n], clearance_fuse, AR)
#println("\nlength of tank=", l_tank)
#println("\nn =", threshold_percent[n])
grav_index = Wfuel_tot / Wtank_total
percent_boiloff = 100 * m_boiloff*gee / (Wfuel * time_flight / 3600)
boiloff_percentage[n] = percent_boiloff
# println("grav index =  ", grav_index)
# println("percent boiloff =  ", percent_boiloff)
Winsulation[n] = Winsul_sum
Wtot_tank[n] = Wtank_total
gravimetric_eff[n] = grav_index
boiloff[n] = m_boiloff
Wmetal[n] = Wtank
Wins_to_tot[n] = Winsulation[n]/Wtot_tank[n]
Wmet_to_tot[n] = Wmetal[n]/Wtot_tank[n]
println("Boil off = $(threshold_percent[n]), η = $grav_index, Wtank  = $(Wtank_total) Total insulation weight = ",Winsul_sum)
println(Vfuel)
#n = n + 1
end
#comments
#gee is gravitational acceleration, 9.81 m/s^2
#rhoFuel is calculated based on 3% quality (energies paper)
#deltap is based on pressure difference across tank wall,



#remove when running (only for test)
#m_boiloff = 0.001 * Wfuel
#thickness_insul = sum(t_cond)
#lshell = 2.38
#r_tank = 1.722
