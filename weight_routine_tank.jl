v_f = 67.3 #specific volume as saturated liquid
v_g = 2.5 #specific volume as saturated vapor
x = 0.03 #quality (suggested by energies)

rhoFuel = 1 / ((1 / v_f) + x * ((1 / v_g) - (1 / v_f)))

gee = 9.81

deltap = 51325#101325 #121590 #based on 6000 ft cabin pressure, 2 bar LH2 pressure (energies)

Rfuse = 3.1#2.07 #radius of fuselage of A320

dRfuse = 0 #extended portion of fuselage at bottom (take as zero due to cylindrical tank)

hconvgas = 0 #no purged gas in foam insulated tank

h_LH2 = 750 #initial estimate based on elsevier pub

Tfuel = 20 #20 K

Tair = 293 #cabin heated temperature

h_v = 447000 #Engineering toolbox, Hydrogen - Thermophysical properties

t_cond = [0.1, 1.524e-5, 0.1, 1.524e-5, 1.57e-2] #assumed

k = [20e-3, 20e-3, 20e-3, 20e-3, 20e-3] #foam thickness

hconvair = 25 #from sciencedirect.com https://www.sciencedirect.com/topics/engineering/convection-heat-transfer-coefficient

time_flight = 18000 #seconds, 5 hr flight

fstring = 0.34

ffadd = 0.2

wfb = 0

nfweb = 0

sigskin = 172.4e6 #AL 2209 Brewer stress for operating conditions (290e6 ultimate operatoin)

rho_insul = [35.24, 14764, 35.24, 14764, 83]

rhoskin = 2825 #Al 2209 Brewer

#Wfuel = 8000
Wfuel = 2000 #8000 * 9.81

threshold_percent = 0.1

mode = 1

clearance_fuse = 0.13 #m , Brewer suggests 5 in clearance

#comments
#gee is gravitational acceleration, 9.81 m/s^2
#rhoFuel is calculated based on 3% quality (energies paper)
#deltap is based on pressure difference across tank wall,



#remove when running (only for test)
m_boiloff = 0.001 * Wfuel
thickness_insul = sum(t_cond)
lshell = 2.38
r_tank = 1.722
