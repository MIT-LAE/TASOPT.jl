using TASOPT
using PyPlot

rhoFuel = 70.0
deltap = 2e5
Rfuse = 2.5
dRfuse = 0.3
hconvgas = 0.0
Tfuel = 20.0
t_cond = [0.45]
material_insul = ["rohacell31"]
iinsuldes = [1]
time_flight = 7*3600.0
fstring = 0.1
ffadd = 0.1
wfb = 0.0
nfweb = 1.0
sigskin = 172.4e6 
rhoskin = 2825.0
Wfuel = 1e5
threshold_percent = 0.15 #use just to oversize tank to account for boiloff
clearance_fuse = 0.1
AR = 2.0
ifuel = 40
qfac = 1.3
z = 11e3
Mair = 0.8
xftank = 15.0

m_boiloff = threshold_percent * Wfuel /(100 * gee) *time_flight/3600

ts = LinRange(0.05,2, 50)
boffs = zeros(length(ts))

for (i,t) in enumerate(ts)
    t_cond = [t]

    outputs_mech = TASOPT.structures.tankWmech(gee, rhoFuel,
                    fstring, ffadd, deltap,
                    Rfuse, dRfuse, wfb, nfweb,
                    sigskin, material_insul, rhoskin,
                    Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

    l_cyl = outputs_mech[2]
    l_tank = outputs_mech[14]
    r_tank = outputs_mech[4]
    Shead = outputs_mech[13]
    (_, mdot_boiloff) = TASOPT.structures.tankWthermal(l_cyl, l_tank, r_tank, Shead, material_insul,
                    hconvgas, 
                    t_cond,
                    Tfuel, z, Mair, xftank,
                    time_flight, ifuel, qfac)
    boffs[i] = mdot_boiloff * gee / Wfuel * 3600
end


p1 = figure()
plt.rc("font", size = 14)
plot(ts*100, boffs*100)
plt.xlabel("Insulation thickness (cm)")
plt.ylabel("Boiloff rate (%/h)")
plt.axis([0,200,0, 1])
plt.grid()

plt.yticks(LinRange(0,1,11))
show()
savefig("Plot.svg")
close(p1)