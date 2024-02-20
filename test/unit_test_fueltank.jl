using TASOPT
using Test

#Sample parameters
rhoFuel = 70.0
deltap = 2e5
Rfuse = 2.5
dRfuse = 0.3
hconvgas = 0.0
Tfuel = 20.0
t_cond = [0.15,0.15,0.15]
material_insul = ["rohacell31", "rohacell31", "polyurethane"]
iinsuldes = [1,2,3]
time_flight = 7*3600.0
fstring = 0.1
ffadd = 0.1
wfb = 0.0
nfweb = 1.0
sigskin = 172.4e6 
rhoskin = 2825.0
Wfuel = 1e5
threshold_percent = 0.15
clearance_fuse = 0.1
AR = 2.0
ifuel = 40
qfac = 1.3
z = 11e3
Mair = 0.8
xftank = 15.0

m_boiloff = threshold_percent * Wfuel /(100 * gee)

@testset "Fuselage tank" begin
    outputs_size = TASOPT.structures.tanksize(gee, rhoFuel, deltap,
                        Rfuse, dRfuse, hconvgas, Tfuel, z, Mair, xftank,
                        t_cond, time_flight, fstring,ffadd,
                        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
                        iinsuldes, ifuel, qfac)

    outputs_size_check = (163057.2735182767, 0.4917517453742062, 12.827138884391198, 0.15, 145.8424348332605, 100150.0, 15.290519877675841, 0.004914365837305676, 0.004908045922391445, 1.9082482546257937, 2138.0744749265928, 21940.159670980705, 31447.703173276, [10055.793235955925, 11105.819793547458, 10286.090143772617], 14.730472773179686, 62907.27351827668)
    for i in 1:length(outputs_size)
        @test outputs_size[i] == outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(gee, rhoFuel,
                  fstring, ffadd, deltap,
                  Rfuse, dRfuse, wfb, nfweb,
                  sigskin, material_insul, rhoskin,
                  Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

    outputs_mech_check = (163057.2735182767, 12.827138884391198, 0.004914365837305676, 1.9082482546257937, 145.8424348332605, 62907.27351827668, 100150.0, 31447.703173276, 0.004908045922391445, 2138.0744749265928, 21940.159670980705, [10055.793235955925, 11105.819793547458, 10286.090143772617], [15.816540761806268, 19.243955325165583, 23.01260555444891, 27.121533644908638], 14.730472773179686)
    for i in 1:length(outputs_mech)
        @test outputs_mech[i] == outputs_mech_check[i]
    end

end