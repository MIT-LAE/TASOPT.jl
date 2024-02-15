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
threshold_percent = 0.1
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

    outputs_size_check = (183216.4510868988, 0.6983467946259261, 16.586193153964686, 0.1, 145.76962283384302, 100100.0, 10.193679918450561, 0.00438231574909625, 0.00437668005497447, 1.7016532053740738, 1516.112316029177, 22690.32935349458, 52249.386304235275, [16003.900534151206, 18479.68715987249, 17765.798610211583], 18.283464043589664, 83116.4510868988)

    for i in 1:length(outputs_size)
        @test outputs_size[i] == outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(gee, rhoFuel,
                  fstring, ffadd, deltap,
                  Rfuse, dRfuse, wfb, nfweb,
                  sigskin, material_insul, rhoskin,
                  Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

    outputs_mech_check = (183216.4510868988, 16.586193153964686, 0.00438231574909625, 1.7016532053740738, 145.76962283384302, 83116.4510868988, 100100.0, 52249.386304235275, 0.00437668005497447, 1516.112316029177, 22690.32935349458, [16003.900534151206, 18479.68715987249, 17765.798610211583], [12.577197176930863, 17.04566248903523, 22.201170720110152, 28.04148242560481], 18.283464043589664)
    for i in 1:length(outputs_mech)
        @test outputs_mech[i] == outputs_mech_check[i]
    end

    l_cyl = outputs_mech_check[2]
    l_tank = outputs_mech_check[14]
    r_tank = outputs_mech_check[4]
    Shead = outputs_mech_check[13]
    outputs_thermal = TASOPT.structures.tankWthermal(l_cyl, l_tank, r_tank, Shead, material_insul,
                      hconvgas, 
                      t_cond,
                      Tfuel, z, Mair, xftank,
                      time_flight, ifuel, qfac)

    outputs_thermal_check = (71.35575943521175, 0.0028315777553655458)

    for i in 1:length(outputs_thermal)
        @test outputs_thermal[i] == outputs_thermal_check[i]
    end

end