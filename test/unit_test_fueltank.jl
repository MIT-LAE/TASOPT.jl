using TASOPT
using Test

#Sample parameters
rhoFuel = 70.0
deltap = 2e5
Rfuse = 2.5
dRfuse = 0.3
hconvgas = 0.0
Tfuel = 20.0
Tair = 220.0
t_cond = [0.175,0.175,0.175]
material_insul = ["rohacell31", "rohacell31", "polyurethane"]
iinsuldes = [1,2,3]
hconvair = 100.0
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

m_boiloff = threshold_percent * Wfuel /(100 * gee)

@testset "Fuselage tank" begin
    outputs_size = TASOPT.structures.tanksize(gee, rhoFuel, deltap,
                        Rfuse, dRfuse, hconvgas, Tfuel, Tair,
                        t_cond, hconvair, time_flight, fstring,ffadd,
                        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
                        iinsuldes, ifuel)

    outputs_size_check = ((183207.07180144673, 0.6982656552257915, 16.58444966535786, 0.1, 145.76962283384302, 100100.0, 10.193679918450561, 0.004382524709694073, 0.0043768887468472436, 1.7017343447742084, 1516.3292033055952, 22690.050286332436, 52239.82136991436, [16001.250923210973, 18476.293116316385, 17762.277330386998], 18.281801485422374, 83107.07180144671))

    for i in 1:length(outputs_size)
        @test outputs_size[i] == outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(gee, rhoFuel,
                  fstring, ffadd, deltap,
                  Rfuse, dRfuse, wfb, nfweb,
                  sigskin, material_insul, rhoskin,
                  Wfuel, m_boiloff, t_cond, clearance_fuse, AR)
    outputs_mech_check = (183207.07180144673, 16.58444966535786, 0.004382524709694073, 1.7017343447742084, 145.76962283384302, 83107.07180144671, 100100.0, 52239.82136991436, 0.0043768887468472436, 1516.3292033055952, 22690.050286332436, [16001.250923210973, 18476.293116316385, 17762.277330386998], [12.578396634683275, 17.04649944176178, 22.201485965411504, 28.0411173527444], 18.281801485422374)
    for i in 1:length(outputs_mech)
        @test outputs_mech[i] == outputs_mech_check[i]
    end

    l_cyl = outputs_mech_check[2]
    l_tank = outputs_mech_check[14]
    r_tank = outputs_mech_check[4]
    Shead = outputs_mech_check[13]
    outputs_thermal = TASOPT.structures.tankWthermal(l_cyl, l_tank, r_tank, Shead, material_insul,
                      hconvgas,  hconvair, 
                      t_cond,
                      Tfuel, Tair, 
                      time_flight, ifuel)

    outputs_thermal_check = (71.35951543110049, 0.002831726802821448)

    for i in 1:length(outputs_thermal)
        @test outputs_thermal[i] == outputs_thermal_check[i]
    end

end