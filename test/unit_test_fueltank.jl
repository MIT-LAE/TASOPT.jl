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

m_boiloff = threshold_percent * Wfuel /(100 * gee) *time_flight/3600

@testset "Fuselage tank" begin
    outputs_size = TASOPT.structures.tanksize(gee, rhoFuel, deltap,
                        Rfuse, dRfuse, hconvgas, Tfuel, z, Mair, xftank,
                        t_cond, time_flight, fstring,ffadd,
                        wfb, nfweb, sigskin, material_insul, rhoskin, Wfuel, threshold_percent, clearance_fuse, AR, 
                        iinsuldes, ifuel, qfac)

    outputs_size_check = (164963.6665454459, 0.4976717044899518, 13.046408143820344, 0.15, 147.15305082277558, 101050.0, 107.03363914373088, 0.004899119998738213, 0.004892819690097861, 1.9023282955100482, 2118.2372960911366, 22180.263341517446, 32213.581025006253, [10290.04791815106, 11376.776229594276, 10546.75687726092], 14.943837319331653, 63913.666545445914)
    for i in 1:length(outputs_size)
        @test outputs_size[i] ≈ outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(gee, rhoFuel,
                  fstring, ffadd, deltap,
                  Rfuse, dRfuse, wfb, nfweb,
                  sigskin, material_insul, rhoskin,
                  Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

    outputs_mech_check = (164963.6665454459, 13.046408143820344, 0.004899119998738213, 1.9023282955100482, 147.15305082277558, 63913.666545445914, 101050.0, 32213.581025006253, 0.004892819690097861, 2118.2372960911366, 22180.263341517446, [10290.04791815106, 11376.776229594276, 10546.75687726092], [15.71855766640118, 19.17908964140118, 22.98910579953093, 27.147618780900714], 14.943837319331653)
    for i in 1:length(outputs_mech)
        @test outputs_mech[i] ≈ outputs_mech_check[i]
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

    outputs_thermal_check = (107.0336391563202, 0.00424736663318731)

    for i in 1:length(outputs_thermal)
        @test outputs_thermal[i] ≈ outputs_thermal_check[i]
    end

end