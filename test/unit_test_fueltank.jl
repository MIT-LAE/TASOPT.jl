#Sample parameters
rhoFuel = 70.0
deltap = 2e5
Rfuse = 2.5
dRfuse = 0.3
hconvgas = 0.0
Tfuel = 20.0
t_cond = [0.15,0.15,0.15]
material_insul = ["rohacell41S", "rohacell41S", "polyurethane"]
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

    outputs_size_check = (164717.55636677786, 0.49077274862420084, 12.790944804923033, 0.004247366632687734, 145.623998835008, 100000.0, 107.03363914373088, 0.004916887075394796, 0.0049105639181476314, 1.9092272513757993, 2141.366879896994, 21900.169398210935, 33298.07257717196, [10956.41232435591, 12098.332463606823, 10243.327789209226], 14.695255169223437, 64717.55636677787)
    for i in 1:length(outputs_size)
        @test outputs_size[i] ≈ outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(gee, rhoFuel,
                  fstring, ffadd, deltap,
                  Rfuse, dRfuse, wfb, nfweb,
                  sigskin, material_insul, rhoskin,
                  Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

    outputs_mech_check = (164717.55636677786, 12.790944804923033, 0.004916887075394796, 1.9092272513757993, 145.623998835008, 64717.55636677787, 100000.0, 33298.07257717196, 0.0049105639181476314, 2141.366879896994, 21900.169398210935, [10956.41232435591, 12098.332463606823, 10243.327789209226], [15.832773779668218, 19.254692964997925, 23.016493195949316, 27.117221510580325], 14.695255169223437)
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