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

    outputs_size_check = (163056.6070900031, 0.4917436948407957, 12.827014595064798, 0.15, 145.8424348332605, 100150.0, 15.290519877675841, 0.0049143865700726356, 0.004908066628495896, 1.9082563051592043, 2138.1015354192245, 21940.127779697446, 31447.010069359996, [10055.587937824652, 11105.57437039221, 10285.847761143132], 14.730356513653929, 62906.607090003075)
    for i in 1:length(outputs_size)
        @test outputs_size[i] == outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(gee, rhoFuel,
                  fstring, ffadd, deltap,
                  Rfuse, dRfuse, wfb, nfweb,
                  sigskin, material_insul, rhoskin,
                  Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

    outputs_mech_check = (163056.6070900031, 12.827014595064798, 0.0049143865700726356, 1.9082563051592043, 145.8424348332605, 62906.607090003075, 100150.0, 31447.010069359996, 0.004908066628495896, 2138.1015354192245, 21940.127779697446, [10055.587937824652, 11105.57437039221, 10285.847761143132], [15.816674215991934, 19.244043611109657, 23.01263752182894, 27.12149818327768], 14.730356513653929)
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

    outputs_thermal_check = (107.03363915621976, 0.004247366633183324)

    for i in 1:length(outputs_thermal)
        @test outputs_thermal[i] == outputs_thermal_check[i]
    end

end