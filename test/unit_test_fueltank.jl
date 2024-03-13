#Sample parameters
ρfuel = 70.0
Rfuse = 2.5
dRfuse = 0.3
Tfuel = 20.0
time_flight = 7*3600.0
wfb = 0.0
nfweb = 1.0
Wfuel = 1e5
ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
hconvgas = 0.0

fuse_tank = TASOPT.fuselage_tank()
fuse_tank.ptank = 2e5
fuse_tank.t_insul = [0.15,0.15,0.15]
fuse_tank.material_insul = ["rohacell41S", "rohacell41S", "polyurethane27"]
fuse_tank.iinsuldes = [1,2,3]
fuse_tank.ftankstiff = 0.1
fuse_tank.ftankadd = 0.1
fuse_tank.sigskin = 476e6 
fuse_tank.rhoskintank = 2825.0
fuse_tank.max_boiloff = 0.15
fuse_tank.clearance_fuse = 0.1
fuse_tank.ARtank = 2.0
fuse_tank.qfac = 1.3

@testset "Fuselage tank" begin
    outputs_size = TASOPT.structures.tanksize(fuse_tank, ρfuel, 
                                    Rfuse, dRfuse, Tfuel, z, Mair, xftank,
                                    time_flight,
                                    wfb, nfweb, Wfuel,
                                    ifuel)
    outputs_size_check = (155918.19958080447, 0.4897700354090134, 12.754661771709657, 0.004247366632687734, 145.623998835008, 100000.0, 107.03363914373088, 0.0035645268979118986, 0.00356120425911817, 1.9102299645909864, 1556.7903219763107, 15850.267823324659, 33161.58142007174, [10913.366527138433, 12048.595810012961, 10199.619082920344], 14.661327209402732, 55918.19958080448)
    for i in 1:length(outputs_size)
        @test outputs_size[i] ≈ outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(fuse_tank, fuse_tank.t_insul, ρfuel,
                                                Rfuse, dRfuse, wfb, nfweb,
                                                Wfuel)

    outputs_mech_check = (155918.19958080447, 12.754661771709657, 0.0035645268979118986, 1.9102299645909864, 145.623998835008, 55918.19958080448, 100000.0, 33161.58142007174, 0.00356120425911817, 1556.7903219763107, 15850.267823324659, [10913.366527138433, 12048.595810012961, 10199.619082920344], [15.844936882988812, 19.260709211262586, 23.014988363379675, 27.106823359265597], 14.661327209402732)
    for i in 1:length(outputs_mech)
        @test outputs_mech[i] ≈ outputs_mech_check[i]
    end

    l_cyl = outputs_mech_check[2]
    l_tank = outputs_mech_check[14]
    r_tank = outputs_mech_check[4]
    Shead = outputs_mech_check[13]
    outputs_thermal = TASOPT.structures.tankWthermal(l_cyl, l_tank, r_tank, Shead, fuse_tank.material_insul,
                      hconvgas, 
                      fuse_tank.t_insul,
                      Tfuel, z, Mair, xftank,
                      time_flight, ifuel, fuse_tank.qfac)

    outputs_thermal_check = (107.0336391563202, 0.00424736663318731)

    for i in 1:length(outputs_thermal)
        @test outputs_thermal[i] ≈ outputs_thermal_check[i]
    end

end