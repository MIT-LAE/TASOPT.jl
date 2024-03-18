#Sample parameters
ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
hconvgas = 0.0
time_flight = 7*3600.0

fuse_tank = TASOPT.fuselage_tank()
fuse_tank.ptank = 2e5
fuse_tank.t_insul = [0.15,0.15,0.15]
fuse_tank.material_insul = ["rohacell41S", "rohacell41S", "polyurethane27"]
fuse_tank.iinsuldes = [1,2,3]
fuse_tank.ftankstiff = 0.1
fuse_tank.ftankadd = 0.1
fuse_tank.UTSinner = 476e6 
fuse_tank.rhoinner = 2825.0
fuse_tank.boiloff_rate = 0.15
fuse_tank.clearance_fuse = 0.1
fuse_tank.ARtank = 2.0
fuse_tank.ullage_frac = 0.1
fuse_tank.qfac = 1.3
fuse_tank.ew = 0.9

fuse_tank.rhofuel = 70.0
fuse_tank.Rfuse = 2.5
fuse_tank.dRfuse = 0.3
fuse_tank.Tfuel = 20.0

fuse_tank.wfb = 0.0
fuse_tank.nfweb = 1.0
fuse_tank.Wfuelintank = 1e5

@testset "Fuselage tank" begin
    outputs_size = TASOPT.structures.tanksize(fuse_tank, z, Mair, xftank,
                                    time_flight,
                                    ifuel)
    outputs_size_check = (155918.19958080447, 0.4897700354090134, 12.754661771709657, 0.004247366632687734, 145.623998835008, 100000.0, 107.03363914373088, 0.0035645268979118986, 0.00356120425911817, 1.9102299645909864, 1556.7903219763107, 15850.267823324659, 33161.58142007174, [10913.366527138433, 12048.595810012961, 10199.619082920344], 14.661327209402732, 55918.19958080448)
    for i in 1:length(outputs_size)
        @test outputs_size[i] ≈ outputs_size_check[i]
    end

    outputs_mech = TASOPT.structures.tankWmech(fuse_tank, fuse_tank.t_insul)

    outputs_mech_check = (155918.19958080447, 12.754661771709657, 0.0035645268979118986, 1.9102299645909864, 145.623998835008, 55918.19958080448, 100000.0, 33161.58142007174, 0.00356120425911817, 1556.7903219763107, 15850.267823324659, [10913.366527138433, 12048.595810012961, 10199.619082920344], 176.2269243480927, [15.844936882988812, 19.260709211262586, 23.014988363379675, 27.106823359265597], 14.661327209402732)
    for i in 1:length(outputs_mech)
        @test outputs_mech[i] ≈ outputs_mech_check[i]
    end

    l_cyl = outputs_mech_check[2]
    l_tank = outputs_mech_check[15]
    r_tank = outputs_mech_check[4]
    Shead = outputs_mech_check[14]
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