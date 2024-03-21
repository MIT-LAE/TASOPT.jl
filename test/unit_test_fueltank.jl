#Sample parameters
ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
time_flight = 7*3600.0

fuse_tank = TASOPT.fuselage_tank()
fuse_tank.Rfuse = 2.5
fuse_tank.dRfuse = 0.3
fuse_tank.wfb = 0.0
fuse_tank.nfweb = 1.0
fuse_tank.clearance_fuse = 0.1

fuse_tank.ptank = 2e5
fuse_tank.ARtank = 2.0
fuse_tank.UTSinner = 476e6 
fuse_tank.rhoinner = 2840.0

fuse_tank.t_insul = [0.15,0.15,0.15]
fuse_tank.material_insul = ["rohacell41S", "rohacell41S", "polyurethane27"]
fuse_tank.size_insulation = true
fuse_tank.iinsuldes = [1,2,3]
fuse_tank.boiloff_rate = 0.15

fuse_tank.ullage_frac = 0.1
fuse_tank.qfac = 1.3
fuse_tank.theta_inner = 1.0
fuse_tank.ftankadd = 0.1
fuse_tank.ew = 0.9

fuse_tank.UTSouter = 476e6 
fuse_tank.rhoouter = 2840.0
fuse_tank.Eouter = 73e9
fuse_tank.poissouter = 0.3
fuse_tank.theta_outer = [1.0, 2.0]
fuse_tank.Ninterm = 1.0

fuse_tank.rhofuel = 70.0
fuse_tank.Tfuel = 20.0
fuse_tank.Wfuelintank = 1e5

@testset "Fuselage tank" begin
    @testset "Foam insulation" begin
        outputs_size = TASOPT.structures.tanksize!(fuse_tank, z, Mair, xftank,
                                        time_flight,
                                        ifuel)
        outputs_size_check = (155809.35798473458, 0.48977003536494124, 12.754661771033133, 0.004247366632687734, 145.623998835008, 100000.0, 107.03363914373088, 0.0035645268979941384, 0.003561204259200333, 1.9102299646350587, 1565.0564653871652, 15934.428537303356, 33161.58141606815, [10913.366525915251, 12048.59580855156, 10199.619081601335], 14.661327208770198, 55809.35798473459)
        for i in 1:length(outputs_size)
            @test outputs_size[i] ≈ outputs_size_check[i]
        end

        outputs_mech = TASOPT.structures.size_inner_tank(fuse_tank, fuse_tank.t_insul)

        outputs_mech_check = (155809.35798473458, 12.754661771033133, 0.0035645268979941384, 1.9102299646350587, 145.623998835008, 55809.35798473459, 100000.0, 33161.58141606815, 0.003561204259200333, 1565.0564653871652, 15934.428537303356, 1524.346321619075, [10913.366525915251, 12048.59580855156, 10199.619081601335], 192.00106755661065, [15.844936883719953, 19.260709211746004, 23.014988363554576, 27.106823359071385], 14.661327208770198)
        for i in 1:length(outputs_mech)
            @test outputs_mech[i] ≈ outputs_mech_check[i]
        end

        l_cyl = outputs_mech_check[2]
        l_tank = outputs_mech_check[16]
        r_tank = outputs_mech_check[4]
        Shead = outputs_mech_check[15]
        outputs_thermal = TASOPT.structures.tankWthermal(fuse_tank, z, Mair, xftank, time_flight, ifuel)

        outputs_thermal_check = (107.03363914390675, 0.004247366632694712)

        for i in 1:length(outputs_thermal)
            @test outputs_thermal[i] ≈ outputs_thermal_check[i]
        end
    end
    
    fuse_tank.t_insul = [0.065]
    fuse_tank.material_insul = ["vacuum"]
    fuse_tank.size_insulation = false

    @testset "Vacuum insulation" begin

        outputs_vac_size = TASOPT.structures.tanksize!(fuse_tank, z, Mair, xftank,
                                            time_flight,
                                            ifuel)
        outputs_vac_size_check = (203190.88704051627, 0.065, 7.833183975356766, 0.0020612680724211444, 145.623998835008, 100000.0, 51.94395542501284, 0.00435715618585557, 0.004353094705443698, 2.335, 2858.4704423225353, 14495.185531216492, 0.0, [0.0], 10.163826819170911, 103190.88704051627)
    
        for i in 1:length(outputs_vac_size)
            @test outputs_vac_size[i] ≈ outputs_vac_size_check[i]
        end

        outputs_vac_mech = TASOPT.structures.size_inner_tank(fuse_tank, fuse_tank.t_insul)
        outputs_vac_mech_check = (124405.4750588923, 7.833183975356766, 0.00435715618585557, 2.335, 145.623998835008, 24405.475058892298, 100000.0, 0.0, 0.004353094705443698, 2858.4704423225353, 14495.185531216492, 1974.6690922223409, [0.0], 166.54673585196167, [23.675162947566548, 25.28189535258484], 10.163826819170911)
        
        for i in 1:length(outputs_vac_mech)
            @test outputs_vac_mech[i] ≈ outputs_vac_mech_check[i]
        end
        Winnertank = outputs_vac_mech_check[1]
        l_cyl = outputs_vac_mech_check[2]

        fuse_tank.Ninterm = 1.0
        Ninterm = TASOPT.structures.optimize_outer_tank(fuse_tank, Winnertank, l_cyl)
        Ninterm_check = 17.3857421875
        @test Ninterm ≈ Ninterm_check

        outputs_vac_outer = TASOPT.structures.size_outer_tank(fuse_tank, Winnertank, l_cyl, Ninterm_check)
        outputs_vac_outer_check = (78785.41198162397, 29309.017167216, 12884.34023240256, 16545.404169455214, 168.1077334143945, 24.993050952321795, 118.12163150975091, 0.008906036506379233, 0.018503571497814696)

        for i in 1:length(outputs_vac_outer)
            @test outputs_vac_outer[i] ≈ outputs_vac_outer_check[i]
        end
    end

    
    
    @testset "Stringer sizing" begin
        θ = 85 * pi /180
        outputs_bendM = TASOPT.structures.stiffeners_bendingM(θ)
        outputs_bendM_check = [1.1726278985544076, 0.09494956354378614]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM[i] ≈ outputs_bendM_check[i]
        end
    end

end