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
        outputs_size_check = (155771.1276307066, 0.48977003536494124, 12.754661771033133, 0.004247366632687734, 145.623998835008, 100000.0, 107.03363914373088, 0.0035645268979941384, 0.003561204259200333, 1.9102299646350587, 1565.0564653871652, 15934.428537303356, 33161.58141606815, [10913.366525915251, 12048.59580855156, 10199.619081601335], 14.661327208770198, 55771.1276307066)
        for i in 1:length(outputs_size)
            @test outputs_size[i] ≈ outputs_size_check[i]
        end

        outputs_mech = TASOPT.structures.size_inner_tank(fuse_tank, fuse_tank.t_insul)

        outputs_mech_check = (155771.1276307066, 12.754661771033133, 0.0035645268979941384, 1.9102299646350587, 145.623998835008, 55771.1276307066, 100000.0, 33161.58141606815, 0.003561204259200333, 1565.0564653871652, 15934.428537303356, 1489.5914543209028, [10913.366525915251, 12048.59580855156, 10199.619081601335], 192.00106755661065, [15.844936883719953, 19.260709211746004, 23.014988363554576, 27.106823359071385], 14.661327208770198)
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
        outputs_vac_size_check = (202582.47739770001, 0.065, 7.833183975356766, 0.0020612680724211444, 145.623998835008, 100000.0, 51.94395542501284, 0.00435715618585557, 0.004353094705443698, 2.335, 2858.4704423225353, 14495.185531216492, 0.0, [0.0], 10.163826819170911, 102582.47739770003) 
        
        for i in 1:length(outputs_vac_size)
            @test outputs_vac_size[i] ≈ outputs_vac_size_check[i]
        end

        outputs_vac_mech = TASOPT.structures.size_inner_tank(fuse_tank, fuse_tank.t_insul)
        outputs_vac_mech_check = (124336.97641373744, 7.833183975356766, 0.00435715618585557, 2.335, 145.623998835008, 24336.976413737444, 100000.0, 0.0, 0.004353094705443698, 2858.4704423225353, 14495.185531216492, 1912.3975966270182, [0.0], 166.54673585196167, [23.675162947566548, 25.28189535258484], 10.163826819170911)

        for i in 1:length(outputs_vac_mech)
            @test outputs_vac_mech[i] ≈ outputs_vac_mech_check[i]
        end
        Winnertank = outputs_vac_mech_check[1]
        l_cyl = outputs_vac_mech_check[2]

        fuse_tank.Ninterm = 1.0
        Ninterm = TASOPT.structures.optimize_outer_tank(fuse_tank, Winnertank, l_cyl)
        Ninterm_check = 17.59375
        @test Ninterm ≈ Ninterm_check

        outputs_vac_outer = TASOPT.structures.size_outer_tank(fuse_tank, Winnertank, l_cyl, Ninterm_check)
        outputs_vac_outer_check = ((78245.50098396259, 29149.49739622203, 12884.34023240256, 16214.09576075701, 168.1077334143945, 24.993050952321795, 118.12163150975091, 0.00885756374812002, 0.018503571497814696))
        for i in 1:length(outputs_vac_outer)
            @test outputs_vac_outer[i] ≈ outputs_vac_outer_check[i]
        end

    end
    
    @testset "Stringer sizing" begin
        θ = 85 * pi /180
        outputs_bendM = TASOPT.structures.stiffeners_bendingM(θ)
        outputs_bendM_check = [1.1693705988362009, 0.09494234539184043 ]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM[i] ≈ outputs_bendM_check[i]
        end

        θ1 = 45 * pi /180
        θ2 = 135 * pi /180
        outputs_bendM_outer = TASOPT.structures.stiffeners_bendingM_outer(θ1, θ2)
        outputs_bendM_outer_check = [pi/2, 0.6166583566771067]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM_outer[i] ≈ outputs_bendM_outer_check[i]
        end

        Wstiff = TASOPT.structures.stiffener_weight("outer", 7e4, fuse_tank.Rfuse, fuse_tank.UTSouter / 4, 
        fuse_tank.rhoouter , θ1, θ2, 10.0, 5.0, fuse_tank.Eouter )
        Wstiff_check = 1302.5839377437662 
        
        @test Wstiff ≈ Wstiff_check
    end

    @testset "Thermal models" begin
        outputs_h = TASOPT.structures.tank_heat_coeffs(21.0, ifuel, fuse_tank.Tfuel, 5.0)
        outputs_h_check = (105.81383176854604, 447000.0)
        for i in 1:length(outputs_h)
            @test outputs_h[i] ≈ outputs_h_check[i]
        end

        outputs_free = TASOPT.structures.freestream_heat_coeff(z, Mair, xftank, 240.0)
        outputs_free_check = (93.08029151543289, 218.06145060705106, 242.96196484061733)
        for i in 1:length(outputs_free)
            @test outputs_free[i] ≈ outputs_free_check[i]
        end

        Taw = outputs_free_check[3]
        Rvac = TASOPT.structures.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0)
        Rvac_check = 0.5500815195658247

        @test Rvac ≈ Rvac_check
    end

end