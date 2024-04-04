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
fuse_tank.inner_material = TASOPT.StructuralAlloy("Al-2219-T87")

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

fuse_tank.outer_material = TASOPT.StructuralAlloy("Al-2219-T87")
fuse_tank.theta_outer = [1.0, 2.0]
fuse_tank.Ninterm = 1.0

Tsat, ρl, ρg, hv = TASOPT.cryo_fuel_properties("LH2", fuse_tank.ptank)

fuse_tank.rhofuel = ρl
fuse_tank.rhofuelgas = ρg
fuse_tank.Tfuel = Tsat
fuse_tank.hvap = hv
fuse_tank.Wfuelintank = 1e5

@testset "Fuselage tank" begin
    @testset "Foam insulation" begin
        outputs_size = TASOPT.structures.tanksize!(fuse_tank, z, Mair, xftank,
                                        time_flight,
                                        ifuel)
        outputs_size_check = (0.004247366632687734, 166.77327116515787, 1.8560115697526485, 39681.4881770504, 17.18053325519887, 62912.821723497895)
        
        for i in 1:length(outputs_size)
            @test outputs_size[i] ≈ outputs_size_check[i]
        end

        outputs_mech = TASOPT.structures.size_inner_tank(fuse_tank, fuse_tank.t_insul)

        outputs_mech_check = (162912.82172349788, 14.233087927314303, 0.0034633543007140293, 1.8560115697526485, 166.77327116515787, 62912.821723497895, 100000.0, 39681.4881770504, 0.003460125968964669, 1435.539306816222, 16809.755775358066, 1438.5597441435732, [12929.372061423137, 14427.140430195976, 12324.97568543129], 203.99444667419192, [14.95824234874573, 18.670635236161818, 22.800445249066772, 27.346434182678916], 17.18053325519887)
        for i in 1:length(outputs_mech)
            @test outputs_mech[i] ≈ outputs_mech_check[i]
        end
        
        l_cyl = outputs_mech_check[2]
        l_tank = outputs_mech_check[16]
        r_tank = outputs_mech_check[4]
        Shead = outputs_mech_check[15]
        outputs_thermal = TASOPT.structures.tankWthermal(fuse_tank, z, Mair, xftank, time_flight, ifuel)

        outputs_thermal_check = (1835.8608916060898, 107.03363914395707, 0.004247366632696709)

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
        outputs_vac_size_check = (0.0022166022747804252, 166.77327116515787, 2.4, 0.0, 10.884763251319967, 114902.51923675802)
        
        for i in 1:length(outputs_vac_size)
            @test outputs_vac_size[i] ≈ outputs_vac_size_check[i]
        end

        outputs_vac_mech = TASOPT.structures.size_inner_tank(fuse_tank, fuse_tank.t_insul)
        outputs_vac_mech_check = (125124.15434149424, 8.219175642639724, 0.00435715618585557, 2.335, 166.77327116515787, 25124.15434149424, 100000.0, 0.0, 0.004353094705443698, 2858.4704423225353, 15209.457128611886, 1913.7422971923515, [0.0], 172.43073912598143, [23.675162947566548, 25.28189535258484], 10.688524675864757)

        for i in 1:length(outputs_vac_mech)
            @test outputs_vac_mech[i] ≈ outputs_vac_mech_check[i]
        end
        Winnertank = outputs_vac_mech_check[1]
        l_cyl = outputs_vac_mech_check[2]

        fuse_tank.Ninterm = 1.0
        Ninterm = TASOPT.structures.optimize_outer_tank(fuse_tank, Winnertank, l_cyl)
        Ninterm_check = 15.38916015625
        @test Ninterm ≈ Ninterm_check

        outputs_vac_outer = TASOPT.structures.size_outer_tank(fuse_tank, Winnertank, l_cyl, Ninterm_check)
        outputs_vac_outer_check = (88153.61722490133, 33033.25329386777, 12808.143626608879, 21490.11147555205, 173.928350628553, 24.993050952321795, 123.9422487239094, 0.009566313929022768, 0.018394143361195915, 10.655963929362116)
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

        Wstiff = TASOPT.structures.stiffener_weight("outer", 7e4, fuse_tank.Rfuse, fuse_tank.outer_material.UTS / 4, 
        fuse_tank.outer_material.ρ, θ1, θ2, 10.0, 5.0, fuse_tank.outer_material.E)
        Wstiff_check = 1653.9008501855653
        
        @test Wstiff ≈ Wstiff_check
    end

    @testset "Thermal models" begin
        outputs_h = TASOPT.structures.tank_heat_coeff(21.0, ifuel, fuse_tank.Tfuel, 5.0)
        outputs_h_check = (133.43555789562762)
        for i in 1:length(outputs_h)
            @test outputs_h[i] ≈ outputs_h_check[i]
        end

        outputs_free = TASOPT.structures.freestream_heat_coeff(z, Mair, xftank, 240.0)
        outputs_free_check = (91.1765275792035, 218.06145060705106, 243.3345022554043)
        for i in 1:length(outputs_free)
            @test outputs_free[i] ≈ outputs_free_check[i]
        end

        outputs_natural = TASOPT.structures.freestream_heat_coeff(0.0, 0.0, xftank, 287.0, fuse_tank.Rfuse)
        outputs_natural_check = (1.6786556204397758, 288.2, 288.2)
        for i in 1:length(outputs_natural)
            @test outputs_natural[i] ≈ outputs_natural_check[i]
        end

        Taw = outputs_free_check[3]
        Rvac = TASOPT.structures.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0)
        Rvac_check = 0.5408710346147904

        @test Rvac ≈ Rvac_check
    end

    @testset "Fuel properties" begin
        ps = [1, 2, 5] * 101325.0 #pressures to test at

        fuel = "LH2"
        outputs_lh2_check = [20.369 70.848 1.3322 448.71e3; 22.965 67.639 2.5133 431.474e3; 27.314 60.711 6.1715 371.36e3]
    
        for (i, p) in enumerate(ps)
            outputs_lh2 = TASOPT.cryo_fuel_properties(fuel, p)
    
            outputs_check = outputs_lh2_check[i, :]
            for j in 1:length(outputs_check)
                @test outputs_lh2[j] ≈ outputs_check[j] rtol = 1e-2 #Only require 1% diff as function uses fits
            end
        end

        fuel = "CH4"
        outputs_ch4_check = [111.67 422.36 1.8164 510.83e3; 120.81 409.78 3.4382 492.921e3; 135.59 384.63 8.1024 457.614e3]
    
        for (i, p) in enumerate(ps)
            outputs_ch4 = TASOPT.cryo_fuel_properties(fuel, p)
    
            outputs_check = outputs_ch4_check[i, :]
            for j in 1:length(outputs_check)
                @test outputs_ch4[j] ≈ outputs_check[j] rtol = 1e-2 #Only require 1% diff as function uses fits
            end
        end
    
    end
end