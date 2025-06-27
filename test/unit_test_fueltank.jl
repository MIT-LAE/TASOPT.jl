#Sample parameters
ac = load_default_model()

ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
ac.options.ifuel = ifuel
ac.para[iaalt, ipcruise1,1] = z
ac.para[iaMach, ipcruise1,1] = Mair
ac.parg[igxftank] = xftank
ac.fuse_tank.tank_count = 1
ac.parg[igWfuel] = 1e5

fuse_tank = ac.fuse_tank
fuse_tank.fueltype = "LH2"
fuse_tank.placement = "front"

fuse_tank.clearance_fuse = 0.1
TSL = 288.2
fuse_tank.TSLtank = [TSL]

fuse_tank.pvent = 2e5
fuse_tank.ARtank = 2.0
fuse_tank.inner_material = TASOPT.StructuralAlloy("Al-2219-T87")

fuse_tank.t_insul = [0.15,0.15,0.15]
fuse_tank.material_insul = [ThermalInsulator("rohacell41s"), ThermalInsulator("rohacell41s"), ThermalInsulator("polyurethane27")]
fuse_tank.sizes_insulation = true
fuse_tank.iinsuldes = [1,2,3]
fuse_tank.boiloff_rate = 0.15

fuse_tank.ullage_frac = 0.1
fuse_tank.qfac = 1.3
fuse_tank.theta_inner = 1.0
fuse_tank.ftankadd = 0.1
fuse_tank.ew = 0.9

fuse_tank.outer_material = StructuralAlloy("Al-2219-T87")
fuse_tank.theta_outer = [1.0, 2.0]
fuse_tank.Ninterm = 1.0

#Calculate fuel temperature and density as a function of pressure
β0 = 1 - fuse_tank.ullage_frac
fuel_mix = TASOPT.SaturatedMixture(fuse_tank.fueltype, fuse_tank.pvent, β0)

Tsat = fuel_mix.liquid.T
ρl = fuel_mix.liquid.ρ
ρg = fuel_mix.gas.ρ
hvap = fuel_mix.hvap

fuse_tank.rhofuel = ρl
fuse_tank.rhofuelgas = ρg
fuse_tank.Tfuel = Tsat
fuse_tank.hvap = hvap

#Fuselage
fuse = ac.fuselage
fuse.layout.cross_section.radius = 2.5
fuse.layout.cross_section.bubble_lower_downward_shift = 0.3

@testset "Fuselage tank" begin
    @testset "Foam insulation" begin

        TASOPT.tanksize!(ac, 1)
        outputs_mech = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)

        outputs_mech_check = (60109.124334566346 , 38147.81559388217, 166.77327116515787, [16.42133035129676, 20.311090752967278, 24.620380536430755, 29.34797626786876], 1.8743852420176998, 15.771807531012154, 12.842701653674803)
        for i in 1:length(outputs_mech)
            @test outputs_mech[i] ≈ outputs_mech_check[i]
        end
        
        outputs_thermal = TASOPT.CryoTank.tankWthermal(fuse, fuse_tank, z, TSL, Mair, xftank, ifuel)

        outputs_thermal_check = (1835.8608916011449,)

        for i in 1:length(outputs_thermal)
            @test outputs_thermal[i] ≈ outputs_thermal_check[i]
        end
    end
    
    fuse_tank.t_insul = [0.065]
    fuse_tank.material_insul = [ThermalInsulator("vacuum")]
    fuse_tank.sizes_insulation = false

    @testset "Vacuum insulation" begin
        TASOPT.tanksize!(ac, 1)

        outputs_vac_mech = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        outputs_vac_mech_check = (24229.845334360776, 0.0, 166.77327116515787, [25.483812169139977, 27.21328988831087], 2.335, 9.994915110929762, 7.525566077704729)

        for i in 1:length(outputs_vac_mech)
            @test outputs_vac_mech[i] ≈ outputs_vac_mech_check[i]
        end
        Winnertank = outputs_vac_mech_check[1]
        l_cyl = outputs_vac_mech_check[7]
        Winner_tot = Winnertank + fuse_tank.Wfuelintank

        fuse_tank.Ninterm = 1.0
        Ninterm = TASOPT.CryoTank.optimize_outer_tank(fuse, fuse_tank, Winner_tot, l_cyl)
        Ninterm_check = 14.025390625
        @test Ninterm ≈ Ninterm_check

        outputs_vac_outer = TASOPT.CryoTank.size_outer_tank(fuse, fuse_tank, Winner_tot, l_cyl, Ninterm_check)
        outputs_vac_outer_check = (87481.40843690369, 31381.827647692957, 13786.613724211644, 20573.498028341648, 171.6223477616751, 26.90237940128641, 117.81758895910225, 0.009560503271690743, 0.018394143361195915, 9.96235436442712)
        for i in 1:length(outputs_vac_outer)
            @test outputs_vac_outer[i] ≈ outputs_vac_outer_check[i]
        end

    end
    
    @testset "Stringer sizing" begin
        θ = 85 * pi /180
        outputs_bendM = TASOPT.CryoTank.stiffeners_bendingM(θ)
        outputs_bendM_check = [1.1749999999999998, 0.09494892091249901]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM[i] ≈ outputs_bendM_check[i]
        end

        θ1 = 45 * pi /180
        θ2 = 135 * pi /180
        outputs_bendM_outer = TASOPT.CryoTank.stiffeners_bendingM_outer(θ1, θ2)
        outputs_bendM_outer_check = [pi/2, 0.6166583566771067]
        for i in 1:length(outputs_bendM)
            @test outputs_bendM_outer[i] ≈ outputs_bendM_outer_check[i]
        end

        perim_vessel = fuse.layout.cross_section.perimeter

        Wstiff = TASOPT.CryoTank.stiffener_weight("outer", 7e4, fuse.layout.radius, perim_vessel, fuse_tank.outer_material.UTS / 4, 
        fuse_tank.outer_material.ρ, θ1, θ2, 10.0, 5.0, fuse_tank.outer_material.E)
        Wstiff_check = 1717.0752091513864
        
        @test Wstiff ≈ Wstiff_check
    end

    @testset "Thermal models" begin
        outputs_h = TASOPT.CryoTank.tank_heat_coeff(21.0, ifuel, fuse_tank.Tfuel, 5.0)
        outputs_h_check = (133.43555789562762)
        for i in 1:length(outputs_h)
            @test outputs_h[i] ≈ outputs_h_check[i]
        end

        outputs_free = TASOPT.CryoTank.freestream_heat_coeff(z, TSL, Mair, xftank, 240.0)
        outputs_free_check = (91.30101055991513, 218.06145060705106, 243.26644796151263)
        for i in 1:length(outputs_free)
            @test outputs_free[i] ≈ outputs_free_check[i]
        end

        outputs_natural = TASOPT.CryoTank.freestream_heat_coeff(0.0, TSL, 0.0, xftank, 287.0, fuse.layout.radius)
        outputs_natural_check = (1.6795624846472137, 288.2, 288.2)
        for i in 1:length(outputs_natural)
            @test outputs_natural[i] ≈ outputs_natural_check[i]
        end

        Taw = outputs_free_check[3]
        Rvac = TASOPT.CryoTank.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0)
        Rvac_check = 0.3565313051181731

        @test Rvac ≈ Rvac_check
    end
end