#Sample parameters
ifuel = 40
z = 11e3
Mair = 0.8
xftank = 15.0
time_flight = 7*3600.0

fuse_tank = TASOPT.fuselage_tank()
fuse_tank.fueltype = "LH2"

fuse_tank.clearance_fuse = 0.1
fuse_tank.TSLtank = 288.2
fuse_tank.TSLtank = 288.2

fuse_tank.pvent = 2e5
fuse_tank.ARtank = 2.0
fuse_tank.inner_material = TASOPT.StructuralAlloy("Al-2219-T87")

fuse_tank.t_insul = [0.15,0.15,0.15]
fuse_tank.material_insul = [ThermalInsulator("rohacell41s"), ThermalInsulator("rohacell41s"), ThermalInsulator("polyurethane27")]
fuse_tank.size_insulation = true
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
fuse_tank.Wfuelintank = 1e5

#Fuselage
fuse = TASOPT.Fuselage()
fuse.layout.cross_section.radius = 2.5
fuse.layout.cross_section.bubble_lower_downward_shift = 0.3

@testset "Fuselage tank" begin
    @testset "Foam insulation" begin
        outputs_size = TASOPT.CryoTank.tanksize!(fuse, fuse_tank, z, Mair, xftank,
                                        time_flight,
                                        ifuel)
        outputs_size_check = (0.004247366632687734, 166.77327116515787, 1.8740119151889265, 38185.08012765128, 15.778042937598197, 60140.568517238586)
        
        for i in 1:length(outputs_size)
            @test outputs_size[i] ≈ outputs_size_check[i]
        end

        outputs_mech = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)

        outputs_mech_check = (160140.5685172386, 12.848564428768595, 0.003496943301341531, 1.8740119151889265, 166.77327116515787, 60140.568517238586, 100000.0, 38185.08012765128, 0.003493683659934613, 1590.6019601067815, 15273.952985313967, 1504.3779940972875, [12469.462987273082, 13880.427717198025, 11835.189423180173], 189.4574311887127, [16.41478963327391, 20.306728650515254, 24.618791939176347, 29.34975388954838], 15.778042937598197)
        for i in 1:length(outputs_mech)
            @test outputs_mech[i] ≈ outputs_mech_check[i]
        end
        
        l_cyl = outputs_mech_check[2]
        l_tank = outputs_mech_check[16]
        r_tank = outputs_mech_check[4]
        Shead = outputs_mech_check[15]
        outputs_thermal = TASOPT.CryoTank.tankWthermal(fuse, fuse_tank, z, Mair, xftank, time_flight, ifuel)

        outputs_thermal_check = (1835.8608916011449, 107.0336391436688, 0.0042473666326852694)

        for i in 1:length(outputs_thermal)
            @test outputs_thermal[i] ≈ outputs_thermal_check[i]
        end
    end
    
    fuse_tank.t_insul = [0.065]
    fuse_tank.material_insul = [ThermalInsulator("vacuum")]
    fuse_tank.size_insulation = false

    @testset "Vacuum insulation" begin

        outputs_vac_size = TASOPT.CryoTank.tanksize!(fuse, fuse_tank, z, Mair, xftank,
                                            time_flight,
                                            ifuel)
        outputs_vac_size_check = (0.0033447002334008723, 166.77327116515787, 2.4, 0.0, 10.17184549721119, 113248.97600936973)
        
        for i in 1:length(outputs_vac_size)
            @test outputs_vac_size[i] ≈ outputs_vac_size_check[i]
        end

        outputs_vac_mech = TASOPT.CryoTank.size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        outputs_vac_mech_check = (124223.11547325406, 7.525566077704729, 0.00435715618585557, 2.335, 166.77327116515787, 24223.115473254056, 100000.0, 0.0, 0.004353094705443698, 3076.8414985998347, 13888.811713075278, 1978.5193563196497, [0.0], 165.15270608329905, [25.483812169139977, 27.21328988831087], 9.994915110929762)

        for i in 1:length(outputs_vac_mech)
            @test outputs_vac_mech[i] ≈ outputs_vac_mech_check[i]
        end
        Winnertank = outputs_vac_mech_check[1]
        l_cyl = outputs_vac_mech_check[2]

        fuse_tank.Ninterm = 1.0
        Ninterm = TASOPT.CryoTank.optimize_outer_tank(fuse, fuse_tank, Winnertank, l_cyl)
        Ninterm_check = 14.025390625
        @test Ninterm ≈ Ninterm_check

        outputs_vac_outer = TASOPT.CryoTank.size_outer_tank(fuse, fuse_tank, Winnertank, l_cyl, Ninterm_check)
        outputs_vac_outer_check = (87481.39236647873, 31381.827647692957, 13786.613724211644, 20573.48341886441, 171.6223477616751, 26.90237940128641, 117.81758895910225, 0.009560503271690743, 0.018394143361195915, 9.96235436442712)
        for i in 1:length(outputs_vac_outer)
            @test outputs_vac_outer[i] ≈ outputs_vac_outer_check[i]
        end

    end
    
    @testset "Stringer sizing" begin
        θ = 85 * pi /180
        outputs_bendM = TASOPT.CryoTank.stiffeners_bendingM(θ)
        outputs_bendM_check = [1.1693705988362009, 0.09494234539184043 ]
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

        outputs_free = TASOPT.CryoTank.freestream_heat_coeff(z, fuse_tank.TSLtank, Mair, xftank, 240.0)
        outputs_free_check = (91.1765275792035, 218.06145060705106, 243.3345022554043)
        for i in 1:length(outputs_free)
            @test outputs_free[i] ≈ outputs_free_check[i]
        end

        outputs_natural = TASOPT.CryoTank.freestream_heat_coeff(0.0, fuse_tank.TSLtank, 0.0, xftank, 287.0, fuse.layout.radius)
        outputs_natural_check = (1.6786556204397758, 288.2, 288.2)
        for i in 1:length(outputs_natural)
            @test outputs_natural[i] ≈ outputs_natural_check[i]
        end

        Taw = outputs_free_check[3]
        Rvac = TASOPT.CryoTank.vacuum_resistance(fuse_tank.Tfuel, Taw, 90.0, 100.0)
        Rvac_check = 0.35637094066735286

        @test Rvac ≈ Rvac_check
    end
end