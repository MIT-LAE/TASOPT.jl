@testset "Materials" verbose=true begin
    database = TASOPT.materials.MaterialProperties
    @testset "Structural Alloys" begin
        Al7075 = StructuralAlloy("Al-7075")
        @test Al7075.ρ == database["Al-7075"]["density"]
        @test Al7075.E == database["Al-7075"]["youngs_modulus"]
        #Check that safety_factors are being applied
        @test Al7075.σmax == database["Al-7075"]["YTS"]/1.1/1.5 
        Al = StructuralAlloy("Al-7075"; max_avg_stress = 1.0, safety_factor = 2.0)
        @test Al.σmax == database["Al-7075"]["YTS"]/1.0/2.0
        @test Al.τmax == database["Al-7075"]["shear_strength"]/1.0/2.0
        @test Al.YTS == database["Al-7075"]["YTS"]
        @test Al.UTS == database["Al-7075"]["UTS"]
        @test Al.USS == database["Al-7075"]["shear_strength"]
    end

    @testset "Conductors" begin
        Cu = Conductor("Cu")
        @test Cu.ρ == database["Cu"]["density"] 
        @test TASOPT.materials.resxden(Cu) == 0.000150528
        @test TASOPT.materials.resistivity(Cu) == database["Cu"]["resistivity"]
        @test TASOPT.materials.resistivity(Cu, 300.0) ≈ 1.7264923200000005e-8
    end

    @testset "Insulators" begin
        PTFE = Insulator("PTFE")
        @test PTFE.ρ == database["PTFE"]["density"]
    end

    @testset "Thermal Insulators" begin
        poly = ThermalInsulator("polyurethane32")
        @test poly.ρ == database["polyurethane32"]["density"]
    end
    
    # Test that unreasonable requests throw errors
    ## Try making something structural with silver, a conductor using
    ## insulator properties or an insulator using conductor properties.
    @testset "Material errors" begin
        @test_throws ErrorException StructuralAlloy("unobtanium")
        @test_throws ErrorException Conductor("unobtanium")
        @test_throws ErrorException Insulator("unobtanium")
        @test_throws ErrorException ThermalInsulator("unobtanium")

        @test_throws ErrorException StructuralAlloy("Ag")
        @test_throws ErrorException Conductor("PTFE")
        @test_throws ErrorException Insulator("Cu")
    end
end
