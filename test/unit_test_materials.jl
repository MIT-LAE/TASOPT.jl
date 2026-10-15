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

    @testset "S1 database: CFRP entries and additivity" begin
        cfrp_keys = [k for k in keys(database) if startswith(k, "CFRP-")]
        @test !isempty(cfrp_keys)  # guard against a vacuous pass
        required = ["density", "E1", "E2", "E_equivalent_buckling", "nu_equivalent_buckling",
            "ply_thickness", "layup", "strain_allowable_tension"]
        for k in cfrp_keys
            @test all(haskey(database[k], r) for r in required)
            # Cured-laminate density must be resin-inclusive, never the dry-tow ~1780/1790
            @test 1400.0 <= database[k]["density"] <= 1700.0
            if haskey(database[k], "fiber_volume_fraction")
                Vf = database[k]["fiber_volume_fraction"]
                @test isapprox(Vf*1770.0 + (1 - Vf)*1301.0, database[k]["density"], rtol = 0.01)
            end
        end
        # Appending the CFRP tables must leave existing alloys unchanged
        Al2219 = StructuralAlloy("Al-2219-T87")
        @test Al2219.ρ == 2840.0
        @test Al2219.E == 73.1e9
        @test Al2219.UTS == 476e6
        @test Al2219.YTS == 393e6
    end

    @testset "S2 T1A: StructuralAlloy is a TankWallMaterial" begin
        @test StructuralAlloy <: TASOPT.materials.TankWallMaterial
    end

    @testset "S2 T1B: CompositeLaminate" begin
        CL = TASOPT.materials.CompositeLaminate
        default = CL("CFRP-IM7-8552-PW-0090")
        @test default isa TASOPT.materials.TankWallMaterial
        @test default.ρ ≈ 1560.0
        @test default.E1 ≈ 76.0e9
        @test !isempty(sprint(show, default))

        # plies per stack is parsed once from the layup string (first bracketed group, also with trailing text)
        @test default.plies_per_stack == 5
        @test CL("CFRP-IM7-8552-PW-0090-VDA").plies_per_stack == 5
        @test CL("CFRP-IM7-8552-PW-QI").plies_per_stack == 4
        @test CL("CFRP-IM7-8552-UD-RING").plies_per_stack == 1
        @test CL("CFRP-IM7-8552-PW-0090"; plies_per_stack = 3).plies_per_stack == 3

        # Alloys, dry fibre, resin and unknown names are not laminates
        @test_throws ErrorException CL("Al-2219-T87")
        @test_throws ErrorException CL("CF-IM7-fiber")
        @test_throws ErrorException CL("CF-8552-resin")
        @test_throws ErrorException CL("unobtanium")
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
