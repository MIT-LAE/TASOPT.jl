# Define a function to check if each value in two structs is equal
function check_struct_equivalence(s1, s2)
    fields_s1 = fieldnames(typeof(s1))
    fields_s2 = fieldnames(typeof(s2))
    
    # Check if both structs have the same fields
    if fields_s1 != fields_s2
        return false
    end
    
    # Check if each field has the same value in both structs
    for field in fields_s1
        val1 = getproperty(s1, field)
        val2 = getproperty(s2, field)
        if typeof(val1) == typeof(val2)
            if typeof(val1) != Float64
                if !check_struct_equivalence(val1, val2)
                    return false
                end
            else
                # println(field)
                @test val1 ≈ val2 
            end
        else
            return false
        end
    end
    
    return true
end

#Simple function to call fly_mission!() and test on- and off-design performance
function test_ac_off_design(ac, PFEI, Wfuel, WTO)
    @testset "Off-design" begin
        TASOPT.fly_mission!(ac, 2; printTO=false)

        @test ac.parm[imPFEI, 2] ≈ PFEI
        @test ac.parm[imWfuel, 2] ≈ Wfuel
        @test ac.parm[imWTO, 2] ≈ WTO
    end
end

@testset "Default sizing" verbose=true begin
    ac = load_default_model()
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius == 1.9558
    
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_sized.jl"))
    # Fuselage
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_structures.jl"))

    size_aircraft!(ac; printiter=false);

    @testset "Fuselage" begin
        @test  check_struct_equivalence(ac_test.fuselage, ac.fuselage)
    end

    @testset "Wing" begin
        @test  check_struct_equivalence(ac_test.wing, ac.wing)
    end

    @testset "Htail" begin
        @test  check_struct_equivalence(ac_test.htail, ac.htail)
    end

    @testset "Vtail" begin
        @test  check_struct_equivalence(ac_test.vtail, ac.vtail)
    end

    @testset "Geometry" begin
        for i in eachindex(parg)
            @test parg[i] ≈ ac.parg[i] 
        end
    end

    @testset "Aero" begin
        for i in eachindex(para)
            @test para[i] ≈ ac.para[i] 
        end
    end

    @testset "Propulsion" begin
        for i in eachindex(pare)
            @test pare[i] ≈ ac.pare[i] 
        end
    end

    test_ac_off_design(ac, 1.0800347958058685, 152152.15193363777, 763293.7221007551)
    
    @test ac.parm[imPFEI] ≈ 0.9402313710034372 rtol=1e-4
end

# Compare imported models with the independently replayed source-study masses.
# The 0.01% tolerance covers rounded published values, not APD fit acceptance.
@testset "Aircraft comparison models" verbose=true begin
    for (class, mass_kg) in (("narrow", 174842.5 * 4.448222 / 9.81),
                             ("regional", 52065.09),
                             ("wide", 775421.06165 * 4.448222 / 9.81))
        @testset "$class" begin
            path = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_$class.toml")
            ac = read_aircraft_model(path)
            result = size_aircraft!(ac; printiter=false)
            @test result.success
            @test ac.is_sized[1]
            @test ac.parg[igWMTO] / 9.81 ≈ mass_kg rtol=1e-4
            @test ac.parg[igWMTO] > ac.parg[igWfuel] + ac.parg[igWpay] > 0
            for mission in 2:size(ac.parm, 2)
                @test TASOPT.fly_mission!(ac, mission; printTO=false)
                @test isfinite(ac.parm[imPFEI, mission])
                @test 0 < ac.parm[imWfuel, mission] < ac.parm[imWTO, mission]
            end
        end
    end
end

@testset "Hydrogen sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 2.54

    size_aircraft!(ac, iter=50; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.0042440592543385 rtol=1e-4

end
