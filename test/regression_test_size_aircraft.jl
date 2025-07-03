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
        TASOPT.fly_mission!(ac, 2)

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

    test_ac_off_design(ac, 1.0869638391729122, 153128.29535348987,  769359.1150444464)
    
    @test ac.parm[imPFEI] ≈ 0.945758611404728 rtol=1e-4
end

@testset "Wide sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 3.0988
    

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.1903760871373523 rtol=1e-4

end

@testset "Regional sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 1.5113

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 0.8483560952994892 rtol=1e-4

end

@testset "Hydrogen sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 2.54

    size_aircraft!(ac, iter=50; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.0076619899926231 rtol=1e-4

end