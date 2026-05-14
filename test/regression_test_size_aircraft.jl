rtol_sizing = 1e-4 #coarse = 10% =1e-1; fine = 0.001% = 1e-5
# Define a function to check if each value in two structs is equal
function check_struct_equivalence(s1, s2; verbose::Bool=false, _path::String="")
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
        fieldpath = _path == "" ? string(field) : _path * "." * string(field)
        if typeof(val1) == typeof(val2)
            if typeof(val1) != Float64
                if !check_struct_equivalence(val1, val2; verbose, _path=fieldpath)
                    return false
                end
            else
                if verbose && !isapprox(val1, val2; rtol=rtol_sizing)
                    println("MISMATCH: $fieldpath  =>  $val1 vs $val2")
                end
                @test val1 ≈ val2 rtol=rtol_sizing
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
        try
            TASOPT.fly_mission!(ac, 2; printTO=false)
            @test_broken ac.parm[imPFEI, 2] ≈ PFEI rtol=rtol_sizing
            @test_broken ac.parm[imWfuel, 2] ≈ Wfuel rtol=rtol_sizing
            @test_broken ac.parm[imWTO, 2] ≈ WTO rtol=rtol_sizing
        catch #indicates broken test if fly_mission!() fails
            @test_broken false
        end
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
    fly_mission!(ac, 1);

    @testset "Fuselage" begin
        @test  check_struct_equivalence(ac_test.fuselage, ac.fuselage, verbose=true)
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
            @test parg[i] ≈ ac.parg[i]  rtol=rtol_sizing
        end
    end

    @testset "Aero" begin
        for i in eachindex(para)
            @test para[i] ≈ ac.para[i]  rtol=rtol_sizing
        end
    end

    @testset "Propulsion" begin
        for i in eachindex(pare)
            @test pare[i] ≈ ac.pare[i]  rtol=rtol_sizing
        end
    end

    test_ac_off_design(ac, 1.0869638391729122, 153128.29535348987,  769359.1150444464)
    
    @test ac.parm[imPFEI] ≈ 0.9456457746362635  rtol=rtol_sizing
end

@testset "Wide sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_wide.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 3.0988 rtol=rtol_sizing
    

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.1903000293875512 rtol=rtol_sizing

end

@testset "Regional sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 1.5113 rtol=rtol_sizing

    size_aircraft!(ac; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 0.8482236162426772 rtol=rtol_sizing

end

@testset "Hydrogen sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    
    include(__TASOPTindices__)

    @test ac.fuselage.layout.radius ≈ 2.54 rtol=rtol_sizing

    size_aircraft!(ac, iter=50; printiter=false);
    
    @test ac.parm[imPFEI] ≈ 1.005873551903147 rtol=rtol_sizing

end
