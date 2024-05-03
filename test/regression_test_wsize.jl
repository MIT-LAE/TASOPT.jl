
@testset "Default sizing" verbose=true begin
    ac = load_default_model()
    
    include(joinpath(TASOPT.__TASOPTroot__, "./misc/index.inc"))

    @test ac.parg[igRfuse] == 1.9558
    
    include(joinpath(TASOPT.__TASOPTroot__, "../test/default_sized.jl"))

    size_aircraft!(ac; printiter=false);

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
            @test pare[i] ≈ ac.pare[i] rtol=1e-6
        end
    end
    
    @test ac.parm[imPFEI] ≈  0.919121257897844

    @testset "LTO" begin
        EIs, mfs = TASOPT.LTO("Default A/C + Engine", ac)
        test_EIs = [35.76509639768726, 26.70084308861548, 10.22948188721507, 6.480767538587702]
        test_mfs = [1.401100172004722, 1.2535426491385873, 0.7015473738995113, 0.4123593208289815]
        for i in eachindex(test_EIs)
            @test test_EIs[i] ≈ EIs[i]
            @test test_mfs[i] ≈ mfs[i]
        end

        EIs, mfs = TASOPT.LTO("Default A/C + Engine", ac, method = "quartic")
        test_EIs = [35.72641687463056, 26.27891743413049, 10.095612208446196, 6.180472562714481]
        for i in eachindex(test_EIs)
            @test test_EIs[i] ≈ EIs[i]
        end

    end

end

@testset "Wide sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/example_widebody.toml"))
    
    include(joinpath(TASOPT.__TASOPTroot__, "./misc/index.inc"))

    @test ac.parg[igRfuse] ≈ 3.0988
    
    include(joinpath(TASOPT.__TASOPTroot__, "../test/wide_sized.jl"))

    size_aircraft!(ac; printiter=false);

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
            @test pare[i] ≈ ac.pare[i] rtol=1e-6
        end
    end
    
    @test ac.parm[imPFEI] ≈ 1.1642157697355595

end

@testset "Regional sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/example_regional.toml"))
    
    include(joinpath(TASOPT.__TASOPTroot__, "./misc/index.inc"))

    @test ac.parg[igRfuse] ≈ 1.5113
    
    include(joinpath(TASOPT.__TASOPTroot__, "../test/regional_sized.jl"))

    size_aircraft!(ac; printiter=false);

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
            @test pare[i] ≈ ac.pare[i] rtol=1e-6
        end
    end
    
    @test ac.parm[imPFEI] ≈ 0.8107247842565991

end