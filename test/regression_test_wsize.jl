
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
        test_EIs = [35.727087484447225, 26.674239512749597, 10.228895832877859, 6.473190938712144]
        test_mfs = [1.4031641016757275, 1.2552506710045581, 0.7022981161085624, 0.4126470600128394]
        for i in eachindex(test_EIs)
            @test test_EIs[i] ≈ EIs[i]
            @test test_mfs[i] ≈ mfs[i]
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