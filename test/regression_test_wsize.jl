
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

@testset "Hydrogen sizing" verbose=true begin
    ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    
    include(joinpath(TASOPT.__TASOPTroot__, "./misc/index.inc"))

    @test ac.parg[igRfuse] ≈ 2.54
    
    include(joinpath(TASOPT.__TASOPTroot__, "../test/hydrogen_sized.jl"))

    size_aircraft!(ac, iter=50; printiter=false);

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
    
    @test ac.parm[imPFEI] ≈ 0.9806661849624643

end