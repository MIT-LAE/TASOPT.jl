using Unicode
@testset "outputs" verbose=true begin

    ac = load_default_model()
    size_aircraft!(ac; printiter=false);

    @testset "output summaries" begin
        f = open(io->TASOPT.weight_buildup(ac,io=io), "temp.txt", "w")
        @test Unicode.normalize(read("weights.txt", String),newline2lf=true) == read("temp.txt", String)
        rm("temp.txt")

        f = open(io->TASOPT.aero(ac,io=io), "temp.txt", "w")
        @test Unicode.normalize(read("aero.txt", String),newline2lf=true) == read("temp.txt", String)
        rm("temp.txt")

        f = open(io->TASOPT.geometry(ac,io=io), "temp.txt", "w")
        @test Unicode.normalize(read("geom.txt", String),newline2lf=true) == read("temp.txt", String)
        rm("temp.txt")
    end

    @testset "output plots" begin
        
        TASOPT.stickfig(ac)
        @test 1 == 1
        TASOPT.high_res_airplane_plot(ac)
        @test 2 == 2
        TASOPT.plot_details(ac)
        @test 3 == 3

    end

end