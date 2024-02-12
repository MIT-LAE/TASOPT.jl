
@testset "outputs" verbose=true begin

    ac = load_default_model()
    size_aircraft!(ac; printiter=false);

    if Sys.iswindows()
        file_weights = "test/weights.txt"
        file_aero = "test/aero.txt"
        file_geom = "test/geom.txt"
    else
        file_weights = "weights.txt"
        file_aero = "aero.txt"
        file_geom = "geom.txt"
    end

    @testset "output summaries" begin
        f = open(io->TASOPT.weight_buildup(ac,io=io), "temp.txt", "w")
        content = read("temp.txt", String)
        content = replace(content, "\n"=> "\r\n") # Convert \n to \r\n
        @test read(file_weights, String) == content
        rm("temp.txt")

        f = open(io->TASOPT.aero(ac,io=io), "temp.txt", "w")
        content = read("temp.txt", String)
        content = replace(content, "\n"=> "\r\n") # Convert \n to \r\n
        @test read(file_aero, String) == content
        rm("temp.txt")

        f = open(io->TASOPT.geometry(ac,io=io), "temp.txt", "w")
        content = read("temp.txt", String)
        content = replace(content, "\n"=> "\r\n") # Convert \n to \r\n
        @test read(file_geom, String) == content
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