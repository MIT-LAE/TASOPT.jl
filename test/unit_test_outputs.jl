@testset "outputs" verbose=true begin

    ac = load_default_model()
    size_aircraft!(ac; printiter=false);

    file_weights = joinpath(TASOPT.__TASOPTroot__, "../test/weights.txt")
    file_aero = joinpath(TASOPT.__TASOPTroot__, "../test/aero.txt")
    file_geom = joinpath(TASOPT.__TASOPTroot__, "../test/geom.txt")
    
    @testset "output summaries" begin
        f = open(io->TASOPT.weight_buildup(ac,io=io), "temp.txt", "w")
        content = read(file_weights, String)
        content = replace(content, "\r\n"=> "\n") # Convert \r\n to \n
        @test content == read("temp.txt", String)
        rm("temp.txt")

        f = open(io->TASOPT.aero(ac,io=io), "temp.txt", "w")
        content = read(file_aero, String)
        content = replace(content, "\r\n"=> "\n") # Convert \r\n to \n
        @test content == read("temp.txt", String)
        rm("temp.txt")

        f = open(io->TASOPT.geometry(ac,io=io), "temp.txt", "w")
        content = read(file_geom, String)
        content = replace(content, "\r\n"=> "\n") # Convert \r\n to \n
        @test content == read("temp.txt", String)
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