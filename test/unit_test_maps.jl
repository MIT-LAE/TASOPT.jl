@testset "Compressor maps" begin
    @testset "Loading tables" begin
        # Test that data has been loaded in correctly for E3 maps;
        # choose the [23, 45] point in Nb_mb matrix
        # NOTE - These values will be changed for different maps
        @test TASOPT.engine.E3fan.Nb_mb[23, 45] ≈ 0.0001317367946127749
        @test TASOPT.engine.E3lpc.Nb_mb[23, 45] ≈ 0.0037887310772633156
        @test TASOPT.engine.E3hpc.Nb_mb[23, 45] ≈ 0.0018892705002216022
    end

    @testset "Accuracy vs NPSS" begin
        # Get sample speed values
        Nfan, _, _ = TASOPT.engine.NcTblMap(1.294, 953.871, 1.537, 1320.644, 1.0, TASOPT.engine.E3fan)
        Nlpc, _, _ = TASOPT.engine.NcTblMap(3.060, 59.010, 4.506, 88.799, 1.0, TASOPT.engine.E3lpc)
        Nhpc, _, _ = TASOPT.engine.NcTblMap(10.811, 23.412, 12.504, 25.546, 1.0, TASOPT.engine.E3hpc)
        
        # Interpolation seems to be accurate at 1-2% relative error
        # Values are taken from the LAE CFM56-7B model developed by Lee
        @test Nfan ≈ 0.769 rtol=1e-2
        @test Nlpc ≈ 0.769 rtol=1e-2
        @test Nhpc ≈ 0.977 rtol=1e-2

        # Get sample efficiency values
        Efan, _, _ = TASOPT.engine.ecTblMap(1.294, 953.871, 1.537, 1320.644, TASOPT.engine.E3fan, 0.9043)
        Elpc, _, _ = TASOPT.engine.ecTblMap(3.060, 59.010, 4.506, 88.799, TASOPT.engine.E3lpc, 0.8901)
        Ehpc, _, _ = TASOPT.engine.ecTblMap(10.811, 23.412, 12.504, 25.546, TASOPT.engine.E3hpc, 0.8803)

        # Interpolation seems to be accurate at 1-2% relative error
        # Values are taken from the LAE CFM56-7B model developed by Lee
        @test Efan ≈ 0.8903 rtol=1e-2
        @test Elpc ≈ 0.8891 rtol=1e-2
        @test Ehpc ≈ 0.8784 rtol=1e-2
    end
end