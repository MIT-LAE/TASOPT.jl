@testset "Compressor maps" begin
    @testset "Gridded tables" begin
        fan, lpc, hpc = TASOPT.engine.E3fan, TASOPT.engine.E3lpc, TASOPT.engine.E3hpc
        
        # Test that data has been loaded in correctly for E3 maps;
        # choose the [23, 45] point in Nb_mb matrix
        # NOTE - These values will be changed for different maps
        @test TASOPT.engine.E3fan.Nb_mb[23, 45] ≈ 0.00014094728533512843
        @test TASOPT.engine.E3lpc.Nb_mb[23, 45] ≈ 0.0037887310772633156
        @test TASOPT.engine.E3hpc.Nb_mb[23, 45] ≈ 0.0018892705002216022
    end
end