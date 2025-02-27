@testset "Weights" begin
    #Test frame
    frame = TASOPT.structures.Frame()
    @test frame.origin == [0.0, 0.0, 0.0]

    #Test Weight constructors
    W = TASOPT.Weight()
    @test W.r == [0.0, 0.0, 0.0]
    @test W.frame === TASOPT.WORLD
    @test W.W == 0.0

    W = TASOPT.Weight(1.0, [1.0, 0.0, 0.0])
    @test W.W == 1.0
    @test W.r == [1.0, 0.0, 0.0]

    W1 = TASOPT.Weight(W = 1.0, frame = frame)
    W2 = TASOPT.Weight(W = 2.0)
    W3 = TASOPT.Weight(W = 2.0, x = 10.0, y = 10.0, z = 10.0)

    # Test exceptions
    @test_throws "Cannot add weights in different frames" W1 + W2
    @test_nowarn W2 + W3

    # Test operations
    W = W2 + W3
    @test W.W == 4.0
    @test W.r == [5.0, 5.0, 5.0]
    @test W.x == 5.0
    @test W.y == 5.0
    @test W.z == 5.0

    W = W2 * 2.0
    @test W.W == 4.0
    @test W.r == [0.0, 0.0, 0.0]

    # Test that * is non-mutating
    @test W2.W == 2.0
    # Test scale! is mutating
    TASOPT.structures.scale!(W2, 3.0)
    @test W2.W == 6.0

    TASOPT.structures.scale!(W2, 1 / 3.0)
    @test W2.W == 2.0

    W_sum = sum([W2, W3])
    W_CoW = TASOPT.structures.center_of_weight([W3, W2])
    @test W_sum.W == W_CoW.W
    @test W_sum.r == W_CoW.r

    M2 = TASOPT.structures.moment(W2)
    @test M2 == [0.0, 0.0, 0.0]
    M3 = TASOPT.structures.moment(W3)
    @test M3 == [-20.0, 20.0, 0.0]
end
