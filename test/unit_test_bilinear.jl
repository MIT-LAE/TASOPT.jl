@testset "Bilinear bounded" begin
    xGrid        = 0.0:0.1:1.0
    yGrid        = 0.0:0.1:1.0
    nominal_vals = zeros(Float64, length(xGrid), length(yGrid))
    derivx_vals  = zeros(Float64, length(xGrid), length(yGrid))
    derivy_vals  = zeros(Float64, length(xGrid), length(yGrid))

    dx = xGrid[2] - xGrid[1]
    Nx = length(xGrid)

    dy = yGrid[2] - yGrid[1]
    Ny = length(yGrid)

    f    = (x, y) -> x*(y^2)
    dfdx = (x, y) -> y^2
    dfdy = (x, y) -> 2*x*y
    
    for i in eachindex(xGrid)
        for j in eachindex(yGrid)
            nominal_vals[i,j] = f(xGrid[i], yGrid[j])
            derivx_vals[i,j]  = dfdx(xGrid[i], yGrid[j])
            derivy_vals[i,j]  = dfdy(xGrid[i], yGrid[j])
        end
    end

    # Test that the correct indices and clamp values are returned for a normal case, an on point case, and either out-of-bounds extreme
    @testset "Index determination" begin
        # -------------------------------------------------------------------  iLow  iHigh  clamp
        @test TASOPT.engine.regGridBoundingIndices(0.25, dx, Nx, xGrid[1]) == (3,    4,     0)
        @test TASOPT.engine.regGridBoundingIndices(0.2,  dx, Nx, xGrid[1]) == (3,    3,     0)
        @test TASOPT.engine.regGridBoundingIndices(-0.1, dx, Nx, xGrid[1]) == (1,    2,    -1)
        @test TASOPT.engine.regGridBoundingIndices(1.1,  dx, Nx, xGrid[1]) == (Nx-1, Nx,    1)
    end

    # Test interpolation in an element, at an element vertex, and along each edge
    # Testing at boundaries is handled by the out-of-bounds index code above such that these cases
    # fully cover the space of possibilities
    @testset "Interpolation" begin
        fapprox_inbounds, _, _ = TASOPT.engine.bilinearBoundedLookup(0.13, 0.12, dx, dy, Nx, Ny, 
                                                            xGrid, yGrid, nominal_vals, 
                                                            derivx_vals, derivy_vals)
        @test fapprox_inbounds ≈ 0.00208

        fapprox_onxy, _, _ = TASOPT.engine.bilinearBoundedLookup(0.1, 0.1, dx, dy, Nx, Ny, 
                                                            xGrid, yGrid, nominal_vals, 
                                                            derivx_vals, derivy_vals)
        @test fapprox_onxy ≈ 0.001

        fapprox_onx, _, _ = TASOPT.engine.bilinearBoundedLookup(0.1, 0.12, dx, dy, Nx, Ny, 
                                                            xGrid, yGrid, nominal_vals, 
                                                            derivx_vals, derivy_vals)
        @test fapprox_onx ≈ 0.0016

        fapprox_ony, _, _ = TASOPT.engine.bilinearBoundedLookup(0.13, 0.1, dx, dy, Nx, Ny, 
                                                            xGrid, yGrid, nominal_vals, 
                                                            derivx_vals, derivy_vals)
        @test fapprox_ony ≈ 0.0013
    end

end