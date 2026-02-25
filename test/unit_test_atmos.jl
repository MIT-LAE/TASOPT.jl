@testset "Atmosphere" begin
    sea_level = TASOPT.atmos(0.0)
    @test sea_level isa TASOPT.AtmosState
    @test sea_level.T ≈ 288.2
    @test sea_level.p ≈ 1.0132e5
    R = (1.4 - 1.0) * 1004.0/ 1.4
    @test sea_level.ρ ≈ sea_level.p/(R*sea_level.T)
    @test sea_level.a ≈ sqrt(1.4*R*sea_level.T)
    @test sea_level.μ ≈ 1.78e-5 

    cruise = TASOPT.atmos(11_000.0)
    @test cruise.T < sea_level.T
    @test cruise.p < sea_level.p
    @test cruise.ρ < sea_level.ρ

    warm_day = TASOPT.atmos(0.0, 10.0)
    @test warm_day.T ≈ sea_level.T + 10.0
    @test warm_day.ρ < sea_level.ρ

    # Struct returns are intentionally non-iterable; tuple unpacking should fail.
    @test_throws MethodError begin
        T, p, ρ, a, μ = TASOPT.atmos(0.0)
    end

    for h in (0.0, 3_000.0, 8_000.0, 12_000.0)
        ρ = TASOPT.atmos(h).ρ
        @test TASOPT.find_altitude_from_density(ρ) ≈ h atol = 1e-5
    end
end
