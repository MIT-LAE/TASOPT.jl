using Random

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

@testset "Atmosphere Property Tests" begin
    # Inspired by https://fsharpforfunandprofit.com/posts/property-based-testing-2/
    # and https://hypothesis.readthedocs.io/en/latest/tutorial/introduction.html#when-to-use-hypothesis-and-property-based-testing
    # after a conversation with I. Ross.
    rng = MersenneTwister(42)

    @testset "Monotonic pressure and density with altitude" begin
        for _ in 1:10
            h1 = rand(rng) * 20_000.0
            h2 = h1 + 100.0 + rand(rng) * 900.0
            s1 = TASOPT.atmos(h1)
            s2 = TASOPT.atmos(h2)
            @test s2.p < s1.p
            @test s2.ρ < s1.ρ
        end
    end

    @testset "Metamorphic behavior under ΔT" begin
        for _ in 1:10
            h = rand(rng) * 20_000.0
            ΔT = -15.0 + 30.0 * rand(rng)
            std = TASOPT.atmos(h, 0.0)
            shifted = TASOPT.atmos(h, ΔT)
            @test shifted.T ≈ std.T + ΔT
            @test shifted.p ≈ std.p
            @test shifted.ρ * shifted.T ≈ std.ρ * std.T rtol = 1e-12
        end
    end

    @testset "Density-altitude round trip" begin
        for _ in 1:10
            h = rand(rng) * 15_000.0
            ΔT = -10.0 + 20.0 * rand(rng)
            ρ = TASOPT.atmos(h, ΔT).ρ
            h_back = TASOPT.find_altitude_from_density(ρ, ΔT)
            @test h_back ≈ h atol = 1e-6
        end
    end

    @testset "Gas constant invariant" begin
        for _ in 1:10
            h = rand(rng) * 15_000.0
            ΔT = -10.0 + 20.0 * rand(rng)
            s = TASOPT.atmos(h, ΔT)
            R_calc = s.p / (s.ρ * s.T)
            @test R_calc ≈ TASOPT.RSL rtol = 1e-12
        end
    end
end
