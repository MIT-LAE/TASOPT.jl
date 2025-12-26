@testset "FlightCondition" verbose=false begin

    @testset "Basic construction from altitude and Mach" begin
        # Test at sea level, Mach 0.3
        fc_sl = FlightCondition(0.0, 0.3)

        @test fc_sl.atmos.alt ≈ 0.0
        @test fc_sl.Mach ≈ 0.3
        @test fc_sl.TAS ≈ fc_sl.Mach * fc_sl.atmos.a
        @test fc_sl.Re_unit ≈ fc_sl.atmos.ρ * fc_sl.TAS / fc_sl.atmos.μ
        @test fc_sl.climb_angle ≈ 0.0

        # Check sea level atmospheric properties (ISA)
        @test fc_sl.atmos.T ≈ 288.2 rtol=0.01  # ~288 K at sea level
        @test fc_sl.atmos.p ≈ 101320.0 rtol=0.01  # ~101 kPa

        # Test at 10 km altitude (typical cruise), Mach 0.8
        fc_cruise = FlightCondition(10000.0, 0.8)

        @test fc_cruise.atmos.alt ≈ 10000.0
        @test fc_cruise.Mach ≈ 0.8
        @test fc_cruise.TAS ≈ fc_cruise.Mach * fc_cruise.atmos.a
        @test altitude(fc_cruise) ≈ 10000.0
        @test altitude_km(fc_cruise) ≈ 10.0

    end

    @testset "Construction with flight path angle" begin
        # Climbing at 3 degrees
        γ = deg2rad(3.0)
        fc_climb = FlightCondition(5000.0, 0.5; climb_angle=γ)

        @test fc_climb.climb_angle ≈ γ
        @test rad2deg(fc_climb.climb_angle) ≈ 3.0
    end

    @testset "Construction with temperature offset" begin
        alt = 10000.0
        Mach = 0.8

        fc_isa = FlightCondition(alt, Mach; ΔT=0.0)
        fc_hot = FlightCondition(alt, Mach; ΔT=15.0)  # ISA + 15K
        fc_cold = FlightCondition(alt, Mach; ΔT=-10.0)  # ISA - 10K

        # Hot day: higher temperature, lower density
        @test fc_hot.atmos.T == fc_isa.atmos.T + 15.0

        # Cold day: lower temperature, higher density
        @test fc_cold.atmos.T == fc_isa.atmos.T - 10.0
    end

    @testset "Derived quantities" begin
        fc = FlightCondition(10000.0, 0.8)

        # Dynamic pressure: q = 0.5 * ρ * V²
        q_expected = 0.5 * fc.atmos.ρ * fc.TAS^2
        @test dynamic_pressure(fc) ≈ q_expected

        # Total temperature: Tt = T * (1 + (γ-1)/2 * M²)
        γ = 1.4
        Tt_expected = fc.atmos.T * (1.0 + (γ - 1.0) / 2.0 * fc.Mach^2)
        @test total_temperature(fc) ≈ Tt_expected

        # Total pressure: pt = p * (1 + (γ-1)/2 * M²)^(γ/(γ-1))
        pt_expected = fc.atmos.p * (1.0 + (γ - 1.0) / 2.0 * fc.Mach^2)^(γ / (γ - 1.0))
        @test total_pressure(fc) ≈ pt_expected
    end

    @testset "Display methods" begin
        fc = FlightCondition(10668.0, 0.785; climb_angle=deg2rad(2.5))

        # Test compact display (single line)
        str_compact = sprint(show, fc)
        @test occursin("FlightCondition", str_compact)
        @test occursin("km", str_compact)
        @test occursin("M=", str_compact)

        str_verbose = sprint(show, MIME"text/plain"(), fc)
        @test occursin("FlightCondition:", str_verbose)
        @test occursin("Altitude:", str_verbose)
        @test occursin("Mach number:", str_verbose)
    end

    @testset "Approximate equality" begin
        fc1 = FlightCondition(10000.0, 0.8)
        fc2 = FlightCondition(10000.0, 0.8)
        fc3 = FlightCondition(10001.0, 0.8)  # Slightly different altitude
        fc4 = FlightCondition(10000.0, 0.81)  # Different Mach

        @test fc1 ≈ fc2
        @test fc1 ≈ fc3 rtol=0.001  # Within 0.1%
        @test !(fc1 ≈ fc4)  # Different Mach
    end
end
