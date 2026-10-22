using TOML

@testset "Model-defined constraint limits" begin
    mktempdir() do dir
        path = joinpath(dir, "limits.toml")
        ac = load_default_model()
        save_aircraft_model(ac, path)
        data = TOML.parsefile(path)
        prop = data["Propulsion"]
        delete!(prop, "max_fan_diameter")
        delete!(prop, "temperature_limits")
        prop["T_max_metal"] = 1275.0
        write_model() = open(io -> TOML.print(io, data), path, "w")
        write_model()
        legacy = read_aircraft_model(path)
        @test legacy.engine.max_fan_diameter == 2.0
        @test legacy.engine.max_compressor_exit_temperature == 900.0
        @test legacy.parg[igTmetal] == 1275.0

        prop["max_fan_diameter"] = "10 ft"
        prop["temperature_limits"] = Dict("compressor_exit" => "950 K", "turbine_metal" => "1350 K")
        write_model()
        custom = read_aircraft_model(path)
        @test custom.engine.max_fan_diameter ≈ 3.048
        @test custom.engine.max_compressor_exit_temperature == 950.0
        @test custom.parg[igTmetal] == 1350.0
        save_aircraft_model(custom, path)
        reread = read_aircraft_model(path)
        @test reread.engine.max_fan_diameter == custom.engine.max_fan_diameter
        @test reread.engine.max_compressor_exit_temperature == custom.engine.max_compressor_exit_temperature
        @test reread.parg[igTmetal] == custom.parg[igTmetal]

        # Relaxing a bound must reduce its residual without changing performance.
        custom.parg[igdfan] = 3.0
        report = aircraft_constraint_report(custom)
        @test report.limits.fan_diameter == custom.engine.max_fan_diameter
        custom.engine.max_fan_diameter *= 2
        @test aircraft_constraint_report(custom).scaled.fan_diameter < report.scaled.fan_diameter
        custom.pare[ieTt3, :, 1] .= 1000.0
        context = MissionConstraintContext(1, 1e6, false)
        report = mission_constraint_report(custom, 1, context)
        @test report.limits.compressor_exit_temperature == 950.0
        @test report.limits.turbine_metal_temperature == 1350.0
        @test report.scaled.compressor_exit_temperature > 0
        custom.engine.max_compressor_exit_temperature = 1100.0
        relaxed = mission_constraint_report(custom, 1, context)
        @test relaxed.values == report.values
        @test relaxed.scaled.compressor_exit_temperature < 0

        for value in (0.0, -1.0, Inf, NaN)
            prop["max_fan_diameter"] = value
            write_model()
            @test_throws ArgumentError read_aircraft_model(path)
        end
        prop["max_fan_diameter"] = 3.0
        for key in ("compressor_exit", "turbine_metal"), value in (0.0, -1.0, Inf, NaN)
            prop["temperature_limits"] = Dict(key => value)
            write_model()
            @test_throws ArgumentError read_aircraft_model(path)
        end
        prop["temperature_limits"] = Dict("typo" => 900.0)
        write_model()
        @test_throws ArgumentError read_aircraft_model(path)
    end
end
