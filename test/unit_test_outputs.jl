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
        p1 = TASOPT.stickfig(ac)
        @test p1 isa Plots.Plot
        p2 = TASOPT.plot_details(ac)
        @test p2 isa Plots.Plot
        p3 = TASOPT.plot737compare(ac)
        @test p3 isa Plots.Plot
        p4 = TASOPT.MomentShear(ac)
        @test p4 isa Plots.Plot
        pr_state = (
            parg = copy(ac.parg),
            parm = copy(ac.parm),
            para = copy(ac.para),
            pare = copy(ac.pare),
        )
        # Altitude order cannot change the minimum among feasible candidates.
        template = deepcopy(ac)
        template.parm = repeat(ac.parm[:, 1:1], 1, 2)
        template.para = repeat(ac.para[:, :, 1:1], 1, 1, 2)
        template.pare = repeat(ac.pare[:, :, 1:1], 1, 1, 2)
        altitudes = [28000, 35000, 42000] .* 0.3048
        mission(h; kwargs...) = TASOPT._payload_range_trial(
            template, 500 * 1852.0, ac.parg[igWpaymax], h; kwargs...)
        selected = mission(altitudes)
        reversed = mission(reverse(altitudes))
        @test selected.feasible == reversed.feasible
        @test selected.fuel ≈ reversed.fuel
        @test selected.altitude == reversed.altitude
        for h in altitudes
            single = mission([h])
            if single.feasible
                @test selected.feasible
                @test selected.fuel <= single.fuel + 1e-6
            end
        end
        @test !selected.feasible || TASOPT.is_feasible(selected.constraints, 1e-6)

        # A restrictive policy must reject the previously selected altitude.
        function exclude_selected(trial, im, context)
            report = TASOPT.mission_constraint_report(trial, im, context)
            residual = trial.para[iaalt, ipcruise1, im] == selected.altitude ? 1.0 : -1.0
            return TASOPT.ConstraintReport(merge(report.raw, (; altitude_band=residual)),
                merge(report.scaled, (; altitude_band=residual)), report.values, report.limits)
        end
        restricted = mission(altitudes; mission_constraints=exclude_selected)
        @test !restricted.feasible || restricted.altitude != selected.altitude
        @test !selected.feasible || !restricted.feasible || restricted.fuel >= selected.fuel - 1e-6
        @test_throws ErrorException mission([selected.altitude]; itermax=0)
        @test_throws ArgumentError TASOPT.payload_range_data(ac; altitude_ft=[])
        @test_throws ArgumentError TASOPT.payload_range_data(ac; range_tol=NaN)
        @test_throws ArgumentError TASOPT.payload_range_data(ac; Rpts=1)
        @test_throws ArgumentError TASOPT.payload_range_data(ac; Ppts=0)

        # Static design violations are fatal; they cannot be fixed by off-design flight.
        blocked(args...) = TASOPT.ConstraintReport((; wing_span=1.0), (; wing_span=1.0), (;), (;))
        @test_throws ErrorException TASOPT.payload_range_data(ac; aircraft_constraints=blocked)

        # Independent analytic boundaries: either fuel or a non-fuel limit can bind.
        for payload in (0.0, 1.0, 3.0), fuel in (20.0, 50.0), operating_limit in (3.0, 100.0)
            analytic(r, p) = (; range=r, payload=p,
                feasible=(1 + p) * r <= fuel && r <= operating_limit,
                constraints=(; raw=(;)))
            point = TASOPT._payload_range_limit(analytic, payload, 1.0, 2.0, 0.01)
            @test 0 <= min(fuel / (1 + payload), operating_limit) - point.range <= 0.01
        end
        impossible(r, p) = (; range=r, feasible=false, constraints=(; raw=(;)))
        @test_throws ErrorException TASOPT._payload_range_limit(impossible, 0.0, 1.0, 2.0, 0.01)

        # Sampling preserves corners and respects an independent curved boundary.
        curved(r, p) = (; range=r, payload=p, feasible=r <= 20 / (1 + p),
                           constraints=(; raw=(;)))
        limit(p) = TASOPT._payload_range_limit(curved, p, 1.0, 2.0, 0.01)
        corners = limit.([3.0, 1.0, 0.0])
        origin_check = curved(1.0, 3.0)
        for nr in (2, 5), np in (2, 4)
            samples = TASOPT._payload_range_samples(corners, origin_check, curved, limit, nr, np)
            @test all(c -> c in samples, corners)
            @test length(samples) == nr + 2(np - 1)
            @test all(p -> p.feasible, samples)
            @test issorted(getproperty.(samples, :range))
            @test issorted(getproperty.(samples, :payload); rev=true)
        end
        @test_throws ErrorException TASOPT._payload_range_samples(
            corners, origin_check, impossible, limit, 3, 2)
        @test isequal(ac.parg, pr_state.parg)
        @test isequal(ac.parm, pr_state.parm)
        @test isequal(ac.para, pr_state.para)
        @test isequal(ac.pare, pr_state.pare)

        # Plotting is one panel and identifies a sizing point outside tolerance.
        plot_data = (; range=[0.0, 1e6, 2e6, 3e6], payload=[2e5, 2e5, 1e5, 0.0],
                      Wempty=4e5, design_feasible=false)
        p5 = TASOPT._payload_range_plot(ac, plot_data)
        @test p5 isa Plots.Plot
        @test length(p5.subplots) == 1
        @test length(p5.series_list) == 2
        @test p5.series_list[2][:label] == "Design point (outside tolerance)"
        with_design = TASOPT._payload_range_plot(ac, merge(plot_data, (; design_feasible=true)))
        @test length(with_design.series_list) == 2
        @test with_design.series_list[2][:label] == "Design point"
        with_OEW = TASOPT._payload_range_plot(ac, plot_data; plots_OEW=true)
        for i in 1:2
            @test with_OEW.series_list[i][:x] == p5.series_list[i][:x]
            @test with_OEW.series_list[i][:y] ≈ p5.series_list[i][:y] .+ plot_data.Wempty / (gee * 1000)
        end
        @test_throws ArgumentError TASOPT.PayloadRange(ac; range_tol=NaN)
        p6 = DragPolar(ac)
        @test p6 isa Plots.Plot
        p7 = plot_airf(ac)
        @test p7 isa Plots.Plot
        p8 = plot_drag_breakdown(ac, show_values=true)
        p9 = plot_drag_breakdown(ac, ip = :, show_fractions=false)
        @test p8 isa Plots.Plot
        @test p9 isa Plots.Plot
    end

end
