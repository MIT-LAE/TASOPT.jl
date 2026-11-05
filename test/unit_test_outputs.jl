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
        # the title shows the governing case; the plot also works without evaluated cases
        function title_text(p)
            join([get(p.attr, :plot_title, ""); [sp[:title] for sp in p.subplots]], " ")
        end
        ac_cases = deepcopy(ac)
        sl = ac_cases.wing.sizing_loads
        sl.results[:pushover] = TASOPT.structures.WingCaseLoads(S_root = 5.0e5, M_root = 6.0e6, S_break = 7.0e4, M_break = 8.0e5)
        sl.results[:taxi] = TASOPT.structures.WingCaseLoads(S_root = 1.0e5, M_root = 2.0e6, S_break = 3.0e4, M_break = 4.0e5)
        sl.governing_root, sl.governing_break = :pushover, :taxi
        p4_cases = TASOPT.MomentShear(ac_cases)
        @test p4_cases isa Plots.Plot
        @test occursin("pushover", title_text(p4_cases))
        @test occursin("taxi", title_text(p4_cases))

        empty!(sl.results)
        @test TASOPT.MomentShear(ac_cases) isa Plots.Plot

        # a cantilever wing keeps the plot as it is: 2 curves, 2 sets of dividing lines, 2 sets of case points, 2 annotations
        annotation_texts(p) = [a[3].str for sp in p.subplots for a in sp[:annotations]]
        @test length(p4.series_list) == 6
        @test length(annotation_texts(p4)) == 2

        # a braced wing: the shear jumps by the strut reaction at ηs, the moment is continuous
        ac_braced = sized_jury_aircraft("strut_braced")
        sl_braced = ac_braced.wing.sizing_loads
        r = sl_braced.results[sl_braced.governing_strut]
        ηs = ac_braced.wing.layout.ηs
        p_braced = TASOPT.MomentShear(ac_braced)
        @test p_braced isa Plots.Plot

        # y-steps of the curves of `subplot` between two consecutive points at the strut attachment
        function steps_at_ηs(subplot)
            steps = Float64[]
            for series in subplot.series_list
                series[:seriestype] == :path || continue
                x, y = collect(series[:x]), collect(series[:y])
                for i in 1:length(x)-1
                    x[i] ≈ ηs && x[i+1] ≈ ηs && push!(steps, y[i+1] - y[i])
                end
            end
            return steps
        end
        shear_steps = steps_at_ηs(p_braced[1])
        moment_steps = steps_at_ηs(p_braced[2])
        @test length(shear_steps) == 1
        @test abs(only(shear_steps)) ≈ r.strut_reaction rtol = 1e-6
        @test all(abs.(moment_steps) .< 1e-6 * abs(r.M_break))
        @test any(t -> occursin(r"\bR\b.*\d.*N", t), annotation_texts(p_braced))

        p5 = TASOPT.PayloadRange(ac)
        @test p5 isa Plots.Plot
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