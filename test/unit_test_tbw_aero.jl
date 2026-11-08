const CoreLogging = Base.CoreLogging
# Tests for braced-wing aerodynamics: N-surface Trefftz plane (wing, tail, strut), strut lift and drag,
# infeasibility flags and aero/geometry output. Cantilever results must match the two-surface implementation.

const trefftz_aero = TASOPT.aerodynamics

"Load an aircraft model from a TOML file without printing log messages."
function load_aircraft_quietly(path = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"))
    return CoreLogging.with_logger(CoreLogging.NullLogger()) do
        read_aircraft_model(path)
    end
end

load_braced_aircraft_quietly() =
    load_aircraft_quietly(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw.toml"))

"Fixed wing/tail geometry and loading (as in `unit_test_trefftz.jl`) on an unsized aircraft; returns the `para` read by `induced_drag!`."
function set_trefftz_test_case!(ac)
    wing, htail = ac.wing, ac.htail
    wing.layout.S = 125.0
    wing.layout.span = 35.486921631434697
    wing.layout.root_span = 3.6068
    wing.layout.ηs = 10.113772664958887 / wing.layout.span
    wing.inboard.λ = 1.0
    wing.outboard.λ = 1.0

    htail.layout.S = 25.0
    htail.layout.span = 15.958117796995291
    htail.layout.root_span = 1.524
    htail.layout.ηs = 1.524 / htail.layout.span
    htail.outboard.λ = 0.25

    para = zeros(size(ac.para, 1))
    para[iaCL] = 1.2
    para[iaCLh] = 0.05 # referenced to the tail area
    para[iarcls] = 0.77
    para[iarclt] = 0.15
    return para
end

"Run `induced_drag!` for the fixed test case and return `(CDi, span efficiency)`; the cached wake geometry is discarded unless `reuse_cached_geometry`."
function trefftz_drag_for_test_case!(ac; reuse_cached_geometry = false)
    para = set_trefftz_test_case!(ac)
    reuse_cached_geometry || (trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0))
    trefftz_aero.induced_drag!(para, ac, ac.options.trefftz_config)
    return para[iaCDi], para[iaspaneff]
end

"Set the strut inputs that determine the strut lift (area, velocity ratio, `lift_fraction_strut_to_total`, `spanload_exponent`)."
function set_strut_lift_inputs!(ac, lift_fraction_strut_to_total; spanload_exponent = 0.5)
    ac.wing.strut.S = 10.0
    ac.wing.strut.local_velocity_ratio = 1.0
    ac.wing.strut.lift_fraction_strut_to_total = lift_fraction_strut_to_total
    ac.wing.strut.spanload_exponent = spanload_exponent
    return ac
end

"Lift coefficient (referenced to the wing area) of the Trefftz points `i_first:i_last`, integrated from the bound circulation."
function trefftz_surface_lift(ac, i_first, i_last)
    geom = trefftz_aero.TREFFTZ_GEOM
    panels = i_first:(i_last - 1)
    return 2 * ac.wing.layout.span / (0.5 * ac.wing.layout.S) *
           sum(geom.gc[i] * (geom.yp[i + 1] - geom.yp[i]) for i in panels)
end

"""
Set the strut and jury strut drag inputs (and what `aircraft_drag!` needs) on an unsized aircraft; returns the `para` slice of
mission 1, point `ip`. The strut carries 20% of the wing and strut lift, which puts its section lift coefficient inside the
airfoil database lift range, up to a share of 30%.
"""
function set_strut_drag_test_case!(ac; ip = 1)
    para = view(ac.para, :, ip, 1)
    lift_inputs = set_trefftz_test_case!(ac)
    for i in (iaCL, iaCLh, iarcls, iarclt)
        para[i] = lift_inputs[i]
    end
    para[iaMach] = 0.70
    para[iaReunit] = 7.0e6

    ac.wing.layout.root_chord = 5.0
    ac.htail.layout.root_chord = 2.0
    ac.vtail.layout.root_chord = 2.0
    ac.vtail.layout.span = 5.0
    ac.vtail.layout.root_span = 1.0
    ac.vtail.outboard.λ = 0.3

    set_strut_lift_inputs!(ac, 0.20)
    ac.wing.strut.S = 60.0
    ac.wing.strut.dihedral = 0.3
    ac.wing.strut.chord = 0.8
    ac.wing.strut.root_chord = 0.8 # no taper
    ac.wing.strut.cos_lambda = 0.9
    ac.wing.strut.local_velocity_ratio = 1.1
    ac.wing.jury.chord = 0.4
    ac.wing.jury.length = 2.5
    return para
end

"Run `aircraft_drag!` for mission 1, point `ip` and return the `para` slice."
function strut_drag_for_test_case!(ac; ip = 1)
    para = view(ac.para, :, ip, 1)
    trefftz_aero.aircraft_drag!(ac, 1, ip, true)
    return para
end

"""
Streamwise section lift coefficient of the main strut: the lift of both struts (`lift_fraction_strut_to_total` of the wing and
strut lift) over the effective area `S cos²(dihedral)`, at the local velocity ratio.
"""
function expected_strut_section_cl(ac, para)
    wing, strut = ac.wing, ac.wing.strut
    CLhtail = para[iaCLh] * ac.htail.layout.S / wing.layout.S
    CL_strut = strut.lift_fraction_strut_to_total * (para[iaCL] - CLhtail)
    return CL_strut * wing.layout.S / (strut.S * cos(strut.dihedral)^2 * strut.local_velocity_ratio^2)
end

"Airfoil database section drag of the main strut (`strut_section_cd` result); the Reynolds number uses the mean of the end chords."
function expected_main_strut_section(ac, para)
    strut = ac.wing.strut
    rV = strut.local_velocity_ratio
    mean_chord = (strut.root_chord + strut.chord) / 2
    return trefftz_aero.strut_section_cd(expected_strut_section_cl(ac, para), strut.thickness_to_chord, para[iaMach] * rV,
        strut.cos_lambda, para[iaReunit] * rV * mean_chord, para[iaaRexp], para[iafexcdw], ac.wing.airsection)
end

"Airfoil database section drag of the jury struts (vertical, unswept, no lift), looked up at the strut local velocity ratio."
function expected_jury_section(ac, para)
    jury = ac.wing.jury
    rV = ac.wing.strut.local_velocity_ratio
    return trefftz_aero.strut_section_cd(0.0, jury.thickness_to_chord, para[iaMach] * rV, 1.0,
        para[iaReunit] * rV * jury.chord, para[iaaRexp], para[iafexcdw], ac.wing.airsection)
end

"Main strut profile drag from its definition: area ratio times the section drag times the cubed local velocity ratio."
function expected_main_strut_drag(ac, para)
    strut = ac.wing.strut
    return strut.S / ac.wing.layout.S * expected_main_strut_section(ac, para).cd * strut.local_velocity_ratio^3
end

"Profile drag of the two jury struts from its definition."
function expected_jury_drag(ac, para)
    jury = ac.wing.jury
    return 2 * jury.chord * jury.length / ac.wing.layout.S * expected_jury_section(ac, para).cd *
           ac.wing.strut.local_velocity_ratio^3
end

@testset "Trefftz surfaces" verbose=true begin
    @testset "Configuration without strut panels" begin
        # MEDIUM: wing 40+12+6 panels, tail 20+0+4 panels, one dummy point after each surface
        config = trefftz_aero.get_trefftz_config("MEDIUM")
        @test trefftz_aero.n_surfaces(config) == 2
        @test (trefftz_aero.i_last_wing(config), trefftz_aero.i_first_tail(config), trefftz_aero.n_points_used(config)) == (59, 60, 84)
    end

    @testset "Configuration with strut panels" begin
        # quality => (strut panels, points without strut, points with strut)
        # The strut has as many panels as the wing has outer panels
        cases = ("COARSE" => (20, 43, 64),
                 "MEDIUM" => (40, 84, 125),
                 "FINE"   => (160, 330, 491))
        for (quality, (n_strut_panels, n_points_cantilever, n_points_braced)) in cases
            @testset "$quality" begin
                config = trefftz_aero.get_trefftz_config(quality; braced = true)
                @test trefftz_aero.n_surfaces(config) == 3
                @test config.strut_panels.n_outer_panels == config.wing_panels.n_outer_panels
                @test config.strut_panels == trefftz_aero.SurfaceDiscretization(n_strut_panels, 0, 0)

                # Surface blocks are contiguous and the strut block ends the point list
                @test trefftz_aero.i_first_tail(config) == trefftz_aero.i_last_wing(config) + 1
                @test trefftz_aero.i_first_strut(config) == trefftz_aero.i_last_tail(config) + 1
                @test trefftz_aero.i_last_strut(config) == trefftz_aero.i_first_strut(config) + n_strut_panels
                @test trefftz_aero.n_points_used(config) == n_points_braced
                @test trefftz_aero.i_last_tail(config) == n_points_cantilever
                @test trefftz_aero.n_strut_points(config.strut_panels) == n_strut_panels + 1
            end
        end
    end

    @testset "Geometry array overflow includes the strut block" begin
        n_slots = length(trefftz_aero.TREFFTZ_GEOM.y)
        n_strut_panels_that_fit = n_slots - 330 - 1 # one more point is the dummy point after the strut
        @test n_slots >= 491 # FINE with a strut
        wing_panels = trefftz_aero.SurfaceDiscretization(160, 48, 24)
        tail_panels = trefftz_aero.SurfaceDiscretization(80, 0, 16)
        config_with_strut(n_strut_panels) = trefftz_aero.TrefftzPlaneConfig(wing_panels, tail_panels;
            strut_panels = trefftz_aero.SurfaceDiscretization(n_strut_panels, 0, 0))

        @test trefftz_aero.n_points_used(config_with_strut(n_strut_panels_that_fit)) == n_slots # just fits
        @test_throws ErrorException config_with_strut(n_strut_panels_that_fit + 1)
    end

    @testset "calculate_array_size includes the strut block" begin
        ac = load_aircraft_quietly()
        htail = ac.htail
        htail.layout.root_span = 1.524
        cantilever_config = trefftz_aero.get_trefftz_config("COARSE")
        braced_config = trefftz_aero.get_trefftz_config("COARSE"; braced = true)

        @test trefftz_aero.calculate_array_size(cantilever_config, htail) == 43
        @test trefftz_aero.calculate_array_size(braced_config, htail) == 64

        # A T-tail drops the tail image panels, not the strut block
        htail.layout.root_span = 0.0
        @test trefftz_aero.calculate_array_size(cantilever_config, htail) == 41
        @test trefftz_aero.calculate_array_size(braced_config, htail) == 62
    end

    @testset "generate_trefftz_points_inclined!" begin
        # NaN-filled arrays detect any entry the function must not touch
        geom = trefftz_aero.TrefftzGeometry{40}()
        for name in (:t, :y, :yp, :z, :zp, :yc, :ycp, :zc, :zcp, :gc)
            fill!(getfield(geom, name), NaN)
        end
        i_start, n_panels = 5, 8
        y_start, z_start, y_end, z_end = 0.5, -1.4, 21.3, 1.7
        i_end = trefftz_aero.generate_trefftz_points_inclined!(i_start, geom, y_start, z_start, y_end, z_end, n_panels)

        @test i_end == i_start + n_panels
        nodes = i_start:i_end
        controls = i_start:(i_end - 1)
        @test (geom.y[i_start], geom.z[i_start]) == (y_start, z_start)
        @test geom.y[i_end] ≈ y_end rtol = 1e-14
        @test geom.z[i_end] ≈ z_end rtol = 1e-14

        # Evenly spaced and collinear, no wake contraction
        @test all(isapprox.(diff(geom.y[nodes]), (y_end - y_start) / n_panels; rtol = 1e-12))
        slope = (z_end - z_start) / (y_end - y_start)
        @test all(isapprox.(geom.z[nodes], z_start .+ slope .* (geom.y[nodes] .- y_start); atol = 1e-12))
        @test geom.yp[nodes] == geom.y[nodes]
        @test geom.zp[nodes] == geom.z[nodes]

        # Control points are the midpoints of consecutive nodes; no circulation yet
        @test geom.yc[controls] ≈ 0.5 .* (geom.y[i_start:(i_end - 1)] .+ geom.y[(i_start + 1):i_end])
        @test geom.zc[controls] ≈ 0.5 .* (geom.z[i_start:(i_end - 1)] .+ geom.z[(i_start + 1):i_end])
        @test all(==(0.0), geom.gc[controls])

        # Nothing outside the strut is written
        @test all(isnan, geom.y[1:(i_start - 1)]) && all(isnan, geom.y[(i_end + 1):end])
        @test all(isnan, geom.yc[i_end:end])
    end

    @testset "Two-element strut: influence coefficients, circulation and induced drag from hand calculation" begin
        # Strut of two panels along z = y from (0, 0) to (2, 2): elements of length sqrt(2), unit normal (-1, 1)/sqrt(2)
        geom = trefftz_aero.TrefftzGeometry{40}()
        i_end = trefftz_aero.generate_trefftz_points_inclined!(1, geom, 0.0, 0.0, 2.0, 2.0, 2)
        @test i_end == 3
        wake_system = trefftz_aero.WakeSystem(geom.yp[1:i_end], geom.zp[1:i_end], geom.ycp[1:(i_end - 1)], geom.zcp[1:(i_end - 1)])
        for element in wake_system.elements
            @test element.length ≈ sqrt(2) rtol = 1e-14
            @test element.unit_normal ≈ [-1.0, 1.0] / sqrt(2) rtol = 1e-14
        end

        # Influence of a wake point on a control point is (r_y + r_z) / (sqrt(2) |r|²) minus the image point (-y, z) term;
        # the wake point (0, 0) coincides with its image and cancels.
        s = sqrt(2)
        influence = [0.0  -6s / 5   -20s / 51;
                     0.0  10s / 13  -28s / 25]
        @test wake_system.influence_matrix ≈ influence rtol = 1e-13

        gc = [1.0, 0.5]
        gw = zeros(3)
        trefftz_aero.calculate_wake_circulation!(gw, [gc; 0.0], [1], [3], 1)
        @test gw == [0.0, 0.5, 0.5]

        # w = bref/(2π) A gw, lift = 2 Σ gc Δy, drag = -Σ gc w Δs, both times bref/(Sref/2)
        bref, Sref = 2.0, 3.0
        w = bref / (2π) * influence * gw
        normalisation = bref / (0.5 * Sref)
        CL_expected = 2 * (gc[1] * 1.0 + gc[2] * 1.0) * normalisation
        CD_expected = -(gc[1] * w[1] + gc[2] * w[2]) * s * normalisation
        CL, CD, span_efficiency = trefftz_aero.compute_induced_drag!(zeros(2), wake_system, gc, gw, bref, Sref)
        @test CL ≈ CL_expected rtol = 1e-13
        @test CD ≈ CD_expected rtol = 1e-13
        @test span_efficiency ≈ CL_expected^2 / (π * bref^2 / Sref) / CD_expected rtol = 1e-13
    end

    @testset "Cantilever induced drag is unchanged and ignores the strut lift input" begin
        ac = load_aircraft_quietly()
        @test trefftz_aero.n_surfaces(ac.options.trefftz_config) == 2
        set_strut_lift_inputs!(ac, 0.5)

        CDi, span_efficiency = trefftz_drag_for_test_case!(ac)

        # Reference values from the two-surface implementation (before strut panels existed)
        @test CDi ≈ 0.05468605653459833 rtol = 1e-12
        @test span_efficiency ≈ 0.8319732645252401 rtol = 1e-12
        @test length(ac.wake_system.points) == 43
        @test ac.wing.strut.CL == 0.0
        @test ac.wing.strut.lift_fraction == 0.0
    end

    @testset "Braced aircraft gets a strut surface in the wake system" begin
        ac_cantilever = load_aircraft_quietly()
        ac_braced = load_braced_aircraft_quietly()
        braced_config = ac_braced.options.trefftz_config
        @test trefftz_aero.n_surfaces(braced_config) == 3

        CDi_cantilever, span_efficiency_cantilever = trefftz_drag_for_test_case!(ac_cantilever)
        CDi_braced, span_efficiency_braced = trefftz_drag_for_test_case!(ac_braced)

        @test length(ac_braced.wake_system.points) == 64

        # A strut without lift has zero circulation and cannot induce drag on the wing or tail
        @test CDi_braced ≈ CDi_cantilever rtol = 1e-8
        @test span_efficiency_braced ≈ span_efficiency_cantilever rtol = 1e-8

        i_strut = trefftz_aero.i_first_strut(braced_config):(trefftz_aero.i_last_strut(braced_config) - 1)
        @test all(==(0.0), trefftz_aero.TREFFTZ_GEOM.gc[i_strut])

        # Switching surface count on unchanged geometry must rebuild the cached wake geometry
        CDi_again, _ = trefftz_drag_for_test_case!(ac_cantilever; reuse_cached_geometry = true)
        @test length(ac_cantilever.wake_system.points) == 43
        @test CDi_again ≈ CDi_cantilever rtol = 1e-12
    end

    @testset "Strut Trefftz line joins the structural strut at the wing root and the wing break" begin
        ac = load_braced_aircraft_quietly()
        config = ac.options.trefftz_config
        wing = ac.wing
        para = set_trefftz_test_case!(ac)
        wing.layout.ηs = 0.625 # set after the test case, which sets the wing break
        wing.strut.z = 3.0
        trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
        trefftz_aero.induced_drag!(para, ac, config)

        geom = trefftz_aero.TREFFTZ_GEOM
        i_first, i_last = trefftz_aero.i_first_strut(config), trefftz_aero.i_last_strut(config)
        y_root, z_root = 0.5 * wing.layout.root_span, wing.layout.z - wing.strut.z
        y_break, z_break = 0.5 * wing.layout.break_span, wing.layout.z
        @test geom.y[i_first] ≈ y_root rtol = 1e-14
        @test geom.z[i_first] ≈ z_root rtol = 1e-14
        @test geom.y[i_last] ≈ y_break rtol = 1e-14
        @test geom.z[i_last] ≈ z_break rtol = 1e-14

        # Same angle as the structural strut: dihedral atan(z / (b/2 (ηs - ηo)))
        half_bay = 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
        @test y_break - y_root ≈ half_bay rtol = 1e-14
        @test atan(geom.z[i_last] - geom.z[i_first], geom.y[i_last] - geom.y[i_first]) ≈ atan(wing.strut.z, half_bay) rtol = 1e-12
    end

    @testset "Strut tip meets the wing break in the Trefftz plane" begin
        ac = load_braced_aircraft_quietly()
        config = ac.options.trefftz_config
        para = set_trefftz_test_case!(ac)
        ac.wing.layout.ηs = 0.625 # set after the test case, which sets the wing break
        trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
        trefftz_aero.induced_drag!(para, ac, config)

        geom = trefftz_aero.TREFFTZ_GEOM
        i_last_strut = trefftz_aero.i_last_strut(config)
        # The wing node at the panel break follows the inboard panels after the image panels
        i_break = trefftz_aero.i_first_wing(config) + config.wing_panels.n_image_panels + config.wing_panels.n_inner_panels
        @test geom.y[i_break] ≈ 0.5 * ac.wing.layout.break_span rtol = 1e-14
        @test geom.yp[i_last_strut] ≈ geom.yp[i_break] rtol = 1e-10
        @test geom.zp[i_last_strut] ≈ geom.zp[i_break] rtol = 1e-10
    end

    @testset "Strut induced drag converges with strut panels" begin
        function braced_cdi_with_strut_panels(n_strut_panels)
            ac = braced_aircraft_with_strut_panels(n_strut_panels)
            set_strut_lift_inputs!(ac, 0.2; spanload_exponent = 0.5)
            para = set_trefftz_test_case!(ac)
            ac.wing.layout.ηs = 0.625 # set after the test case, which sets the wing break
            ac.wing.strut.z = 3.0
            trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
            trefftz_aero.induced_drag!(para, ac, ac.options.trefftz_config)
            return para[iaCDi]
        end
        CDi = Dict(n => braced_cdi_with_strut_panels(n) for n in (40, 80, 160))
        @test abs(CDi[80] / CDi[160] - 1) < 5e-3
        @test abs(CDi[40] / CDi[160] - 1) < 2e-2
    end

    @testset "Wake geometry is rebuilt only when the geometry changes" begin
        wake_points(ac, i_points) = [ac.wake_system.points[i] for i in i_points]
        all_at_height(points, z) = all(point -> isapprox(point[2], z; atol = 1e-12), points)

        @testset "Braced aircraft: strut and wing height" begin
            ac = load_braced_aircraft_quietly()
            config = ac.options.trefftz_config
            wing = ac.wing
            para = set_trefftz_test_case!(ac)
            i_wing = trefftz_aero.i_first_wing(config):trefftz_aero.i_last_wing(config)
            i_strut = trefftz_aero.i_first_strut(config):trefftz_aero.i_last_strut(config)

            trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
            trefftz_aero.induced_drag!(para, ac, config)

            # Unchanged geometry keeps the wake system
            wake_system = ac.wake_system
            trefftz_aero.induced_drag!(para, ac, config)
            @test ac.wake_system === wake_system

            changes = (
                "strut.z"   => () -> wing.strut.z += 0.1,
                "layout.z"  => () -> wing.layout.z += 0.3,
                "layout.ηs" => () -> wing.layout.ηs += 0.05)
            for (name, change!) in changes
                @testset "$name" begin
                    wake_system = ac.wake_system
                    strut_before = wake_points(ac, i_strut)
                    wing_before = wake_points(ac, i_wing)

                    change!()
                    trefftz_aero.induced_drag!(para, ac, config)
                    strut_after = wake_points(ac, i_strut)
                    wing_after = wake_points(ac, i_wing)

                    @test ac.wake_system !== wake_system
                    yo = 0.5 * wing.layout.root_span
                    yop = yo * config.wing_root_contraction
                    y_break_wake = trefftz_aero.get_wake_contraction(trefftz_aero.WINGWAKE, 0.5 * wing.layout.break_span, yo, yop)
                    @test all(isapprox.(strut_after[1], (yop, wing.layout.z - wing.strut.z); rtol = 1e-14))
                    @test all(isapprox.(strut_after[end], (y_break_wake, wing.layout.z); rtol = 1e-14))
                    @test all_at_height(wing_after, wing.layout.z)
                    @test strut_after != strut_before
                    name == "layout.z" && @test wing_after != wing_before
                end
            end
        end

        @testset "Cantilever aircraft: tail height" begin
            ac = load_aircraft_quietly()
            config = ac.options.trefftz_config
            para = set_trefftz_test_case!(ac)
            i_tail = trefftz_aero.i_first_tail(config):trefftz_aero.i_last_tail(config)

            trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
            trefftz_aero.induced_drag!(para, ac, config)
            wake_system, CDi_before = ac.wake_system, para[iaCDi]

            ac.htail.layout.z += 0.4
            trefftz_aero.induced_drag!(para, ac, config)
            @test ac.wake_system !== wake_system
            @test all_at_height(wake_points(ac, i_tail), ac.htail.layout.z)
            @test !isapprox(para[iaCDi], CDi_before; rtol = 1e-4)

            # The updated wake must give the result of a wake built from scratch
            CDi_updated = para[iaCDi]
            trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
            trefftz_aero.induced_drag!(para, ac, config)
            @test CDi_updated ≈ para[iaCDi] rtol = 1e-12

            # Moving wing and tail by the same height leaves the induced drag unchanged
            ac.htail.layout.z -= 0.4
            CDi, _ = trefftz_drag_for_test_case!(ac)
            ac.wing.layout.z += 0.4
            ac.htail.layout.z += 0.4
            CDi_shifted, _ = trefftz_drag_for_test_case!(ac; reuse_cached_geometry = true)
            @test CDi_shifted ≈ CDi rtol = 1e-12
        end
    end

    @testset "Strut lift inputs are read from [Wing.Strut]" begin
        path = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw.toml")
        toml = replace(read(path, String),
            r"(?m)^(\s*lift_fraction_strut_to_total\s*=\s*)[^\s#]+" => s"\g<1>0.4",
            r"(?m)^(\s*spanload_exponent\s*=\s*)[^\s#]+" => s"\g<1>0.8")
        mktempdir() do dir
            modified_path = joinpath(dir, "strut_lift.toml")
            write(modified_path, toml)
            strut = load_aircraft_quietly(modified_path).wing.strut
            @test strut.lift_fraction_strut_to_total == 0.4
            @test strut.spanload_exponent == 0.8
        end
    end

    @testset "Strut carries lift in addition to the wing" begin
        CL, CLhtail = 1.2, 0.05 * 25.0 / 125.0 # as in set_trefftz_test_case!, referenced to the wing area
        lift_fraction_strut_to_total, spanload_exponent = 0.1, 0.5

        ac_no_lift = load_braced_aircraft_quietly()
        set_strut_lift_inputs!(ac_no_lift, 0.0)
        CDi_no_lift, _ = trefftz_drag_for_test_case!(ac_no_lift)

        ac = load_braced_aircraft_quietly()
        set_strut_lift_inputs!(ac, lift_fraction_strut_to_total; spanload_exponent)
        CDi, span_efficiency = trefftz_drag_for_test_case!(ac)
        config = ac.options.trefftz_config
        strut = ac.wing.strut
        geom = trefftz_aero.TREFFTZ_GEOM
        i_first, i_last = trefftz_aero.i_first_strut(config), trefftz_aero.i_last_strut(config)

        @testset "Strut lift coefficient" begin
            CL_strut = lift_fraction_strut_to_total * (CL - CLhtail)
            @test CL_strut > 0.0
            @test strut.CL ≈ CL_strut rtol = 1e-12
            @test strut.lift_fraction ≈ CL_strut / CL rtol = 1e-12

            # The strut takes lift from the wing; the total stays CL
            strut_lift = trefftz_surface_lift(ac, i_first, i_last)
            wing_lift = trefftz_surface_lift(ac, trefftz_aero.i_first_wing(config), trefftz_aero.i_last_wing(config))
            tail_lift = trefftz_surface_lift(ac, trefftz_aero.i_first_tail(config), trefftz_aero.i_last_tail(config))
            @test strut_lift ≈ CL_strut rtol = 1e-10
            @test tail_lift ≈ CLhtail rtol = 1e-10
            @test wing_lift ≈ CL - CLhtail - CL_strut rtol = 1e-10
            @test wing_lift + tail_lift + strut_lift ≈ CL rtol = 1e-10

            # Stored in para, and cleared at zero lift
            para = set_trefftz_test_case!(ac)
            trefftz_aero.induced_drag!(para, ac, config)
            @test para[iaCLstrut] == strut.CL > 0.0
            para[iaCL] = 0.0
            trefftz_aero.induced_drag!(para, ac, config)
            @test (strut.CL, strut.lift_fraction, para[iaCLstrut]) == (0.0, 0.0, 0.0)
        end

        @testset "Strut circulation shape" begin
            n = trefftz_aero.n_total(config.strut_panels)
            gc = geom.gc[i_first:(i_last - 1)]
            # ξ = 0 at the fuselage, ξ = 1 at the wing attachment
            ξ = [(k - 0.5) / n for k in 1:n]
            shape = (4 .* ξ .* (1 .- ξ)) .^ spanload_exponent

            # Loading follows (4ξ(1-ξ))^p over every panel, pinned to zero at both ends of the strut
            @test all(>(0.0), gc)
            @test gc ./ gc[1] ≈ shape ./ shape[1] rtol = 1e-10

            # The strut sheds its root circulation at the fuselage end and its tip circulation at the wing end;
            # the wake circulation of the strut sums to zero
            gw = trefftz_aero.gw[i_first:i_last]
            @test gw[1] ≈ -gc[1] rtol = 1e-10
            @test gw[end] ≈ gc[end] rtol = 1e-10
            @test sum(gw) ≈ 0.0 atol = 1e-12 * gc[1]
        end

        @testset "Induced drag" begin
            @test isfinite(CDi) && CDi > 0.0 && span_efficiency > 0.0
            # The strut takes about 10 % of the lift, which changes the induced drag by several percent
            @test abs(CDi - CDi_no_lift) / CDi_no_lift > 1e-3
        end
    end

    @testset "Strut circulation shape on every strut panel" begin
        geom = trefftz_aero.TrefftzGeometry{40}()
        i_first, n_panels, p = 3, 6, 0.5
        i_last = i_first + n_panels
        trefftz_aero.generate_trefftz_points_inclined!(i_first, geom, 1.0, 1.0, 5.0, 2.0, n_panels)
        trefftz_aero.set_strut_circulation_shape!(geom, i_first, i_last, p)

        ξ = [(k - 0.5) / n_panels for k in 1:n_panels]
        gc = geom.gc[i_first:(i_last - 1)]
        @test all(>(0.0), gc)
        @test gc ≈ (4 .* ξ .* (1 .- ξ)) .^ p rtol = 1e-14
    end

    @testset "Strut circulation is rebuilt when its inputs change on cached geometry" begin
        function strut_circulation(ac)
            config = ac.options.trefftz_config
            return copy(trefftz_aero.TREFFTZ_GEOM.gc[trefftz_aero.i_first_strut(config):(trefftz_aero.i_last_strut(config) - 1)])
        end

        ac = load_braced_aircraft_quietly()
        set_strut_lift_inputs!(ac, 0.5; spanload_exponent = 0.5)
        trefftz_drag_for_test_case!(ac)
        gc_before = strut_circulation(ac)

        ac.wing.strut.spanload_exponent = 1.0
        CDi_updated, _ = trefftz_drag_for_test_case!(ac; reuse_cached_geometry = true)
        gc_updated = strut_circulation(ac)
        @test !isapprox(gc_updated, gc_before; rtol = 1e-3)

        # Same as a wake built from scratch
        CDi_fresh, _ = trefftz_drag_for_test_case!(ac)
        @test CDi_updated ≈ CDi_fresh rtol = 1e-12
        @test gc_updated ≈ strut_circulation(ac) rtol = 1e-12
    end

    @testset "Strut and jury drag from the airfoil database" begin
        function strut_drag_case(config; ip = 1)
            ac = load_braced_aircraft_quietly()
            ac.options.opt_wing_config = config
            para = set_strut_drag_test_case!(ac; ip)
            return ac, para
        end

        @testset "Main strut drag is the airfoil database section drag at the strut lift, thickness, sweep and Reynolds number" begin
            ac, para = strut_drag_case(WingConfig.StrutBraced)
            strut, air = ac.wing.strut, ac.wing.airsection
            rV = strut.local_velocity_ratio

            # Both struts carry 20% of the wing and strut lift, over the area S cos²(dihedral), at the local velocity
            cl_strut = 0.20 * (1.2 - 0.05 * 25.0 / 125.0) * 125.0 / (60.0 * cos(0.3)^2 * rV^2)
            @test expected_strut_section_cl(ac, para) ≈ cl_strut rtol = 1e-14
            section = trefftz_aero.strut_section_cd(cl_strut, strut.thickness_to_chord, para[iaMach] * rV, strut.cos_lambda,
                para[iaReunit] * rV * strut.chord, para[iaaRexp], para[iafexcdw], air)
            @test isempty(section.exceeded)
            @test air.cl[1] <= cl_strut / strut.cos_lambda^2 <= air.cl[end]
            @test para[iafexcdw] != 1.0

            CDstrut_untapered = strut_drag_for_test_case!(ac)[iaCDstrut]
            @test CDstrut_untapered ≈ strut.S / ac.wing.layout.S * section.cd * rV^3 rtol = 1e-12
            @test CDstrut_untapered ≈ expected_main_strut_drag(ac, para) rtol = 1e-12

            # 0.8 m at the wing end, 1.6 m at the fuselage end: the Reynolds number uses the mean chord 1.2 m
            strut.taper_ratio = 0.5
            strut.root_chord = strut.chord / strut.taper_ratio
            mean_chord = 1.2
            section_tapered = trefftz_aero.strut_section_cd(cl_strut, strut.thickness_to_chord, para[iaMach] * rV,
                strut.cos_lambda, para[iaReunit] * rV * mean_chord, para[iaaRexp], para[iafexcdw], air)
            CDstrut = strut_drag_for_test_case!(ac)[iaCDstrut]
            @test CDstrut ≈ strut.S / ac.wing.layout.S * section_tapered.cd * rV^3 rtol = 1e-12
            @test CDstrut ≈ CDstrut_untapered * (mean_chord / 0.8)^para[iaaRexp] rtol = 1e-12
        end

        @testset "Truss-braced wing adds the jury strut drag to the strut drag; a strut-braced wing has no jury term" begin
            ac, para = strut_drag_case(WingConfig.TrussBraced)
            expected_strut = expected_main_strut_drag(ac, para)
            expected_jury = expected_jury_drag(ac, para)
            @test expected_jury > 0

            # The jury carries no lift: it is looked up at the lowest database lift coefficient, without the lift penalty
            jury, air = ac.wing.jury, ac.wing.airsection
            jury_section = expected_jury_section(ac, para)
            @test jury_section.exceeded == [:cl]
            at_floor = trefftz_aero.strut_section_cd(air.cl[1], jury.thickness_to_chord, para[iaMach] * ac.wing.strut.local_velocity_ratio,
                1.0, para[iaReunit] * ac.wing.strut.local_velocity_ratio * jury.chord, para[iaaRexp], para[iafexcdw], air)
            @test jury_section.cd == at_floor.cd

            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] ≈ expected_strut + expected_jury rtol = 1e-12

            # The jury drag is part of the total drag
            ac_strut, para_strut = strut_drag_case(WingConfig.StrutBraced)
            strut_drag_for_test_case!(ac_strut)
            @test para_strut[iaCDstrut] ≈ expected_strut rtol = 1e-12
            @test para[iaCD] - para_strut[iaCD] ≈ expected_jury rtol = 1e-9
        end

        @testset "Strut drag is part of the total drag" begin
            ac, _ = strut_drag_case(WingConfig.TrussBraced)
            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] > 0
            components = (iaCDi, iaCDfuse, iaCDwing, iaCDover, iaCDhtail, iaCDvtail, iaCDstrut, iaCDnace, iadCDBLIf, iadCDBLIw)
            @test para[iaCD] ≈ sum(para[i] for i in components) rtol = 1e-14
        end

        @testset "Strut drag with stored wing profile drag (computes_wing_direct = false), as in the mission iteration" begin
            ac, para = strut_drag_case(WingConfig.StrutBraced)
            CDstrut_direct = strut_drag_for_test_case!(ac)[iaCDstrut]

            trefftz_aero.aircraft_drag!(ac, 1, 1, false)
            @test para[iaCDstrut] > 0
            @test para[iaCDstrut] ≈ expected_main_strut_drag(ac, para) rtol = 1e-12
            @test para[iaCDstrut] ≈ CDstrut_direct rtol = 1e-12
        end

        @testset "Lift dependence enters through the database lookup" begin
            ac, para = strut_drag_case(WingConfig.StrutBraced)
            air = ac.wing.airsection
            CDstrut_by_share = Float64[]
            for lift_fraction_strut_to_total in (0.20, 0.30)
                ac.wing.strut.lift_fraction_strut_to_total = lift_fraction_strut_to_total
                cl_perp = expected_strut_section_cl(ac, para) / ac.wing.strut.cos_lambda^2
                @test air.cl[1] <= cl_perp <= air.cl[end]
                @test isempty(expected_main_strut_section(ac, para).exceeded)

                CDstrut = strut_drag_for_test_case!(ac)[iaCDstrut]
                @test CDstrut ≈ expected_main_strut_drag(ac, para) rtol = 1e-12
                push!(CDstrut_by_share, CDstrut)
            end
            @test CDstrut_by_share[2] > CDstrut_by_share[1]
        end

        @testset "Unsized strut or jury gives zero drag, not NaN" begin
            ac, para = strut_drag_case(WingConfig.TrussBraced)
            expected_strut = expected_main_strut_drag(ac, para)
            expected_jury = expected_jury_drag(ac, para)

            ac.wing.jury.chord = 0.0
            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] ≈ expected_strut rtol = 1e-12
            @test isfinite(para[iaCD])

            ac.wing.jury.chord = 0.4
            ac.wing.jury.length = 0.0
            @test strut_drag_for_test_case!(ac)[iaCDstrut] ≈ expected_strut rtol = 1e-12

            ac.wing.jury.length = 2.5
            ac.wing.strut.S = 0.0
            @test strut_drag_for_test_case!(ac)[iaCDstrut] ≈ expected_jury rtol = 1e-12

            ac.wing.strut.S = 60.0
            ac.wing.strut.chord = 0.0
            ac.wing.strut.root_chord = 0.0
            @test strut_drag_for_test_case!(ac)[iaCDstrut] ≈ expected_jury rtol = 1e-12

            ac.wing.jury.chord = 0.0
            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] == 0.0
            @test isfinite(para[iaCD])
        end

        @testset "Cantilever wing has no strut drag" begin
            ac = load_aircraft_quietly()
            set_strut_drag_test_case!(ac)
            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] == 0.0
            CD = para[iaCD]
            @test isfinite(CD)

            # Strut and jury inputs do not affect the cantilever drag
            ac.wing.strut.S *= 2
            ac.wing.jury.chord *= 2
            @test strut_drag_for_test_case!(ac)[iaCD] ≈ CD rtol = 1e-12
        end
    end

    @testset "Strut section drag from the airfoil database" begin
        air = load_braced_aircraft_quietly().wing.airsection
        strut_section_cd = trefftz_aero.strut_section_cd
        airfoil_lookup = trefftz_aero.airfun

        # Database lookup at the unswept lift and perpendicular Mach number, scaled to Re_chord, with the sweep factor on cdp
        function expected_section_cd(cl, τ, Mach, cosΛ, Re_chord, aRexp, fexcd)
            cdf1, cdp1, _, _ = airfoil_lookup(max(cl / cosΛ^2, air.cl[1]), τ, Mach * cosΛ, air)
            Refac = (Re_chord / air.Re)^aRexp
            return (cdf = cdf1 * Refac * fexcd, cdp = cdp1 * Refac * fexcd * cosΛ)
        end

        # Section drag at the database reference Reynolds number unless `Re_chord` is given, no excrescences unless `fexcd` is given
        section_cd(cl, τ, Mach, cosΛ = 1.0; Re_chord = air.Re, aRexp = -0.15, fexcd = 1.0) =
            strut_section_cd(cl, τ, Mach, cosΛ, Re_chord, aRexp, fexcd, air)

        @testset "inside the database, matches the hand calculation from airfun" begin
            cl, τ, Mach, cosΛ, Re_chord, aRexp, fexcd = 0.5, 0.12, 0.75, 0.94, 5.0e6, -0.15, 1.03
            result = strut_section_cd(cl, τ, Mach, cosΛ, Re_chord, aRexp, fexcd, air)
            expected = expected_section_cd(cl, τ, Mach, cosΛ, Re_chord, aRexp, fexcd)
            @test result.cdf ≈ expected.cdf rtol = 1e-12
            @test result.cdp ≈ expected.cdp rtol = 1e-12
            @test result.cd ≈ expected.cdf + expected.cdp rtol = 1e-12
            @test result.cd ≈ result.cdf + result.cdp rtol = 1e-14
            @test result.cd > 0
            @test isempty(result.exceeded)
        end

        @testset "lift enters through the database lookup at cl / cosΛ²" begin
            cosΛ = 0.94
            low, high = section_cd(0.5, 0.12, 0.75, cosΛ), section_cd(0.7, 0.12, 0.75, cosΛ)
            for (cl, result) in ((0.5, low), (0.7, high))
                _, cdp1, _, _ = airfoil_lookup(cl / cosΛ^2, 0.12, 0.75 * cosΛ, air)
                @test result.cdp ≈ cdp1 * cosΛ rtol = 1e-12
                @test isempty(result.exceeded)
            end
            @test high.cdp != low.cdp
        end

        @testset "thickness: drag grows with t/c, and the penalty applies outside the database" begin
            @test section_cd(0.5, 0.145, 0.80).cd > section_cd(0.5, 0.10, 0.80).cd
            @test isempty(section_cd(0.5, 0.145, 0.80).exceeded)

            thick = section_cd(0.5, 0.20, 0.80)
            expected = expected_section_cd(0.5, 0.20, 0.80, 1.0, air.Re, -0.15, 1.0)
            @test thick.cdp ≈ expected.cdp rtol = 1e-12
            @test thick.cd > section_cd(0.5, 0.145, 0.80).cd
            @test thick.exceeded == [:thickness]
        end

        @testset "Mach drag rise, reduced by sweep" begin
            @test section_cd(0.5, 0.145, 0.80).cd > section_cd(0.5, 0.145, 0.60).cd
            @test section_cd(0.5, 0.145, 0.80, 0.9).cd < section_cd(0.5, 0.145, 0.80, 1.0).cd
        end

        @testset "Reynolds number scaling with the exponent aRexp" begin
            aRexp = -0.15
            Re1, Re2 = 4.0e6, 1.0e7
            result1 = section_cd(0.5, 0.12, 0.75, 0.94; Re_chord = Re1, aRexp)
            result2 = section_cd(0.5, 0.12, 0.75, 0.94; Re_chord = Re2, aRexp)
            @test result1.cd / result2.cd ≈ (Re1 / Re2)^aRexp rtol = 1e-12
            @test result1.cdf / result2.cdf ≈ (Re1 / Re2)^aRexp rtol = 1e-12
            @test result1.cdp / result2.cdp ≈ (Re1 / Re2)^aRexp rtol = 1e-12
        end

        @testset "excrescence factor scales both drag components" begin
            base = section_cd(0.5, 0.12, 0.75, 0.94)
            scaled = section_cd(0.5, 0.12, 0.75, 0.94; fexcd = 1.2)
            @test scaled.cdf ≈ 1.2 * base.cdf rtol = 1e-12
            @test scaled.cdp ≈ 1.2 * base.cdp rtol = 1e-12
        end
    end
end

"""
Wing, tail, strut geometry and pitching-moment inputs of `balance_aircraft!` on an unsized aircraft: sweep 20°, span 35 m,
strut root 1 m and wing attachment 10.5 m from the centreline.
"""
function set_strut_moment_test_case!(ac, lift_fraction_strut_to_total; spanload_exponent = 0.5)
    wing = ac.wing
    wing.layout.S = 125.0
    wing.layout.span = 35.0
    wing.layout.ηs = 0.6
    wing.layout.sweep = 20.0
    wing.layout.root_chord = 5.0
    wing.layout.root_span = 2.0
    ac.htail.layout.S = 25.0
    ac.htail.layout.root_chord = 2.0
    ac.htail.layout.box_x = 30.0
    set_strut_lift_inputs!(ac, lift_fraction_strut_to_total; spanload_exponent)
    ac.wing.strut.local_velocity_ratio = 1.1

    para = view(ac.para, :, ipcruise1, 1)
    para[iaCL] = 1.2
    para[iaCLh] = 0.05
    para[iaCMw1] = -0.1
    para[iaCMh1] = -0.3
    ac.parg[igdCLhdCL] = 0.2
    return ac
end

"Neutral point and centre of pressure from `balance_aircraft!` at the cruise point without changing the geometry."
function balance_neutral_point_and_cp!(ac)
    TASOPT.balance_aircraft!(ac, 1, ipcruise1, 0.5, 1.0, 0.0, TrimVar.None)
    return ac.para[iaxNP, ipcruise1, 1], ac.para[iaxCP, ipcruise1, 1]
end

const sized_aero_output_aircraft_cache = Dict{Any,Any}()

"Aircraft of wing configuration `config` (\"cantilever\", \"strut_braced\" or \"truss_braced\"), sized once and cached; `deepcopy` before changing it."
function sized_aero_output_aircraft(config; lift_fraction_strut_to_total = 0.25)
    return get!(sized_aero_output_aircraft_cache, (config, lift_fraction_strut_to_total)) do
        CoreLogging.with_logger(CoreLogging.NullLogger()) do
            txt = replace(read(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"), String),
                          r"^[ \t]*(wing_config|has_strut)[ \t]*=.*\n"m => "")
            txt = replace(txt, "[Options]" => "[Options]\n    wing_config = \"$config\""; count = 1)
            ac = mktempdir() do dir
                path = joinpath(dir, "aero_output_test.toml")
                write(path, txt)
                read_aircraft_model(path)
            end
            config == "cantilever" || (ac.wing.layout.ηs = 0.7)
            ac.wing.strut.lift_fraction_strut_to_total = lift_fraction_strut_to_total
            size_aircraft!(ac; printiter = false)
            ac
        end
    end
end

@testset "Strut lift share" verbose=true begin
    function wing_with_sloped_strut(lift_fraction_strut_to_total)
        ac = load_braced_aircraft_quietly()
        set_strut_lift_inputs!(ac, lift_fraction_strut_to_total)
        ac.wing.layout.S = 125.0
        ac.wing.strut.local_velocity_ratio = 1.1
        ac.wing.strut.dihedral = 0.3
        return ac
    end

    @testset "strut_lift_share from the strut definitions" begin
        ac = wing_with_sloped_strut(0.5)
        wing, config = ac.wing, ac.options.trefftz_config

        # The share is the input; strut velocity ratio and dihedral do not enter
        @test TASOPT.strut_lift_share(wing, config) == 0.5
        @test TASOPT.strut_lift_share(wing) == 0.5

        # No lift without a lift fraction, without a strut area or wing area, or in a cantilever Trefftz configuration
        wing.strut.lift_fraction_strut_to_total = 0.0
        @test TASOPT.strut_lift_share(wing, config) == 0.0
        wing.strut.lift_fraction_strut_to_total = 0.5
        wing.strut.S = 0.0
        @test TASOPT.strut_lift_share(wing, config) == 0.0
        @test TASOPT.strut_lift_share(wing) == 0.0
        wing.strut.S = 10.0
        wing.layout.S = 0.0
        @test TASOPT.strut_lift_share(wing) == 0.0
        wing.layout.S = 125.0
        @test TASOPT.strut_lift_share(wing, load_aircraft_quietly().options.trefftz_config) == 0.0
    end

    @testset "sized aircraft keeps the strut lift of its first cruise point" begin
        ac = sized_aero_output_aircraft("strut_braced")
        para = ac.para[:, ipcruise1, 1]
        CLhtail = para[iaCLh] * ac.htail.layout.S / ac.wing.layout.S
        k = TASOPT.strut_lift_share(ac.wing, ac.options.trefftz_config)

        @test 0.0 < k < 0.3
        @test ac.wing.strut.CL ≈ k * (para[iaCL] - CLhtail) rtol = 1e-6
        @test ac.wing.strut.lift_fraction ≈ ac.wing.strut.CL / para[iaCL] rtol = 1e-6
        @test para[iaCLstrut] == ac.wing.strut.CL
    end
end

@testset "Strut lift in the pitching moment" verbose=true begin
    @testset "Strut lift arm and blended arm" begin
        ac = set_strut_moment_test_case!(load_braced_aircraft_quietly(), 0.5; spanload_exponent = 0.0)
        wing, config = ac.wing, ac.options.trefftz_config
        co_CMw1 = -0.5

        k = 0.5
        cantilever_config = load_aircraft_quietly().options.trefftz_config
        wing.strut.lift_fraction_strut_to_total = 0.0

        # Lift centroid between strut root (1 m) and wing attachment (10.5 m): the centroid of the symmetric loading is at 1/2 for every p > 0
        for p in (0.5, 1.0, 2.0)
            wing.strut.spanload_exponent = p
            @test TASOPT.strut_lift_arm(wing) ≈ -tand(20.0) * (1.0 + 9.5 * 0.5) rtol = 1e-12
        end
        @test TASOPT.strut_lift_centroid_fraction(0.5) == 0.5

        # The blend keeps the cantilever arm bit-for-bit without strut lift
        @test TASOPT.strut_blended_lift_arm(co_CMw1, wing, config) === co_CMw1
        wing.strut.lift_fraction_strut_to_total = 0.5
        @test TASOPT.strut_blended_lift_arm(co_CMw1, wing, config) ≈
              (1 - k) * co_CMw1 + k * (-tand(20.0) * (1.0 + 9.5 * 0.5)) rtol = 1e-12
        @test TASOPT.strut_blended_lift_arm(co_CMw1, wing, cantilever_config) === co_CMw1
    end

    @testset "balance_aircraft! neutral point and centre of pressure" begin
        ac_no_strut_lift = set_strut_moment_test_case!(load_braced_aircraft_quietly(), 0.0)
        ac_strut_lift = set_strut_moment_test_case!(load_braced_aircraft_quietly(), 0.5)
        ac_cantilever = set_strut_moment_test_case!(load_aircraft_quietly(), 0.5)
        ac_cantilever_no_strut_lift = set_strut_moment_test_case!(load_aircraft_quietly(), 0.0)

        xNP_no_strut_lift, xCP_no_strut_lift = balance_neutral_point_and_cp!(ac_no_strut_lift)
        xNP, xCP = balance_neutral_point_and_cp!(ac_strut_lift)

        # Without strut lift a braced and a cantilever aircraft balance identically; a cantilever ignores the strut lift input
        @test (xNP_no_strut_lift, xCP_no_strut_lift) == balance_neutral_point_and_cp!(ac_cantilever_no_strut_lift)
        @test (xNP_no_strut_lift, xCP_no_strut_lift) == balance_neutral_point_and_cp!(ac_cantilever)

        # The blended arm replaces the wing lift arm co*CMw1 in both moment terms
        wing, htail, para = ac_strut_lift.wing, ac_strut_lift.htail, ac_strut_lift.para[:, ipcruise1, 1]
        co_CMw1 = wing.layout.root_chord * para[iaCMw1]
        arm_shift = TASOPT.strut_blended_lift_arm(co_CMw1, wing, ac_strut_lift.options.trefftz_config) - co_CMw1
        Sh_S = htail.layout.S / wing.layout.S
        @test arm_shift != 0.0
        @test xNP - xNP_no_strut_lift ≈ arm_shift * (ac_strut_lift.parg[igdCLhdCL] * Sh_S - 1.0) rtol = 1e-10
        @test xCP - xCP_no_strut_lift ≈ -arm_shift * (para[iaCL] - para[iaCLh] * Sh_S) / para[iaCL] rtol = 1e-10
    end
end

@testset "Infeasibility flags" verbose=true begin
    @testset "record, clear and report" begin
        ac = load_aircraft_quietly()

        # Recording returns the aircraft, never warns, and keeps each flag once in order
        @test_logs min_level=CoreLogging.Warn begin
            @test TASOPT.record_infeasibility!(ac, :flag_b) === ac
            TASOPT.record_infeasibility!(ac, :flag_a)
            TASOPT.record_infeasibility!(ac, :flag_b)
        end
        @test ac.infeasibilities == [:flag_b, :flag_a]

        # Reporting emits a single warning that names every flag, and none when empty
        logs, _ = Test.collect_test_logs(min_level = CoreLogging.Warn) do
            TASOPT.report_infeasibilities(ac)
        end
        @test length(logs) == 1
        @test occursin("flag_a", logs[1].message) && occursin("flag_b", logs[1].message)
        empty!(ac.infeasibilities)
        @test_logs min_level=CoreLogging.Warn TASOPT.report_infeasibilities(ac)
    end

    @testset "flags do not accumulate and stale flags are cleared; baseline carries the span and fuel-volume flags" begin
        ac = load_aircraft_quietly()
        for _ in 1:2
            logs, _ = Test.collect_test_logs(min_level = CoreLogging.Warn) do
                size_aircraft!(ac, printiter = false)
            end
            # The default aircraft is wider than its span limit and carries 2.2% more fuel than its wing tanks hold
            @test ac.infeasibilities == [:wing_span_exceeded, :wing_fuel_volume_exceeded]
            warning = only(filter(l -> occursin("nfeasib", string(l.message)), logs))
            @test occursin("wing_span_exceeded", warning.message) && occursin("wing_fuel_volume_exceeded", warning.message)
        end

        # Flags left over from an earlier sizing do not survive the next one
        TASOPT.record_infeasibility!(ac, :stale_flag)
        @test ac.infeasibilities == [:wing_span_exceeded, :wing_fuel_volume_exceeded, :stale_flag]
        Test.collect_test_logs(min_level = CoreLogging.Warn) do
            size_aircraft!(ac, printiter = false)
        end
        @test ac.infeasibilities == [:wing_span_exceeded, :wing_fuel_volume_exceeded]
    end

    @testset "wing span limit" begin
        ac = load_aircraft_quietly()
        @test ac.wing.layout.max_span > 0

        # Flagged only when the span strictly exceeds a positive limit
        for (span, max_span, flagged) in ((40.0, 35.0, true), (35.0, 35.0, false), (30.0, 35.0, false), (40.0, 0.0, false))
            empty!(ac.infeasibilities)
            ac.wing.layout.span = span
            ac.wing.layout.max_span = max_span
            @test TASOPT.check_wing_span!(ac) === ac
            @test (:wing_span_exceeded in ac.infeasibilities) == flagged
        end
    end

    @testset "wing fuel volume" begin
        W_fuel_capacity = 2.0e5
        ac = load_aircraft_quietly()
        @test ac.options.has_wing_fuel

        # Flagged only when the fuel strictly exceeds the tank capacity
        for (ratio, flagged) in ((1.01, true), (1.0, false), (0.9, false))
            empty!(ac.infeasibilities)
            @test TASOPT.check_wing_fuel_volume!(ac, ratio * W_fuel_capacity, W_fuel_capacity) === ac
            @test (:wing_fuel_volume_exceeded in ac.infeasibilities) == flagged
        end

        # Fuel outside the wing is not limited by the wing tank capacity
        empty!(ac.infeasibilities)
        ac.options.has_wing_fuel = false
        @test TASOPT.check_wing_fuel_volume!(ac, 1.5 * W_fuel_capacity, W_fuel_capacity) === ac
        @test isempty(ac.infeasibilities)
    end

    @testset "minimum gauge governing the wing cap thickness" begin
        ac = load_aircraft_quietly()
        @test ac.wing.inboard.min_cap_thickness == 0.0 && ac.wing.outboard.min_cap_thickness == 0.0

        ac.wing.inboard.min_cap_thickness = 0.05
        ac.wing.outboard.min_cap_thickness = 0.05
        Test.collect_test_logs(min_level = CoreLogging.Warn) do
            size_aircraft!(ac, printiter = false)
        end
        @test :wing_cap_min_gauge in ac.infeasibilities
    end

    @testset "wing outside the airfoil database" begin
        ac = load_aircraft_quietly()
        τ_min = ac.wing.airsection.τ[1]

        t_over_c = τ_min - 0.002
        ac.wing.inboard.cross_section.thickness_to_chord = t_over_c
        ac.wing.outboard.cross_section.thickness_to_chord = t_over_c
        Test.collect_test_logs(min_level = CoreLogging.Warn) do
            size_aircraft!(ac, printiter = false)
        end
        @test :airfoil_database_exceeded in ac.infeasibilities
    end

    @testset "strut attachment outboard limit" begin
        @test TASOPT.STRUT_ETA_S_MAX == 0.9

        ac = load_braced_aircraft_quietly()
        @test ac.options.has_strut
        for ηs in (0.95, 0.9, 0.6)
            empty!(ac.infeasibilities)
            ac.wing.layout.ηs = ηs
            @test TASOPT.check_strut_attachment!(ac) === ac
            @test (:strut_attachment_outboard_limit in ac.infeasibilities) == (ηs > TASOPT.STRUT_ETA_S_MAX)
        end

        # A cantilever wing has no strut
        ac = load_aircraft_quietly()
        @test !ac.options.has_strut
        ac.wing.layout.ηs = 0.95
        @test TASOPT.check_strut_attachment!(ac) === ac
        @test isempty(ac.infeasibilities)
    end
end

# Aerodynamics output

"Text printed by `aero` for `ac`."
aero_text(ac) = sprint(io -> TASOPT.aero(ac; io = io))

"Number printed after the label matched by `pattern` (a regex with one capture group) in `text`, or `nothing` if there is no such line."
function printed_aero_value(text, pattern)
    m = match(pattern, text)
    return m === nothing ? nothing : parse(Float64, m.captures[1])
end

@testset "Aerodynamics output" verbose=true begin
    # Printed values are rounded to 5 decimals
    printing_tolerance = 5e-6

    for (name, has_strut) in (("cantilever", false), ("strut_braced", true), ("truss_braced", true))
        @testset "$name" begin
            ac = sized_aero_output_aircraft(name)
            text = aero_text(ac)
            para = ac.para[:, ipcruise1, 1]
            Sref = ac.wing.layout.S

            @test ac.options.has_strut == has_strut

            @testset "strut drag and lift lines are printed only with a strut" begin
                CDstrut = printed_aero_value(text, r"^CDstrut\s*=\s*([-0-9.]+)"m)
                CDstrut_area = printed_aero_value(text, r"^CDstrut\s*×\s*Sref\s*=\s*([-0-9.]+) m²"m)
                CLstrut = printed_aero_value(text, r"^CLstrut\s*=\s*([-0-9.]+)"m)
                lift_fraction = printed_aero_value(text, r"^Strut lift fraction\s*=\s*([-0-9.]+)"m)
                if has_strut
                    @test para[iaCDstrut] > 0
                    @test something(CDstrut, NaN) ≈ para[iaCDstrut] atol = printing_tolerance
                    @test something(CDstrut_area, NaN) ≈ Sref * para[iaCDstrut] atol = printing_tolerance
                    @test 0 < ac.wing.strut.CL < para[iaCL]
                    @test something(CLstrut, NaN) ≈ ac.wing.strut.CL atol = printing_tolerance
                    @test something(lift_fraction, NaN) ≈ ac.wing.strut.lift_fraction atol = printing_tolerance
                else
                    @test CDstrut === nothing
                    @test CDstrut_area === nothing
                    @test ac.wing.strut.CL == 0.0 # a cantilever ignores the strut lift input
                    @test CLstrut === nothing
                    @test lift_fraction === nothing
                end
            end

            @testset "printed drag components add up to the printed drag" begin
                printed = Dict(m.captures[1] => parse(Float64, m.captures[2])
                               for m in eachmatch(r"^(CD\w*)\s*=\s*([-0-9.]+)"m, text))
                components = [value for (label, value) in printed if label != "CD"]
                @test printed["CD"] ≈ para[iaCD] atol = printing_tolerance
                @test haskey(printed, "CDstrut") == has_strut
                @test sum(components) ≈ printed["CD"] atol = 1e-4
            end
        end
    end
end

# Geometry output

"Text printed by `geometry` for `ac`."
geometry_text(ac) = sprint(io -> TASOPT.geometry(ac; io = io))

"Word printed after the label matched by `pattern` (a regex with one capture group) in `text`, or `nothing` if there is no such line."
function printed_geometry_word(text, pattern)
    m = match(pattern, text)
    return m === nothing ? nothing : String(m.captures[1])
end

@testset "Geometry output" verbose=true begin
    # Lengths, areas and angles are printed with at least two decimals, so the printed value is within 0.01 of the model value
    printing_tolerance = 0.01

    for (name, has_strut, has_jury) in (("cantilever", false, false), ("strut_braced", true, false), ("truss_braced", true, true))
        @testset "$name" begin
            ac = sized_aero_output_aircraft(name)
            text = geometry_text(ac)
            strut = ac.wing.strut
            jury = ac.wing.jury
            sizing_loads = ac.wing.sizing_loads
            printed(label) = printed_aero_value(text, Regex("^" * label * raw"\s*=\s*([-0-9.]+)", "m"))
            printed_word(label) = printed_geometry_word(text, Regex("^" * label * raw"\s*=\s*(\w+)", "m"))

            @test ac.options.has_strut == has_strut

            @testset "strut and jury layout is printed only when present" begin
                @test occursin("Strut Layout:", text) == has_strut
                strut_labels = ("l_sp", "c_strut", "S_strut", "z_strut", "y_strut", "dihedral", "t/c")
                if has_strut
                    layout = ac.wing.layout
                    bay = 0.5 * layout.span * (layout.ηs - layout.ηo)
                    @test strut.dihedral ≈ atan(strut.z / bay) > 0
                    @test something(printed("l_sp"), NaN) ≈ strut.length atol = printing_tolerance
                    @test something(printed("c_strut"), NaN) ≈ strut.chord atol = printing_tolerance
                    @test something(printed("S_strut"), NaN) ≈ strut.S atol = printing_tolerance
                    @test something(printed("z_strut"), NaN) ≈ strut.z atol = printing_tolerance
                    @test something(printed("y_strut"), NaN) ≈ 0.5 * layout.root_span atol = printing_tolerance
                    @test something(printed("dihedral"), NaN) ≈ rad2deg(strut.dihedral) atol = printing_tolerance
                    @test something(printed("t/c"), NaN) ≈ strut.thickness_to_chord atol = printing_tolerance
                    @test printed_word("sizing case") == string(sizing_loads.governing_strut)
                    @test printed_word("compression case") == string(sizing_loads.governing_strut_compression)
                    @test printed_word("sizing criterion") == string(strut.sizing_criterion)
                else
                    @test strut.dihedral == 0.0
                    @test all(isnothing, printed.(strut_labels))
                    @test isnothing(printed_word("sizing case"))
                end
                if has_jury
                    @test something(printed("l_jury"), NaN) ≈ jury.length atol = printing_tolerance
                    @test something(printed("c_jury"), NaN) ≈ jury.chord atol = printing_tolerance
                    @test something(printed("eta_jury"), NaN) ≈ jury.eta_wing atol = printing_tolerance
                    @test something(printed("t/c_jury"), NaN) ≈ jury.thickness_to_chord atol = printing_tolerance
                else
                    @test all(isnothing, printed.(("l_jury", "c_jury", "eta_jury", "t/c_jury")))
                end
            end
        end
    end
end

# Drag breakdown plots

"""
Drag component bars of the plot `p` as `label => values` pairs in plot order, one value per bar. A stacked bar is stored as
a `:shape` series whose polygon runs from the top of the segment down to its bottom, so the value is the height of the segment.
"""
function drag_bar_segments(p)
    segments = Pair{String,Vector{Float64}}[]
    for series in p.series_list
        label = string(series[:label])
        if series[:seriestype] == :shape && startswith(label, "CD")
            corners = reshape(filter(!isnan, series[:y]), 5, :)
            push!(segments, label => abs.(corners[1, :] .- corners[2, :]))
        end
    end
    return segments
end

@testset "Drag breakdown plots" verbose=true begin
    breakdown_points = [ipclimb1, ipcruise1, ipcruisen, ipdescent1, ipdescentn]
    base_labels = ["CDi", "CDnace", "CDvtail", "CDhtail", "CDwing", "CDfuse"]

    @testset "cantilever has no strut series" begin
        ac = sized_aero_output_aircraft("cantilever")
        @test first.(drag_bar_segments(plot_drag_breakdown(ac))) == base_labels
        @test first.(drag_bar_segments(plot_details(ac))) == base_labels
    end

    for name in ("strut_braced", "truss_braced")
        @testset "$name" begin
            ac = sized_aero_output_aircraft(name)
            para = ac.para[:, :, 1]
            CDstrut = [para[iaCDstrut, i] for i in breakdown_points]
            CD = [para[iaCD, i] for i in breakdown_points]

            @testset "plot_drag_breakdown fractions" begin
                segments = drag_bar_segments(plot_drag_breakdown(ac))
                @test first.(segments) == [base_labels; "CDstrut"]
                @test sum(last.(segments)) ≈ ones(length(breakdown_points)) atol = 1e-6
                @test last(last(segments)) ≈ CDstrut ./ CD atol = 1e-12
            end

            @testset "plot_details" begin
                p = plot_details(ac)
                @test p isa Plots.Plot
                segments = drag_bar_segments(p)
                @test first.(segments) == [base_labels; "CDstrut"]
                @test only(last(last(segments))) ≈ para[iaCDstrut, ipcruise1] / para[iaCD, ipcruise1] atol = 1e-12
            end
        end
    end
end
