const CoreLogging = Base.CoreLogging
# Tests for the Trefftz plane generalisation to N lifting surfaces (wing, tail and, for braced
# wings, a strut). The wing-only-plus-tail (cantilever) results must stay identical to the
# two-surface implementation.

const trefftz_aero = TASOPT.aerodynamics

"Load an aircraft model from a TOML file without printing log messages."
function load_aircraft_quietly(path = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml"))
    return CoreLogging.with_logger(CoreLogging.NullLogger()) do
        read_aircraft_model(path)
    end
end

load_braced_aircraft_quietly() =
    load_aircraft_quietly(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw.toml"))

"""
Set a fixed wing and tail geometry and loading on an unsized aircraft, and return the
`para` vector that `induced_drag!` reads. Geometry values are the ones used in the Trefftz
integration test of `unit_test_trefftz.jl`; the lift split gives a wing CL of 1.19 and a tail CL of 0.01.
"""
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
    para[iaCLh] = 0.05 # tail CL referenced to the tail area; normalised by the wing area it is 0.01
    para[iarcls] = 0.77
    para[iarclt] = 0.15
    return para
end

"""
Run `induced_drag!` on `ac` for the fixed test case; return `(CDi, span efficiency)`.
The cached wake geometry is discarded first unless `reuse_cached_geometry` is true, so that each
test starts from a freshly built wake system.
"""
function trefftz_drag_for_test_case!(ac; reuse_cached_geometry = false)
    para = set_trefftz_test_case!(ac)
    reuse_cached_geometry || (trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0))
    trefftz_aero.induced_drag!(para, ac, ac.options.trefftz_config)
    return para[iaCDi], para[iaspaneff]
end

"""
Set the strut inputs that determine the strut lift on an unsized aircraft: an area, the free-stream velocity at the
strut, `cl_ratio` and `spanload_exponent`. `set_trefftz_test_case!` leaves the strut inputs alone, so this can be
called before `trefftz_drag_for_test_case!`.
"""
function set_strut_lift_inputs!(ac, cl_ratio; spanload_exponent = 0.5)
    ac.wing.strut.S = 10.0
    ac.wing.strut.local_velocity_ratio = 1.0
    ac.wing.strut.cl_ratio = cl_ratio
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
Set the inputs of the strut profile drag, and everything else `aircraft_drag!` needs, on an unsized aircraft; return
the `para` slice of mission 1, point 1 that `aircraft_drag!` writes to. The values are round numbers with no
particular meaning, chosen so that the strut and the jury strut each have a distinct, non-zero drag.
"""
function set_strut_drag_test_case!(ac)
    para = view(ac.para, :, 1, 1)
    lift_inputs = set_trefftz_test_case!(ac)
    for i in (iaCL, iaCLh, iarcls, iarclt)
        para[i] = lift_inputs[i]
    end
    para[iaMach] = 0.78
    para[iaReunit] = 7.0e6

    ac.wing.layout.root_chord = 5.0
    ac.htail.layout.root_chord = 2.0
    ac.vtail.layout.root_chord = 2.0
    ac.vtail.layout.span = 5.0
    ac.vtail.layout.root_span = 1.0
    ac.vtail.outboard.λ = 0.3

    set_strut_lift_inputs!(ac, 0.0)
    ac.wing.strut.chord = 0.8
    ac.wing.strut.cos_lambda = 0.9
    ac.wing.strut.local_velocity_ratio = 1.1
    ac.wing.jury.chord = 0.4
    ac.wing.jury.length = 2.5
    return para
end

"Run `aircraft_drag!` for mission 1, point 1 and return the `para` slice, which holds `iaCDstrut` and `iaCD`."
function strut_drag_for_test_case!(ac)
    para = view(ac.para, :, 1, 1)
    trefftz_aero.aircraft_drag!(ac, 1, 1, true)
    return para
end

"Strut profile drag coefficient of the main strut for the inputs set on `ac` and `para`, from its definition."
function expected_main_strut_drag(ac, para)
    strut, Sw = ac.wing.strut, ac.wing.layout.S
    Re_factor = (para[iaReunit] * strut.chord / para[iaRerefs])^para[iaaRexp]
    return strut.S / Sw * (para[iacdfs] + para[iacdps] * strut.cos_lambda^3) * Re_factor * strut.local_velocity_ratio^3
end

"Strut profile drag coefficient of the two jury struts: vertical, unswept and without lift, from its definition."
function expected_jury_drag(ac, para)
    jury, Sw = ac.wing.jury, ac.wing.layout.S
    Re_factor = (para[iaReunit] * jury.chord / para[iaRerefs])^para[iaaRexp]
    return 2 * jury.chord * jury.length / Sw * (para[iacdfs] + para[iacdps]) * Re_factor *
           ac.wing.strut.local_velocity_ratio^3
end

@testset "Trefftz surfaces" verbose=true begin

    @testset "Configuration without strut panels" begin
        config = trefftz_aero.get_trefftz_config("MEDIUM")
        @test trefftz_aero.n_surfaces(config) == 2
        @test trefftz_aero.n_total(config.strut_panels) == 0

        # MEDIUM: wing 40+12+6 panels, tail 20+0+4 panels, one dummy point after each surface
        @test trefftz_aero.i_first_wing(config) == 1
        @test trefftz_aero.i_last_wing(config) == 59
        @test trefftz_aero.i_first_tail(config) == 60
        @test trefftz_aero.i_last_tail(config) == 84
        @test trefftz_aero.n_points_used(config) == 84

        # braced = false is the same as not asking for a braced layout
        unbraced = trefftz_aero.get_trefftz_config("MEDIUM"; braced = false)
        @test trefftz_aero.n_surfaces(unbraced) == 2
        @test trefftz_aero.n_points_used(unbraced) == 84
    end

    @testset "Configuration with strut panels" begin
        # quality => (strut panels, points without strut, points with strut)
        # The strut has halfway between the wing's and the tail's total panel count, rounded down
        cases = ("COARSE" => (20, 43, 64),
                 "MEDIUM" => (41, 84, 126),
                 "FINE"   => (164, 330, 495))
        for (quality, (n_strut_panels, n_points_cantilever, n_points_braced)) in cases
            @testset "$quality" begin
                config = trefftz_aero.get_trefftz_config(quality; braced = true)
                @test trefftz_aero.n_surfaces(config) == 3
                @test config.strut_panels == trefftz_aero.SurfaceDiscretization(n_strut_panels, 0, 0)

                # Each surface block follows the previous one without a gap or an overlap
                @test trefftz_aero.i_first_tail(config) == trefftz_aero.i_last_wing(config) + 1
                @test trefftz_aero.i_first_strut(config) == trefftz_aero.i_last_tail(config) + 1
                @test trefftz_aero.i_last_strut(config) ==
                      trefftz_aero.i_first_strut(config) + n_strut_panels

                # The strut block ends the point list
                @test trefftz_aero.n_points_used(config) == n_points_braced
                @test trefftz_aero.i_last_strut(config) == n_points_braced
                @test trefftz_aero.i_last_tail(config) ==
                      trefftz_aero.n_points_used(trefftz_aero.get_trefftz_config(quality))
                @test trefftz_aero.i_last_tail(config) == n_points_cantilever

                # The strut has its panels plus one dummy point after it
                @test trefftz_aero.n_strut_points(config.strut_panels) == n_strut_panels + 1
            end
        end
    end

    @testset "Geometry array overflow includes the strut block" begin
        # FINE wing and tail use 330 of the points in the geometry arrays
        n_slots = length(trefftz_aero.TREFFTZ_GEOM.y)
        n_strut_panels_that_fit = n_slots - 330 - 1 # one more point is the dummy point after the strut
        @test n_slots >= 495 # FINE with a strut
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

        # With a T-tail the tail image panels (2 in COARSE) drop out, the strut block does not
        htail.layout.root_span = 0.0
        @test trefftz_aero.calculate_array_size(cantilever_config, htail) == 41
        @test trefftz_aero.calculate_array_size(braced_config, htail) == 62
    end

    @testset "generate_trefftz_points_inclined!" begin
        # Fill the arrays with NaN so that any entry the function must not touch is detected
        function nan_filled_geometry()
            geom = trefftz_aero.TrefftzGeometry{40}()
            for name in (:t, :y, :yp, :z, :zp, :yc, :ycp, :zc, :zcp, :gc)
                fill!(getfield(geom, name), NaN)
            end
            return geom
        end

        @testset "Rising segment starting at an offset index" begin
            geom = nan_filled_geometry()
            i_start, n_panels = 5, 8
            y_start, z_start, y_end, z_end = 0.5, -1.4, 21.3, 1.7
            i_end = trefftz_aero.generate_trefftz_points_inclined!(
                i_start, geom, y_start, z_start, y_end, z_end, n_panels)

            @test i_end == i_start + n_panels
            nodes = i_start:i_end
            controls = i_start:(i_end - 1)

            @test geom.y[i_start] == y_start
            @test geom.z[i_start] == z_start
            @test geom.y[i_end] ≈ y_end rtol = 1e-14
            @test geom.z[i_end] ≈ z_end rtol = 1e-14

            # Evenly spaced and collinear
            @test all(isapprox.(diff(geom.y[nodes]), (y_end - y_start) / n_panels; rtol = 1e-12))
            slope = (z_end - z_start) / (y_end - y_start)
            @test all(isapprox.(diff(geom.z[nodes]) ./ diff(geom.y[nodes]), slope; rtol = 1e-12))
            @test all(isapprox.(geom.z[nodes], z_start .+ slope .* (geom.y[nodes] .- y_start); atol = 1e-12))

            # No wake contraction
            @test geom.yp[nodes] == geom.y[nodes]
            @test geom.zp[nodes] == geom.z[nodes]

            # Control points are the midpoints of consecutive nodes, and the strut carries no circulation
            @test geom.yc[controls] ≈ 0.5 .* (geom.y[i_start:(i_end - 1)] .+ geom.y[(i_start + 1):i_end])
            @test geom.zc[controls] ≈ 0.5 .* (geom.z[i_start:(i_end - 1)] .+ geom.z[(i_start + 1):i_end])
            @test geom.ycp[controls] == geom.yc[controls]
            @test geom.zcp[controls] == geom.zc[controls]
            @test all(==(0.0), geom.gc[controls])

            # Nothing outside the strut nodes (and the control points between them) is written
            @test all(isnan, geom.y[1:(i_start - 1)])
            @test all(isnan, geom.y[(i_end + 1):end])
            @test all(isnan, geom.yc[1:(i_start - 1)])
            @test all(isnan, geom.yc[i_end:end])
            @test all(isnan, geom.gc[i_end:end])
        end

        @testset "Falling segment" begin
            geom = nan_filled_geometry()
            i_end = trefftz_aero.generate_trefftz_points_inclined!(1, geom, 2.0, 3.0, 10.0, -1.0, 4)
            @test i_end == 5
            @test geom.z[1:5] ≈ [3.0, 2.0, 1.0, 0.0, -1.0]
            @test geom.y[1:5] ≈ [2.0, 4.0, 6.0, 8.0, 10.0]
        end

        @testset "Single panel" begin
            geom = nan_filled_geometry()
            i_end = trefftz_aero.generate_trefftz_points_inclined!(7, geom, 1.0, 0.0, 5.0, 2.0, 1)
            @test i_end == 8
            @test (geom.y[7], geom.z[7]) == (1.0, 0.0)
            @test (geom.y[8], geom.z[8]) == (5.0, 2.0)
            @test (geom.yc[7], geom.zc[7]) == (3.0, 1.0)
            @test (geom.ycp[7], geom.zcp[7]) == (3.0, 1.0)
            @test geom.gc[7] == 0.0
            @test isnan(geom.yc[8]) && isnan(geom.gc[8])
        end
    end

    @testset "Cantilever induced drag is unchanged" begin
        ac = load_aircraft_quietly()
        @test trefftz_aero.n_surfaces(ac.options.trefftz_config) == 2

        CDi, span_efficiency = trefftz_drag_for_test_case!(ac)

        # Reference values were computed with the two-surface implementation (before strut panels
        # existed) by running this same test case. The tolerance only allows for platform-dependent
        # rounding in the linear algebra.
        @test CDi ≈ 0.05468605653459833 rtol = 1e-12
        @test span_efficiency ≈ 0.8319732645252401 rtol = 1e-12
        @test length(ac.wake_system.points) == 43
        @test length(ac.wake_system.elements) == 42
    end

    @testset "Braced aircraft gets a strut surface in the wake system" begin
        ac_cantilever = load_aircraft_quietly()
        ac_braced = load_braced_aircraft_quietly()
        braced_config = ac_braced.options.trefftz_config
        @test trefftz_aero.n_surfaces(braced_config) == 3
        @test trefftz_aero.n_surfaces(ac_cantilever.options.trefftz_config) == 2

        # Both aircraft use COARSE, so the strut is the only difference between the two configurations
        @test trefftz_aero.n_points_used(braced_config) ==
              trefftz_aero.n_points_used(ac_cantilever.options.trefftz_config) + 21

        CDi_cantilever, span_efficiency_cantilever = trefftz_drag_for_test_case!(ac_cantilever)
        CDi_braced, span_efficiency_braced = trefftz_drag_for_test_case!(ac_braced)

        @test length(ac_braced.wake_system.points) == 64
        @test length(ac_braced.wake_system.elements) == 63
        @test isfinite(CDi_braced) && CDi_braced > 0.0
        @test isfinite(span_efficiency_braced) && 0.0 < span_efficiency_braced <= 1.0

        # The strut carries no lift, so its wake has zero circulation and cannot induce drag on the
        # wing or tail; only rounding differences remain (rtol = 1e-8 is far above them).
        @test CDi_braced ≈ CDi_cantilever rtol = 1e-8
        @test span_efficiency_braced ≈ span_efficiency_cantilever rtol = 1e-8

        i_strut = trefftz_aero.i_first_strut(braced_config):(trefftz_aero.i_last_strut(braced_config) - 1)
        @test all(==(0.0), trefftz_aero.TREFFTZ_GEOM.gc[i_strut])

    end

    @testset "Strut is an inclined straight line from the fuselage to the wing" begin
        ac = load_braced_aircraft_quietly()
        config = ac.options.trefftz_config
        trefftz_drag_for_test_case!(ac)
        wing = ac.wing
        @test wing.strut.z > 0.0

        i_nodes = trefftz_aero.i_first_strut(config):trefftz_aero.i_last_strut(config)
        strut_y = [ac.wake_system.points[i][1] for i in i_nodes]
        strut_z = [ac.wake_system.points[i][2] for i in i_nodes]
        @test length(i_nodes) == trefftz_aero.n_total(config.strut_panels) + 1

        # Fuselage attachment is strut.z below the wing box; the wing attachment is at the wing height
        @test strut_y[1] == wing.strut.y_root
        @test strut_z[1] ≈ wing.layout.z - wing.strut.z rtol = 1e-14
        @test strut_y[end] ≈ 0.5 * wing.layout.break_span rtol = 1e-14
        @test strut_z[end] ≈ wing.layout.z rtol = 1e-14

        # The strut rises towards the wing along a line of constant slope
        @test all(diff(strut_y) .> 0.0)
        @test all(diff(strut_z) .> 0.0)
        slope = (strut_z[end] - strut_z[1]) / (strut_y[end] - strut_y[1])
        @test all(isapprox.(diff(strut_z) ./ diff(strut_y), slope; rtol = 1e-10))
    end

    @testset "Inclined strut elements have finite influence coefficients" begin
        ac = load_braced_aircraft_quietly()
        config = ac.options.trefftz_config
        trefftz_drag_for_test_case!(ac)
        influence_matrix = ac.wake_system.influence_matrix

        n_points = trefftz_aero.n_points_used(config)
        @test size(influence_matrix) == (n_points - 1, n_points)
        @test all(isfinite, influence_matrix)

        # Strut control points are rows of the strut block; the strut wake points are columns of it
        i_strut_elements = trefftz_aero.i_first_strut(config):(trefftz_aero.i_last_strut(config) - 1)
        i_strut_points = trefftz_aero.i_first_strut(config):trefftz_aero.i_last_strut(config)
        @test all(isfinite, influence_matrix[i_strut_elements, :])
        @test all(isfinite, influence_matrix[:, i_strut_points])
        @test any(!=(0.0), influence_matrix[i_strut_elements, :])
    end

    @testset "Switching between two and three surfaces on unchanged geometry" begin
        ac_cantilever = load_aircraft_quietly()
        ac_braced = load_braced_aircraft_quietly()

        # The wake geometry is cached on the geometry: the surface count must also trigger a rebuild
        CDi_cantilever, _ = trefftz_drag_for_test_case!(ac_cantilever)
        trefftz_drag_for_test_case!(ac_braced; reuse_cached_geometry = true)
        @test length(ac_braced.wake_system.points) == 64
        CDi_cantilever_again, _ = trefftz_drag_for_test_case!(ac_cantilever; reuse_cached_geometry = true)
        @test length(ac_cantilever.wake_system.points) == 43
        @test CDi_cantilever_again ≈ CDi_cantilever rtol = 1e-12
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

            # Each change is made on top of the previous state, without resetting the cached geometry
            changes = (
                "strut.y_root" => () -> wing.strut.y_root += 0.2,
                "strut.z"      => () -> wing.strut.z += 0.1,
                "layout.z"     => () -> wing.layout.z += 0.3,
                "layout.ηs"    => () -> wing.layout.ηs += 0.05)
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
                    @test strut_after[1][1] == wing.strut.y_root
                    @test strut_after[1][2] ≈ wing.layout.z - wing.strut.z rtol = 1e-14
                    @test strut_after[end][1] ≈ 0.5 * wing.layout.break_span rtol = 1e-14
                    @test strut_after[end][2] ≈ wing.layout.z rtol = 1e-14
                    @test all_at_height(wing_after, wing.layout.z)
                    @test strut_after != strut_before
                    name == "layout.z" && @test wing_after != wing_before
                end
            end
        end

        @testset "Cantilever aircraft: wing and tail height" begin
            for (name, surface_of) in ("wing" => ac -> ac.wing, "htail" => ac -> ac.htail)
                @testset "$name.layout.z" begin
                    ac = load_aircraft_quietly()
                    config = ac.options.trefftz_config
                    surface = surface_of(ac)
                    para = set_trefftz_test_case!(ac)
                    i_wing = trefftz_aero.i_first_wing(config):trefftz_aero.i_last_wing(config)
                    i_tail = trefftz_aero.i_first_tail(config):trefftz_aero.i_last_tail(config)

                    trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
                    trefftz_aero.induced_drag!(para, ac, config)
                    wake_system = ac.wake_system
                    CDi_before = para[iaCDi]

                    surface.layout.z += 0.4
                    trefftz_aero.induced_drag!(para, ac, config)

                    @test ac.wake_system !== wake_system
                    @test all_at_height(wake_points(ac, i_wing), ac.wing.layout.z)
                    @test all_at_height(wake_points(ac, i_tail), ac.htail.layout.z)

                    # The height difference between wing and tail changes the induced drag. The result
                    # with the updated wake must be the one of a wake built from scratch.
                    CDi_updated, span_efficiency_updated = para[iaCDi], para[iaspaneff]
                    @test !isapprox(CDi_updated, CDi_before; rtol = 1e-4)
                    trefftz_aero.TREFFTZ_GEOMETRY_HASH[] = UInt64(0)
                    trefftz_aero.induced_drag!(para, ac, config)
                    @test CDi_updated ≈ para[iaCDi] rtol = 1e-12
                    @test span_efficiency_updated ≈ para[iaspaneff] rtol = 1e-12
                end
            end

            # Moving the wing and the tail by the same height leaves the induced drag unchanged
            ac = load_aircraft_quietly()
            CDi, span_efficiency = trefftz_drag_for_test_case!(ac)
            ac.wing.layout.z += 0.4
            ac.htail.layout.z += 0.4
            CDi_shifted, span_efficiency_shifted = trefftz_drag_for_test_case!(ac; reuse_cached_geometry = true)
            @test CDi_shifted ≈ CDi rtol = 1e-12
            @test span_efficiency_shifted ≈ span_efficiency rtol = 1e-12
        end
    end

    @testset "Cantilever induced drag ignores the strut lift input" begin
        ac = load_aircraft_quietly()
        set_strut_lift_inputs!(ac, 1.0)

        CDi, span_efficiency = trefftz_drag_for_test_case!(ac)

        @test CDi ≈ 0.05468605653459833 rtol = 1e-12
        @test span_efficiency ≈ 0.8319732645252401 rtol = 1e-12
        @test ac.wing.strut.CL == 0.0
        @test ac.wing.strut.lift_fraction == 0.0
    end

    @testset "Strut lift inputs" begin
        @testset "Defaults" begin
            for ac in (load_braced_aircraft_quietly(), load_aircraft_quietly())
                @test ac.wing.strut.cl_ratio == 0.0
                @test ac.wing.strut.spanload_exponent == 0.0
            end
        end

        @testset "Read from [Wing.Strut]" begin
            path = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_tbw.toml")
            toml = replace(read(path, String),
                r"(?m)^(\s*cl_ratio\s*=\s*)[^\s#]+" => s"\g<1>0.7",
                r"(?m)^(\s*spanload_exponent\s*=\s*)[^\s#]+" => s"\g<1>0.8")
            mktempdir() do dir
                modified_path = joinpath(dir, "strut_lift.toml")
                write(modified_path, toml)
                strut = load_aircraft_quietly(modified_path).wing.strut
                @test strut.cl_ratio == 0.7
                @test strut.spanload_exponent == 0.8
            end
        end

        @testset "No strut lift when cl_ratio is zero" begin
            ac = load_braced_aircraft_quietly()
            set_strut_lift_inputs!(ac, 0.0)
            trefftz_drag_for_test_case!(ac)
            @test ac.wing.strut.CL == 0.0
            @test ac.wing.strut.lift_fraction == 0.0
            i_strut = trefftz_aero.i_first_strut(ac.options.trefftz_config):trefftz_aero.i_last_strut(ac.options.trefftz_config)
            @test all(==(0.0), trefftz_aero.gw[i_strut])
        end
    end

    @testset "Strut carries lift in addition to the wing" begin
        CL, CLhtail = 1.2, 0.05 * 25.0 / 125.0 # CL and tail CL referenced to the wing area, as in set_trefftz_test_case!
        cl_ratio, spanload_exponent = 1.0, 0.5

        ac_no_lift = load_braced_aircraft_quietly()
        set_strut_lift_inputs!(ac_no_lift, 0.0)
        CDi_no_lift, _ = trefftz_drag_for_test_case!(ac_no_lift)

        ac = load_braced_aircraft_quietly()
        set_strut_lift_inputs!(ac, cl_ratio; spanload_exponent)
        CDi, span_efficiency = trefftz_drag_for_test_case!(ac)
        config = ac.options.trefftz_config
        strut = ac.wing.strut
        geom = trefftz_aero.TREFFTZ_GEOM
        i_first, i_last = trefftz_aero.i_first_strut(config), trefftz_aero.i_last_strut(config)

        @testset "Strut lift coefficient" begin
            CL_strut = cl_ratio * (strut.S / ac.wing.layout.S) * strut.local_velocity_ratio^2 * (CL - CLhtail)
            @test CL_strut > 0.0
            @test strut.CL ≈ CL_strut rtol = 1e-12
            @test strut.lift_fraction ≈ CL_strut / CL rtol = 1e-12

            # The strut lift is added to the wing's own lift: the wing target is reduced and the total stays CL
            strut_lift = trefftz_surface_lift(ac, i_first, i_last)
            wing_lift = trefftz_surface_lift(ac, trefftz_aero.i_first_wing(config), trefftz_aero.i_last_wing(config))
            tail_lift = trefftz_surface_lift(ac, trefftz_aero.i_first_tail(config), trefftz_aero.i_last_tail(config))
            @test strut_lift ≈ CL_strut rtol = 1e-10
            @test tail_lift ≈ CLhtail rtol = 1e-10
            @test wing_lift ≈ CL - CLhtail - CL_strut rtol = 1e-10
            @test wing_lift + tail_lift + strut_lift ≈ CL rtol = 1e-10
        end

        @testset "Strut circulation shape" begin
            n = trefftz_aero.n_total(config.strut_panels)
            gc = geom.gc[i_first:(i_last - 1)]
            # Each strut is half a wing: ξ = 0 at the fuselage (like a wing root) and ξ = 1 at the wing attachment
            ξ = [(k - 0.5) / n for k in 1:n]
            shape = (1 .- ξ .^ 2) .^ spanload_exponent

            # The loading is largest at the fuselage end; no tip vortex where the strut is attached to the wing
            @test gc[1] > 0.0
            @test gc[end] == 0.0
            @test all(>(0.0), gc[1:(end - 1)])

            # Up to the last panel the circulation follows (1-ξ²)^p and decreases from the fuselage to the wing
            @test gc[1:(end - 1)] ./ gc[1] ≈ shape[1:(end - 1)] ./ shape[1] rtol = 1e-10
            @test all(diff(gc[1:(end - 1)]) .< 0.0)

            # The wake circulation vanishes at both strut ends (no vortex at the fuselage centreline, where the left
            # and right strut cancel, and none at the wing attachment); the interior values telescope, so the net
            # strength equals the loading of the first panel
            gw = trefftz_aero.gw[i_first:i_last]
            @test gw[1] == 0.0
            @test gw[end] == 0.0
            @test sum(gw) ≈ gc[1] rtol = 1e-10
        end

        @testset "Induced drag" begin
            @test isfinite(CDi) && CDi > 0.0
            @test isfinite(span_efficiency) && span_efficiency > 0.0
            # The strut takes about 8 % of the lift from the wing (CL_strut / CL), which changes the induced drag by
            # several percent; 1e-3 is far above rounding and well below that change.
            @test abs(CDi - CDi_no_lift) / CDi_no_lift > 1e-3
        end
    end

    @testset "Strut-braced induced drag differs from the same-span cantilever" begin
        # Same wing span, area and tail: the strut surface is the only difference. With strut lift, its wake interferes
        # with the wing and tail wakes, so the Trefftz-plane drag must change relative to the cantilever.
        ac_cantilever = load_aircraft_quietly()
        CDi_cantilever, span_efficiency_cantilever = trefftz_drag_for_test_case!(ac_cantilever)

        ac_braced = load_braced_aircraft_quietly()
        set_strut_lift_inputs!(ac_braced, 1.0)
        CDi_braced, span_efficiency_braced = trefftz_drag_for_test_case!(ac_braced)

        @test ac_braced.wing.layout.span == ac_cantilever.wing.layout.span
        @test ac_braced.wing.layout.S == ac_cantilever.wing.layout.S
        @test ac_braced.wing.strut.lift_fraction > 0.0
        @test isfinite(CDi_braced) && CDi_braced > 0.0

        # The strut takes about 8 % of the lift, which raises CDi by about 14 % (span efficiency drops from 0.83 to 0.73)
        # because the strut's own, short-span loading is inefficient; 5 % is far above rounding and below that change.
        @test abs(CDi_braced - CDi_cantilever) / CDi_cantilever > 0.05
        @test abs(span_efficiency_braced - span_efficiency_cantilever) / span_efficiency_cantilever > 0.05
    end

    @testset "Strut circulation is rebuilt when its inputs change on cached geometry" begin
        function strut_circulation(ac)
            config = ac.options.trefftz_config
            return copy(trefftz_aero.TREFFTZ_GEOM.gc[trefftz_aero.i_first_strut(config):(trefftz_aero.i_last_strut(config) - 1)])
        end

        @testset "spanload_exponent" begin
            ac = load_braced_aircraft_quietly()
            set_strut_lift_inputs!(ac, 1.0; spanload_exponent = 0.5)
            trefftz_drag_for_test_case!(ac)
            gc_elliptic = strut_circulation(ac)

            ac.wing.strut.spanload_exponent = 1.0
            CDi_updated, _ = trefftz_drag_for_test_case!(ac; reuse_cached_geometry = true)
            gc_updated = strut_circulation(ac)
            @test !isapprox(gc_updated, gc_elliptic; rtol = 1e-3)

            # The result is the one of a wake built from scratch with the new exponent
            CDi_fresh, _ = trefftz_drag_for_test_case!(ac)
            @test CDi_updated ≈ CDi_fresh rtol = 1e-12
            @test gc_updated ≈ strut_circulation(ac) rtol = 1e-12
        end

        @testset "cl_ratio from zero" begin
            ac = load_braced_aircraft_quietly()
            set_strut_lift_inputs!(ac, 0.0)
            trefftz_drag_for_test_case!(ac)

            ac.wing.strut.cl_ratio = 1.0
            CDi_updated, _ = trefftz_drag_for_test_case!(ac; reuse_cached_geometry = true)
            @test ac.wing.strut.CL > 0.0
            CDi_fresh, _ = trefftz_drag_for_test_case!(ac)
            @test CDi_updated ≈ CDi_fresh rtol = 1e-12
        end
    end

    @testset "Strut profile drag" begin
        function strut_drag_case(config)
            ac = load_braced_aircraft_quietly()
            ac.options.opt_wing_config = config
            para = set_strut_drag_test_case!(ac)
            return ac, para
        end

        @testset "Strut-braced wing: strut only, with Reynolds number scaling" begin
            ac, para = strut_drag_case(WingConfig.StrutBraced)
            expected = expected_main_strut_drag(ac, para)
            para = strut_drag_for_test_case!(ac)
            @test expected > 0
            @test para[iaCDstrut] ≈ expected rtol = 1e-12
        end

        @testset "Reynolds number scaling follows Reunit, Rerefs and the strut chord" begin
            ac, para = strut_drag_case(WingConfig.StrutBraced)
            CDstrut = strut_drag_for_test_case!(ac)[iaCDstrut]
            aRexp = para[iaaRexp]

            para[iaRerefs] *= 2
            @test strut_drag_for_test_case!(ac)[iaCDstrut] ≈ CDstrut * 2.0^(-aRexp) rtol = 1e-12
            para[iaRerefs] /= 2

            para[iaReunit] *= 2
            @test strut_drag_for_test_case!(ac)[iaCDstrut] ≈ CDstrut * 2.0^aRexp rtol = 1e-12
            para[iaReunit] /= 2

            ac.wing.strut.chord *= 2
            @test strut_drag_for_test_case!(ac)[iaCDstrut] ≈ CDstrut * 2.0^aRexp rtol = 1e-12
        end

        @testset "Truss-braced wing adds the jury strut drag to the strut drag" begin
            ac, para = strut_drag_case(WingConfig.TrussBraced)
            expected_strut = expected_main_strut_drag(ac, para)
            expected_jury = expected_jury_drag(ac, para)
            para = strut_drag_for_test_case!(ac)
            @test expected_jury > 0
            @test para[iaCDstrut] ≈ expected_strut + expected_jury rtol = 1e-12

            # The jury drag is part of the total drag
            ac_strut, _ = strut_drag_case(WingConfig.StrutBraced)
            para_strut = strut_drag_for_test_case!(ac_strut)
            @test para[iaCD] - para_strut[iaCD] ≈ expected_jury rtol = 1e-9
        end

        @testset "Jury strut drag uses the jury chord, length and the strut Reynolds number scaling" begin
            ac, para = strut_drag_case(WingConfig.TrussBraced)
            CDjury = strut_drag_for_test_case!(ac)[iaCDstrut] - expected_main_strut_drag(ac, para)
            aRexp = para[iaaRexp]

            ac.wing.jury.length *= 2
            CDjury_long = strut_drag_for_test_case!(ac)[iaCDstrut] - expected_main_strut_drag(ac, para)
            @test CDjury_long ≈ 2 * CDjury rtol = 1e-9

            ac.wing.jury.length /= 2
            para[iaRerefs] *= 2
            CDjury_rerefs = strut_drag_for_test_case!(ac)[iaCDstrut] - expected_main_strut_drag(ac, para)
            @test CDjury_rerefs ≈ CDjury * 2.0^(-aRexp) rtol = 1e-9
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
            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] ≈ expected_jury rtol = 1e-12
            @test isfinite(para[iaCD])

            ac.wing.strut.S = 10.0
            ac.wing.strut.chord = 0.0
            @test strut_drag_for_test_case!(ac)[iaCDstrut] ≈ expected_jury rtol = 1e-12

            ac.wing.jury.chord = 0.0
            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] == 0.0
            @test isfinite(para[iaCD])
        end

        @testset "Cantilever wing has no strut drag" begin
            ac = load_aircraft_quietly()
            para = set_strut_drag_test_case!(ac)
            para = strut_drag_for_test_case!(ac)
            @test para[iaCDstrut] == 0.0
            CD = para[iaCD]
            @test isfinite(CD)

            # Strut and jury inputs and the strut reference Reynolds number do not affect the cantilever drag
            ac.wing.strut.S *= 2
            ac.wing.strut.chord *= 2
            ac.wing.jury.chord *= 2
            para[iaRerefs] *= 2
            @test strut_drag_for_test_case!(ac)[iaCD] ≈ CD rtol = 1e-12
        end
    end
end

"""
Set a wing, tail and strut geometry, and the pitching-moment inputs of `balance_aircraft!`, on an unsized aircraft.
The strut lift inputs come from `set_strut_lift_inputs!`. The wing has a sweep of 20°, a span of 35 m, the break at
60 % of the span and a strut root 1 m from the centreline, so the strut attaches to the wing 10.5 m from the centreline.
"""
function set_strut_moment_test_case!(ac, cl_ratio; spanload_exponent = 0.5)
    wing = ac.wing
    wing.layout.S = 125.0
    wing.layout.span = 35.0
    wing.layout.ηs = 0.6
    wing.layout.sweep = 20.0
    wing.layout.root_chord = 5.0
    wing.strut.y_root = 1.0
    ac.htail.layout.S = 25.0
    ac.htail.layout.root_chord = 2.0
    ac.htail.layout.box_x = 30.0
    set_strut_lift_inputs!(ac, cl_ratio; spanload_exponent)
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

@testset "Strut lift in the pitching moment" verbose=true begin
    @testset "Strut lift share, lift arm and blended arm" begin
        ac = set_strut_moment_test_case!(load_braced_aircraft_quietly(), 0.5; spanload_exponent = 0.0)
        wing, config = ac.wing, ac.options.trefftz_config
        co_CMw1 = -0.5

        # Strut share of the lift: cl_ratio * strut area / wing area * velocity ratio squared
        k = 0.5 * (10.0 / 125.0) * 1.1^2
        @test TASOPT.strut_lift_share(wing, config) ≈ k rtol = 1e-14
        wing.strut.cl_ratio = 0.0
        @test TASOPT.strut_lift_share(wing, config) == 0.0
        cantilever_config = load_aircraft_quietly().options.trefftz_config
        wing.strut.cl_ratio = 0.5
        @test TASOPT.strut_lift_share(wing, cantilever_config) == 0.0
        wing.strut.cl_ratio = 0.0

        # Strut lift centroid between the strut root (1 m) and the wing attachment (10.5 m) for
        # loadings (1-ξ²)^p, with ξ = 0 at the fuselage end: 1/2 (p = 0), 3/8 (p = 1), 4/(3π) (p = 1/2)
        for (p, ξ_cp) in ((0.0, 0.5), (1.0, 3 / 8), (0.5, 4 / (3π)))
            wing.strut.spanload_exponent = p
            @test TASOPT.strut_lift_arm(wing) ≈ -tand(20.0) * (1.0 + 9.5 * ξ_cp) rtol = 1e-12
        end

        # The blend keeps the cantilever arm bit-for-bit without strut lift
        @test TASOPT.strut_blended_lift_arm(co_CMw1, wing, config) === co_CMw1
        wing.strut.cl_ratio = 0.5
        @test TASOPT.strut_blended_lift_arm(co_CMw1, wing, config) ≈
              (1 - k) * co_CMw1 + k * (-tand(20.0) * (1.0 + 9.5 * 4 / (3π))) rtol = 1e-12
        @test TASOPT.strut_blended_lift_arm(co_CMw1, wing, cantilever_config) === co_CMw1
    end

    @testset "balance_aircraft! neutral point and centre of pressure" begin
        ac_no_strut_lift = set_strut_moment_test_case!(load_braced_aircraft_quietly(), 0.0)
        ac_strut_lift = set_strut_moment_test_case!(load_braced_aircraft_quietly(), 0.5)
        ac_cantilever = set_strut_moment_test_case!(load_aircraft_quietly(), 0.5)
        ac_cantilever_no_strut_lift = set_strut_moment_test_case!(load_aircraft_quietly(), 0.0)

        xNP_no_strut_lift, xCP_no_strut_lift = balance_neutral_point_and_cp!(ac_no_strut_lift)
        xNP, xCP = balance_neutral_point_and_cp!(ac_strut_lift)

        # Without strut lift, a braced and a cantilever aircraft balance identically; a cantilever
        # aircraft ignores the strut lift input
        @test (xNP_no_strut_lift, xCP_no_strut_lift) == balance_neutral_point_and_cp!(ac_cantilever_no_strut_lift)
        @test (xNP_no_strut_lift, xCP_no_strut_lift) == balance_neutral_point_and_cp!(ac_cantilever)

        # The strut lift replaces the wing lift arm co*CMw1 by the blended arm in both moment terms
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
        @test ac.infeasibilities == Symbol[]

        # Recording returns the aircraft, never warns, and keeps each flag once in first-occurrence order
        @test_logs min_level=CoreLogging.Warn begin
            @test TASOPT.record_infeasibility!(ac, :flag_b) === ac
            TASOPT.record_infeasibility!(ac, :flag_a)
            TASOPT.record_infeasibility!(ac, :flag_b)
        end
        @test ac.infeasibilities == [:flag_b, :flag_a]

        # Reporting emits a single warning that names every flag
        logs, _ = Test.collect_test_logs(min_level = CoreLogging.Warn) do
            TASOPT.report_infeasibilities(ac)
        end
        @test length(logs) == 1
        @test logs[1].level == CoreLogging.Warn
        @test occursin("flag_a", logs[1].message) && occursin("flag_b", logs[1].message)
        @test ac.infeasibilities == [:flag_b, :flag_a]

        TASOPT.clear_infeasibilities!(ac)
        @test isempty(ac.infeasibilities)
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
            @test count(l -> occursin("nfeasib", string(l.message)), logs) == 1
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

        # The flag is recorded only when the span strictly exceeds a positive limit
        for (span, max_span, flagged) in ((40.0, 35.0, true), (35.0, 35.0, false), (30.0, 35.0, false), (40.0, 0.0, false))
            TASOPT.clear_infeasibilities!(ac)
            ac.wing.layout.span = span
            ac.wing.layout.max_span = max_span
            @test TASOPT.check_wing_span!(ac) === ac
            @test (:wing_span_exceeded in ac.infeasibilities) == flagged
        end
    end

    @testset "sized regional aircraft is within its span limit" begin
        ac = load_aircraft_quietly(joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_regional.toml"))
        Test.collect_test_logs(min_level = CoreLogging.Warn) do
            size_aircraft!(ac, printiter = false)
        end
        @test 0 < ac.wing.layout.span < ac.wing.layout.max_span
        @test !(:wing_span_exceeded in ac.infeasibilities)
    end

    @testset "wing fuel volume" begin
        W_fuel_capacity = 2.0e5
        ac = load_aircraft_quietly()
        @test ac.options.has_wing_fuel

        # The flag is recorded only when the fuel strictly exceeds the tank capacity
        for (ratio, flagged) in ((1.01, true), (1.0, false), (0.9, false))
            TASOPT.clear_infeasibilities!(ac)
            @test TASOPT.check_wing_fuel_volume!(ac, ratio * W_fuel_capacity, W_fuel_capacity) === ac
            @test (:wing_fuel_volume_exceeded in ac.infeasibilities) == flagged
        end

        # Fuel stored outside the wing is not limited by the wing tank capacity
        TASOPT.clear_infeasibilities!(ac)
        ac.options.has_wing_fuel = false
        @test TASOPT.check_wing_fuel_volume!(ac, 1.5 * W_fuel_capacity, W_fuel_capacity) === ac
        @test isempty(ac.infeasibilities)
    end

    @testset "sized cryogenic aircraft has no wing fuel volume flag" begin
        ac = load_aircraft_quietly(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
        @test !ac.options.has_wing_fuel
        Test.collect_test_logs(min_level = CoreLogging.Warn) do
            size_aircraft!(ac, iter = 50, printiter = false)
        end
        @test !(:wing_fuel_volume_exceeded in ac.infeasibilities)

        # The tank volume accounts for the ullage: the fuel weight is below the full-liquid capacity
        fuel_fill_ratio = ac.parg[igWfuel] / ac.parg[igWfmax]
        @test fuel_fill_ratio ≈ ac.parg[igrhofuel] / ac.fuse_tank.rhofuel rtol = 1e-6
        @test fuel_fill_ratio < 1.0
        @test isapprox(fuel_fill_ratio, 1 - ac.fuse_tank.ullage_frac, atol = 0.01)
    end

    @testset "minimum gauge governing the wing cap thickness" begin
        ac = load_aircraft_quietly()
        @test ac.wing.inboard.min_cap_thickness == 0.0 && ac.wing.outboard.min_cap_thickness == 0.0

        # A minimum gauge far above the strength-required cap thickness governs the wing cap thickness
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

        # A wing thinner than the database lower bound is penalised by the extrapolation penalty
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
            TASOPT.clear_infeasibilities!(ac)
            ac.wing.layout.ηs = ηs
            @test TASOPT.check_strut_attachment!(ac) === ac
            @test (:strut_attachment_outboard_limit in ac.infeasibilities) == (ηs > TASOPT.STRUT_ETA_S_MAX)
        end

        # A cantilever wing has no strut to attach
        ac = load_aircraft_quietly()
        @test !ac.options.has_strut
        ac.wing.layout.ηs = 0.95
        @test TASOPT.check_strut_attachment!(ac) === ac
        @test isempty(ac.infeasibilities)
    end

    @testset "wing beam-column ratio" begin
        @test TASOPT.WING_BEAM_COLUMN_RATIO_MAX == 0.5

        # Axial load over the cantilever Euler load π² EI / (4 L²), for either sign of the axial load
        EI, L, P_axial = 1e7, 5.0, -1e5
        ratio = abs(P_axial) * 4 * L^2 / (π^2 * EI)
        @test TASOPT.wing_beam_column_ratio(P_axial, EI, L) ≈ ratio rtol = 1e-14
        @test TASOPT.wing_beam_column_ratio(-P_axial, EI, L) ≈ ratio rtol = 1e-14
        @test TASOPT.wing_beam_column_ratio(0.0, EI, L) == 0.0

        # Doubling the length quarters the Euler load
        @test TASOPT.wing_beam_column_ratio(P_axial, EI, 2L) ≈ 4 * ratio rtol = 1e-14
    end
end

# ============================================================================
# Aerodynamics output: strut drag and lift lines, and the drag components closing with the total
# ============================================================================

"Text printed by `aero` for `ac`."
aero_text(ac) = sprint(io -> TASOPT.aero(ac; io = io))

"Number printed after the label matched by `pattern` (a regex with one capture group) in `text`, or `nothing` if there is no such line."
function printed_aero_value(text, pattern)
    m = match(pattern, text)
    return m === nothing ? nothing : parse(Float64, m.captures[1])
end

const sized_aero_output_aircraft_cache = Dict{String,Any}()

"""
Aircraft of wing configuration `config` ("cantilever", "strut_braced" or "truss_braced"), sized once and shared by all
later requests. A braced wing has its strut attached at ηs = 0.7, which sizes with a strut. Tests that change the aircraft
must `deepcopy` it first.
"""
function sized_aero_output_aircraft(config)
    return get!(sized_aero_output_aircraft_cache, config) do
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
            size_aircraft!(ac; printiter = false)
            ac
        end
    end
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
                    # The default strut carries no lift
                    @test something(CLstrut, NaN) == 0.0
                    @test something(lift_fraction, NaN) == 0.0
                else
                    @test CDstrut === nothing
                    @test CDstrut_area === nothing
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
                # Eight to ten components, each rounded to 5e-6
                @test 8 <= length(components) <= 10
                @test sum(components) ≈ printed["CD"] atol = 1e-4
            end
        end
    end

    @testset "strut lift is printed for a strut with a lift coefficient ratio" begin
        # `wing.strut.CL` and `wing.strut.lift_fraction` hold the values of the last drag evaluation of a sizing, which is
        # at zero lift, so the expected values come from a drag evaluation at the first cruise point
        for name in ("strut_braced", "truss_braced")
            @testset "$name" begin
                ac = deepcopy(sized_aero_output_aircraft(name))
                ac.wing.strut.cl_ratio = 0.5
                trefftz_aero.aircraft_drag!(ac, 1, ipcruise1, true)
                text = aero_text(ac)

                @test 0 < ac.wing.strut.CL < ac.para[iaCL, ipcruise1, 1]
                @test something(printed_aero_value(text, r"^CLstrut\s*=\s*([-0-9.]+)"m), NaN) ≈ ac.wing.strut.CL atol = 5e-6
                @test something(printed_aero_value(text, r"^Strut lift fraction\s*=\s*([-0-9.]+)"m), NaN) ≈
                      ac.wing.strut.lift_fraction atol = 5e-6
            end
        end
    end
end

# ============================================================================
# Geometry output
# ============================================================================

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

            @testset "strut layout is printed only with a strut" begin
                @test occursin("Strut Layout:", text) == has_strut
                if has_strut
                    @test strut.length > 0
                    @test something(printed("l_sp"), NaN) ≈ strut.length atol = printing_tolerance
                    @test something(printed("c_strut"), NaN) ≈ strut.chord atol = printing_tolerance
                    @test something(printed("S_strut"), NaN) ≈ strut.S atol = printing_tolerance
                    @test something(printed("z_strut"), NaN) ≈ strut.z atol = printing_tolerance
                    @test something(printed("y_strut"), NaN) ≈ strut.y_root atol = printing_tolerance
                    @test something(printed("dihedral"), NaN) ≈ rad2deg(strut.dihedral) atol = printing_tolerance
                    @test something(printed("t/c"), NaN) ≈ strut.thickness_to_chord atol = printing_tolerance
                    @test printed_word("sizing case") == string(sizing_loads.governing_strut)
                    @test printed_word("compression case") == string(sizing_loads.governing_strut_compression)
                    @test printed_word("sizing criterion") == string(strut.sizing_criterion)
                else
                    for label in ("l_sp", "c_strut", "S_strut", "z_strut", "y_strut", "dihedral", "t/c")
                        @test printed(label) === nothing
                    end
                    for label in ("sizing case", "compression case", "sizing criterion")
                        @test printed_word(label) === nothing
                    end
                end
            end

            @testset "jury lines are printed only for a truss-braced wing" begin
                if has_jury
                    @test jury.length > 0
                    @test something(printed("l_jury"), NaN) ≈ jury.length atol = printing_tolerance
                    @test something(printed("c_jury"), NaN) ≈ jury.chord atol = printing_tolerance
                    @test something(printed("eta_jury"), NaN) ≈ jury.eta_wing atol = printing_tolerance
                else
                    @test printed("l_jury") === nothing
                    @test printed("c_jury") === nothing
                    @test printed("eta_jury") === nothing
                end
            end

            @testset "strut dihedral is the angle of the in-plane strut with the horizontal" begin
                if has_strut
                    layout = ac.wing.layout
                    bay = 0.5 * layout.span * (layout.ηs - layout.ηo)
                    @test strut.dihedral > 0
                    @test strut.dihedral ≈ atan(strut.z / bay)
                else
                    @test strut.dihedral == 0.0
                end
            end
        end
    end
end

# ============================================================================
# Drag breakdown plots
# ============================================================================

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
        @test first.(drag_bar_segments(plot_drag_breakdown(ac, show_fractions = false))) == base_labels
        @test plot_details(ac) isa Plots.Plot
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

            @testset "plot_drag_breakdown values" begin
                segments = drag_bar_segments(plot_drag_breakdown(ac, show_fractions = false))
                @test first.(segments) == [base_labels; "CDstrut"]
                @test last(last(segments)) ≈ CDstrut atol = 1e-12
                @test sum(last.(segments)) ≈ CD atol = 1e-6
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
