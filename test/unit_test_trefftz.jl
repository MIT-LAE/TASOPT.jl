"""
Unit tests for Trefftz plane induced drag calculations

This test file covers the refactored functions for Trefftz plane analysis:
- WakeGeometry.jl: Geometric structs and influence matrix
- induced_drag.jl: Panel point generation and circulation calculations
"""

using Test
using TASOPT
using LinearAlgebra
using StaticArrays
const aero = TASOPT.aerodynamics

@testset "Trefftz calculations" verbose=true begin

    @testset "Wake geometry primitives" begin
        # Test Point2D creation and basic operations
        @testset "Point2D" begin
            p1 = aero.Point2D(1.0, 2.0)
            p2 = aero.Point2D(4.0, 6.0)

            @test p1[1] == 1.0
            @test p1[2] == 2.0
            @test norm(p2 - p1) ≈ 5.0  # 3-4-5 triangle
        end

        # Test WakeElement construction and properties
        @testset "WakeElement" begin
            p1 = aero.Point2D(0.0, 0.0)
            p2 = aero.Point2D(3.0, 4.0)
            elem = aero.WakeElement(p1, p2)

            @test elem.length ≈ 5.0
            @test elem.Δy ≈ 3.0
            @test elem.Δz ≈ 4.0

            # Control point should be midpoint by default
            @test elem.control_point[1] ≈ 1.5
            @test elem.control_point[2] ≈ 2.0

            # Unit normal should be perpendicular and unit length
            @test norm(elem.unit_normal) ≈ 1.0
            @test dot(elem.unit_normal, p2 - p1) ≈ 0.0  atol=1e-12# Perpendicular
            # HA! There's a gotcha here ^ don't try to just compare to zero, need to use atol. 
            # See `isapprox` documentation for details.

            custom_cp = aero.Point2D(1.0, 1.0)
            elem_custom = aero.WakeElement(p1, p2; control_point=custom_cp)
            @test elem_custom.control_point == custom_cp

            @test_throws ArgumentError aero.WakeElement(p1, p1)
        end

        # Test mirroring function about z-axis
        @testset "mirror_point" begin
            p = aero.Point2D(3.0, 5.0)
            mirrored = aero.mirror_point(p)
            @test mirrored[1] ≈ -3.0
            @test mirrored[2] ≈ 5.0  # z unchanged
        end
    end

    @testset "Bunching transformations" begin
        @testset "bunch_transform" begin
            # Test endpoints
            @test aero.bunch_transform(0.0, 0.5) ≈ 0.0
            @test aero.bunch_transform(1.0, 0.5) ≈ 1.0

            # Test no bunching with bunch=0
            @test aero.bunch_transform(0.3, 0.0) ≈ 0.3
        end

        @testset "inv_bunch_transform" begin
            # Test round-trip consistency
            bunch = 0.5
            for t in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
                t_bunched = aero.bunch_transform(t, bunch)
                t_recovered = aero.inv_bunch_transform(t_bunched, bunch)
                @test t_recovered ≈ t rtol=1e-12
            end

            # Test with different bunch factors
            for bunch in [0.25, 0.5, 0.75]
                t = 0.3
                t_bunched = aero.bunch_transform(t, bunch)
                t_recovered = aero.inv_bunch_transform(t_bunched, bunch)
                @test t_recovered ≈ t rtol=1e-12
            end
        end
    end

    @testset "Wake contraction" begin
        # Test fuselage wake contraction
        @testset "FUSEWAKE contraction" begin
            yo = 2.0   # Half-root span
            yop = 0.4  # Contracted span (20% of yo)
            yexp = (yo / yop)^2

            # At centerline
            y_wake = aero.get_wake_contraction(
                aero.FUSEWAKE, 0.0, yo, yop)
            @test y_wake ≈ 0.0

            # At fuselage edge (should map to contracted edge)
            y_wake = aero.get_wake_contraction(
                aero.FUSEWAKE, yo, yo, yop)
            @test y_wake ≈ yop

            # Intermediate point (should use power law)
            y_wake = aero.get_wake_contraction(
                aero.FUSEWAKE, 0.5 * yo, yo, yop)
            @test y_wake ≈ yop * (0.5)^yexp 

            # Test error for y outside fuselage region
            @test_throws ErrorException aero.get_wake_contraction(
                aero.FUSEWAKE, yo * 1.1, yo, yop)
        end

        # Test wing wake contraction (mass conservation)
        @testset "WINGWAKE contraction" begin
            yo = 2.0
            yop = 0.4

            # Just outside fuselage edge
            y = yo * 1.01
            y_wake = aero.get_wake_contraction(
                aero.WINGWAKE, y, yo, yop)
            @test y_wake ≈ sqrt(y^2 - yo^2 + yop^2)

            y = 10.0
            y_wake = aero.get_wake_contraction(
                aero.WINGWAKE, y, yo, yop)
            @test y_wake ≈ sqrt(y^2 - yo^2 + yop^2)

            # Test error for y inside fuselage region
            @test_throws ErrorException aero.get_wake_contraction(
                aero.WINGWAKE, yo * 0.9, yo, yop)
        end
    end

    @testset "Circulation calculations" begin
        @testset "calculate_wake_circulation!" begin
            # Simple test with 2 surfaces (wing + tail)
            nsurf = 2

            # Wing: 5 panels, Tail: 3 panels
            n_wing_pts = 6  # 5 panels + 1
            n_tail_pts = 4  # 3 panels + 1
            n_total = n_wing_pts + n_tail_pts

            gc = [1.0, 0.9, 0.7, 0.4, 0.1, 0.0,  # Wing bound circulation (decreasing)
                  0.5, 0.3, 0.1, 0.0]              # Tail bound circulation
            gw = zeros(n_total)

            i_first = [1, n_wing_pts + 1]
            i_last = [n_wing_pts, n_total]

            aero.calculate_wake_circulation!(gw, gc, i_first, i_last, nsurf)
            # # Check wing wake circulation
            @test gw[1] ≈ 0.0  # Centerline
            @test gw[2] ≈ gc[1] - gc[2]  # 1.0 - 0.9 = 0.1
            @test gw[3] ≈ gc[2] - gc[3]  # 0.9 - 0.7 = 0.2
            @test gw[n_wing_pts] ≈ gc[n_wing_pts-1]  # Tip: 0.1

            # # Check tail wake circulation
            @test gw[n_wing_pts + 1] ≈ 0.0  # Centerline
            @test gw[n_total] ≈ gc[n_total-1]  # Tip
        end

        @testset "scale_circulation!" begin
            # Create simple test case
            gc = [1.0, 1.0, 1.0]  # Constant circulation
            gw = [0.0, 0.0, 0.0, 1.0]  # Corresponding wake
            yp = [0.0, 1.0, 2.0, 3.0]  # Even spacing

            nsurf = 1
            i_first = [1]
            i_last = [4]
            CLsurfsp = [0.5]  # Target CL
            bref = 6.0
            Sref = 10.0

            aero.scale_circulation!(gc, gw, yp, i_first, i_last,
                CLsurfsp, bref, Sref, nsurf)

            # Calculate actual CL from scaled circulation
            cl_test = 0.0
            for i = 1:3
                dy = yp[i+1] - yp[i]
                cl_test += gc[i] * dy
            end
            cl_test = cl_test * 2.0 * bref / (0.5 * Sref)
            @test cl_test ≈ CLsurfsp[1] rtol=1e-16

            # Verify gw was also scaled
            @test gw[4] ≈ gc[3]  # Tip condition maintained
        end
    end

    @testset "calculate_influence_coefficient" begin
        # Test basic influence calculation
        normal = aero.Point2D(0.0, 1.0)  # Normal pointing up

        # Point to the right (+y)
        r = aero.Point2D(1.0, 0.0)
        influence = aero.calculate_influence_coefficient(r, normal)
        @test influence ≈ 1.0  # y*nz / r² = 1*1 / 1 = 1

        # Point above (+z)
        r = aero.Point2D(0.0, 1.0)
        influence = aero.calculate_influence_coefficient(r, normal)
        @test influence ≈ 0.0  # y*nz - z*ny / r² = 0*1 - 1*0 = 0

        # Test singularity handling (very close point)
        r_small = aero.Point2D(1e-20, 1e-20)
        influence = aero.calculate_influence_coefficient(r_small, normal)
        @test influence ≈ 0.0  # Should be clipped to avoid infinity
    end

    @testset "Generate WakeSystem" begin
        # Test creation from coordinate vectors
        yp = [0.0, 1.0, 2.0, 3.0]
        zp = [0.0, 0.0, 0.0, 0.0]
        ycp = [0.5, 1.5, 2.5]
        zcp = [0.0, 0.0, 0.0]

        ws = aero.WakeSystem(yp, zp, ycp, zcp)

        @test length(ws.points) == 4
        @test length(ws.elements) == 3
        @test size(ws.influence_matrix) == (3, 4)

        # Test control points match input
        @test all(aero.ctrl_ys(ws) .≈ ycp)
        @test all(aero.ctrl_zs(ws) .≈ zcp)

        # Test error on mismatched dimensions
        @test_throws AssertionError aero.WakeSystem(
            yp, zp[1:2], ycp, zcp)  # Wrong zp length
    end

    @testset "Full Trefftz analysis (integration test)" begin

        ac = load_default_model()
        nsurf = 2
        trefftz_config = aero.TrefftzPlaneConfig(
            aero.SurfaceDiscretization(20, 6, 3),  # Wing
            aero.SurfaceDiscretization(10, 0, 2),  # Tail
            k_tip = 16.0,
            bunch = 0.5,
            wing_root_contraction = 0.2,
            tail_root_contraction = 1.0
        )

        # Create wing and htail structures with test geometry
        wing = ac.wing
        wing.layout.span = 35.486921631434697
        wing.layout.root_span = 3.6067999999999998
        wing.layout.ηs = 10.113772664958887 / wing.layout.span
        wing.layout.z = -1.6764000000000001
        wing.inboard.λ = 1.0
        wing.outboard.λ = 1.0

        htail = ac.htail
        htail.layout.span = 15.958117796995291
        htail.layout.root_span = 1.5240000000000000
        htail.layout.ηs = 1.5240000000000000 / htail.layout.span
        htail.layout.z = 0.0
        htail.inboard.λ = 1.0
        htail.outboard.λ = 1.0

        wing.layout.S = Sref = 124.68530761144433
        wing.layout.span = bref = 35.486921631434697
        po = [1.0, 1.0]
        gammat = [0.15, 0.25]
        gammas = [0.77, 1.0]
        fLo = -0.3
        specifies_CL = true
        CLsurfsp = [1.2502595055642693, 0.011976022033901442]

        # Expected results from Fortran code
        fort_CLsurf = [1.2502595055642693, 0.011976022033901442]
        fort_CLtp = 1.2622355275981709
        fort_CDtp = 0.060382619569389735
        fort_sefftp = 0.83156768339673048

        # Use module-level geometry and work arrays
        geom = aero.TREFFTZ_GEOM

        aero.ensure_trefftz_current!(ac, po, gammat, gammas, trefftz_config)
        second_hash = aero.TREFFTZ_GEOMETRY_HASH[]   
        aero.ensure_trefftz_current!(ac, po, gammat, gammas, trefftz_config)
        @test aero.TREFFTZ_GEOMETRY_HASH[] == second_hash  # No change on same params
        
        ac.wake_system = aero._build_trefftz_geometry!(wing, htail, po, gammat, gammas, trefftz_config)
        
        ifrst = [aero.i_first_wing(trefftz_config), aero.i_first_tail(trefftz_config)]
        ilast = [aero.i_last_wing(trefftz_config), aero.i_last_tail(trefftz_config)]

        aero.calculate_wake_circulation!(aero.gw, aero.TREFFTZ_GEOM.gc, ifrst, ilast, nsurf)
       
        aero.scale_circulation!(aero.TREFFTZ_GEOM.gc, aero.gw, aero.TREFFTZ_GEOM.yp, ifrst, ilast,
                CLsurfsp, bref, Sref, nsurf)


        CLtp, CDtp, sefftp = aero.compute_induced_drag!(aero.vnc, ac.wake_system, 
        aero.TREFFTZ_GEOM.gc, aero.gw, bref, Sref)

        @test fort_CLtp ≈ CLtp
        @test fort_CDtp ≈ CDtp
        @test fort_sefftp ≈ sefftp

    end

end  # @testset "Trefftz calculations"
