using StaticArrays
using Random

@testset "Cryogenic tank" begin
    p_atm = TASOPT.p_atm
    Δp = 0.001 * p_atm
    @testset "Thermodynamic properties" begin
        #Hydrogen
        species = "H2"

        #Data from NIST for saturated vapor
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [18.232, 20.369, 22.965, 27.314]
        ρs = [0.71653, 1.3322, 2.5133, 6.1715]
        us = [364.55, 372.65, 379.15, 378.1] * 1e3
        hs = [435.25, 448.71, 459.78, 460.2] * 1e3
        
        @testset "Hydrogen gas" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.gas_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.5
                @test props.ρ ≈ ρs[i] atol = 1e-1
                @test props.u ≈ us[i] rtol = 1e-2
                @test props.h ≈ hs[i] rtol = 1e-2

                #Check derivative by finite differences
                props2 = TASOPT.CryoTank.gas_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Data from NIST for saturated liquid
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [18.232, 20.369, 22.965, 27.314]
        ρs = [73.142, 70.848, 67.639, 60.711]
        us = [-20.876, -1.4302, 25.31, 80.495] * 1e3
        hs = [-20.183, 0, 28.306, 88.84] * 1e3

        @testset "Hydrogen liquid" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.liquid_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.5
                @test props.ρ ≈ ρs[i] rtol = 1e-2
                @test props.u ≈ us[i] atol = 2e3
                @test props.h ≈ hs[i] atol = 2e3

                props2 = TASOPT.CryoTank.liquid_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Methane
        species = "CH4"

        #Data from NIST for saturated vapor
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [103.87, 111.67, 120.81, 135.59]
        ρs = [0.96233, 1.8164, 3.4382, 8.1024]
        us = [444.65, 455.05, 466.28, 481.52] * 1e3
        hs = [497.3, 510.83, 525.22, 544.05] * 1e3

        @testset "Methane gas" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.gas_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.6
                @test props.ρ ≈ ρs[i] atol = 1e-1
                @test props.u ≈ us[i] rtol = 1e-2
                @test props.h ≈ hs[i] rtol = 1e-2

                props2 = TASOPT.CryoTank.gas_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end

        #Data from NIST for saturated liquid
        ps = [0.5, 1.0, 2.0, 5.0] * p_atm 
        Tsats = [103.87, 111.67, 120.81, 135.59]
        ρs = [433.5, 422.36, 408.66, 384.63]
        us = [-27.112, 0, 31.803, 85.119] * 1e3
        hs = [-26.995, 0, 32.299, 86.436] * 1e3

        @testset "Methane liquid" begin
            for (i,p) in enumerate(ps)
                props = TASOPT.CryoTank.liquid_properties(species, p)
                @test props.Tsat ≈ Tsats[i] atol = 0.6
                @test props.ρ ≈ ρs[i] rtol = 1e-2
                @test props.u ≈ us[i] atol = 2e3
                @test props.h ≈ hs[i] atol = 2e3

                props2 = TASOPT.CryoTank.liquid_properties(species, p + Δp)
                ρ_pFD = (props2.ρ - props.ρ) / Δp
                u_pFD = (props2.u - props.u) / Δp
                @test props.ρ_p ≈ ρ_pFD rtol = 1e-2
                @test props.u_p ≈ u_pFD rtol = 1e-2
            end
        end
    end

    @testset "Saturated mixtures" begin
        species = "H2"
        
        p = p_atm
        β = 0.5

        @testset "SaturatedMixture assembly" begin
            # gas_properties / liquid_properties → SaturatedPhaseProps stored directly
            # in mix.gas and mix.liquid; verify the constructor assembles them correctly.
            gas = TASOPT.CryoTank.gas_properties(species, p)
            liq = TASOPT.CryoTank.liquid_properties(species, p)
            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)

            @test mix.gas === mix.gas       # sanity: field exists
            @test mix.p ≈ p
            @test mix.T ≈ gas.Tsat         # temperature comes from saturation curve

            # Verify mix.gas and mix.liquid hold the phase properties
            @test mix.gas.Tsat ≈ gas.Tsat
            @test mix.gas.ρ    ≈ gas.ρ
            @test mix.gas.h    ≈ gas.h
            @test mix.gas.u    ≈ gas.u
            @test mix.gas.ρ_p  ≈ gas.ρ_p
            @test mix.gas.u_p  ≈ gas.u_p
            @test mix.liquid.ρ    ≈ liq.ρ
            @test mix.liquid.h    ≈ liq.h
            @test mix.liquid.u    ≈ liq.u
            @test mix.liquid.ρ_p  ≈ liq.ρ_p
            @test mix.liquid.u_p  ≈ liq.u_p

            # Verify the mixture-level scalars against the mixing rule
            x = 1 / (1 + (liq.ρ/gas.ρ) * (β / (1 - β)))
            ρ = 1 / (x / gas.ρ + (1-x) / liq.ρ)
            u = x * gas.u + (1-x) * liq.u
            h = x * gas.h + (1-x) * liq.h
            x_p = (-x/gas.ρ^2 * gas.ρ_p - (1-x)/liq.ρ^2 * liq.ρ_p) / (1/liq.ρ - 1/gas.ρ)
            u_p = x * gas.u_p + (1 - x) * liq.u_p + x_p * (gas.u - liq.u)
            ϕ      = 1 / (u_p * ρ)
            hvap   = gas.h - liq.h
            ρ_star = gas.ρ / (liq.ρ - gas.ρ)

            @test mix.ρ     ≈ ρ
            @test mix.x     ≈ x
            @test mix.h     ≈ h
            @test mix.u     ≈ u
            @test mix.β     ≈ β
            @test mix.u_p   ≈ u_p
            @test mix.ϕ     ≈ ϕ
            @test mix.hvap  ≈ hvap
            @test mix.ρ_star ≈ ρ_star
        end
        @testset "Mixture functions" begin     

            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)

            #Test update to mixture by comparing against new object
            p2 = 2 * p_atm 
            β2 = 0.3
            TASOPT.CryoTank.update_pβ!(mix, p2, β2)

            mix2 = TASOPT.CryoTank.SaturatedMixture(species, p2, β2)

            @test mix2.p ≈ mix.p
            @test mix2.T ≈ mix.T
            @test mix2.ρ ≈ mix.ρ
            @test mix2.x ≈ mix.x
            @test mix2.h ≈ mix.h
            @test mix2.u ≈ mix.u
            @test mix2.β ≈ mix.β
            @test mix2.u_p ≈ mix.u_p
            @test mix2.ϕ ≈ mix.ϕ
            @test mix2.hvap ≈ mix.hvap
            @test mix2.ρ_star ≈ mix.ρ_star

            #Check that density is the same during pressure change with convert_β_same_ρ
            mix = TASOPT.CryoTank.SaturatedMixture(species, p, β)
            ρorig = mix.ρ

            β3 = TASOPT.CryoTank.convert_β_same_ρ(species, p2, p, β)
            TASOPT.CryoTank.update_pβ!(mix, p2, β3)
            @test mix.ρ ≈ ρorig

        end
    end

    species = "h2"
    p = 2.0 * p_atm
    β = 0.5 
    mix = TASOPT.CryoTank.SaturatedMixture("h2", p, β)

    Q = 2e4
    mdot = 1.0
    W = 0.0
    mdot_vent = 0.0

    xout = 0.0 #Taking liquid out
    xvent = 1.0 #Venting gas
    α = 2.0
    V = 100.0
    pmax = 5*p_atm

    @testset "Model derivatives" begin
        #Check function for dp/dt
        dp_dt = TASOPT.CryoTank.dpdt(mix, Q, W, mdot, xout, mdot_vent, xvent, V, α)
        dp_dt_validate = α * mix.ϕ / V * (Q + W - mdot * mix.hvap * (xout + mix.ρ_star) - mdot_vent * mix.hvap * (xvent + mix.ρ_star))
        dp_dt_regression = 6.995703982056127

        @test dp_dt ≈ dp_dt_validate
        @test dp_dt ≈ dp_dt_regression

        #Check function for dβ/dt
        dβ_dt = TASOPT.CryoTank.dβdt(mix, dp_dt_regression, mdot + mdot_vent, V)
        dβ_dt_validate = (-(mdot + mdot_vent) / V - dp_dt_validate * (β * mix.liquid.ρ_p + (1-β)*mix.gas.ρ_p)) / (mix.liquid.ρ - mix.gas.ρ)
        dβ_dt_regression = -0.00015294953663597238

        @test dβ_dt ≈ dβ_dt_validate
        @test dβ_dt ≈ dβ_dt_regression

        #Check function for venting mass flow rate
        mdot_vent_out = TASOPT.CryoTank.venting_mass_flow(mix, Q, W, mdot, xout, xvent) #This venting should give zero dp/dt
        vent_validate = (Q - mdot * mix.hvap * (xout + mix.ρ_star)) / (mix.hvap * (xvent + mix.ρ_star))
        vent_regression = 0.007274844258979718

        @test mdot_vent_out ≈ vent_validate
        @test mdot_vent_out ≈ vent_regression

        dpdt_venting = TASOPT.CryoTank.dpdt(mix, Q, W, mdot, xout, mdot_vent_out, xvent, V, α) #Check pressure derivative with venting; it should be 0
        dpdt_venting_sanity = 0.0

        @test dpdt_venting ≈ dpdt_venting_sanity

        #Check function for boiloff rate
        mdot_liq = (1 - xout) * mdot + (1 - xvent) * mdot_vent
        mdotboil = TASOPT.CryoTank.mdot_boiloff(mix, dβ_dt_regression, dp_dt_regression, mdot_liq, V)
        mdotboil_validate = - (dβ_dt_regression * V * mix.liquid.ρ + β * V * mix.liquid.ρ_p * dp_dt_regression + mdot)
        mdotboil_regression = 0.04266684331352266

        @test mdotboil ≈ mdotboil_validate
        @test mdotboil ≈ mdotboil_regression

        #Check function that computes all derivatives for inputs
        Q_calc(t) = Q
        W_calc(t) = W
        mdot_calc(t) = mdot

        u = TASOPT.CryoTank.tank_inputs(Q_calc, W_calc, mdot_calc)
        params = TASOPT.CryoTank.tank_params(mix, V, pmax, xout, xvent, α)

        y = @SVector [p, β, V * mix.ρ, 0.0, 0.0, 0.0]
        t = 0.0

        dydt = TASOPT.CryoTank.TankDerivatives(t, y, u, params)
        @test dydt[1] ≈ dp_dt_regression
        @test dydt[2] ≈ dβ_dt_regression
        @test dydt[3] ≈ -mdot
        @test dydt[4] ≈ mdot
        @test dydt[5] ≈ 0.0
        @test dydt[6] ≈ mdotboil_regression
    end
end

# ── Property-based tests ────────────────────────────────────────────────────
# Inspired by https://fsharpforfunandprofit.com/posts/property-based-testing-2/
# Strategy: test *relationships* and *invariants* that must hold for all valid
# inputs.
@testset "SaturatedPhaseProps property tests" begin
    rng = MersenneTwister(42)
    p_atm_val = TASOPT.p_atm
    p_lo = 0.1 * p_atm_val   # lower bound of polynomial fit domain
    p_hi = 10.0 * p_atm_val  # upper bound

    for species in ("H2", "CH4")
        @testset "$(species): saturation temperature" begin
            # gas_properties and liquid_properties are both fits to the same
            # saturation curve - Tsat must be identical at any given pressure.
            for _ in 1:10
                p = p_lo + rand(rng) * (p_hi - p_lo)
                gp = TASOPT.CryoTank.gas_properties(species, p)
                lp = TASOPT.CryoTank.liquid_properties(species, p)
                @test gp.Tsat ≈ lp.Tsat rtol = 1e-12
            end
        end

        @testset "$(species): Tsat monotone with pressure (Clausius-Clapeyron)" begin
            for _ in 1:10
                p1 = p_lo + rand(rng) * (p_hi - p_lo - 0.1 * p_atm_val)
                p2 = p1 + 0.05 * p_atm_val + rand(rng) * 0.5 * p_atm_val
                gp1 = TASOPT.CryoTank.gas_properties(species, p1)
                gp2 = TASOPT.CryoTank.gas_properties(species, p2)
                @test gp2.Tsat > gp1.Tsat
            end
        end

        @testset "$(species): thermodynamic identity h = u + p/ρ" begin
            # Follows from the definition h ≡ u + pv, v = 1/ρ.
            # h and u are independent NIST polynomial fits, so the identity holds
            # to ~1% inside the fit range (0.5–9 atm) but degrades at the edges.
            # We stay away from the edges.
            p_lo_inner = 0.5 * p_atm_val
            p_hi_inner = 9.0 * p_atm_val
            for _ in 1:10
                p = p_lo_inner + rand(rng) * (p_hi_inner - p_lo_inner)
                for props in (TASOPT.CryoTank.gas_properties(species, p),
                              TASOPT.CryoTank.liquid_properties(species, p))
                    @test props.h ≈ props.u + p / props.ρ rtol = 1e-2
                end
            end
        end

        @testset "$(species): derivative consistency" begin
            # ρ_p and u_p are analytic derivatives of the same polynomials.
            # Try one-sided finite difference to within 0.05 %.
            δ = 10.0  # 10 Pa — small relative to ~1e5 Pa but large enough to avoid cancellation
            for _ in 1:10
                p = p_lo + rand(rng) * (p_hi - p_lo - δ)
                for (f, label) in ((TASOPT.CryoTank.gas_properties, "gas"),
                                   (TASOPT.CryoTank.liquid_properties, "liquid"))
                    props  = f(species, p)
                    props2 = f(species, p + δ)
                    @test props.ρ_p ≈ (props2.ρ - props.ρ) / δ rtol = 5e-4
                    @test props.u_p ≈ (props2.u - props.u) / δ rtol = 5e-4
                end
            end
        end
    end

    
end

# ── Scaled fuselage cross-section ───────────────────────────────────────────
# `scaled_cross_section(cs, R)` returns (perimeter, area) of a cross-section that
# is geometrically similar to `cs` but has radius `R`. It is used by the tank
# mechanical and thermal sizing (size_inner_tank, size_outer_tank, residuals_Q).
# Strategy: (i) compare against closed-form geometry written independently of
# `area`/`get_perimeter`, and (ii) test invariants of a similarity transform.
@testset "Scaled cross-section" begin
    scs  = TASOPT.scaled_cross_section
    SB   = TASOPT.SingleBubble
    MB   = TASOPT.MultiBubble
    get_perimeter = TASOPT.structures.get_perimeter
    cs_area       = TASOPT.structures.area

    # Independent closed-form references -------------------------------------
    # Single bubble: circle of radius R whose lower half is moved down by ΔR,
    # joined by a rectangle of size 2R x ΔR.
    sb_perimeter(R, ΔR) = 2π * R + 2ΔR
    sb_area(R, ΔR)      = π * R^2 + 2R * ΔR

    # Double bubble (one web): union of two circles of radius R whose centres are
    # at y = ±dy, plus the same 2R x ΔR block. The union area is built from the
    # circle-circle "lens" overlap, and the outer perimeter from the arc that
    # lies outside the other circle, so this does not reuse TASOPT's web_geometry.
    function db_reference(R, ΔR, dy)
        d     = 2dy                                  # centre distance
        lens  = 2R^2 * acos(d / (2R)) - (d / 2) * sqrt(4R^2 - d^2)
        area  = 2π * R^2 - lens + 2R * ΔR
        φ     = acos(dy / R)                         # half-angle buried inside the other circle
        perim = 2 * R * (2π - 2φ) + 2ΔR
        return perim, area
    end

    sb  = SB(radius = 2.5, bubble_lower_downward_shift = 0.3)
    mb  = MB(radius = 2.0, bubble_lower_downward_shift = 0.4,
             bubble_center_y_offset = 0.5, n_webs = 1)

    @testset "SingleBubble: closed form" begin
        for R in (0.5, 1.0, 2.5, 3.7, 6.0)
            p, A = scs(sb, R)
            k    = sb.bubble_lower_downward_shift / sb.radius   # similarity ratio ΔR/R
            @test p ≈ sb_perimeter(R, k * R)
            @test A ≈ sb_area(R, k * R)
        end
        # Hand-computed value: R = 2.5, ΔR = 0.3
        p, A = scs(sb, 2.5)
        @test p ≈ 5π + 0.6
        @test A ≈ 6.25π + 1.5
    end

    @testset "SingleBubble: identity and ratio" begin
        # Scaling to its own radius must return the unscaled geometry
        p, A = scs(sb, sb.radius)
        @test p ≈ get_perimeter(sb)
        @test A ≈ cs_area(sb)

        # A circle (ΔR = 0) must reduce to the circle formulas
        circle = SB(radius = 1.7, bubble_lower_downward_shift = 0.0)
        for R in (0.3, 1.7, 4.0)
            p, A = scs(circle, R)
            @test p ≈ 2π * R
            @test A ≈ π * R^2
        end
    end

    @testset "MultiBubble: closed form" begin
        for (R0, ΔR0, dy0) in ((2.0, 0.4, 0.5), (3.0, 0.2, 1.0), (1.5, 0.0, 0.3))
            cs = MB(radius = R0, bubble_lower_downward_shift = ΔR0,
                    bubble_center_y_offset = dy0, n_webs = 1)
            for s in (0.4, 1.0, 1.9)                 # scale factors
                R = s * R0
                p, A = scs(cs, R)
                p_ref, A_ref = db_reference(R, s * ΔR0, s * dy0)
                @test p ≈ p_ref
                @test A ≈ A_ref
            end
        end
        # Special geometry with a hand-computed answer: dy = R/2 -> θ = 30°,
        # so the web sits at sin(θ) = 1/2 and sin(2θ) = √3/2.
        cs = MB(radius = 2.0, bubble_lower_downward_shift = 0.0,
                bubble_center_y_offset = 1.0, n_webs = 1)
        p, A = scs(cs, 2.0)
        @test p ≈ (2π + 4 * π / 6) * 2.0
        @test A ≈ (π + (2π / 6 + √3 / 2)) * 2.0^2
    end

    @testset "MultiBubble: identity and limits" begin
        p, A = scs(mb, mb.radius)
        @test p ≈ get_perimeter(mb)
        @test A ≈ cs_area(mb)

        # With the bubble centres on top of each other (dy -> 0) the double
        # bubble degenerates to a single bubble.
        mb_dy0 = MB(radius = 2.0, bubble_lower_downward_shift = 0.4,
                    bubble_center_y_offset = 1e-12, n_webs = 1)
        sb_eq  = SB(radius = 2.0, bubble_lower_downward_shift = 0.4)
        for R in (1.0, 2.0, 3.5)
            @test collect(scs(mb_dy0, R)) ≈ collect(scs(sb_eq, R)) rtol = 1e-9
        end

        # A MultiBubble with no webs is also a single bubble
        mb_noweb = MB(radius = 2.0, bubble_lower_downward_shift = 0.4,
                      bubble_center_y_offset = 0.5, n_webs = 0)
        for R in (1.0, 2.0, 3.5)
            @test collect(scs(mb_noweb, R)) ≈ collect(scs(sb_eq, R))
        end

        # A double bubble is larger than a single bubble of the same radius, so it
        # encloses more area and has a longer perimeter.
        p_mb, A_mb = scs(mb, 2.0)
        p_sb, A_sb = scs(sb_eq, 2.0)
        @test A_mb > A_sb
        @test p_mb > p_sb
    end

    cross_sections = (
        "SingleBubble"            => sb,
        "SingleBubble (ΔR = 0)"   => SB(radius = 1.9, bubble_lower_downward_shift = 0.0),
        "MultiBubble, 1 web"      => mb,
        "MultiBubble, 2 webs"     => MB(radius = 2.2, bubble_lower_downward_shift = 0.3,
                                        bubble_center_y_offset = 0.6, n_webs = 2),
    )

    @testset "Similarity invariants: $name" for (name, cs) in cross_sections
        p0, A0 = scs(cs, cs.radius)
        radii  = (0.25, 0.8, 1.3, 2.0, 5.0) .* cs.radius

        @testset "Analytic sensitivities: dp/dR = p/R, dA/dR = 2A/R" begin
            R  = 1.4 * cs.radius
            h  = 1e-6 * R
            p, A   = scs(cs, R)
            pp, Ap = scs(cs, R + h)
            pm, Am = scs(cs, R - h)
            @test (pp - pm) / (2h) ≈ p / R rtol = 1e-6
            @test (Ap - Am) / (2h) ≈ 2A / R rtol = 1e-6
        end

        @testset "Returns two Float64 values" begin
            out = scs(cs, 1.1 * cs.radius)
            @test out isa Tuple{Float64, Float64}
        end
    end

    @testset "Fuselage cross-section in the tank sizing path" begin
        # The default aircraft has a SingleBubble fuselage; the tank code passes
        # `fuse.layout.cross_section` straight to scaled_cross_section.
        ac_cs = load_default_model()
        cs    = ac_cs.fuselage.layout.cross_section
        @test cs isa SB
        p, A = scs(cs, cs.radius)
        @test p ≈ get_perimeter(cs)
        @test A ≈ cs_area(cs)
        # Inner tank radius must always yield smaller area/perimeter than the fuselage
        p_in, A_in = scs(cs, cs.radius - 0.1)
        @test p_in < p
        @test A_in < A
    end
end
