# Contains functions for calculating wing geometry, loading, and pitching moment.

"""
    set_wing_geometry!(W, CL, qinf, wing)

Sizes wing area, span, root chord from `q`, `CL`, `W`, `AR` at given point (taken as start-of-cruise in `size_aircraft!`).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `W::Float64`: Aircraft weight.
    - `CL::Float64`: Lift coefficient.
    - `qinf::Float64`: Freestream dynamic head.
    - `wing::TASOPT.structures.Wing`: Wing structure 

See Sections 2.5 and 3.4.1 of the [TASOPT Technical Desc](@ref dreladocs).
"""
function set_wing_geometry!(W, CL, q_inf, wing)
    wing.layout.S = W / (q_inf * CL)
    wing.layout.span = sqrt(wing.layout.S * wing.layout.AR)

    wing.layout.ηs = max(wing.layout.ηs, wing.layout.ηo)

    ηo = wing.layout.ηo
    ηs = wing.layout.ηs

    Kc =
        ηo +
        0.5 * (1.0 + wing.inboard.λ) * (ηs - ηo) +
        0.5 * (wing.inboard.λ + wing.outboard.λ) * (1.0 - ηs)

    wing.layout.root_chord = wing.layout.S / (Kc * wing.layout.span)
    wing.inboard.co = wing.layout.root_chord
    wing.outboard.co = wing.inboard.co * wing.inboard.λ
end # set_wing_geometry


"""
    wing_section_cls(wing, gammat, gammas,
            CL, CLhtail,
	        fduo, fdus, fdut)

Calculates section cl at  eta = ηo,ηs,1. Formerly, `wingcl()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Wing::TASOPT.Wing`: Wing Structure
    - `γt::Float64`, `γs::Float64`: Wing lift distribution "taper" ratios for outer and inner wing sections, respectively.
    - `CL::Float64`, `CLhtail::Float64`: Overall lift coefficient of wing and horizontal tail, respectively.
    - `duo::Float64`, `dus::Float64`, `dut::Float64`: Velocity-change fractions at wing root, break ("snag"), and tip due to fuselage flow.

    **Outputs:**
    - `clo::Float64`, `cls::Float64`, `clt::Float64`: Section lift coefficient at root, wing break ("snag"), and tip.

See Sections 2.5 and 2.6 of the [TASOPT Technical Desc](@ref dreladocs). Called by [`aircraft_drag!`](@ref).
"""
function wing_section_cls(wing, γt, γs, CL, CLhtail, duo, dus, dut)

    cosL = cosd(wing.layout.sweep)

    ηo, ηs = wing.layout.ηo, wing.layout.ηs

    Kc =
        ηo +
        0.5 * (1.0 + wing.inboard.λ) * (ηs - ηo) +
        0.5 * (wing.inboard.λ + wing.outboard.λ) * (1.0 - ηs) #surface area factor S = co*bo*K

    Ko = 1.0 / (wing.layout.AR * Kc) #root chord def'n factor

    Kp =
        ηo + #planform loading factor
        0.5 * (1.0 + γs) * (ηs - ηo) +
        0.5 * (γs + γt) * (1.0 - ηs) +
        wing.fuse_lift_carryover * ηo +
        2.0 * wing.tip_lift_loss * Ko * γt * wing.outboard.λ

    cl1 = (CL - CLhtail) / cosL^2 * (Kc / Kp)

    clo = cl1 / (1.0 + duo)^2
    cls = cl1 * γs / wing.inboard.λ / (1.0 + dus)^2
    clt = cl1 * γt / wing.outboard.λ / (1.0 + dut)^2

    return clo, cls, clt

end # wing_section_cls

"""
    wing_lift_span_factor(wing, gammat, gammas)

Returns the factor ``K_p`` that relates the wing root loading ``p_o`` to the net lift carried by the wing, ``p_o K_p b``
(``p_o K_p b / 2`` per wing), for the spanload of [`wing_loading`](@ref): the root, inboard and outboard trapezoids, the fuselage
lift carryover and the tip lift loss.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.wing`: Wing structure.
    - `gammat::Float64`: Wing outer section taper ratio, `wing.outboard.λ * rclt`.
    - `gammas::Float64`: Wing inner section taper ratio, `wing.inboard.λ * rcls`.

    **Outputs:**
    - `Kp::Float64`: Spanload factor [-].
"""
function wing_lift_span_factor(wing, gammat, gammas)
    ηo, ηs = wing.layout.ηo, wing.layout.ηs

    Kc =
        ηo +
        0.5 * (1.0 + wing.inboard.λ) * (ηs - ηo) +
        0.5 * (wing.inboard.λ + wing.outboard.λ) * (1.0 - ηs)

    Ko = 1.0 / (wing.layout.AR * Kc)

    Kp =
        ηo +
        0.5 * (1.0 + gammas) * (ηs - ηo) +
        0.5 * (gammas + gammat) * (1.0 - ηs) +
        wing.fuse_lift_carryover * ηo +
        2.0 * wing.tip_lift_loss * Ko * gammat * wing.outboard.λ

    return Kp
end # wing_lift_span_factor

"""
    wing_loading(wing, rclt, rcls, N, W, Lhtail)

Computes wing root ("center") loading ``p_o`` to balance the net load. Formerly, `wingpo()`.

```math
N*W - L_{h tail} \times 2*∫p(η) dy + 2ΔL₀ + 2ΔLₜ = N*W - (L_{htail}).
```

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.wing`: Wing structure.
    - `rclt::Float64`: tip  /root cl ratio (clt/clo)
    - `rcls::Float64`: break/root cl ratio (cls/clo)
    - `N::Float64`: Max vertical load factor for wing bending loads
    - `W::Float64`: Aircraft Weight
    - `Lhtail::Float64`: Worst-case (most negative) tail lift expected in the critical sizing case

    **Outputs:**
    - `po::Float64`: Wing root loading magnitude.

See Section 2.6.2 of the [TASOPT Technical Desc](@ref dreladocs).
"""
function wing_loading(wing, rclt, rcls, N, W, Lhtail)

    γt, γs = wing.outboard.λ * rclt, wing.inboard.λ * rcls

    Kp = wing_lift_span_factor(wing, γt, γs)

    po = (N * W - Lhtail) / (Kp * wing.layout.span)

    return po
end # wing_loading

"""
    tail_loading!(tail,S,qne; t_fac = 1.0)

Calculates stabilizer span, root chord, and root loading based on the 
never-exceed dynamic pressure, maximum CL, sweep, and aspect ratio. Formerly, `tailpo!()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `tail::TASOPT.structures.tail`: Tail structure.
    - `S::Float64`: Stabilizer area.
    - `qne::Float64`: Never-exceed dynamic pressure.
    - `t_fac::Float64`: Tail Factor (1 for Htail/Wing, 2 for Vtail).
    
    **Outputs:**
    - `po::Float64`: Stabilizer root loading.
    - `b::Float64`: Stabilizer wingspan.

See [Geometry](@ref geometry) or Section 2.3.2 and 2.9.6 of the [TASOPT Technical Description](@ref dreladocs).
"""
function tail_loading!(tail, S, qne; t_fac = 1.0)
    S < 0 && throw(ArgumentError("tail_loading!: negative tail area S = $S; the sizing iterate has left the feasible region"))
    b  = sqrt(S*tail.layout.AR*t_fac)
    tail.layout.root_chord = S/(0.5*b*(1.0+tail.outboard.λ))
    po = qne*S*tail.CL_max/b * 2.0/(1.0 + tail.outboard.λ)
    tail.outboard.co = tail.layout.root_chord*tail.inboard.λ
    tail.inboard.co = tail.layout.root_chord
    return po,b
end

"""
    wing_CM(b,bs,bo, sweep, Xaxis,
                       λt, λs, γt, γs,
                       AR, fLo, fLt, cmpo, cmps, cmpt)

Calculates components of wing pitching moment (``C_M``) about wing root axis:

``C_M = C_{M,0} + C_{M,1} (C_L - C_{L,surf})``

``ΔC_{m, surf} = ΔC_{m, 0} + dCₘ/dCL × (C_L - C_{L,h})``

Formerly, `surfcm()`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `b::Float64`: Span.
      - `bs::Float64`: Outer panel break span.
      - `bo::Float64`: Root (fuselage) span.
      - `sweep::Float64`: Sweep, degrees.
      -	`Xaxis::Float64`: Surface axis position.
      - `λt::Float64`: Outer-panel chord taper ratio  ct/co.
      - `λs::Float64`: Inner-panel chord taper ratio  cs/co.
      - `γt::Float64`: Outer-panel load  taper ratio  pt/po.
      - `γs::Float64`: Inner-panel load  taper ratio  ps/po.
      - `AR::Float64`: Surface aspect ratio.
      - `fLo::Float64`, `fLt::Float64` : Wing root and tip load adjustment factors.
      - `cmpo::Float64`,`cmps::Float64`,`cmpt::Float64`: Perpendicular sectional lift coefficient at wing root, break ("snag"), and tip.

      **Outputs:**
      - `CM0::Float64`: Zero-lift surface pitching moment.
      - `CM1::Float64`: Surface pitching moment including lift contribution.


See Section 2.6.3 of the [TASOPT Technical Desc](@ref dreladocs).
See also [`wing_profiledrag_scaled`](@ref) and [`wing_profiledrag_direct`](@ref).
"""
function wing_CM(b, bs, bo, sweep, Xaxis, λt, λs, γt, γs, AR, fLo, fLt, cmpo, cmps, cmpt)

    cosL = cosd(sweep)
    tanL = tand(sweep)

    ηo = bo / b
    ηs = bs / b

    Kc = ηo + 0.5 * (1.0 + λs) * (ηs - ηo) + 0.5 * (λs + λt) * (1.0 - ηs)

    Ko = 1.0 / (AR * Kc)

    Kp =
        ηo +
        0.5 * (1.0 + γs) * (ηs - ηo) +
        0.5 * (γs + γt) * (1.0 - ηs) +
        fLo * ηo +
        2.0 * fLt * Ko * γt * λt

    C1 =
        (1.0 + 0.5 * (λs + γs) + λs * γs) * (ηs - ηo) +
        (λs * γs + 0.5 * (λs * γt + γs * λt) + λt * γt) * (1.0 - ηs)

    C2 =
        (1.0 + 2.0 * γs) * (ηs - ηo)^2 +
        (γs + 2.0 * γt) * (1.0 - ηs)^2 +
        3.0 * (γs + γt) * (ηs - ηo) * (1.0 - ηs)

    C3 =
        (cmpo * (3.0 + 2.0 * λs + λs^2) + cmps * (3.0 * λs^2 + 2.0 * λs + 1.0)) *
        (ηs - ηo) +
        (
            cmps * (3.0 * λs^2 + 2.0 * λs * λt + λt^2) +
            cmpt * (3.0 * λt^2 + 2.0 * λs * λt + λs^2)
        ) * (1.0 - ηs)

    CM1 =
        (1.0 / Kp) * (
            ηo * (1.0 + fLo) * (Xaxis - 0.25) + (Xaxis - 0.25) * cosL^2 * C1 / 3.0 -
            (tanL / Ko) * C2 / 12.0 +
            2.0 *
            fLt *
            λt *
            γt *
            (Ko * λt * (Xaxis - 0.25) * cosL^2 - 0.5 * (1.0 - ηo) * tanL)
        )

    CM0 = (cosL^4 / Kc) * C3 / 12.0


    return CM0, CM1
end # wing_CM
