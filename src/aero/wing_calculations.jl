"""
    wingpo(wing, rclt, rcls, N, W, Lhtail)

Computes wing root ("center") loading ``p_o`` to balance the net load.

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
function wingpo(wing, rclt, rcls, N, W, Lhtail)

    γt, γs = wing.outboard.λ * rclt, wing.inboard.λ * rcls

    ηo, ηs = wing.layout.ηo, wing.layout.ηs

    Kc =
        ηo +
        0.5 * (1.0 + wing.inboard.λ) * (ηs - ηo) +
        0.5 * (wing.inboard.λ + wing.outboard.λ) * (1.0 - ηs)

    Ko = 1.0 / (wing.layout.AR * Kc)

    Kp =
        ηo +
        0.5 * (1.0 + γs) * (ηs - ηo) +
        0.5 * (γs + γt) * (1.0 - ηs) +
        wing.fuse_lift_carryover * ηo +
        2.0 * wing.tip_lift_loss * Ko * γt * wing.outboard.λ

    po = (N * W - Lhtail) / (Kp * wing.layout.span)

    return po
end # wingpo

"""
    tailpo!(tail,S,qne; t_fac = 1.0)

Calculates stabilizer span, root chord, and root loading based on the 
never-exceed dynamic pressure, maximum CL, sweep, and aspect ratio.

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
function tailpo!(tail, S, qne; t_fac = 1.0)

    b  = sqrt(S*tail.layout.AR*t_fac)
    tail.layout.root_chord = S/(0.5*b*(1.0+tail.outboard.λ))
    po = qne*S*tail.CL_max/b * 2.0/(1.0 + tail.outboard.λ)
    tail.outboard.co = tail.layout.root_chord*tail.inboard.λ
    tail.inboard.co = tail.layout.root_chord
    return po,b
end

"""
    wingcl(wing, gammat, gammas,
            CL, CLhtail,
	        fduo, fdus, fdut)

Calculates section cl at  eta = ηo,ηs,1

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
function wingcl(wing, γt, γs, CL, CLhtail, duo, dus, dut)

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

end # wingcl

"""
    set_wing_geometry!(W,CL,qinf,wing)

Sizes wing area, span, root chord from `q`, `CL`, `W`, `AR` at given point (taken as start-of-cruise in `wsize`).

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
end # wingsc

"""
      wing_profiledrag_direct(wing, γt, γs,
            Mach, CL, CLhtail, 
            Reco, aRexp, kSuns, fexcd,
            fduo, fdus, fdut)

Calculates wing or tail surface profile `CD` by calculating the performance of wing segments explicitly via airfoil data (found in [`./src/air/C.air`] and accessed by [`airfun`], [`airtable`]).
Called by [`aircraft_drag!`](@ref) if `computes_wing_direct` flag set to true.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `Wing::TASOPT.Wing`: Wing Structure.
      - `γt::Float64`: Outer-panel load  taper ratio  pt/po.
      - `γs::Float64`: Inner-panel load  taper ratio  ps/po.
      - `Mach::Float64`: Mach number.
      - `CL::Float64`: Wing CL.
      - `CLhtail::Float64`: Htail CL.
      - `Reco::Float64`: Reynolds number for co.
      - `aRexp::Float64`: Re-scaling exponent.
      - `kSuns::Float64`: Shock-unsweep area constant.
      - `fexcd::Float64`: Excrescence cd scaling factor.
      - `fduo::Float64`, `fdus::Float64`, `fdut::Float64`: Velocity-change fractions at wing root, break ("snag"), and tip due to fuselage flow.

      **Outputs:**
      - `clpo::Float64`,`clps::Float64`,`clpt::Float64`: Perpendicular sectional lift coefficient at wing root, break ("snag"), and tip.
      - `CDfwing::Float64`: Friction profile cd in perp. plane.
      - `CDpwing::Float64`: Pressure profile cd in perp. plane.
      - `CDwing::Float64`: Overall profile CD.
      - `CDover::Float64`: Fuselage added CD due to lift carryover.

See Sections 2.14.3 and 3.8.3 of TASOPT Technical Desc.
See also [`aircraft_drag!`](@ref), [`wing_profiledrag_scaled`](@ref), [`surfcm`], and [`airfun`].
"""
function wing_profiledrag_direct(
      wing, γt, γs,
      Mach, CL, CLhtail, 
      Reco, aRexp, kSuns, fexcd,
      fduo, fdus, fdut)

      #     n = 4     #  ~0.30% error
      #     n = 6     #  ~0.13% error
      n::Int = 8     #  ~0.07% error
      #     n = 12    #  ~0.07% error

      cosL = cosd(wing.layout.sweep)

      AR = wing.layout.span^2 / wing.layout.S

      ηo = wing.layout.ηo
      ηs = wing.layout.ηs

      #---- set center clpo
      Kc = ηo +
           0.5 * (1.0 + wing.inboard.λ) * (ηs - ηo) +
           0.5 * (wing.inboard.λ + wing.outboard.λ) * (1.0 - ηs)

      Kp0 = ηo +
            0.5 * (1.0 + γs) * (ηs - ηo) +
            0.5 * (γs + γt) * (1.0 - ηs)

      Ko = 1.0 / (Kc * AR)
      Kp = Kp0 + wing.fuse_lift_carryover * ηo + 2.0 * wing.tip_lift_loss * Ko * γt * wing.outboard.λ
      
      clp1 = (CL - CLhtail) / cosL^2 * wing.layout.S / (Kp * wing.layout.span * wing.layout.root_chord)

      #---- set break and tip clp for passing back
      clpo = clp1 / (1.0 + fduo)^2
      clps = clp1 * γs / wing.inboard.λ / (1.0 + fdus)^2
      clpt = clp1 * γt / wing.outboard.λ / (1.0 + fdut)^2

      #c---- area of exposed wing
      #      Swing = co*b * ( 0.5*(1.0    +λs)*(ηs-ηo) + 0.5*(λs+λt)*(1.0 -ηs) )

      #---- fraction of exposed wing subject to wing-root shock unsweep
      CDfwing = 0.0
      CDpwing = 0.0
      CDwing = 0.0
      Snorm = 0.0
      ARe = wing.airsection.Re

      #remember \ will give you a float which makes the iterator
      # also work on floats, which can make some cases type unstable.
      @inbounds for i = 1:n÷2 
            frac = (i - 0.5) / (n / 2)
            η = ηo * (1.0 - frac) + ηs * frac
            dη = (ηs - ηo) / (n / 2)
            P = 1.0*(1.0 - frac) + γs * frac
            C = 1.0*(1.0 - frac) + wing.inboard.λ * frac
            toc = wing.inboard.cross_section.thickness_to_chord * (1.0 - frac) + wing.outboard.cross_section.thickness_to_chord * frac
            fdu = fduo * (1.0 - frac) + fdus * frac

            #wing root shock "unsweep" function
            fSuns = exp(-(η - ηo) * wing.layout.span / (kSuns * C * 2.0 * wing.layout.root_chord))

            clp = clp1 * (P / C) / (1.0 + fdu)^2
            Rec = Reco * C * (1.0 + fdu)
            Mperp = Mach * cosL * (1.0 + fdu)

            cdf1, cdp1, cdwbar, cm = airfun(clp, toc, Mperp, wing.airsection)
            #println(cdf1, " ", cdp1, " ", cm)

            Refac = (Rec / ARe)^aRexp
            cdf = cdf1 * Refac * fexcd * (1.0 + fdu)^3
            cdp = cdp1 * Refac * fexcd * (1.0 + fdu)^3 *
                  (fSuns + (1.0 - fSuns) * cosL^2) * cosL

            cd = cdf + cdp

            Snorm = Snorm + C * dη
            CDwing = CDwing + C * dη * cd

            CDfwing = CDfwing + C * dη * cdf
            CDpwing = CDpwing + C * dη * cdp
      end

      @inbounds for i = n/2+1:n
            frac = (i - n / 2 - 0.5) / (n / 2)
            η = ηs * (1.0 - frac) + 1.0 * frac
            dη = (1.0 - ηs) / (n / 2)
            P = γs * (1.0 - frac) + γt * frac
            C = wing.inboard.λ * (1.0 - frac) + wing.outboard.λ * frac
            toc = wing.outboard.cross_section.thickness_to_chord * (1.0 - frac) + wing.outboard.cross_section.thickness_to_chord * frac
            fdu = fdus * (1.0 - frac) + fdut * frac

            fSuns = exp(-(η - ηo) * wing.layout.span / (kSuns * C * 2.0 * wing.layout.root_chord))

            clp = clp1 * (P / C) / (1.0 + fdu)^2
            Rec = Reco * C * (1.0 + fdu)
            Mperp = Mach * cosL * (1.0 + fdu)

            cdf1, cdp1, cdwbar, cm = airfun(clp, toc, Mperp, wing.airsection)

            Refac = (Rec / ARe)^aRexp
            cdf = cdf1 * Refac * fexcd * (1.0 + fdu)^3
            cdp = cdp1 * Refac * fexcd * (1.0 + fdu)^3 *
                  (fSuns + (1.0 - fSuns) * cosL^2) * cosL

            cd = cdf + cdp

            Snorm = Snorm + C * dη
            CDwing = CDwing + C * dη * cd

            CDfwing = CDfwing + C * dη * cdf
            CDpwing = CDpwing + C * dη * cdp
      end

      CDwing = CDwing / Snorm
      CDfwing = CDfwing / Snorm
      CDpwing = CDpwing / Snorm

      #c---- fuselage carryover CD
      #      Cdiss = 0.001
      #      dCp = (Kc/Kp)*(CL-CLhtail)*(1.0+fLo)
      #      CDover = 2.0*Cdiss*(ηo/Kc)*0.5*dCp*(3.0 + 0.0625*dCp^2)
      CDover = 0.0

      return clpo, clps, clpt, CDfwing, CDpwing, CDwing, CDover

end # wing_profiledrag_direct

"""
    wing_profiledrag_scaled(S, 
    b, bs, bo, 
    λt, λs, sweep, 
    co, cdf, cdp, 
    Reco, Reref, 
    aRexp, kSuns, fCDcen)

Computes wing or tail surface profile CD from pre-computed chord quantities and corrections.
Called by [`aircraft_drag!`](@ref) if `computes_wing_direct` flag set to false.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `S::Float64`: reference area.
      - `b::Float64`: span.
      - `bs::Float64`: outer panel break span.
      - `bo::Float64`: wing-root (fuselage) span.
      - `λt::Float64`: outer-panel taper ratio  ``ct/co``.
      - `λs::Float64`: inner-panel taper ratio  ``cs/co``.
      - `sweep::Float64`: wing sweep, degrees.
      - `co::Float64`: wing root chord.
      - `cdf::Float64`: friction profile cd.
      - `cdp::Float64`: pressure profile cd.
      - `Reco::Float64`: Reynolds number for co.
      - `Reref::Float64`: reference Reynolds number for cd scaling.
      - `aRexp::Float64`: Re-scaling exponent.
      - `kSuns::Float64`: shock-unsweep area constant.
      - `fCDcen::Float64` : related to fraction of wing BLI (see Eqns. 619 - 621).

      **Outputs:**
      - `CDsurf`: overall profile CD.
      - `CDover`: fuselage added CD due to lift carryover.

See Sections 2.14.3 and 3.8.3 of the [TASOPT Technical Desc](@ref dreladocs).

"""
function wing_profiledrag_scaled(S,
      b, bs, bo, λt, λs, sweep, co,
      cdf, cdp, Reco, Reref, aRexp, kSuns,
      fCDcen)


      cosL = cosd(sweep)

      ηo = bo / b
      ηs = bs / b

      #---- area of exposed surf
      Ssurf = co * 0.5 * b * ((1.0 + λs) * (ηs - ηo) +
                              (λs + λt) * (1.0 - ηs))

      #---- fraction of exposed surf subject to root shock unsweep
      fSuns = kSuns * co^2 / (0.5 * Ssurf)

      #---- overall profile cd with Re and sweep effects, referenced to root chord co
      Refac = (Reco / Reref)^aRexp
      cdo = (cdf + cdp * (fSuns + (1.0 - fSuns) * cosL^2) * cosL) * Refac
      #
      #---- integrate cd over tapered surf, assuming  cd ~ Re^aRexp  dependence
      lsxa = λs^aRexp
      lxa = λt^aRexp
      #
      if (abs(λs - 1.0) < 0.02)
            #----- asymptotic form for λs ~ 1
            lsfac = 1.0 - aRexp * (1.0 - λs) / (1.0 + λs)
      else
            #----- exact form, but singular if λs = 1
            lsfac = 2.0 / (2.0 + aRexp) *
                    (1.0 - lsxa * λs^2) /
                    (1.0 - λs^2)
      end
      #
      if (abs(λt - λs) < 0.02)
            #----- asymptotic form for λt ~ λs
            lfac = 1.0 - aRexp * (λs - λt) / (λs + λt)
      else
            #----- exact form, but singular if λt = λs
            lfac = 2.0 / (2.0 + aRexp) *
                   (lsxa * λs^2 - lxa * λt^2) /
                   (λs^2 - λt^2)
      end

      CDsurf = (co * 0.5 * b / S) * cdo *
               (2.0 * ηo * fCDcen +
                (1.0 + λs) * (ηs - ηo) * lsfac +
                (λs + λt) * (1.0 - ηs) * lfac)

      CDover = 0.0

      return CDsurf, CDover
end # wing_profiledrag_scaled


"""
    surfcm(b,bs,bo, sweep, Xaxis,
                       λt, λs, γt, γs,
                       AR, fLo, fLt, cmpo, cmps, cmpt)

Calculates components of wing pitching moment (``C_M``) about wing root axis:

``C_M = C_{M,0} + C_{M,1} (C_L - C_{L,surf})``

``ΔC_{m, surf} = ΔC_{m, 0} + dCₘ/dCL × (C_L - C_{L,h})``

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
function surfcm(b, bs, bo, sweep, Xaxis, λt, λs, γt, γs, AR, fLo, fLt, cmpo, cmps, cmpt)

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
end # surfcm
