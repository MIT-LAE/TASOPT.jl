"""
      surfcd2(wing, γt, γs,
            Mach, CL, CLhtail, 
            Reco, aRexp, kSuns, fexcd,
            fduo, fdus, fdut)

Calculates wing or tail surface profile `CD` by calculating the performance of wing segments explicitly via airfoil data (found in [`./src/air/C.air`] and accessed by [`airfun`], [`airtable`]).
Called by [`cdsum!`](@ref) if `icdfun` flag set to 1.

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
See also [`cdsum!`](@ref), [`surfcd`](@ref), [`surfcm`], and [`airfun`].

!!! compat "Future Changes" 
      This function will be renamed for clarity of use.
"""
function surfcd2(
      wing, γt, γs,
      Mach, CL, CLhtail, 
      Reco, aRexp, kSuns, fexcd,
      fduo, fdus, fdut)

      #     n = 4     #  ~0.30% error
      #     n = 6     #  ~0.13% error
      n::Int = 8     #  ~0.07% error
      #     n = 12    #  ~0.07% error

      #     call tset(time0)

      cosL = cosd(wing.layout.sweep)

      AR = wing.layout.b^2 / wing.layout.S

      ηo = wing.outboard.layout.b / wing.layout.b
      ηs = wing.inboard.layout.b / wing.layout.b

      #---- set center clpo
      Kc = ηo +
           0.5 * (1.0 + wing.inboard.layout.λ) * (ηs - ηo) +
           0.5 * (wing.inboard.layout.λ + wing.outboard.layout.λ) * (1.0 - ηs)

      Kp0 = ηo +
            0.5 * (1.0 + γs) * (ηs - ηo) +
            0.5 * (γs + γt) * (1.0 - ηs)

      Ko = 1.0 / (Kc * AR)
      Kp = Kp0 + wing.inboard.lift_rolloff * ηo + 2.0 * wing.outboard.lift_rolloff * Ko * γt * wing.outboard.layout.λ
      
      clp1 = (CL - CLhtail) / cosL^2 * wing.layout.S / (Kp * wing.layout.b * wing.layout.chord)

      #---- set break and tip clp for passing back
      clpo = clp1 / (1.0 + fduo)^2
      clps = clp1 * γs / wing.inboard.layout.λ / (1.0 + fdus)^2
      clpt = clp1 * γt / wing.outboard.layout.λ / (1.0 + fdut)^2

      #c---- area of exposed wing
      #      Swing = co*b * ( 0.5*(1.0    +λs)*(ηs-ηo)
      #     &               + 0.5*(λs+λt)*(1.0 -ηs) )

      #      write(*,*)
      #      write(*,*) CL, CLhtail, S/(Kp*b*co)

      #---- fraction of exposed wing subject to wing-root shock unsweep
      CDfwing = 0.0
      CDpwing = 0.0
      CDwing = 0.0
      Snorm = 0.0
      ARe = airsection.Re

      for i = 1:n/2
            frac = (float(i) - 0.5) / float(n / 2)
            η = ηo * (1.0 - frac) + ηs * frac
            dη = (ηs - ηo) / float(n / 2)
            P = 1.0 - frac + γs * frac
            C = 1.0 - frac + wing.inboard.layout.λ * frac
            toc = wing.inboard.layout.chord_thickness * (1.0 - frac) + wing.outboard.layout.chord_thickness * frac
            fdu = fduo * (1.0 - frac) + fdus * frac

            #wing root shock "unsweep" function
            fSuns = exp(-(η - ηo) * wing.layout.b / (kSuns * C * 2.0 * wing.layout.chord))

            clp = clp1 * (P / C) / (1.0 + fdu)^2
            Rec = Reco * C * (1.0 + fdu)
            Mperp = Mach * cosL * (1.0 + fdu)

            cdf1, cdp1, cdwbar, cm = airfun(clp, toc, Mperp, airsection)
            #println(cdf1, " ", cdp1, " ", cm)

            Refac = (Rec / ARe)^aRexp
            cdf = cdf1 * Refac * fexcd * (1.0 + fdu)^3
            cdp = cdp1 * Refac * fexcd * (1.0 + fdu)^3 *
                  (fSuns + (1.0 - fSuns) * cosL^2) * cosL

            cd = cdf + cdp

            #        write(*,'(1x,8f12.6)') η, fSuns, clp, cdf+cdp*cosL^3, cd

            Snorm = Snorm + C * dη
            CDwing = CDwing + C * dη * cd

            CDfwing = CDfwing + C * dη * cdf
            CDpwing = CDpwing + C * dη * cdp
      end

      for i = n/2+1:n
            frac = (float(i - n / 2) - 0.5) / float(n / 2)
            η = ηs * (1.0 - frac) + 1.0 * frac
            dη = (1.0 - ηs) / float(n / 2)
            P = γs * (1.0 - frac) + γt * frac
            C = wing.inboard.layout.λ * (1.0 - frac) + wing.outboard.layout.λ * frac
            toc = wing.outboard.layout.chord_thickness * (1.0 - frac) + wing.outboard.layout.chord_thickness * frac
            fdu = fdus * (1.0 - frac) + fdut * frac

            fSuns = exp(-(η - ηo) * wing.layout.b / (kSuns * C * 2.0 * wing.layout.chord))

            clp = clp1 * (P / C) / (1.0 + fdu)^2
            Rec = Reco * C * (1.0 + fdu)
            Mperp = Mach * cosL * (1.0 + fdu)

            cdf1, cdp1, cdwbar, cm = airfun(clp, toc, Mperp, airsection)

            Refac = (Rec / ARe)^aRexp
            cdf = cdf1 * Refac * fexcd * (1.0 + fdu)^3
            cdp = cdp1 * Refac * fexcd * (1.0 + fdu)^3 *
                  (fSuns + (1.0 - fSuns) * cosL^2) * cosL

            cd = cdf + cdp

            #        write(*,'(1x,8f12.6)') η, fSuns, clp, cdf+cdp*cosL^3, cd

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

      #     call tadd(time0,t_surfcd2)

      return clpo, clps, clpt, CDfwing, CDpwing, CDwing, CDover

end # surfcd2

"""
    surfcd(S, 
    b, bs, bo, 
    λt, λs, sweep, 
    co, cdf, cdp, 
    Reco, Reref, 
    aRexp, kSuns, fCDcen)

Computes wing or tail surface profile CD from pre-computed chord quantities and corrections.
Called by [`cdsum!`](@ref) if `icdfun` flag set to 0.

!!! compat "Future Changes" 
      This function may be renamed for clarity of use.

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
function surfcd(S,
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

      #      if(λt==0.25) 
      #      write(*,*)
      #      write(*,*) fCDcen, lsfac, lfac
      #      write(*,*) ηo, ηs-ηo, 1.0-ηs
      #      write(*,*) ηo, ηs-ηo, 1.0-ηs
      #      write(*,*) 2.0 * ηo,
      #     &          (1.0    +λs)*(ηs-ηo),
      #     &          (λs+λt)*(1.0 -ηs), 
      #     &          2.0 * ηo
      #     &         +(1.0    +λs)*(ηs-ηo)
      #     &         +(λs+λt)*(1.0 -ηs)
      #        write(*,*) 2.0 * ηo * fCDcen,
      #     &          (1.0    +λs)*(ηs-ηo)*lsfac,
      #     &          (λs+λt)*(1.0 -ηs)*lfac,
      #     &           2.0 * ηo * fCDcen
      #     &         +(1.0    +λs)*(ηs-ηo)*lsfac
      #     &         +(λs+λt)*(1.0 -ηs)*lfac
      #      write(*,*) Sh/0.3048^2 , CDsurf
      #      end

      CDover = 0.0

      #     call tadd(time0,t_surfcd)

      return CDsurf, CDover
end # surfcd

