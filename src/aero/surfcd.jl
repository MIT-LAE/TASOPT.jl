"""
    surfcd2(S,b, bs, bo,λt, λs, γt, γs,
    toco, tocs, toct,Mach, sweep, co,
    CL, CLhtail, fLo, fLt,Reco, 
    aRexp, kSuns, fexcd,fduo, fdus, fdut)

Calculates wing or tail surface profile `CD`.

# Inputs
- `S::Float64`: reference area.
- `b::Float64`: span.
- `bs::Float64`: outer panel break span.
- `bo::Float64`: wing-root (fuselage) span.
- `λt::Float64`: outer-panel chord taper ratio  ct/co.
- `λs::Float64`: inner-panel chord taper ratio  cs/co.
- `γt::Float64`: outer-panel load  taper ratio  pt/po.
- `γs::Float64`: inner-panel load  taper ratio  ps/po.
- `toco::Float64`: root  airfoil t/c.
- `tocs::Float64`: break airfoil t/c.
- `toct::Float64`: tip   airfoil t/c.
- `sweep::Float64`: wing sweep, degrees.
- `co::Float64`: wing root chord.
- `Reco::Float64`: Reynolds number for co.
- `aRexp::Float64`: Re-scaling exponent.
- `kSuns::Float64`: shock-unsweep area constant.
- `fexcd::Float64`: excrescence cd scaling factor.
- `fduo::Float64`: fractional overspeeds at wing stations.
- `fdus::Float64`:
- `fdut::Float64`:

# Outputs
- `clpo::Float64`: cl at wing stations.
- `clps::Float64`: 
- `clpt::Float64`:
- `CDfwing::Float64`: friction profile cd in perp. plane.
- `CDpwing::Float64`: pressure profile cd in perp. plane.
- `CDwing::Float64`: overall profile CD.
- `CDover::Float64`: fuselage added CD due to lift carryover.

See A. 2. 13 Wing Profile Drag of TASOPT docs.
See also [`surfcm`](@ref) and [`airtable`](@ref).
"""
function surfcd2(
      S,
      b, bs, bo,
      λt, λs, γt, γs,
      toco, tocs, toct,
      Mach, sweep, co,
      CL, CLhtail, fLo, fLt,
      Reco, aRexp, kSuns, fexcd,
      fduo, fdus, fdut)

      #     n = 4     #  ~0.30% error
      #     n = 6     #  ~0.13% error
      n::Int = 8     #  ~0.07% error
      #     n = 12    #  ~0.07% error

      #     call tset(time0)

      cosL = cos(sweep * π / 180.0)

      AR = b^2 / S

      ηo = bo / b
      ηs = bs / b

      #---- set center clpo
      Kc = ηo +
           0.5 * (1.0 + λs) * (ηs - ηo) +
           0.5 * (λs + λt) * (1.0 - ηs)

      Kp0 = ηo +
            0.5 * (1.0 + γs) * (ηs - ηo) +
            0.5 * (γs + γt) * (1.0 - ηs)

      Ko = 1.0 / (Kc * AR)
      Kp = Kp0 + fLo * ηo + 2.0 * fLt * Ko * γt * λt

      clp1 = (CL - CLhtail) / cosL^2 * S / (Kp * b * co)

      #---- set break and tip clp for passing back
      clpo = clp1 / (1.0 + fduo)^2
      clps = clp1 * γs / λs / (1.0 + fdus)^2
      clpt = clp1 * γt / λt / (1.0 + fdut)^2

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
            C = 1.0 - frac + λs * frac
            toc = toco * (1.0 - frac) + tocs * frac
            fdu = fduo * (1.0 - frac) + fdus * frac

            fSuns = exp(-(η - ηo) * b / (kSuns * C * 2.0 * co))

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
            C = λs * (1.0 - frac) + λt * frac
            toc = tocs * (1.0 - frac) + toct * frac
            fdu = fdus * (1.0 - frac) + fdut * frac

            fSuns = exp(-(η - ηo) * b / (kSuns * C * 2.0 * co))

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
    surfcd(S, b, bs, bo, λt, λs, 
    sweep, co, cdf, cdp, 
    Reco, Reref, aRexp, kSuns, fCDcen)

Calculates wing or tail surface profile CD

# Inputs
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

# Outputs
- `CDsurf`: overall profile CD.
- `CDover`: fuselage added CD due to lift carryover.
"""
function surfcd(S,
      b, bs, bo, λt, λs, sweep, co,
      cdf, cdp, Reco, Reref, aRexp, kSuns,
      fCDcen)


      cosL = cos(sweep * π / 180.0)

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

