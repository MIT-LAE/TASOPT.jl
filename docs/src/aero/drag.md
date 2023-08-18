# Drag calculations

The key drag contributions are assumed to come from the fuselage, wing and tail surfaces and the lift induced drag calcualted at the trefftz plane


## Axisymmetric fuselage drag 
The fuselage profile drag is determined by a quasi-axisymmetric coupled viscous-inviscid calculation. See "Simplified Viscous/Inviscid Calculation for Nearly-Axisymmetric Bodies" by M. Drela.

This method does not require reliance on any wetted area approximations or fineness-ratio correlations.

![ADfuse](../assets/ADfuse.png)

```@docs

aerodynamics.axisol!(xnose,xend,xblend1,xblend2, Amax, 
	anose, btail, iclose,
	Mach, nc, nldim,
      xl, zl, sl, dyl, ql)

aerodynamics.blsys(simi,lami,wake,direct, Mach, uinv,
                      hksep, x,b,rn,th,ds,ue,
                      h , h_th, h_ds,
                      hk, hk_th, hk_ds, hk_ue,
                      hc, hc_th, hc_ds, hc_ue,
                      hs, hs_th, hs_ds, hs_ue,
                      cf, cf_th, cf_ds, cf_ue,
                      di, di_th, di_ds, di_ue,
                      xm,bm,rnm,thm,dsm,uem, 
                      hm , hm_thm, hm_dsm,
                      hkm, hkm_thm, hkm_dsm, hkm_uem,
                      hcm, hcm_thm, hcm_dsm, hcm_uem,
                      hsm, hsm_thm, hsm_dsm, hsm_uem,
                      cfm, cfm_thm, cfm_dsm, cfm_uem,
                      dim, dim_thm, dim_dsm, dim_uem)

aerodynamics.blax(ndim, n,ite, xi, bi, rni, uinv, Reyn, Mach, fexcr)

aerodynamics.blvar(simi,lami,wake, Reyn,Mach, fexcr,
                      x, θ ,δs ,ue )

fusebl!(pari, parg, para, ip)
```

## Total drag calculation
```@docs
cdsum!(pari,parg,para,pare, icdfun)
```

## Wing and tail surfaces
```@docs
aerodynamics.airtable(fname)

aerodynamics.surfcd2(
      S,
      b, bs, bo,
      λt, λs, γt, γs,
      toco, tocs, toct,
      Mach, sweep, co,
      CL, CLhtail, fLo, fLt,
      Reco, aRexp, kSuns, fexcd,
      AMa, Acl, Atau, ARe, A,
      fduo, fdus, fdut)

aerodynamics.surfcd(S,
      b, bs, bo, λt, λs, sweep, co,
      cdf, cdp, Reco, Reref, aRexp, kSuns,
      fCDcen)
```

## Treftz plane drag calculation
```@docs
aerodynamics.trefftz1(nsurf, npout, npinn, npimg,
	Sref, bref,
	b,bs,bo,bop, zcent,
	po,gammat,gammas, fLo,ktip,
	Lspec,CLsurfsp,t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)
```

## Other utilites

```@docs
aerodynamics.cfturb
```
```@setup cfturb
include("../../../src/aero/cdsum.jl")

```
For example, the turbulent flat plate ``C_f`` for a ``Re`` of ``10e6`` can be calculated as follows:

```@example cfturb
Re = 10e6
cfturb(Re)
```
