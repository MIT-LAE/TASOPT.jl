# Drag

The key drag contributions are assumed to come from the fuselage, wing and tail surfaces, and the lift-induced drag calculated at the Trefftz plane. Wave drag is not explicitly modelled.


## [Axisymmetric fuselage drag](@id axi)
The fuselage profile drag is determined by a quasi-axisymmetric coupled viscous-inviscid calculation. See [Simplified Viscous/Inviscid Analysis for Nearly-Axisymmetric Bodies](../assets/drela_TASOPT_2p16/axibl.pdf) by M. Drela.

This method does not require any wetted area approximations or fineness-ratio correlations, but does require the geometry to be specified in the form of a
cross-sectional area distribution $A{\scriptstyle (x)}$ and a
perimeter distribution $b_0{\scriptstyle (x)}$, shown in the
Figure below. For a round cross-section these are, of course, related. To allow treating
more general fuselage cross-sections, they are assumed to be specified
separately. The cross section sizes and shapes can vary along the body,
provided the variation is reasonably smooth.

![ADfuse](../assets/ADfuse.png)

```@eval
using Markdown
Markdown.parse_file(joinpath("../..", "src/aero","theory_fuse_profile_drag.md"))
```

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

aerodynamics.fusebl!(fuse, parm, para, ip)
```

---

## [Trefftz plane drag calculation](@id trefftz)

Trefftz plane analysis computes the induced drag from lifting surfaces. The lift distributions are propagated downstream, accounting for streamline contraction from fuselage thickness variation as shown in the Figure below. 

![](../assets/trefftz.png)
Two shaded streamtubes are shown. Wake center radius $y'_o$ is nonzero due to the fuselage viscous wake displacement area.

The vorticity in the wake is numerially integrated at collocation points to determine the overall induced drag.

![T](../assets/tpvort.png)

Trefftz Plane vortices $i,i\!+\!1 \ldots$ and collocation points
$i\!+\!1/2$ used for velocity, impulse, and kinetic energy calculations.
Left/right symmetry is exploited.  

```@eval
using Markdown
Markdown.parse_file(joinpath("../..", "src/aero","theory_trefftz_plane.md"))
```

```@docs
aerodynamics.cditrp(para, wing, htail)

aerodynamics.trefftz1(nsurf, npout, npinn, npimg,
	Sref, bref,
	b,bs,bo,bop, zcent,
	po,gammat,gammas, fLo,ktip,
	Lspec,CLsurfsp,t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)
```
---

## Wing and tail surfaces

Lifting surface drag is determined via `surfcd` (when constant airfoil section `cdf` and `cdp` are already determined), and `surfcd2` (when an explicit modelling and integration is desired). Airfoil performance is accessed via a lookup of precomputed airfoil data, `airfun`.

```@docs
aerodynamics.surfcd2(wing, γt, γs,
            Mach, CL, CLhtail, 
            Reco, aRexp, kSuns, fexcd,
            fduo, fdus, fdut)

aerodynamics.surfcd(S,
      b, bs, bo, λt, λs, sweep, co,
      cdf, cdp, Reco, Reref, aRexp, kSuns,
      fCDcen)

aerodynamics.airtable(fname)

aerodynamics.airfun(cl, τ, Mach, air::aerodynamics.airfoil)

```

---

## Total drag calculation
```@docs
aerodynamics.cdsum!(parg, para, pare, wing, htail, vtail, icdfun)
```
---

## Other utilities

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
