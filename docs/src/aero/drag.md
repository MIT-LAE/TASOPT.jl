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
aerodynamics._axisymm_flow(xnose,xend,xblend1,xblend2, Amax, 
	anose, btail, iclose,
	Mach, nc, nldim,
      xl, zl, sl, dyl, ql)

aerodynamics._BL_station_system(is_selfsimilar, is_laminar, is_wake, solves_direct, Mach, uinv,
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

aerodynamics._axisymm_BL(ndim, n,ite, xi, bi, rni, uinv, Reyn, Mach, fexcr)

aerodynamics._BL_station_vars(is_selfsimilar, is_laminar, is_wake, Reyn,Mach, fexcr,
                      x, θ ,δs ,ue )

aerodynamics.fuselage_drag!(fuse, parm, para, ip)
```

---

## [Trefftz plane drag calculation](@id trefftz)

Trefftz plane analysis computes the induced drag from lifting surfaces. The lift distributions are propagated downstream, accounting for streamline contraction from fuselage thickness variation as shown in the Figure below. 

![](../assets/trefftz.png)
Two shaded streamtubes are shown. Wake center radius $y'_o$ is nonzero due to the fuselage viscous wake displacement area.

The vorticity in the wake is numerically integrated at collocation points to determine the overall induced drag.

![T](../assets/tpvort.png)

Trefftz Plane vortices $i,i\!+\!1 \ldots$ and collocation points
$i\!+\!1/2$ used for velocity, impulse, and kinetic energy calculations.
Left/right symmetry is exploited.  

```@eval
using Markdown
Markdown.parse_file(joinpath("../..", "src/aero","theory_trefftz_plane.md"))
```

```@docs
aerodynamics.induced_drag!(para, wing, htail)

aerodynamics._trefftz_analysis(nsurf, npout, npinn, npimg,
	Sref, bref,
	b,bs,bo,bop, zcent,
	po,gammat,gammas, fLo,ktip,
	specifies_CL,CLsurfsp,t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)
```
---

## Wing and tail surfaces

Lifting surface drag is determined via [`wing_profiledrag_scaled`](@ref aerodynamics.wing_profiledrag_scaled) (when constant airfoil section `cdf` and `cdp` are already determined), and [`wing_profiledrag_direct`](@ref aerodynamics.wing_profiledrag_direct) (when an explicit modelling and integration is desired). Airfoil performance is accessed via a lookup of precomputed airfoil data, `airfun`.

```@docs
aerodynamics.wing_profiledrag_direct(wing, γt, γs,
            Mach, CL, CLhtail, 
            Reco, aRexp, kSuns, fexcd,
            fduo, fdus, fdut)

aerodynamics.wing_profiledrag_scaled(S,
      b, bs, bo, λt, λs, sweep, co,
      cdf, cdp, Reco, Reref, aRexp, kSuns,
      fCDcen)
```
### Airfoil section data

The wing airfoil performance is represented by a parameterized transonic
airfoil family spanning a range of thicknesses, whose performance is
determined by 2D viscous/inviscid CFD calculation for a range of lift
coefficients and Mach numbers. Together with suitable sweep corrections,
this gives reliable profile+wave drag of the wing in cruise and high
climb and high descent.

The pressure and friction drag coefficients for the wing are obtained from the 2D drag database via cubic interpolation. For the databases provided, each airfoil has been designed (by Mark Drela) for a well-defined transonic drag rise, so that the database returns $c_{d_f}$ and $c_{d_p}$ values
representative of the best transonic airfoil technology.

Specifically, the perpendicular-plane friction and
pressure drag coefficients are then obtained from the 
database having the form 

$$\begin{aligned}
c_{d_f} & = & 
f_{\rm w_{excr}} \:
\bar{c}_{d_f}(c_{\ell_{\scriptscriptstyle \perp}},M_{\scriptscriptstyle \perp}, {\textstyle \frac{t}{c}} ) 
             \left( \frac{R\!e_c}{R\!e_{\rm ref}} \right)^{\! a_{\scriptscriptstyle Re}} \\
c_{d_p} & = & 
f_{\rm w_{excr}} \:
\bar{c}_{d_p}(c_{\ell_{\scriptscriptstyle \perp}},M_{\scriptscriptstyle \perp}, {\textstyle \frac{t}{c}} ) 
             \left( \frac{R\!e_c}{R\!e_{\rm ref}} \right)^{\! a_{\scriptscriptstyle Re}} \\
\mathrm{where} 
\hspace{5ex}
M_{\scriptscriptstyle \perp}& = & M_{{\scriptscriptstyle \infty}}\, \cos \Lambda \\
\textstyle \frac{t}{c} & = & \bar{h} \\
R\!e_c& = & \frac{\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}\,c}{\mu_{\scriptscriptstyle \infty}} \\
a_{\scriptscriptstyle Re}& \simeq & -0.15
\end{aligned}$$ 

and $f_{\rm w_{excr}} \geq 1$ is an empirical specified
factor to account for wing excrescence drag sources, and
$R\!e_{\rm ref}$ is a reference Reynolds number at which the database
functions $\bar{c}_{d_f}, \bar{c}_{d_p}$ were computed. The chord
Reynolds number $R\!e_c$ could of course be treated as an additional
parameter in the database, but at a considerable increase in the size of
the database and the computational effort needed to construct it. The
value of the Re-scaling exponent $a_{\scriptscriptstyle Re}\simeq -0.15$
is appropriate for fully-turbulent flow. See the theory 📖 block below for more details.

```@eval
using Markdown
Markdown.parse_file(joinpath("../..", "src/aero","theory_airfun_and_splines.md"))
```
```@docs
aerodynamics.airtable(fname)

aerodynamics.airfun(cl, τ, Mach, air::aerodynamics.airfoil)

```

---

## Total drag calculation
```@docs
aerodynamics.aircraft_drag!(ac, imission, ip, computes_wing_direct)
```
---

## Other utilities

```@docs
aerodynamics.cfturb
```
```@setup cfturb
include("../../../src/aero/drag.jl")

```
For example, the turbulent flat plate ``C_f`` for a ``Re`` of ``10e6`` can be calculated as follows:

```@example cfturb
Re = 10e6
cfturb(Re)
```
