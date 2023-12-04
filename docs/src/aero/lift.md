# Lift

Surface lift distributions $\tilde{p}$ are defined in terms of a
baseline piecewise-linear distribution $p{\scriptstyle (\eta)}$ defined
like the chord planform, but with its own taper ratios $\gamma_s$ and
$\gamma_t$. 

![](../assets/pload.png)
Piecewise-linear aerodynamic load ``\tilde{p} {\scriptstyle (\eta)}``,
with modifications at center and tip.

The segments are integrated to compute the lift contribution to the wing root load, ``p_o``, as required by structural computations. The sectional lift distribution is treated as an input to induced drag and surface drag computations.

```@eval
using Markdown
Markdown.parse_file(joinpath("../..", "src/aero","theory_lifting.md"))
```

```@docs
aerodynamics.wingpo(b, bs, bo,
               λt, λs, γt, γs,
               AR, N, W, Lhtail, fLo, fLt)

aerodynamics.wingcl(b,bs,bo,
            λt,λs,γt,γs,
            sweep,AR,CL,CLhtail,fLo,fLt,
            duo,dus,dut)
```
