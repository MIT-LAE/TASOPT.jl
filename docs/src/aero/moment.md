# Pitching moment
The pitching moments of lifting surfaces are computed by integration of the wing loading with reference to a prescribable wing axis.

```@eval
using Markdown
Markdown.parse_file(joinpath("../..", "src/aero","theory_pitching.md"))
```

```@docs
aerodynamics.wing_CM(b,bs,bo, sweep, Xaxis,
                       λt, λs, γt, γs,
                       AR, fLo, fLt, cmpo, cmps, cmpt)
```