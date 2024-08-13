# Pitching moment
The pitching moments of lifting surfaces are computed by integration of the wing loading with reference to a prescribable wing axis.

```@eval
using Markdown
Markdown.parse_file(joinpath("../..", "src/aero","theory_pitching.md"))
```

```@docs
aerodynamics.surfcm(wing, γt,γs,cmpo,cmps,cmpt)
```