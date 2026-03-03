# Atmospheric properties

```@docs
TASOPT.atmos
```
This function can be used to return the atmospheric properties at a given altitude as follows:
```@setup atmos
using TASOPT.atmosphere
```
```@example atmos
h = 10_000.0 # m
atmos_state = atmos(h)
(atmos_state.T, atmos_state.p, atmos_state.ρ, atmos_state.a, atmos_state.μ)
```
