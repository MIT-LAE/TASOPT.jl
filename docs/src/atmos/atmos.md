# Atmospheric properties

```@docs
TASOPT.atmos
```
This function can be used to return the atmospheric properties at a given altitude as follows:
```@setup atmos
include("../../../src/atmos/atmos.jl")
using .atmosphere
```
```@example atmos
h = 10.0 # km
T,p,ρ,a,μ = atmos(h)
T,p,ρ,a,μ
```