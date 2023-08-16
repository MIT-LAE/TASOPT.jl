# TAESOPT.jl Documentation

`TAESOPT.jl` is a multi-disciplinary aircraft design and optimization code implemented in Julia based on [`TASOPT`](http://web.mit.edu/drela/Public/web/tasopt/) by Mark Drela.

It can currently model tube-and-wing aircraft using 2D viscous-invisicd CFD to calculate aerodynamic performance, simple beam bending theory to size the wings, and thermodynamic cycle modeling to calculate engine performance.


```@setup bench
include("../../src/aero/trefftz.jl")
cd("../../test")
#using BenchmarkTools, Pkg
#println("Current location $(pwd())")
#Pkg.activate("../")
#include("../tasopt.jl")
#include("input.jl")
#include("benchmark.jl")

```

```@example bench
benchmark_drag()
```