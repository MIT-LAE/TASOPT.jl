# TAESOPT.jl Documentation

`TAESOPT.jl` is a multi-disciplinary aircraft design and optimization code implemented in Julia based on [`TASOPT`](http://web.mit.edu/drela/Public/web/tasopt/) by Mark Drela.

It can currently model tube-and-wing aircraft using 2D viscous-invisicd CFD to calculate aerodynamic performance, simple beam bending theory to size the wings, and thermodynamic cycle modeling to calculate engine performance.


```@docs
wsize(pari, parg, parm, para, pare,
            itermax, wrlx1, wrlx2, wrlx3,
            initwgt, initeng, iairf, Ldebug, printiter, saveODperf)
```