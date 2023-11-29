# TASOPT.jl Documentation

`TASOPT.jl` is a multi-disciplinary aircraft design and optimization code implemented in Julia based on [`TASOPT`](http://web.mit.edu/drela/Public/web/tasopt/) by Mark Drela.

It can currently model tube-and-wing aircraft using 2D viscous-inviscid CFD to calculate aerodynamic performance, simple beam bending theory to size the wings, and thermodynamic cycle modeling to calculate engine performance.

## Getting started

There are several workflows that are possible to use `TASOPT.jl`. We outline here the most common few.

### Simple install

The easiest way to run `TASOPT.jl` would be to add the package using the julia package manager using the github repository.

You can do this by starting a Julia session and then activating the package manager by typing `]` and then entering:
```julia-repl
pkg> add "git@github.mit.edu:LAE/TAESOPT.jl.git"
```

You can then import `TASOPT` as you would with any Julia package:
```julia-repl
julia> using TASOPT
```
### Local development

If you are going to develop the source code of `TASOPT.jl` you might benefit from a local clone of the git repository which
can then fit into a workflow using [`Revise.jl`](https://timholy.github.io/Revise.jl/stable/) for example.

Step 1: Clone the git repo locally
```bash
git clone git@github.mit.edu:LAE/TAESOPT.jl.git
```

Step 2: `cd` to the folder where TASOPT is cloned

Step 3: Use `Pkg` to install/ develop the package

```julia
pkg> dev .
```

You should now be able to import TASOPT from within any Julia script in your base environment.

!!! tip "Tip"
    If you are using `Revise.jl` be sure to first import `Revise` before importing `TASOPT`

    ```julia
    using Revise
    using TASOPT
    ```

## Loading and running a TASOPT aircraft model

`TASOPT.jl` uses TOML files to define aircraft models. You can find an example input file here `/src/IO/default_input.toml`.

The following section gives an example of using `TASOPT.jl`:

Start by importing `TASOPT.jl` and then loading the default aircraft model.
```julia-repl
julia> using TASOPT
julia> example_ac = load_default_model()
Loading default aircraft model
Name: Default TASOPT Model;
Wpay = 172.0 kN
Des. Range  = 5.56e6 km
Cruise Mach = 0.8
```
Alternatively you can load your desired input file by using
```julia-repl
julia> example_ac = read_aircraft_model("../src/IO/input.toml") # MODIFY <path> appropriately
```

`example_ac` is an instance of an `aircraft` type, that is a thin wrapper for 
a couple of arrays that store, for example, the geometric `parg`,
 aerodynamic (`para`), engine (`pare`).

You can size this aircraft by running
```julia-repl
julia> size_aircraft!(example_ac)
Max payload weight was not set, setting Wpaymax = Wpay
Wfuel initial = 132502.37055588452
iterw             errW            errW1            WMTO           Wfuel          Wftank          Wtesys            Wmot            Wgen         Wtshaft           Wwing            span            area          HTarea           xwbox 
    1  +1.00000000e+00  +1.00000000e+00  6.51600802e+05  1.32502371e+05  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  8.60730957e+04  0.00000000e+00  0.00000000e+00  1.35396163e+01  1.73736000e+01
    2  +1.00000000e+00  +1.66001430e-02  7.70269325e+05  2.16922911e+05  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.10711244e+05  3.25978109e+01  1.05209631e+02  3.71505481e+01  1.73736000e+01
 ...
   15  +7.45835713e-09  -7.45835713e-09  7.76336715e+05  2.12378504e+05  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  1.04209402e+05  3.54371733e+01  1.24335966e+02  4.20539315e+01  1.62730055e+01

Takeoff:
 #        lTO         l1        lBF       dmax
 1   6474.360   5179.488   8416.667    355.380
 2   6474.360   5534.868   8485.441      3.718
 3   6474.360   5538.586   8485.689      0.000
```

```@example
using Dates # hide
println("Documentation built $(Dates.now()) with Julia $(VERSION)") # hide
```