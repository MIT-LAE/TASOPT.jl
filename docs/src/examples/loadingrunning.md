# [Loading and running a model] (@id firstexample)

`TASOPT.jl` uses TOML files to define aircraft models. You can find an example input file at `/src/IO/default_input.toml`. The majority of aircraft parameters and assumptions are defined here, and it's a useful resource for understanding the parameters and typical values.

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
