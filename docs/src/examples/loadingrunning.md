# [Loading and running a model] (@id firstexample)

`TASOPT.jl` uses TOML files to define aircraft models. You can find example input files at `/example/defaults/`. The majority of aircraft parameters and assumptions are defined here, and it's a useful resource for understanding the parameters and typical values.

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
Alternatively, you can load your desired input file (perhaps a modified version of a default file) by using:
```julia-repl
julia> example_ac = read_aircraft_model("../input.toml") # MODIFY <path> appropriately
```

`example_ac` is an instance of an `aircraft` `struct` (what Julia calls composite types); it's a thin wrapper for other `structs` representing aircraft components, and additional "`par`" arrays that store design and performance quantities. Refer to the sections on data structures for an [introduction](@ref datastructs_basics) and [details](@ref datastructs).

You can size this aircraft by running:
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

Once sized, an aircraft's performance can be determined for other missions by modifying or appending to `Mission` fields within an `input.toml` and running [`fly_mission!()`](@ref).

```julia-repl
julia> include(__TASOPTindices__)    #provides vars to access par array parameters (here, imPFEI)

julia> example_ac.parm[imPFEI,:]     #after sizing, the design mission (the first) is evaluated and the mission fuel weight is known
2-element Vector{Float64}:
 0.9443825885056822
 0.0

julia> fly_mission!(example_ac, 2)   #evaluate the second mission

julia> example_ac.parm[imPFEI,:]     #the second mission's fuel weight is known
2-element Vector{Float64}:
 0.9443825885056822
 1.080316809817944
```
