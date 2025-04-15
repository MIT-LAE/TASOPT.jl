# [Data structure basics](@id datastructs_basics) 

Performance and design data is held largely in the `par` arrays (a holdover from FORTRAN TASOPT) along with a growing body of `struct`s to represent cohesive components and systems. An `aircraft` `struct` wraps these arrays, `struct`s, and auxiliary information.

!!! compat "Future Changes"
    We don't like this hybrid approach either. It's the legacy of Fortran.

    In a future major revision, we aim to completely replace the `par` array system with the `struct`-oriented approach.

## `par` arrays

Four arrays contain both prescribed inputs *and* computed outputs, some of which are multi-dimensional:

  - **`parg`**`::AbstractVector{Float64}`:  **geometric** quantities, weights, structural values, and other values inherent to an instantiated design (not operation).
  - **`parm`**`::AbstractArray{Float64, 2}`:  **mission**-prescribing parameters. The second dimension allows the specification of multiple mission profiles.
  - **`para`**`::AbstractArray{Float64, 3}`:  **aerodynamic** performance quantities. The second dimension captures the variation over a mission. The third dimension allows the specification of multiple mission profiles.
  - **`pare`**`::AbstractArray{Float64, 3}`:  **engine** perfomance quantities. As for `para`, the second and third dimensions capture flight-point and mission dependencies, respectively.

Data in the `par` arrays are accessed via Integer indices defined at `src/data_structs/index.inc`. These indices can be added to a namespace via `include(__TASOPTindices__)`:

```julia
using TASOPT
#using __TASOPTroot__, which fetches the src directory
include(joinpath(__TASOPTroot__, "data_structs/index.inc"))

#or more concisely
include(__TASOPTindices__)
```

The variable names of these indices indicate which `par` array they should access and hint at the quantity in question. For example, `ieTfuel` evaluates to `2` and retrieves the model's fuel temperature via `pare[ieTfuel]`. 

Note that for the multi-dimensional `par` arrays, indexing with a single Integer only retrieves the value for the first flight point of the first mission (namely, the design mission). Additional indexing is required to access data from different flight points or missions. Indices for specific flight points are defined in `index.inc` and should be used when indexing `pare` or `para`, e.g., `ipstatic` for static ground condition or `ipcruise1` for the start of cruise.


```@example dataaccess
using TASOPT
include(__TASOPTindices__)
ac = load_default_model()

println("Single element: ", size(ac.pare[ieTfuel]))
println(ac.pare[ieTfuel])

println("Full slices: ", size(ac.pare[ieTfuel,:,:]))
println(ac.pare[ieTfuel,:,:])

println("All missions at cruise start: ", size(ac.pare[ieTfuel,ipcruise1,:]))
println(ac.pare[ieTfuel,ipcruise1,:])

```



can be included via the convenience variable `__TASOPTindices__`






## `aircraft` `struct`

An `aircraft` is composed of `par` array fields, title and description fields, and a `is_sized` flag to indicate its status. An optional `fuse_tank` field is present as a trial for future `struct`-based development. All fields are dot-accessible and array elements can be changed (e.g., `ac.parg[igS] = 20`), though the `struct` itself is not mutable.

Refer to the [`struct` reference page](@ref datastructs) for add'l details.
