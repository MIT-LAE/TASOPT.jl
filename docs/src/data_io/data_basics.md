# Data structure basics

Data is held solely in the `par` arrays (a holdover from FORTRAN TASOPT). A shallow `aircraft` struct wraps these arrays along with auxiliary information.

## `par` arrays

All data, both prescribed inputs *and* computed outputs,  are stored in 5 arrays, some of which are multi-dimensional:

  - **`pari`**`::AbstractVector{Int64}`:  **integer** flags indicating categorical design choices and analysis types
  - **`parg`**`::AbstractVector{Float64}`:  **geometric** quantities, weights, structural values, and other values inherent to an instantiated design (not operation).
  - **`parm`**`::AbstractArray{Float64, 2}`:  **mission**-prescribing parameters. The second dimension allows the specification of multiple mission profiles.
  - **`para`**`::AbstractArray{Float64, 3}`:  **aerodynamic** performance quantities. The second dimension captures the variation over a mission. The third dimension allows the specification of multiple mission profiles.
  - **`pare`**`::AbstractArray{Float64, 3}`:  **engine** perfomance quantities. As for `para`, the second and third dimensions capture flight-point and mission dependencies, respectively.


## `aircraft` struct

An aircraft is composed of `par` array fields, title and description fields, and a `sized` flag to indicate its status. An optional `fuse_tank` field is present as a trial for future struct-based development. All fields are dot-accessible and array elements can be changed (e.g., `ac.parg[igS] = 20`), though the struct itself is not mutable.

Refer to the [struct reference page](@ref datastructs) for add'l details.


!!! compat "Future Changes"
    We don't like this either.

    In a future major revision, we aim to replace the `par` array system with a struct-oriented approach.
