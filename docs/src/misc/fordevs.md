# [Notes for devs](@id fordevs)

!!! info
    `TASOPT.jl` is very much a WIP. Get at us on github with ideas, contributions, and bugs, but search the issues first.  
    Thanks! 🙂

!!! tip "Tips"
    - Refer to the [data structures](@ref datastructs) to see where input file parameters end up.
    - Look out for `!!! compat` admonishments marking where things will likely change in the future.
    - Any remaining references to NPSS are currently non-functional. We're working on replacing detailed (but efficient) engine performance modelling.

## Benchmarks
- Examples of aircraft-sizing benchmarking and profiling files are provided in `test/benchmark_sizing.jl`. These can be run after making changes to the code to check if there has been a speed-up or slow-down.
- Some individual functions are additionally benchmarked in `test/benchmark_elements.jl`.
- Developers can similarly create other benchmarks for new features, models, and functions.


## [Performance tips](@id perftips)

Below are a few performance tips based on lessons learnt during TASOPT.jl development. The Julia docs has a section on performance that you can find [here](https://docs.julialang.org/en/v1/manual/performance-tips/), so the goal is not to repeat everything but expand on sections that are rather terse in there. Two key things to keep in mind to get fast Julia code are:
1. _Writing type stable code_ for which the compiler can generate performant code.
2. _Reducing unnecessary allocations_.


## Reducing allocations and profiling

The `test/benchmark_elements.jl` file shows some examples of using `BenchmarkTools.jl` to benchmark functions in Julia. 

Sometimes you need more than just the number of allocations from benchmarking to actually go eliminate or reduce the allocations. _Where_ the allocations are being made can be non-obvious in some cases. [`Coverage.jl`](https://github.com/JuliaCI/Coverage.jl) is a useful package for getting a better sense of where in your code the allocations are coming from.

Specifically you want to follow the steps [here](https://github.com/JuliaCI/Coverage.jl?tab=readme-ov-file#memory-allocation). Reproduced here for convenience:

!!! details "How to track down allocations"
    Start julia with tracking allocations enabled:
    ```bash
    julia --track-allocation=user
    ```
    Then run all the commands you want to profile (this is to ensure they compile first), then clear the memory allocation tracking by running `Profile.clear_malloc_data()`; run your commands again and then quit julia. For example:

    ```julia-repl
    using TASOPT, Profile
    julia> Re = 10e6
    1.0e7

    julia> TASOPT.aerodynamics.cfturb(Re)
    0.002954557862895432

    julia> Profile.clear_malloc_data()

    julia> TASOPT.aerodynamics.cfturb(Re)
    0.002954557862895432

    julia> exit()
    ```

    Then look at the directory where these files live (i.e., the source code) and you should see some additional files with annotations showing where the allocations were made. 

    You can do all the above steps without needing `Coverage.jl`. Where `Coverage.jl` becomes useful is to analyze large parts of the code by doing:

    ```julia-repl
    using Coverage
    analyze_malloc(dirnames)  # could be "." for the current directory, or "src", etc.
    ```

## Custom types and type inference 

While defining new types (i.e., `structs`) you need to think about type inference and how the compiler can or cannot learn the types of the downstream calculations. See this section [here in the Julia manual that has some examples](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-fields-with-abstract-type). I'll list a more TASOPT.jl relevant example here to emphasize the point. 

Let's take the `airfoil` example. Consider an airfoil database as follows:
```julia
struct airfoil
	Re::AbstractFloat
	Ma::AbstractVector{Float64}
	cl::AbstractVector{Float64}
	thickness::AbstractVector{Float64}
end

# Create an instance
air_unstable = airfoil(10e6, Ma_array, cl_array, toc_array)
```

For the above structure the Julia compiler will *not* be able to generate high performance code. This is fundamentally because the type of `air.Re` cannot be determined by the type of `a`. For example the compiler can't know from the type of `a` if `cl_array` was a `Vector{Float64}` or not and won't be able to create type stable code.

```julia-repl
julia> typeof(a.cl), typeof(a.cl) <: AbstractVector{Float64}, typeof(a.cl) <: Vector{Float64}
(LinRange{Float64, Int64}, true, false)
```

We can do better by declaring the `struct` in such a way that the type of `cl` is inferred from the type of the wrapper object. Like,
```julia
struct airfoil{T<:AbstractFloat, V<:AbstractVector{Float64}}
	Re::T
	Ma::V
	cl::V
	thickness::V
end

# Create an instance
air_stable = airfoil(10e6, Ma_array, cl_array, toc_array)
```

Now if `Ma_array`, `cl_array`, and `toc_array` were all of type `Vector{Float64}` then we'd end up with something like this:
```julia-repl
julia> typeof(a_stable)
airfoil{Float64, Vector{Float64}}
```
 
 But if they were  `<:AbstractRange` you'd get:
 ```
julia> typeof(a_stable)
airfoil{Float64, LinRange{Float64, Int64}}
```

In this case given the type, **not the value**, of `a` the compiler can correctly infer 
the type of `a.cl`, and generate appropriate code. 

## Type instability from using variables in parent scopes

Here let's look at a common pattern that can be tempting to write and what the performance penalty is.

```julia-repl
ia = 3
function do_something1(A)
    a = A[ia]
    if a == 0.0
        return 1.0
    else
        return 2.0
    end
end
```

Consider the code above - it takes in some Array, loads a particular value and then compares it and does something to it. 

```julia-repl
julia> A = rand(3,3)
3×3 Matrix{Float64}:
 0.355865  0.81659    0.529105
 0.210353  0.978487   0.671198
 0.734191  0.0497119  0.72487
 
julia> @code_warntype do_something1(A)
MethodInstance for do_something1(::Matrix{Float64})
  from do_something1(A) @ Main REPL[6]:1
Arguments
  #self#::Core.Const(Main.do_something1)
  A::Matrix{Float64}
Locals
  a::Any
Body::Float64
1 ─      (a = Base.getindex(A, Main.ia))
│   %2 = Main.:(==)::Core.Const(==)
│   %3 = a::Any
│   %4 = (%2)(%3, 0.0)::Any
└──      goto #3 if not %4
2 ─      return 1.0
3 ─      return 2.0
```

The compiler can't determine the type of a from the arguments alone (**even though** it correctly identifies that the input argument type is `Matrix{Float64}`). This is because the type of `ia` is not specified.

You can fix that by doing something like this:
```julia-repl
const ia2 = 3
function do_something2(A)
    a = A[ia2]
    if a == 0.0
        return 1.0
    else
        return 2.0
    end
end

ia3::Int = 3
function do_something3(A)
    a = A[ia3]
    if a == 0.0
        return 1.0
    else
        return 2.0
    end
end
```

Both of the above now return type stable code:

```julia-repl
julia> @code_warntype do_something2(A)
MethodInstance for do_something2(::Matrix{Float64})
  from do_something2(A) @ Main REPL[3]:1
Arguments
  #self#::Core.Const(Main.do_something2)
  A::Matrix{Float64}
Locals
  a::Float64
Body::Float64
1 ─      (a = Base.getindex(A, Main.ia2))
│   %2 = Main.:(==)::Core.Const(==)
│   %3 = a::Float64
│   %4 = (%2)(%3, 0.0)::Bool
└──      goto #3 if not %4
2 ─      return 1.0
3 ─      return 2.0

julia> @code_warntype do_something3(A)
MethodInstance for do_something3(::Matrix{Float64})
  from do_something3(A) @ Main REPL[11]:1
Arguments
  #self#::Core.Const(Main.do_something3)
  A::Matrix{Float64}
Locals
  a::Float64
Body::Float64
1 ─      (a = Base.getindex(A, Main.ia3))
│   %2 = a::Float64
│   %3 = (%2 == 0.0)::Bool
└──      goto #3 if not %3
2 ─      return 1.0
3 ─      return 2.0
```

The relevant sections in the performance docs are [here](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-untyped-global-variables).

## Static arrays and performance

[`StaticArrays`](https://github.com/JuliaArrays/StaticArrays.jl) is a package that provides functionality for *statically sized* (i.e., the size is determined from the *type*, doesn't **have** to be immutable) arrays.

Consider the `Weight` types defined in TASOPT.

```julia
@kwdef struct Weight <: AbstractLoad
	"""Weight [N]"""
	W::Float64
	"""Location {x,y,z} [m]"""
	r::SVector{3, Float64} = SA[0.0,0.0,0.0]
	"""Coordinate Frame"""
	frame::Frame = WORLD

end
```

Then extending `Base.+` and adding a new function to do center of mass calculations:
```julia
import Base.+

function +(W1::T, W2::T) where T<:Weight
	total_W = W1.W + W2.W
	Weight(total_W, (W1.r*W1.W + W2.r*W2.W)/total_W)
end # function +

"""
    center_of_weight(W_array::AbstractArray{Weight})

Calculates the coordinates of the center of mass/weight and returns a `Weight`
type of the equivalent weight and at the center of mass.
"""
function center_of_weight(W_array::AbstractArray{Weight})
    total_weight = 0.0
    r̄ = SVector{3,Float64}(zeros(3))
    for weight in W_array
        total_weight += weight.W
        r̄ = r̄ + weight.W * weight.r
    end
    return Weight(W = total_weight, r = r̄./total_weight)
end

```

Now let's look at performance:

```julia
julia> Ws = [W1, W2, W3, W4] #Assume these are already defined
4-element Vector{Weight}:
 Weight(10.0, [0.0, 0.0, 0.0], Frame(0, [0.0, 0.0, 0.0]))
 Weight(10.0, [10.0, 0.0, 0.0], Frame(0, [0.0, 0.0, 0.0]))
 Weight(10.0, [10.0, 10.0, 0.0], Frame(0, [0.0, 0.0, 0.0]))
 Weight(10.0, [0.0, 10.0, 0.0], Frame(0, [0.0, 0.0, 0.0]))

julia> center_of_weight(Ws)
Weight(40.0, [5.0, 5.0, 0.0], Frame(0, [0.0, 0.0, 0.0]))

julia> sum(Ws)
Weight(40.0, [5.0, 5.0, 0.0], Frame(0, [0.0, 0.0, 0.0]))

julia> @benchmark center_of_weight($Ws)
BenchmarkTools.Trial: 10000 samples with 997 evaluations.
 Range (min … max):  19.475 ns … 571.464 ns  ┊ GC (min … max): 0.00% … 92.64%
 Time  (median):     20.186 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   21.935 ns ±  18.992 ns  ┊ GC (mean ± σ):  3.39% ±  3.80%

  ▃▃▇█▆▃▂▃▄▂▁▁ ▁▁▁▁▁▁ ▂▂▁▁▁▁▂▂▂▁ ▁                             ▂
  █████████████████████████████████▇█▇▇█▇▆▇▇▇▆▇▆▇▆█▆▇▇▆▅▅▅▆▆▁▅ █
  19.5 ns       Histogram: log(frequency) by time      31.5 ns <

 Memory estimate: 80 bytes, allocs estimate: 1.

julia> @benchmark sum($Ws)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  7.458 ns … 21.667 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     7.708 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   7.791 ns ±  0.836 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

    █▁▃                                                       
  ▂▄███▄▃▂▂▂▂▁▁▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▂▂▂▂▂▂▂▂ ▂
 7.46 ns        Histogram: frequency by time        11.5 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.
```

By extending `Base.+` we got the `sum` function for free cause it just knows how to add things together. But you see that the static array approach seems to take much longer, that's because the static array definition here isn't done correctly. This is an easy to make mistake, look at the following comparison: 

```julia
function center_of_weight(W_array::AbstractArray{Weight})
    total_weight = 0.0
    r̄ = SVector{3}([0.0, 0.0, 0.0])
    for weight in W_array
        total_weight += weight.W
        r̄ = r̄ + weight.W * weight.r
    end
    return Weight(W = total_weight, r = r̄./total_weight)
end
julia> @benchmark center_of_weight($Ws)
BenchmarkTools.Trial: 10000 samples with 998 evaluations.
 Range (min … max):  17.995 ns … 541.165 ns  ┊ GC (min … max): 0.00% … 94.78%
 Time  (median):     18.745 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   19.894 ns ±  19.000 ns  ┊ GC (mean ± σ):  3.84% ±  3.87%

  ▄▆▆▅▇██▇▄▂                    ▁▃▃▃▂▂▁▁                       ▂
  ███████████▇▆▆▅▆▅▃▅▁▅▅▅▅▅▇█▇██████████▇█▇▇▆▆▆▆▆▄▅▄▅▄▅▃▅▅▄▄▅▇ █
  18 ns         Histogram: log(frequency) by time      25.8 ns <

 Memory estimate: 80 bytes, allocs estimate: 1.
```

 VERSUS:
```julia
function center_of_weight(W_array::AbstractArray{Weight})
    total_weight = 0.0
    r̄ = SVector{3}(0.0, 0.0, 0.0)
    for weight in W_array
        total_weight += weight.W
        r̄ = r̄ + weight.W * weight.r
    end
    return Weight(W = total_weight, r = r̄./total_weight)
end


julia> @benchmark center_of_weight($Ws)
BenchmarkTools.Trial: 10000 samples with 1000 evaluations.
 Range (min … max):  3.917 ns … 17.292 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     4.000 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   4.046 ns ±  0.431 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

          ▃       █       ▆       ▁        ▃       ▁         ▁
  ▃▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁▁▁▁▁██▁▁▁▁▁▁▁█▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▆ █
  3.92 ns      Histogram: log(frequency) by time     4.21 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.

```

**All that really changed was the little square brackets!**
`SVector{3}(0.0, 0.0, 0.0)` vs `SVector{3}([0.0, 0.0, 0.0])`
The latter results in 1 allocation, which, for such a small calculation, is a significant increase in the time required!
