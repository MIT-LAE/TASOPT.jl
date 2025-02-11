# [Notes for devs](@id fordevs)

!!! info
    `TASOPT.jl` is very much a WIP. Get at us on github, but search the issues first.  
    Thanks! 🙂

!!! tip "Tips"
    - Refer to the [data structures](@ref datastructs) to see where input file parameters end up.
    - Look out for `!!! compat` admonishments marking where things will likely change in the future.
    - References to NPSS are currently non-functional. We're working on replacing this functionality efficiently.

!!! tip "Benchmarks"
    - Examples of aircraft-sizing benchmarking and profiling files are provided in `test/benchmark_sizing.jl`. These can be run after making changes to the code to check if there has been a speed-up or slow-down.
    - Some individual functions are additionally benchmarked in `test/benchmark_elements.jl`.
    - Developers can similarly create other benchmarks for new features, models, and functions.


## [Performance tips](@id perftips)

Below are a few performance tips based on lessons learnt during TASOPT.jl development. The Julia docs has a section on performance that you can find [here](https://docs.julialang.org/en/v1/manual/performance-tips/), so the goal is not to repeat everything but expand on sections that are rather terse in there. 

### Custom types and type inference 

While defining new types you need to think about type inference and how the compiler can or cannot learn the types of the downstream calculations. See this section [here in the Julia manual that has some examples](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-fields-with-abstract-type). I'll list a more TASOPT.jl relevant example here to emphasize the point. 

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

We can do better by declaring the struct in such a way that the type of `cl` is inferred from the type of the wrapper object. Like,
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