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