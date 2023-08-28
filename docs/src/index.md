# TASOPT.jl Documentation

`TASOPT.jl` is a multi-disciplinary aircraft design and optimization code implemented in Julia based on [`TASOPT`](http://web.mit.edu/drela/Public/web/tasopt/) by Mark Drela.

It can currently model tube-and-wing aircraft using 2D viscous-invisicd CFD to calculate aerodynamic performance, simple beam bending theory to size the wings, and thermodynamic cycle modeling to calculate engine performance.

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

!!! tip If you are using `Revise.jl` be sure to first import `Revise` before importing `TASOPT`

    ```julia
    using Revise
    using TASOPT
    ```