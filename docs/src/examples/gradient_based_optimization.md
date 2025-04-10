# Example for a Gradient Based Optimization

This example shows how to perform a gradient based optimization on a TASOPT model. It uses the IPOPT optimization module using a central, relative finite difference sensitivity module.

**_NOTE:_**  We are currently working on a autodiff based sensitivity module.

For this example, the design variables are:

1. Aspect Ratio: `AR`
2. Lift Coefficient: `Cl`  
3. Wing Sweep: `Λ`
4. Tt4 at Cruise: `Tt4`
5. Fan Pressure ratio: `pif`

## Initialiation and loading models

Start the script importing the necessary modules.

```julia
# Import modules
using TASOPT
using JuMP
using Ipopt
using Test
include(__TASOPTindices__)
```

Set the relative tolerance used in the Finite difference method

```julia
epsilon = 1e-5
```

Load the default `aircraft` model and size it to get initial values:

```julia
default_model = load_default_model()
size_aircraft!(default_model)
```

## Setting Optimization Parameters

### Set the input parameters

```julia
input_params = [
    :(ac.parg[igAR]), 
    :(ac.para[iaCL,ipcruise1:ipcruise2,1]),
    :(ac.parg[igsweep]),
    :(ac.pare[ieTt4, ipcruise1:ipcruise2, 1]),
    :(ac.pare[iepif, ipcruise1, 1]) ,
]

# Formatting params for usage 
params = map(p -> TASOPT.format_params(TASOPT.expr_to_string(p)), input_params)
```

### Set the Upper and Lower limits for all design variables as well as initial values

```julia
lower      = [9.0 , 0.53, 25.0, 1400.0, 1.25]
upper      = [11.0, 0.60, 30.0, 1650.0, 2.0 ] 
initial    = [10.5, 0.57, 26.0, 1580.0, 1.685]

```

## Objective and Constraint Functions

In this optimization we will have a constraint on max fuel weight, max span, min climb gradient, and max Tt3. Set functions that return the objective value and constraints of your optimization

```julia
function pfei_fn(ac)
    return ac.parm[imPFEI]
end

function wfuel_fn(ac)
    return ac.parg[igWfuel]/ac.parg[igWfmax]
end

function span_fn(ac)
    return ac.parg[igb]/ac.parg[igbmax]
end

function gtoc_fn(ac)
    return ac.para[iagamV, ipclimbn,1]/ac.parg[iggtocmin]
end

function tt3_fn(ac)
    maxtt3 = 900
    return maximum(ac.pare[ieTt3, :, 1])/maxtt3
end

# Make an array that stores all these functions
con_f_arr = [
    pfei_fn, wfuel_fn, span_fn, tt3_fn, gtoc_fn
]
```

## Objective Function

```julia
# Function that sizes aircraft and returns both objective value and constraint values
function sizing_ac(x::T...) where {T<:Real}
    ac = deepcopy(default_model)
    # Set params
    for (i,x_i) in enumerate(x)
        field_path, index = params[i]
        TASOPT.setNestedProp!(ac, field_path, x_i, index)
    end

    try
        size_aircraft!(ac,printiter=false)
        return [con_f_arr[i](ac) for i in 1:length(con_f_arr)]
    catch
        println("size_aircraft! FAILED")
        return [Inf for i in 1:length(con_f_arr)]
    end
    
end
```

## Gradient Function

```julia
function gradient_all(x::T...) where {T<:Real}
    ac = deepcopy(default_model)
    # Set params
    for (i,x_i) in enumerate(x)
        field_path, index = params[i]
        TASOPT.setNestedProp!(ac, field_path, x_i, index)
    end
    return TASOPT.get_sensitivity(input_params; model_state=ac, eps=epsilon, optimizer=true, f_out_fn=con_f_arr)
end
```

## Memoization

In order to reduce the function calls to `size_aircraft` we shall be using memoization to chache objective and constraint values for a specific point in the design space

```julia
function memoize(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    function foo_i(i, x::T...) where {T<:Real}
        if x !== last_x
            last_x, last_f = x, foo(x...)
        end
        return last_f[i]::T
    end
    return [(x...) -> foo_i(i, x...) for i in 1:n_outputs]
end

function memoize_sensitivity(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    function foo_s!(i, g::AbstractVector{T}, x::T...) where {T<:Real}
        if x !== last_x
            last_x, last_f = x, foo(x...)
        end
        if (size(g)[1] >0)
            for (k,grads) in enumerate(g)
                g[k] = last_f[i][k]
            end
        end
    end
    return [(x...) -> foo_s!(i, x...) for i in 1:n_outputs]
end

memoized_size_ac = memoize(sizing_ac, length(con_f_arr))
memoized_fd_ac = memoize_sensitivity(gradient_all, length(con_f_arr))
```

## Running the optimization

```julia
model = Model(Ipopt.Optimizer)
# set_silent(model)
set_optimizer_attribute(model, "tol", 1e-3)
set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
set_optimizer_attribute(model, "obj_scaling_factor", 1.0)
set_optimizer_attribute(model, "nlp_scaling_max_gradient", Float64(1))
set_optimizer_attribute(model, "nlp_scaling_min_value", 1e-8)

@variable(model, lower[i] <= x[i=1:length(initial)] <= upper[i], start=initial[i])

@operator(model, f_size_ac, length(initial), memoized_size_ac[1],memoized_fd_ac[1])
@operator(model, f_con_wfuel, length(initial), memoized_size_ac[2],memoized_fd_ac[2])
@operator(model, f_con_b, length(initial), memoized_size_ac[3],memoized_fd_ac[3])
@operator(model, f_con_tt3, length(initial), memoized_size_ac[4],memoized_fd_ac[4])
@operator(model, f_con_gtoc, length(initial), memoized_size_ac[5],memoized_fd_ac[5])

@objective(model, Min, f_size_ac(x...))
@constraint(model, f_con_wfuel(x...) <= 1)
@constraint(model, f_con_b(x...) <= 1)
@constraint(model, f_con_tt3(x...) <= 1)
@constraint(model, f_con_gtoc(x...) >= 1)

optimize!(model)

println("Finished optimization with $(value.(x))")

```
