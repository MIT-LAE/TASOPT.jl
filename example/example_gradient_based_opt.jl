using JuMP
using Ipopt
using Test
using TASOPT
include(__TASOPTindices__)

# Set relative tolerance for Finite difference method
epsilon = 1e-5

# Set design variables
input_params = [
    :(ac.parg[igAR]), 
    :(ac.para[iaCL,ipcruise1:ipcruise2,1]),
    :(ac.parg[igsweep]),
    :(ac.pare[ieTt4, ipcruise1:ipcruise2, 1]),
    :(ac.pare[iepif, ipcruise1, 1]) ,
]
params = map(p -> TASOPT.format_params(TASOPT.expr_to_string(p)), input_params)

# Set box constraints
lower      = [9.0 , 0.53, 25.0, 1400.0, 1.25]
upper      = [11.0, 0.60, 30.0, 1650.0, 2.0 ] 
initial    = [10.5, 0.57, 26.0, 1580.0, 1.685]

# Load default model
default_model = load_default_model() 
size_aircraft!(default_model)

# Objective and Constraint Functions
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

# Create an array that stores above functions
# This is used when Memoizing the results 
con_f_arr = [
    pfei_fn, wfuel_fn, span_fn, tt3_fn, gtoc_fn
]

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

# Function that uses TASOPT Sensitivity module to get Finite difference of both objective and constraints
function gradient_all(x::T...) where {T<:Real}
    ac = deepcopy(default_model)
    # Set params
    for (i,x_i) in enumerate(x)
        field_path, index = params[i]
        TASOPT.setNestedProp!(ac, field_path, x_i, index)
    end
    return TASOPT.get_sensitivity(input_params; model_state=ac, eps=epsilon, optimizer=true, f_out_fn=con_f_arr)
end


# Memoization Functions
# The below functions Take a function `foo` and return a vector of length `n_outputs`, where element
# `i` is a function that returns the equivalent of `foo(x...)[i]`.
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

# Set memoization functions
memoized_size_ac = memoize(sizing_ac, length(con_f_arr))
memoized_fd_ac = memoize_sensitivity(gradient_all, length(con_f_arr))

# Start optimization setup
model = Model(Ipopt.Optimizer)

# Scaling parameters for faster convergence
set_optimizer_attribute(model, "tol", 1e-3)
set_optimizer_attribute(model, "nlp_scaling_method", "gradient-based")
set_optimizer_attribute(model, "obj_scaling_factor", 1.0)
set_optimizer_attribute(model, "nlp_scaling_max_gradient", Float64(1))
set_optimizer_attribute(model, "nlp_scaling_min_value", 1e-8)

# Set design variables
@variable(model, lower[i] <= x[i=1:length(initial)] <= upper[i], start=initial[i])

# Set objective value and constraint functions
@operator(model, f_size_ac, length(initial), memoized_size_ac[1],memoized_fd_ac[1])
@operator(model, f_con_wfuel, length(initial), memoized_size_ac[2],memoized_fd_ac[2])
@operator(model, f_con_b, length(initial), memoized_size_ac[3],memoized_fd_ac[3])
@operator(model, f_con_tt3, length(initial), memoized_size_ac[4],memoized_fd_ac[4])
@operator(model, f_con_gtoc, length(initial), memoized_size_ac[5],memoized_fd_ac[5])

# Objective function
@objective(model, Min, f_size_ac(x...))

# Constraints
@constraint(model, f_con_wfuel(x...) <= 1)
@constraint(model, f_con_b(x...) <= 1)
@constraint(model, f_con_tt3(x...) <= 1)
@constraint(model, f_con_gtoc(x...) >= 1)

optimize!(model)

println("Finished optimization with $(value.(x))")