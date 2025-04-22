"""
`sensitivity` is a module that calculates the gradients of different TASOPT params
using central finite difference (relative change)

Example Usage:

using TASOPT
include(__TASOPTindices__)
# List of the parameters you want to update as symbols
    params = [
    :(ac.parg[igetas]),
    :(ac.pare[ieepolf,:,:]),
    :(ac.pare[iepihc,ipclimbn:ipcruise1,1]),
    :(ac.fuselage.layout.cross_section.radius)
]
epsilon = 1e-5
default_model = load_default_model()
size_aircraft!(default_model)
TASOPT.get_sensitivity(params, model_state = default_model, eps = epsilon)

If you want the default model as the model state and epsilon as 1e-5 
you can also call the function directly with just the params:

TASOPT.get_sensitivity(params)
"""


"""
    expr_to_string(expr)

`expr_to_string` converts a Julia expression representing a `struct` field access or indexed array access into a string format.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `expr`: A Julia expression representing a `struct` field access or indexed array access.

    **Outputs:**
    - A string representation of the input expression.

"""
function expr_to_string(expr)
    if typeof(expr) == Symbol
        # If it's just a symbol (like ac), return it as a string
        return string(expr)
    elseif expr.head == :. 
        # For nested fields (like ac.fuselage.floor.material.ρ)
        return string(expr)
    elseif expr.head == :ref
        # For indexed access (like ac.parg[igetas])
        base = string(expr.args[1])
        indices = join(map(arg_to_string, expr.args[2:end]), ", ")
        return "$base[$indices]"
    else
        error("Unsupported expression type")
    end
end

"""
    arg_to_string(arg)

`arg_to_string` converts various types of Julia expressions or values into their string representations, particularly useful for handling array indices and ranges.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `arg`:  Expression, symbol, integer, or other value representing an array index or range.

    **Outputs:**
    - A string representation of the input argument.

"""
function arg_to_string(arg)
    if typeof(arg) == Symbol
        return string(arg)
    elseif typeof(arg) == Expr && arg.head == :call && arg.args[1] == :(:)
        # For range expressions like ipcruise1:ipcruise2
        return "$(arg_to_string(arg.args[2])):$(arg_to_string(arg.args[3]))"
    elseif typeof(arg) == Int
        return string(arg)
    else
        return string(arg)
    end
end


"""
    format_params(param_str)

`format_params` processes a string representation of a `struct` field access or indexed array access and converts it into a tuple format suitable for sensitivity analysis.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `param_str`: A string representing a `struct` field access or indexed array access (e.g., "parg[igetas]" or "fuselage.floor.material.ρ").

    **Outputs:**
    - A tuple containing:
      1. A vector of symbols representing the field path.
      2. Either `nothing` for nested `struct` fields, or a vector of processed indices for array access.

!!! note "Processing Details"
    - Removes 'ac.' prefix if present.
    - Handles nested `struct` fields (e.g., "fuselage.floor.material.ρ").
    - Processes array indices, including:
      - Single indices
      - Colon (:) for full range
      - Ranges (e.g., ipcruise1:ipcruise2)
    - Pads index vector with ones if fewer than 3 indices are provided.

"""
function format_params(param_str)
    # Remove 'ac.' prefix if present
    param_str = replace(param_str, r"^ac\." => "")
    
    # Split the string into field path and index parts
    parts = split(param_str, '[')
    
    if length(parts) == 1  # No brackets, it's a nested struct field
        field_path = Symbol.(split(parts[1], '.'))
        return (field_path, nothing)
    else
        # Extract field path
        field_path = Symbol.(split(parts[1], '.'))
        
        # Extract and process index part
        index_str = strip(parts[2], ']')
        indices = split(index_str, ',')
        
        # Process each index
        processed_indices = []
        for idx in indices
            if strip(idx) == ":"
                push!(processed_indices, Colon())
            elseif contains(idx, ':')  # Range
                range_parts = split(idx, ':')
                idx_1 = eval(Symbol(strip(range_parts[1])))
                idx_2 = eval(Symbol(strip(range_parts[2])))
                push!(processed_indices, UnitRange{Int}(idx_1,idx_2))
            elseif tryparse(Float64, idx) !== nothing
                idx_1 = parse(Int64, idx)
                push!(processed_indices, idx_1)
            else
                push!(processed_indices, eval(Symbol(strip(idx))))
            end
        end
        
        # Pad with ones if less than 3 indices
        while length(processed_indices) < 3
            push!(processed_indices, 1)
        end
        
        return (field_path, processed_indices)
    end
end

"""
    getNestedProp(ac, field_path::Vector, index=nothing)

`getNestedProp` retrieves the value of a nested property from an object, optionally accessing specific indices of an array property.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac`: The root object to traverse.
    - `field_path::Vector`: A vector of symbols representing the path to the desired property.
    - `index=nothing`: Optional. If provided, should be a vector of up to 3 indices or `Colon()` objects for array access.

    **Outputs:**
    - The value of the specified nested property, or a specific element/slice if indices are provided.

!!! note "Behavior"
    - Traverses the object structure using the provided `field_path`.
    - If `index` is `nothing`, returns the entire property value.
    - If `index` is provided:
      - For 1D or 2D array access (when `index[2]` or `index[3]` is `Colon()`), returns `val[index[1]]`.
      - For 3D array access, returns `val[index[1],index[2],index[3]]`.

"""
function getNestedProp(ac, field_path::Vector, index=nothing)
    val = ac
    for field in field_path
        val = getfield(val, field)
    end
    if isnothing(index)
        return val 
    elseif index[2] == Colon() || index[3] == Colon()
        return val[index[1]]
    else
        return val[index[1],index[2],index[3]]
    end
end

"""
    setNestedProp!(ac, field_path::Vector, x, index=nothing)

`setNestedProp!` sets the value of a nested property in an object, optionally updating specific indices of an array property.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac`: The root object to traverse and modify.
    - `field_path::Vector`: A vector of symbols representing the path to the property to be modified.
    - `x`: The new value to set.
    - `index=nothing`: Optional. If provided, should be a vector of up to 3 indices, `Colon()` objects, or `UnitRange{Int64}` for array access.

    **Outputs:**
    - The function modifies the object in-place and does not return a value.

!!! note "Behavior"
    - Traverses the object structure using the provided `field_path`.
    - If `index` is `nothing`, sets the entire property value.
    - If `index` is provided:
      - For array access with `Colon()` or `UnitRange{Int64}`, uses broadcasting (`.=`) to set values.
      - For specific index access, sets the value directly.
"""
function setNestedProp!(ac, field_path::Vector, x, index=nothing)
    val = ac
    for (i, field) in enumerate(field_path[1:end-1])
        val = getfield(val, field)
    end
    # Modify the deepest field
    if isnothing(index)
        setfield!(val, field_path[end], x)
    elseif index[2] == Colon() || index[3] == Colon() || typeof(index[2]) == UnitRange{Int64} || typeof(index[3]) == UnitRange{Int64}
        getfield(val, field_path[end])[index[1],index[2],index[3]] .= x
    else
        getfield(val, field_path[end])[index[1],index[2],index[3]] = x
    end
end

"""
    check_if_multiple_vals(x,i)

`check_if_multiple_vals` Is a helper function to transpose the sensitivities array in case 
    there is an array of output functions and range of input parameters
"""
function check_if_multiple_vals(x,i)
    if isa(x[1], Float64)
        return x[i]
    else
        return map(y -> y[i], x)
    end
end

"""
    central_diff_run(eps, par, model_state)

`central_diff_run` calculates the finite differences using relative central difference
    **Inputs:**
    - `eps`: Epsilon, the relative difference to change x 
    - `par`: The specific parameter to be modified.
    - `model_state`: The model state to be taken as "default" over which the FD is calculated.

    **Outputs:**
    - The function returns the finite difference
"""
function central_diff_run(eps, par, model_state; optimizer=false, f_out_fn=nothing)
    field_path, index = par
    # Get property from default model
    x = getNestedProp(model_state, field_path, index)
    x_both = [x*(1+eps), x*(1-eps)]
    
    # Initialize results array for each perturbation
    f_results = []
    
    for x_i in x_both
        # Reset model
        ac = deepcopy(model_state)
        setNestedProp!(ac, field_path, x_i, index)
        
        # Initialize outputs for this perturbation
        if isnothing(f_out_fn)
            f_out = Inf
            try
                TASOPT.size_aircraft!(ac, printiter=false)
                f_out = ac.parm[imPFEI]
            catch
                println("size_aircraft! FAILED")
            end
            push!(f_results, f_out)
        else
            # Handle single function or array of functions
            f_out = []
            try
                TASOPT.size_aircraft!(ac, printiter=false)
                if isa(f_out_fn, Array)
                    # Multiple output functions
                    for fn in f_out_fn
                        push!(f_out, fn(ac))
                    end
                else
                    # Single output function
                    push!(f_out, f_out_fn(ac))
                end
            catch
                println("size_aircraft! FAILED")
                if isa(f_out_fn, Array)
                    f_out = fill(Inf, length(f_out_fn))
                else
                    f_out = [Inf]
                end
            end
            push!(f_results, f_out)
        end
    end
    
    # Calculate finite differences
    if typeof(x_both[1]) == Vector{Float64} || typeof(x_both[2]) == Vector{Float64}
        vecSens = []
        for (idx, each_x) in enumerate(x_both[1])
            if isa(f_results[1], Array)
                # Multiple output functions
                fd_i = [(f1 - f2) / (x_both[1][idx] - x_both[2][idx]) 
                       for (f1, f2) in zip(f_results[1], f_results[2])]
            else
                # Single output
                fd_i = (f_results[1] - f_results[2]) / (x_both[1][idx] - x_both[2][idx])
            end
            push!(vecSens, fd_i)
        end
        
        if optimizer
            if isa(f_results[1], Array)
                # Return average sensitivity for each output function
                return [sum(getindex.(vecSens, i))/length(vecSens) for i in 1:length(f_results[1])]
            else
                return sum(vecSens)/length(vecSens)
            end
        else
            @warn "Specific range of values given ... returning vector of sensitivities"
            return vecSens
        end
    else
        if isa(f_results[1], Array)
            # Multiple output functions
            return [(f1 - f2) / (x_both[1] - x_both[2]) 
                   for (f1, f2) in zip(f_results[1], f_results[2])]
        else
            # Single output
            return (f_results[1] - f_results[2]) / (x_both[1] - x_both[2])
        end
    end
end

"""
    plot_sensitivities(sensitivities::Vector, output_file::String="sensitivity_plot.png")

`plot_sensitivities` plots the finite differences
    **Inputs:**
    - `sensitivities`: List of caclulated sensitivities
    - `output_file`: File location of the plot figure
"""
function plot_sensitivities(sensitivities::Vector, output_file::String="sensitivity_plot.png")

    #TODO: Sensitivity plotter should list the parameter name as yticks

    # Initialize plot
    p = plot(
        legend=false,
        xlabel="Sensitivity",
        ylabel="Element Number",
        title="Sensitivity Plot",
        dpi=300,
        size=(800, 600),
        yticks=(1:length(sensitivities), string.(1:length(sensitivities)))
    )

    # Store the y-axis positions
    y_positions = []

    # Loop through the sensitivities to plot each value
    for i in 1:length(sensitivities)
        if isa(sensitivities[i], Number)
            bar!(p, [i], [sensitivities[i]], color=:blue, bar_width=0.5, orientation=:horizontal)
            push!(y_positions, i)
        elseif isa(sensitivities[i], Vector)
            # Plot each element in the sub-list as a separate bar, side by side
            for j in 1:length(sensitivities[i])
                offset = (j - 1) * 0.25  # Slight offset for side-by-side bars
                bar!(p, [i + offset], [sensitivities[i][j]], color=:green, bar_width=0.25, orientation=:horizontal)
            end
            push!(y_positions, i)  # Append only once to avoid extra ticks
        end
    end

    # Draw a vertical line at x = 0
    vline!(p, [0], linestyle=:dash, color=:black, linewidth=1)

    # Invert y-axis to match the order of elements
    yflip!(p, true)

    # Save the figure to an image file
    # savefig(p, output_file)
    # println("Plot saved as $output_file")

    return p
end

"""
    get_sensitivity(input_params; model_state=nothing, eps=1e-5)

`get_sensitivity` starts the finite differences calculation
    **Inputs:**
    - `input_params`: List of Input Params as symbols.
    - `model_state`: The model state to be taken as "default" over which the FD is calculated.
    - `eps`: Epsilon, the relative difference to change x 
    
    
    **Outputs:**
    - The function returns the finite difference of each input param in a vector.  
"""
function get_sensitivity(input_params; model_state=nothing, eps=1e-5, optimizer=false, f_out_fn=nothing)
    if isnothing(model_state)
        @info "Aircraft model not provided. Using Default sized model"
        model_state = load_default_model()
        TASOPT.size_aircraft!(model_state, printiter=false)
    end
    
    if isnothing(f_out_fn)
        @info "Output function not defined, using PFEI"
    elseif isa(f_out_fn, Array) && !optimizer
        @info "Multiple output functions provided, will return sensitivities for each"
    end
    
    params = map(p -> format_params(expr_to_string(p)), input_params)
    fd_array = []
    
    for param in params
        finite_diff = central_diff_run(eps, param, model_state; 
                                     optimizer=optimizer, f_out_fn=f_out_fn)
        push!(fd_array, finite_diff)
    end
    
    if isa(f_out_fn, Array)
        # Multiple output function, transpose array
        fd_array = [map(x -> check_if_multiple_vals(x,i), fd_array) for i in 1:length(fd_array[1])]
    end
    return fd_array
end