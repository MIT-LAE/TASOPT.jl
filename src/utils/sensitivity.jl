# SCRIPT USAGE
# using TASOPT
# include(__TASOPTindices__)
# # List of the parameters you want to update
# params = [
#     :(ac.parg[igetas]),
#     :(ac.para[iaalt, ipcruise1:ipcruise2, 1]),
#     :(ac.pare[ieepolf,:,:]),
#     :(ac.fuselage.floor.material.ρ)
# ]
# epsilon = 1e-5
# default_model = load_default_model()
# size_aircraft!(default_model)
# TASOPT.get_sensitivity(params, model_state = default_model, eps = epsilon)

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
                println("d ",idx)
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

# Function to set nested fields dynamically
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

# Central difference calculation
function central_diff_run(eps, par, model_state)
    field_path, index = par
    # Get property from default model
    x = getNestedProp(model_state, field_path, index)
    x_both = [x*(1+eps), x*(1-eps)]
    f = []
    for x_i in x_both
        # Reset model
        ac = deepcopy(model_state)
        setNestedProp!(ac, field_path, x_i, index)
        size_aircraft!(ac)
        f_out = ac.parm[imPFEI]
        push!(f, f_out)
    end
    if typeof(x_both[1]) == Vector{Float64} || typeof(x_both[2]) == Vector{Float64}
        @warn  "Specific range of values given ... returning vector of sensitivities "
        vecSens = []
        for (idx, each_x) in enumerate(x_both[1])
            fd_i =  (f[1] - f[2]) / (x_both[1][idx] - x_both[2][idx])
            push!(vecSens, fd_i)
        end
        return vecSens
    else
        return (f[1] - f[2]) / (x_both[1] - x_both[2])
    end
end

"""
get_sensitivity(model_state, input_params, eps=1e-5)

what format to give params?


"""
function get_sensitivity(input_params; model_state=nothing, eps=1e-5)
    if isnothing(model_state)
        @info "Aircraft model not provided. Using Default sized model"
        model_state = load_default_model()
        TASOPT.size_aircraft!(model_state)
    end
    params = map(p -> format_params(expr_to_string(p)), input_params)
    fd_array = []
    for param in params
        finite_diff = central_diff_run(eps, param, model_state)
        push!(fd_array, finite_diff)
    end
    return fd_array
end
