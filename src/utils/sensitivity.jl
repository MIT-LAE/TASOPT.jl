using TASOPT
include(__TASOPTindices__)
# List of the parameters you want to update
params = [
    # AERO params
    ([:para],[iaalt, ipcruise1:ipcruise2, 1]),

    # Geometry params
    # ([:parg], [igetas,1,1]),
    # ([:parg], [igfLt,1,1]), 
    # ([:parg], [igfLo,1,1]), 
#-1401.57092229
#-1401.54289115
    # ENGINE params
    ([:pare], [ieepolf,:,:]), #-0.439
    # ([:pare], [ieTt4,ipcruise1,1]), 
    # (:pare, ieepollc), #-0.256
    # (:pare, ieepolhc), #-0.725
    # (:pare, ieepolht), #-0.614
    # (:pare, ieepollt), #-0.6
    # (:parg, igTmetal), #-4.177e-5

    # # MATERIAL params
    # ([:fuselage, :material, :E], nothing),
    # ([:fuselage, :floor, :material, :τmax], nothing),
    # ([:fuselage, :floor, :material, :σmax], nothing),
    # ([:fuselage, :floor, :material, :ρ], nothing),
    # ([:fuselage, :skin, :material, :ρ], nothing),
    # ([:fuselage, :cone, :material, :ρ], nothing),
    # ([:fuselage, :cone, :material, :G], nothing),
    # ([:fuselage, :bendingmaterial_h, :material, :ρ], nothing),
    # ([:fuselage, :bendingmaterial_h, :material, :σmax], nothing),
    # ([:fuselage, :bendingmaterial_v, :material, :ρ], nothing),
    # ([:fuselage, :bendingmaterial_v, :material, :σmax], nothing),
    # ([:parg], igsigcap),
    # ([:parg], igtauweb),
    # ([:parg], igEcap),
    # ([:parg], igrhocap),
    # ([:parg], igrhoweb),


]

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

function calculate_fd(params, model_state, eps=1e-5)
    fd_array = []
    for param in params
        finite_diff = central_diff_run(eps, param, model_state)
        push!(fd_array, finite_diff)
    end
    return fd_array
end

epsilon = 1e-5
default_model = load_default_model()
size_aircraft!(default_model)
calculate_fd(params,default_model,epsilon)


