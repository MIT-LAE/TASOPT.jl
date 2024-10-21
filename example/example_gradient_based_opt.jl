using TASOPT
using JuMP
using Ipopt
using PythonPlot
using Test

include(__TASOPTindices__)

epsilon = 1e-5
default_model = load_default_model()
size_aircraft!(default_model)

design_variables = []
objective_array = []

input_params = [
    :(ac.parg[igAR]), 
    :(ac.para[iaCL,ipcruise1:ipcruise2,1]),
    :(ac.parg[igsweep]),
    :(ac.pare[ieTt4, ipcruise1:ipcruise2, 1]),
    :(ac.pare[iepif, ipcruise1, 1]) ,
]

    
params = map(p -> TASOPT.format_params(TASOPT.expr_to_string(p)), input_params)


lower      = [9.0 , 0.53, 25.0, 1400.0, 1.25]
upper      = [11.0, 0.60, 30.0, 1650.0, 2.0 ] 
initial    = [10.5, 0.57, 26.0, 1580.0, 1.685]


function size_ac(x...)
    # println("State: ",x)
    ac = deepcopy(default_model)
    # Set params
    for (i,x_i) in enumerate(x)
        field_path, index = params[i]
        TASOPT.setNestedProp!(ac, field_path, x_i, index)
    end
    # Size aircraft
    push!(design_variables, x)
    try
        size_aircraft!(ac,printiter=false)
        push!(objective_array, ac.parm[imPFEI])
        return ac.parm[imPFEI]
    catch
        println("wsize FAILED")
        push!(objective_array, NaN)
        return 1.0e12
    end
end

function fdiff_all!(g, x...)
    # Reset model
    ac = deepcopy(default_model)
    # Set current state
    for (i,x_i) in enumerate(x)
        field_path, index = params[i]
        TASOPT.setNestedProp!(ac, field_path, x_i, index)
    end
    # Calculate gradient
    if (size(g)[1] >0)
        gradients = TASOPT.get_sensitivity(input_params; model_state=ac, eps=epsilon, optimizer=true)
        for (k,grads) in enumerate(gradients)
            g[k] = grads
        end
    end
end

function exampleOPT()
    model = Model(Ipopt.Optimizer)

    # For almost localy solved (Takes much shorter time, almost solves):
    # set_optimizer_attribute(model, "acceptable_tol", 1e-2)
    # set_optimizer_attribute(model, "acceptable_iter", 1)
    
    @variable(model, lower[i] <= x[i=1:size(initial)[1]] <= upper[i], start=initial[i])

    function eval_f(x...)
        return size_ac(x...)
    end

    register(model, :eval_f, size(initial)[1], eval_f, fdiff_all!; autodiff = false)

    @NLobjective(model, Min, eval_f(x...))

    optimize!(model)
    println("""
    termination_status = $(termination_status(model))
    primal_status      = $(primal_status(model))
    objective_value    = $(objective_value(model))
    """)

    return [value(x[i]) for i in 1:2]
end


opt_time = @elapsed (x) = exampleOPT()

println("Optimization completed in $opt_time with x = $x")


# To plot results:
figure()
savedir = joinpath(__TASOPTroot__,"../example/optimization/")
if !isdir(savedir)
    # If it doesn't exist, create the "optimization" directory
    mkdir(savedir)
    println("The 'optimization' directory has been created.")
end

fig, ax = subplots(size(initial)[1]+1,1, figsize = (12,14), dpi=600)
ax[0].plot(objective_array)
# ax[0].set_xlabel("Iterations")
ax[0].set_ylabel("PFEI (J/Nm)")

for (idx, xi) in enumerate(initial)
    ax[idx].plot([x[idx] for x in design_variables])
    ax[idx].set_ylabel("Variable $idx")
end
ax[size(initial)[1]].set_xlabel("Iterations")
suptitle("Optimization outputs")
figname2 = "Gradient_Optimizer_iterations"
fig.savefig(savedir*figname2*".svg")



