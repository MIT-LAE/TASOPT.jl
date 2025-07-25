# TASOPT.jl Optimization Example
# 1. Import modules
using TASOPT
using Printf
using Plots
using NLopt

# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# Import indices for calling aircraft parameters
include(__TASOPTindices__)

# Initialize arrays for tracking optimization progress
xarray = []
farray = []
PFEIarray = []
CDarray = []
OPRarray = []
plot_obj = nothing

# Load default aircraft model
ac = load_default_model()
# ac.engine.model.weight_model_name = "fitzgerald_adv"
size_aircraft!(ac)

# For constraint logging
struct Constraint
    name::String
    value::Float64
    limit::Float64
    penalty::Float64
end
iter = 0
# Objective function
function obj(x, grad)
    # Update wing parameters
    wing = ac.wing
    wing.layout.AR = x[1]                                    # Aspect Ratio 
    wing.layout.sweep = x[3]                                 # Sweep angle [deg]
    wing.inboard.λ = x[5]                                    # Inner panel taper ratio
    wing.outboard.λ = x[6]                                   # Outer panel taper ratio
    wing.inboard.cross_section.thickness_to_chord = x[7]    # Root thickness-to-chord ratio
    wing.outboard.cross_section.thickness_to_chord = x[8]   # Spanbreak thickness-to-chord ratio

    # Update flight condition parameters
    ac.para[iaCL, ipclimb1+1:ipdescentn-1, 1] .= x[2]           # Cruise lift coefficient
    ac.para[iaalt, ipclimbn:ipcruise1, 1] .= x[4]           # Cruise altitude [ft]
    ac.para[iarcls, ipclimb2:ipdescent4, 1] .= x[9]         # Break/root CL ratio = cls/clo
    ac.para[iarclt, ipclimb2:ipdescent4, 1] .= x[10]        # Tip/root CL ratio = clt/clo
    
    # Update engine parameters
    ac.pare[ieTt4, ipcruise1:ipcruise2, 1] .= x[11]         # Turbine inlet temperature [K]
    ac.pare[iepihc, ipcruise1, 1] = x[12]        # High pressure compressor pressure ratio
    ac.pare[iepif, ipcruise1, 1] = x[13]                   # Fan pressure ratio
    ac.pare[iepilc, ipcruise1, 1] = 3.0          # Low pressure compressor pressure ratio (fixed)
    ac.pare[ieBPR, ipcruise1, 1] = x[14] # Bypass ratio

    # Size aircraft with new parameters
    try
        TASOPT.size_aircraft!(ac, iter=50, printiter=false)
    catch e
        println("Aircraft sizing failed: $e")
        return 1e10  # Return large penalty value
    end

    # Extract objective function value (PFEI)
    f = ac.parm[imPFEI]
    
    # Store optimization history
    push!(PFEIarray, f)
    push!(xarray, copy(x))
    push!(CDarray, ac.para[iaCD, ipcruise1, 1])
    push!(OPRarray, ac.pare[iept3, ipcruise1, 1] / ac.pare[iept2, ipcruise1, 1])

    # Apply constraints with penalty functions
    total_penalty = 0.0
    violated_constraints = [] 
    # 1. Maximum span constraint
    bmax = wing.layout.max_span
    b = wing.span
    if b > bmax
        constraint = b/bmax - 1.0
        penalty = 25.0 * ac.parg[igWpay] * constraint^2
        total_penalty += penalty
        push!(violated_constraints, Constraint("Wing span", b, bmax, penalty))
    end

    # 2. Minimum climb gradient constraint
    gtocmin = ac.parg[iggtocmin]
    gtoc = ac.para[iagamV, ipclimbn, 1]
    if gtoc < gtocmin
        constraint = 1.0 - gtoc/gtocmin
        penalty = 1.0 * ac.parg[igWpay] * constraint^2
        total_penalty += penalty
        push!(violated_constraints, Constraint("Climb gradient", gtoc, gtocmin, penalty))
    end

    # 3. Maximum turbine temperature constraint
    Tt3max = 900.0  # [K]
    Tt3 = maximum(ac.pare[ieTt3, :, 1])
    if Tt3 > Tt3max
        constraint = Tt3/Tt3max - 1.0
        penalty = 5.0 * ac.parg[igWpay] * constraint^2
        total_penalty += penalty
        push!(violated_constraints, Constraint("Turbine Tt3", Tt3, Tt3max, penalty))

    end

    # 4. Fuel volume constraint
    Wfmax = ac.parg[igWfmax]
    Wf = ac.parg[igWfuel]
    if Wf > Wfmax
        constraint = Wf/Wfmax - 1.0
        penalty = 10.0 * ac.parg[igWpay] * constraint^2
        total_penalty += penalty
        push!(violated_constraints, Constraint("Fuel volume", Wf, Wfmax, penalty))
    end

    # Optional additional constraints (commented out for flexibility)
    # 5. Maximum metal temperature constraint
    # Tvanemax = 1333.33  # [K]
    # Tvane = maximum(ac.pare[ieTmet1, :, 1])
    # if Tvane > Tvanemax
    #     constraint = Tvane/Tvanemax - 1.0
    #     penalty = 5.0 * ac.parg[igWpay] * constraint^2
    #     total_penalty += penalty
    # end

    # 6. Maximum takeoff weight constraint
    # WTOmax = ac.parg[igWMTO]
    # WTO = ac.parm[imWTO, 1]
    # if WTO > WTOmax
    #     constraint = WTO/WTOmax - 1.0
    #     penalty = 10.0 * ac.parg[igWpay] * constraint^2
    #     total_penalty += penalty
    # end

    # 7. Maximum fan diameter constraint
    dfanmax = 2.0  # [m]
    dfan = ac.parg[igdfan]
    if dfan > dfanmax
        constraint = dfan/dfanmax - 1.0
        penalty = ac.parg[igWpay] * constraint^2
        total_penalty += penalty
        push!(violated_constraints, Constraint("Fan diameter", dfan, dfanmax, penalty))
    end

    lBF = ac.parm[imlBF] # [m]
    lBFmax = 2.4e3
    if lBF > lBFmax
        constraint = lBF/lBFmax - 1.0
        penalty = ac.parg[igWpay] * constraint^2
        total_penalty += penalty
        push!(violated_constraints, Constraint("Balanced Field Length", lBF, lBFmax, penalty))
    end
    

    # Final objective function value
    f_total = f + total_penalty
    
    # Print progress
    OPR = ac.pare[iept3, ipcruise1, 1] / ac.pare[iept2, ipcruise1, 1]
    FPR = ac.pare[iept21, ipcruise1, 1]/ac.pare[iept2, ipcruise1, 1]
    global iter += 1
    
    diags_to_print = Dict(
        "PFEI" => f,
        "MTOM" => ac.parg[igWMTO]/1e3/gee,
        "Wfuse" => ac.fuselage.weight/1e3/gee,
        "Wwing" => ac.wing.weight/1e3/gee,
        "Weng" => ac.parg[igWebare]/1e3/gee,
        "Wnace" => ac.parg[igWnace]/1e3/gee,
        "FPR" => FPR,
        "OPR" => OPR,
        "Dfan" => ac.parg[igdfan],
        "TSFC[g/kN]" => ac.pare[ieTSFC, ipcruise1, 1]/gee *1e3 *1e3,
    )

    

    if iter == 1 || iter % 10 == 0
        @printf("%-5s", "Iter")
        for key in keys(diags_to_print)
            @printf("│ %-10s", key)
        end
        println()
    end
    @printf("%-5d", iter)
    for (key, val) in diags_to_print
        @printf("│ %10.3f", val)
    end
    for c in violated_constraints
        @printf("│ %10s", c.name*"⚠️")
    end
    println()

    # println("Iteration: X̄ = $(round.(x, digits=3))  ⇨  PFEI = $(round(f, digits=4)), " *
    #         "f_total = $(round(f_total, digits=3)), OPR = $(round(OPR, digits=3)), FPR = $(round(FPR, digits=3))")
    # if !isempty(violated_constraints)
    # println("⚠️ Violated Constraints:")
    #     for c in violated_constraints
    #         @printf("  • %-20s Value: %8.4f | Limit: %8.4f | Penalty: %.1f\n",
    #                 c.name, c.value, c.limit, c.penalty)
    #     end
    # end
    push!(farray, f_total)
    
    return f_total
end

# Design variable bounds and initial values
# Variables: [AR, CL, sweep, altitude, λ_in, λ_out, t/c_root, t/c_span, rcls, rclt, Tt4, π_hc, π_f, BPR]

lower = [6.0,  0.45, 25.0, 10000.0, 0.65, 0.1,  0.125, 0.125, 0.9,  0.7,  1400.0, 10.0, 1.25, 1.0]
upper = [18.0, 0.75, 30.0, 20000.0, 0.85, 0.4,  0.15,  0.15,  1.3,  1.0,  1650.0, 15.0, 2.0, 20.0]

# Set initial values based on default aircraft or user specification
initial = [
    ac.wing.layout.AR,                                      # Current aspect ratio
    ac.para[iaCL, ipcruise1, 1],                           # Current cruise CL
    ac.wing.layout.sweep,                                   # Current sweep angle
    ac.para[iaalt, ipcruise1, 1],                          # Current cruise altitude
    ac.wing.inboard.λ,                                     # Current inner taper ratio
    ac.wing.outboard.λ,                                    # Current outer taper ratio
    ac.wing.inboard.cross_section.thickness_to_chord,     # Current root t/c
    ac.wing.outboard.cross_section.thickness_to_chord,    # Current span t/c
    ac.para[iarcls, ipcruise1, 1],                         # Current rcls
    ac.para[iarclt, ipcruise1, 1],                         # Current rclt
    ac.pare[ieTt4, ipcruise1, 1],                          # Current Tt4
    ac.pare[iepihc, ipcruise1, 1],                         # Current π_hc
    ac.pare[iepif, ipcruise1, 1],                          # Current π_f
    ac.pare[ieBPR, ipcruise1, 1],                          # Current BPR
]

# Ensure initial values are within bounds
for i in 1:length(initial)
    initial[i] = max(lower[i], min(upper[i], initial[i]))
end

# Initial step sizes for optimization
initial_dx = [0.5, 0.05, 0.1, 200.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 100.0, 0.5, 0.05, 1.0]

println(length(initial), length(initial_dx), length(upper), length(lower))
# Optimization settings
f_tol_rel = 1e-6
maxeval = 1000  # Maximum number of function evaluations

# Set up NLopt optimizer
opt = NLopt.Opt(:LN_NELDERMEAD, length(initial))
# Alternative optimizers:
# opt = NLopt.Opt(:LN_BOBYQA, length(initial))    # Good for smooth functions
# opt = NLopt.Opt(:LN_COBYLA, length(initial))    # Handles constraints well
# opt = NLopt.Opt(:LN_SBPLX, length(initial))     # Subplex algorithm

# Configure optimizer
opt.lower_bounds = lower
opt.upper_bounds = upper
opt.min_objective = obj
opt.initial_step = initial_dx
opt.ftol_rel = f_tol_rel
opt.maxeval = maxeval

# Print optimization setup
println("="^60)
println("TASOPT.jl Optimization Setup")
println("="^60)
println("Optimizer: $(opt.algorithm)")
println("Number of design variables: $(length(initial))")
println("Relative tolerance: $f_tol_rel")
println("Maximum evaluations: $maxeval")
println("Initial design variables:")
for (i, val) in enumerate(initial)
    println("  x[$i] = $(round(val, digits=3)) ∈ [$(lower[i]), $(upper[i])]")
end
println("="^60)

# Run optimization
println("Starting optimization...")
(optf, optx, ret) = NLopt.optimize(opt, initial)

# opt_time = @elapsed begin
#     try
#         (optf, optx, ret) = NLopt.optimize(opt, initial)
#     catch e
#         println("Optimization failed: $e")
#         optf, optx, ret = Inf, initial, :FAILURE
#     end
# end

# numevals = opt.numevals

# # Print results
# println("="^60)
# println("Optimization Results")
# println("="^60)
# println("Final objective value: $(round(optf, digits=6))")
# println("Optimization return code: $ret")
# println("Number of function evaluations: $numevals")
# println("Total optimization time: $(round(opt_time/60, digits=2)) minutes")
# println("Optimal design variables:")
# for (i, val) in enumerate(optx)
#     println("  x[$i] = $(round(val, digits=4))")
# end
# println("="^60)

# # Create output directory
# savedir = "./example/optimization/"
# if !isdir(savedir)
#     mkpath(savedir)  # Create directory and any necessary parent directories
#     println("Created output directory: $savedir")
# end

# # Generate and save aircraft details plot
# println("Generating aircraft details plot...")
# try
#     figname = "Opt_tutorial_ac_details"
#     summplot = TASOPT.plot_details(ac, plot_obj=plot_obj)
#     savefig(summplot, joinpath(savedir, figname * ".png"))
#     println("Saved: $(joinpath(savedir, figname * ".png"))")
# catch e
#     println("Warning: Could not generate aircraft details plot: $e")
# end

# # Generate optimization history plots
# println("Generating optimization history plots...")
# try
#     # Create 2x2 layout for optimization history
#     layout = @layout [a b; c d]

#     # Individual plots
#     p1 = plot(1:length(PFEIarray), PFEIarray, 
#               xlabel="Iteration", ylabel="PFEI (J/Nm)", 
#               title="Payload-Fuel Energy Intensity",
#               linewidth=2, marker=:circle, markersize=3)

#     p2 = plot(1:length(farray), farray, 
#               xlabel="Iteration", ylabel="Objective Function", 
#               title="Total Objective (with penalties)",
#               linewidth=2, marker=:circle, markersize=3,
#               yscale=:log10)

#     p3 = plot(1:length(CDarray), CDarray, 
#               xlabel="Iteration", ylabel="Drag Coefficient", 
#               title="Cruise Drag Coefficient",
#               linewidth=2, marker=:circle, markersize=3)

#     p4 = plot(1:length(OPRarray), OPRarray, 
#               xlabel="Iteration", ylabel="Overall Pressure Ratio", 
#               title="Engine Overall Pressure Ratio",
#               linewidth=2, marker=:circle, markersize=3)

#     # Combine plots
#     combined_plot = plot(p1, p2, p3, p4,    
#                         layout=layout,
#                         size=(1200, 800),
#                         plot_title="Optimization History",
#                         titlefontsize=16)

#     # Save optimization history plot
#     figname2 = "Opt_tutorial_iterations"
#     savefig(combined_plot, joinpath(savedir, figname2 * ".png"))
#     println("Saved: $(joinpath(savedir, figname2 * ".png"))")
# catch e
#     println("Warning: Could not generate optimization history plots: $e")
# end

# println("\nOptimization complete!")
# println("Output files saved in: $savedir")