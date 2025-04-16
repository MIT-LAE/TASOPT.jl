# 1. Import modules
using TASOpt
using Printf
using Plots
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand
include(__TASOPTindices__)
# import indices for calling ac.parameters
using NLopt

# Declare dataframes used for plotting later
xarray = []
farray = []
PFEIarray = []
CDarray = []
OPRarray = []
plot_obj = nothing
ft_to_m = 0.3048

# Load default model
ac = load_default_model() #read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/opt_input.toml"))
size_aircraft!(ac)

# Objective function
function obj(x, grad)
	ac.parg[igAR] = x[1] # Aspect Ratio 
    ac.para[iaCL,ipcruise1:ipcruise2,1] .=  x[2] #* ft_to_m # Cruise Altitude
    ac.parg[igsweep] = x[3] # CL
    ac.para[iaalt, ipclimb2:ipcruise1, 1] .= x[4] # Wing sweep 
    ac.parg[iglambdas] = x[5] #inner_panel_taper_ratio
    ac.parg[iglambdat] = x[6] #outer_panel_taper_ratio
    ac.parg[ighboxo] = x[7] #root_thickness_to_chord
    ac.parg[ighboxs] = x[8] #spanbreak_thickness_to_chord
    ac.para[iarcls, ipclimb2 : ipdescent4, 1] .= x[9]   #  rcls    break/root cl ratio = cls/clo
    ac.para[iarclt, ipclimb2 : ipdescent4, 1] .= x[10]   #  rclt    tip  /root cl ratio = clt/clo
    ac.pare[ieTt4, ipcruise1:ipcruise2, 1] .= x[11] # Tt4
    ac.pare[iepihc, ipclimb1 : ipdescentn, 1] .= x[12] # High Pressure Compressor Pressure Ratio
    ac.pare[iepif, ipcruise1, 1] .= x[13] #Fan PR 
    ac.pare[iepilc, ipclimb1 : ipdescentn, 1] .= 3 # Low Pressure Compressure Pressure Ratio set to 3

    # Sizing aircraft with new ac.parameters
    TASOPT.size_aircraft!(ac, iter =50, printiter=false)
    f = ac.parm[imPFEI]
    push!(PFEIarray, ac.parm[imPFEI])
    push!(xarray, x)
    push!(CDarray, ac.para[iaCD, ipcruise1, 1])
    push!(OPRarray, ac.pare[iept3]/ac.pare[iept2])

    # Max span constriant
    bmax = ac.parg[igbmax]
    b    = ac.parg[igb]
    constraint  = b/bmax - 1.0
    penfac  = 25.0* ac.parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # # Min climb gradient
    gtocmin = ac.parg[iggtocmin]
    gtoc    = ac.para[iagamV, ipclimbn,1]
    constraint = 1.0 - gtoc/gtocmin
    penfac = 1.0*ac.parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # # Max Tt3 at TOC 
    Tt3max = 900 
    Tt3    = maximum(ac.pare[ieTt3, :, 1])
    constraint = Tt3/Tt3max - 1
    penfac = 5.0*ac.parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # # Max Tmetal at Rotation
    # Tvanemax = 1333.33 
    # Tvane    = maximum(ac.pare[ieTmet1, :, 1]) #Rankine
    # constraint = Tvane/Tvanemax - 1
    # penfac = 5.0*ac.parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2
    
    # Ensure aircraft weight makes sense
    # WTOmax = ac.parg[igWMTO]
    # WTO = ac.parm[imWTO,1]
    # constraint = WTO/WTOmax - 1.0
    # penfac = 10*ac.parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # Ensure fuel volume makes sense
    Wfmax = ac.parg[igWfmax]
    Wf    = ac.parg[igWfuel]
    constraint = Wf/Wfmax - 1.0
    penfac = 10*ac.parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # # Max Fan diameters
    # dfanmax = 2.0
    # dfan = ac.parg[igdfan]
    # constraint = dfan/dfanmax - 1.0
    # penfac = ac.parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # # Ensure fans will fit within fuselage
    # daftfanmax = min(3.0, 0.9*ac.fuselage.layout.radius)
    # daftfan = ac.parg[igdaftfan]
    # constraint = daftfan/daftfanmax - 1.0
    # penfac = ac.parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # # Ensure Fans will fit on wing
    # bo = ac.parg[igbo]
    # b  = ac.parg[igb]
    # lmax = b/2*4/5 - (bo/2+2*ac.parg[igdfan])
    # lfans = ac.parg[igneng]/2*ac.parg[igdfan]*1.25
    # constraint = lfans/lmax - 1.0
    # penfac = ac.parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2
    
    println("X̄ = $x  ⇨  PFEI = $(ac.parm[imPFEI]) f = $f, OPR = $(ac.pare[iept3]/ac.pare[iept2]),")
    push!(farray, f)
    
    return f
end

# Set lower and upper limits
# DESIGN VARIABLES

lower      = [9.0 , 0.53, 25.0, 10000.0, 0.65, 0.1,  0.125,    0.125,    0.9,   0.7, 1400.0, 10.0,  1.25, 2.98]
upper      = [11.0, 0.60, 30.0, 10900.0, 0.85, 0.4,  0.15,   0.15,   1.3,   1.0, 1650.0, 15.0, 2.0, 3.02] 
initial    = [10.5, 0.57, 26.0, 10668.0, 0.7, 0.25, 0.1268, 0.1266, 1.238, 0.9, 1580.0, 12.0, 1.685, 3.0 ]

# Set initial changes
initial_dx = [ 0.5,  0.05, 0.1, 100.0, 0.01,0.01,0.01,   0.01,   0.01,   0.01, 100, 0.5,0.05, 0.01]

# # Set initial values
initial =[
        ac.parg[igAR], 33000.0, 0.57, ac.parg[igsweep], 
        ac.parg[iglambdas], ac.parg[iglambdat], ac.parg[ighboxo], 
        ac.parg[ighboxs], ac.para[iarcls, ipcruise1,1], ac.para[iarclt, ipcruise1,1], 1587, 11.46, 1.66, 3.0
]

# Set FTOL
f_tol_rel = 1e-5

# Set Optimization module
opt = NLopt.Opt(:LN_NELDERMEAD, length(initial))
# Other Optimization algorithms are also possible:
# # opt = NLopt.Opt(:LN_BOBYQA, length(initial))
# # opt = NLopt.Opt(:LN_COBYLA, length(initial))

# Set Optimization ac.parameters
opt.lower_bounds = lower
opt.upper_bounds = upper
opt.min_objective = obj
opt.initial_step = initial_dx
opt.ftol_rel = f_tol_rel

opt_time = @elapsed (optf, optx, ret) = NLopt.optimize(opt, initial)
numevals = opt.numevals # the number of function evaluations

println("got $optf at $optx after $numevals iterations which took $(opt_time/60) min (returned $ret)")

savedir = "./example/optimization/"
if !isdir(savedir)
    # If it doesn't exist, create the "optimization" directory
    mkdir(savedir)
    println("The 'optimization' directory has been created.")
end

figname = "Opt_tutorial_ac_details"
summplot = TASOPT.plot_details(ac, plot_obj=plot_obj)
savefig(summplot, savedir*figname*".png")


## Second figure
# Create a 2x2 layout
layout = @layout [a b; c d]

p1 = plot(PFEIarray, xlabel="Iterations", ylabel="PFEI (J/Nm)", title="")
p2 = plot(farray, yscale=:log10, xlabel="Iterations", ylabel="Objective f", title="")
p3 = plot(CDarray, xlabel="Iterations", ylabel="CD", title="")
p4 = plot(OPRarray, xlabel="Iterations", ylabel="OPR", title="")

# Create the plot
p = plot(p1, p2, p3, p4,    
    layout = layout,
    size=(1200, 800),
    plot_title="Optimization outputs"
)

# Save the plot
figname2 = "Opt_tutorial_iterations"
savefig(p, savedir * figname2 * ".png")