# This scipt optimize an aircraft model for dual fuel injection for the same payload (All fuel in the wing)[Modified from Aditeya]

# 1) Load TASOPT
using TASOPT
using NLopt
using PyPlot

# using Tables
include("../src/misc/index.inc")
ft_to_m = 0.3048
track_fig = nothing
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
#  load default model
# example_ac = load_default_model() # simply a synonym to read_aircraft_model()
# Alternatively you can load your desired input file 
nameAircraftModel = "../src/IO/experiment_input_3000.toml"
ac = read_aircraft_model(nameAircraftModel) # MODIFY <path> appropriately
saveName = "Etha3000nmi"
# 2.5) Change fuel type
ac.pari[iifuel] = 32 #(JetA:24 Ethanol:32 JetAEtha31%Blend: 322431)
ac.parg[igrhofuel] = 789.0 #(JetA:817.0 Ethanol:789.0 JetAEtha31%Blend: 805)
rhoFuelShell = 789.0 #Fuel density inside the cargo space [kg/m3]
# Objective function
xarray = []
farray = []
PFEIarray = []
CDarray = []
OPRarray = []
function obj(x, grad)
	ac.parg[igAR] = x[1] # Aspect Ratio 
    ac.para[iaalt, ipcruise1, :] .=  x[2] * ft_to_m # Cruise Altitude
    ac.para[iaCL, ipclimb1+1:ipdescentn-1, :] .= x[3] # CL
    ac.parg[igsweep] = x[4] # Wing sweep 
    ac.parg[iglambdas] = x[5] #inner_panel_taper_ratio
    ac.parg[iglambdat] = x[6] #outer_panel_taper_ratio
    ac.parg[ighboxo] = x[7] #root_thickness_to_chord
    ac.parg[ighboxs] = x[8] #spanbreak_thickness_to_chord
    ac.para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= x[9]   #  rcls    break/root cl ratio = cls/clo
    ac.para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= x[10]   #  rclt    tip  /root cl ratio = clt/clo
    ac.pare[ieTt4, ipcruise1:ipcruise2, :] .= x[11] # Tt4
    ac.pare[iepihc, ipclimb1+1 : ipdescentn-1, :] .= x[12] # High Pressure Compressor Pressure Ratio
    ac.pare[iepif, ipclimbn, :] .= x[13] #Fan PR 
    ac.pare[iepilc, :, :] .= 3 # Low Pressure Compressure Pressure Ratio set to 3

    # Sizing aircraft with new ac.parameters
    try
        TASOPT.size_aircraft!(ac, iter =500, printiter=false)
        f = ac.parm[imPFEI]
    catch
        println("sizing fails")
        f = 1000.
    end
    

    # Max span constriant
    # bmax = ac.parg[igbmax]
    # b    = ac.parg[igb]
    # constraint  = b/bmax - 1.0
    # penfac  = 0.01* ac.parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2
    
    # Ensure aircraft weight makes sense
    WTOmax = ac.parg[igWMTO]
    WTO = ac.parm[imWTO,1]
    constraint = WTO/WTOmax - 1.0
    penfac = 0.01*ac.parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # Ensure fuel volume makes sense
    Wfmax = ac.parg[igWfmax]
    Wf    = ac.parg[igWfuel]
    ## Compute the additional fuel tank volume available from the cargo bay
    AFuse = ac.parg[igAfuse] #fuselage crosssection area [m2]
    lShell = ac.parg[igxshell2]-ac.parg[igxshell1] #length of the cylindrical sector [m]
    WCargo = AFuse*0.45*lShell*rhoFuelShell*gee #[N]
    ## Finish additional fuel tank calculation
    constraint = Wf/(Wfmax+WCargo) - 1.0
    println("Fuel Weight = $(Wf) N")
    println("Wing Tank Weight = $(Wfmax) N")
    println("Cargo Tank Max Weight = $(WCargo) N")
    penfac = 0.01*ac.parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2
    
    println("X̄ = $x  ⇨  PFEI = $(ac.parm[imPFEI]) f = $f, OPR = $(ac.pare[iept3]/ac.pare[iept2]),")
    push!(xarray, x)
    push!(farray, f)
    push!(PFEIarray, ac.parm[imPFEI])
    push!(CDarray, ac.para[iaCD, ipcruise1, 1])
    push!(OPRarray, ac.pare[iept3]/ac.pare[iept2])
    return f
end

# Set lower and upper limits
#             AR    Alt(ft)  Cl     Λ     λs  λt   hboxo   hboxs   rcls    rclt     Tt4CR   iepihc iepif
lower      = [6.0 , 20000.0, 0.40, 10.0, 0.1, 0.1, 0.10,   0.10,   0.1,    0.1,     700.0,  6.,      0.]
upper      = [12.0, 60000.0, 0.65, 40.0, 1.0, 1.0, 0.15,   0.15,   1.4,    1.0,     2000.0, 15.,     10.] 

# Set initial changes
initial_dx = [0.5,  1000.0,  0.05, 0.1, 0.01, 0.01,0.01,   0.01,   0.01,   0.01,     100.0, 0.5,     0.2]

# Set initial values
initial =[
        ac.parg[igAR], 33000.0, 0.57, ac.parg[igsweep], 
        ac.parg[iglambdas], ac.parg[iglambdat], ac.parg[ighboxo], 
        ac.parg[ighboxs], ac.para[iarcls, ipcruise1,1], ac.para[iarclt, ipcruise1,1], 1587, 11.46, 1.66
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

figure()
savedir = "optimization"
if !isdir(savedir)
    # If it doesn't exist, create the "optimization" directory
    mkdir(savedir)
    println("The 'optimization' directory has been created.")
end
figname = saveName*"_Opt_details"
global track_fig = TASOPT.plot_details(ac; ax = track_fig)
plt.savefig(savedir*figname*".png")

fig, ax = plt.subplots(2,2, figsize = (12,8))
ax[1].plot(PFEIarray)
ax[1].set_xlabel("Iterations")
ax[1].set_ylabel("PFEI (J/Nm)")
ax[2].semilogy(farray)
ax[2].set_xlabel("Iterations")
ax[2].set_ylabel("Objective f")
ax[3].plot(CDarray)
ax[3].set_xlabel("Iterations")
ax[3].set_ylabel("CD")
ax[4].plot(OPRarray)
ax[4].set_xlabel("Iterations")
ax[4].set_ylabel("OPR")
plt.suptitle("Optimization outputs")
figname2 = saveName*"_Opt_iterations"
fig.savefig(savedir*figname2*".png")