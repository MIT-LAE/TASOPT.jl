# 1. Import modules
using PyPlot
using TASOPT
using Printf
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand
include("../src/misc/index.inc")
# import indices for calling parameters
using NLopt

# Declare dataframes used for plotting later
xarray = []
farray = []
PFEIarray = []
CDarray = []
OPRarray = []
track_fig = nothing
ft_to_m = 0.3048


# DESIGN VARIABLES
# AR Alt Cl  Λ λs  λt  hboxo   hboxs   rcls    rclt 

# Load default model
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/opt_input.toml"))

# Size aircraft once to get initial values
#time_wsize = @elapsed size_aircraft!(ac)

parg = ac.parg 
para = ac.para 
pari = ac.pari
parm = ac.parm
pare = ac.pare

# Objective function
function obj(x, grad)
	parg[igAR] = x[1] # Aspect Ratio 
    para[iaalt, ipcruise1, :] .=  x[2] * ft_to_m # Cruise Altitude
    para[iaCL, ipclimb1+1:ipdescentn-1, :] .= x[3] # CL
    parg[igsweep] = x[4] # Wing sweep 
    parg[iglambdas] = x[5] #inner_panel_taper_ratio
    parg[iglambdat] = x[6] #outer_panel_taper_ratio
    parg[ighboxo] = x[7] #root_thickness_to_chord
    parg[ighboxs] = x[8] #spanbreak_thickness_to_chord
    para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= x[9]   #  rcls    break/root cl ratio = cls/clo
    para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= x[10]   #  rclt    tip  /root cl ratio = clt/clo
    pare[ieTt4, ipcruise1:ipcruise2, :] .= x[11]
    pare[iepilc] = 3
    pare[iepihc] = x[12] #lpc PR fixed  = 3 hpc design var
    pare[iepif] = x[13] #Fan PR
    # Sizing aircraft with new parameters
    TASOPT.size_aircraft!(ac, iter =50, printiter=false)
    f = parm[imPFEI]
    push!(PFEIarray, parm[imPFEI])
    push!(xarray, x)
    push!(CDarray, para[iaCD, ipcruise1, 1])
    push!(OPRarray, pare[iept3]/pare[iept2])

    # Max span constriant
    # bmax = parg[igbmax]
    # b    = parg[igb]
    # constraint  = b/bmax - 1.0
    # penfac  = 25.0* parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # # Min climb gradient
    # gtocmin = parg[iggtocmin]
    # gtoc    = para[iagamV, ipclimbn,1]
    # constraint = 1.0 - gtoc/gtocmin
    # penfac = 1.0*parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # # Max Tt3 at TOC 
    # Tt3max = 900 
    # Tt3    = maximum(pare[ieTt3, :, 1])
    # constraint = Tt3/Tt3max - 1
    # penfac = 5.0*parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # # Max Tmetal at Rotation
    # Tvanemax = 1333.33 
    # Tvane    = maximum(pare[ieTmet1, :, 1]) #Rankine
    # constraint = Tvane/Tvanemax - 1
    # penfac = 5.0*parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2
    
    # Ensure aircraft weight makes sense
    WTOmax = ac.parg[igWMTO]
    WTO = ac.parm[imWTO,1]
    constraint = WTO/WTOmax - 1.0
    penfac = 10*parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # Ensure fuel volume makes sense
    Wfmax = parg[igWfmax]
    Wf    = parg[igWfuel]
    constraint = Wf/Wfmax - 1.0
    penfac = 10*parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # # Max Fan diameters
    # dfanmax = 2.0
    # dfan = parg[igdfan]
    # constraint = dfan/dfanmax - 1.0
    # penfac = parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # # Ensure fans will fit within fuselage
    # daftfanmax = min(3.0, 0.9*parg[igRfuse])
    # daftfan = parg[igdaftfan]
    # constraint = daftfan/daftfanmax - 1.0
    # penfac = parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2

    # # Ensure Fans will fit on wing
    # bo = parg[igbo]
    # b  = parg[igb]
    # lmax = b/2*4/5 - (bo/2+2*parg[igdfan])
    # lfans = parg[igneng]/2*parg[igdfan]*1.25
    # constraint = lfans/lmax - 1.0
    # penfac = parg[igWpay]
    # f = f + penfac*max(0.0, constraint)^2
    
    println("X̄ = $x  ⇨  PFEI = $(parm[imPFEI]) f = $f, OPR = $(pare[iept3]/pare[iept2]),")
    push!(farray, f)
    
    return f
end

# Set lower and upper limits
# DESIGN VARIABLES
#             AR    Alt(ft)  Cl     Λ     λs  λt   hboxo   hboxs   rcls    rclt     Tt4CR   iepihc iepif
lower      = [7.0 , 20000.0, 0.40, 10.0, 0.1, 0.1, 0.10,   0.10,   0.1,    0.1,     700.0,  6,      0]
upper      = [12.0, 60000.0, 0.65, 40.0, 1.0, 1.0, 0.15,   0.15,   1.4,    1.0,     2000.0, 15,     10] 

# Set initial changes
initial_dx = [ 0.5, 1000.0,  0.05, 0.1,  0.01,0.01,0.01,   0.01,   0.01,   0.01, 100, 0.5,0.5]

# Set initial values
initial =[
        parg[igAR], 33000.0, 0.57, parg[igsweep], 
        parg[iglambdas], parg[iglambdat], parg[ighboxo], 
        parg[ighboxs], para[iarcls, ipcruise1,1], para[iarclt, ipcruise1,1], 1587, 11.46, 1.66
]

# Set FTOL
f_tol_rel = 1e-5

# Set Optimization module
opt = NLopt.Opt(:LN_NELDERMEAD, length(initial))
# Other Optimization modules:
# # opt = NLopt.Opt(:LN_BOBYQA, length(initial))
# # opt = NLopt.Opt(:LN_COBYLA, length(initial))
# # opt = NLopt.Opt(:LN_SBPLX, length(initial))
# # opt = NLopt.Opt(:GN_DIRECT_L, length(initial))

# Set Optimization parameters
opt.lower_bounds = lower
opt.upper_bounds = upper
opt.min_objective = obj
opt.initial_step = initial_dx
opt.ftol_rel = f_tol_rel

opt_time = @elapsed (optf, optx, ret) = NLopt.optimize(opt, initial)
numevals = opt.numevals # the number of function evaluations

println("got $optf at $optx after $numevals iterations which took $(opt_time/60) min (returned $ret)")

savedir = "./example/optimization/"
figure()
if !isdir(savedir)
    # If it doesn't exist, create the "optimization" directory
    mkdir(savedir)
    println("The 'optimization' directory has been created.")
end

figname = "Opt_tutorial_ac_details"
global track_fig = TASOPT.plot_details(ac; ax = track_fig)
plt.savefig(savedir*figname*".png")
# savemodel("./Models/"*figname*".mdl", pari, parg, parm, para, pare, parpt, parmot, pargen)

fig, ax = plt.subplots(2,2, figsize = (12,8), dpi = 300)
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
figname2 = "./example/Optimization/local_test/Opt_tutorial_iterations"
fig.savefig(figname2*".png")