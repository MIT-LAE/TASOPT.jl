# Example for a Multi-variable Optimization

![Optimization Iteration Plot](../assets/Opt_tutorial_iterations.png)

To run a multi variable optimization run on an aircraft model first determing your design variables. For this example the design variables are:

1. Aspect Ratio: `AR`
2. Cruise Altitude: `Alt` 
3. Lift Coefficient: `Cl`  
4. Wing Sweep: `Λ` 
5. Inner panel taper ratio: `λs`  
6. Outer panel taper ratio: `λt`  
7. Root thickness to chord: `hboxo`   
8. Spanbreak thickness to chord: `hboxs`   
9. Break/root cl ratio = cls/clo: `rcls`    
10. Tip/root cl ratio = clt/clo: `rclt` 

## Initialiation and loading models
Start the script importing `TASOPT.jl`, `PyPlot`, `index.inc`, `NLopt`.
```julia
# Import modules
using PyPlot
using TASOPT
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand
include("../src/misc/index.inc")
# import indices for calling parameters
using NLopt
```
Initialize arrays used for plotting and variables for simplifying usage
```julia
xarray = []
farray = []
PFEIarray = []
CDarray = []
WMTOarray = []
track_fig = nothing
ft_to_m = 0.3048
parg = ac.parg 
para = ac.para 
pari = ac.pari
parm = ac.parm
pare = ac.pare
```
Load aircraft model and size it to get initial values:
```julia
# Load default model
ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/opt_input.toml"))
#     datafile
# Size aircraft once to get initial values
size_aircraft!(ac)
```
## Setting Optimization Parameters
This example uses a Nedler Mead optimization aimed towards optimizing for passenger fuel emission index (PFEI) while checking for other constraints.
#### Set the Upper and Lower limits for all design variables:
```julia
# DESIGN VARIABLES
#             AR    Alt(ft)  Cl     Λ     λs  λt   hboxo   hboxs   rcls    rclt 
lower      = [7.0 , 20000.0, 0.40, 10.0, 0.1, 0.1, 0.10,   0.10,   0.1,    0.1]
upper      = [12.0, 60000.0, 0.65, 40.0, 1.0, 1.0, 0.15,   0.15,   1.4,    1.0] 

```
#### Set the initial values for all design variables:
```julia
initial =[
        parg[igAR], 33000.0, 0.57, parg[igsweep], 
        parg[iglambdas], parg[iglambdat], parg[ighboxo], 
        parg[ighboxs], para[iarcls, ipcruise1], para[iarclt, ipcruise1], 
]
```
#### Set initial dx values for all design variables:
```julia
initial_dx = [ 0.5, 1000.0,  0.05, 0.1,  0.01,0.01,0.01,   0.01,   0.01,   0.01]
```
#### Set other optimization factors:
```julia
# Set FTOL
f_tol_rel = 1e-5

# Set Optimization module
opt = NLopt.Opt(:LN_NELDERMEAD, length(initial))
# Other Optimization algorithms are also possible:
# # opt = NLopt.Opt(:LN_BOBYQA, length(initial))
# # opt = NLopt.Opt(:LN_COBYLA, length(initial))

# Set Optimization parameters
opt.lower_bounds = lower
opt.upper_bounds = upper
opt.min_objective = obj
opt.initial_step = initial_dx
opt.ftol_rel = f_tol_rel
```
## Objective Function
```julia
function obj(x, grad)
    parg[igAR] = x[1] # Aspect Ratio 
    para[iaalt, ipcruise1, :] .=  x[2] * ft_to_m # Cruise Altitude
    para[iaCL, ipclimb1+1:ipdescentn-1, :] .= x[3] # CL
    parg[igsweep] = x[4] # Wing sweep 
    parg[iglambdas]  = x[5] #inner_panel_taper_ratio
    parg[iglambdat]  = x[6] #outer_panel_taper_ratio
    parg[ighboxo] = x[7] #root_thickness_to_chord
    parg[ighboxs] = x[8] #spanbreak_thickness_to_chord
    para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= x[9]   #  rcls    break/root cl ratio = cls/clo
    para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= x[10]   #  rclt    tip  /root cl ratio = clt/clo

    # Sizing aircraft with new parameters
    TASOPT.size_aircraft!(ac, iter = 50)
    f = parm[imPFEI]
    push!(PFEIarray, parm[imPFEI])
    push!(xarray, x)
    push!(CDarray, para[iaCD, ipcruise1, 1])
    push!(WMTOarray, parg[igWMTO])

    # Max span constriant
    bmax = parg[igbmax]
    b    = parg[igb]
    constraint  = b/bmax - 1.0
    penfac  = 25.0* parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # Min climb gradient
    gtocmin = parg[iggtocmin]
    gtoc    = para[iagamV, ipclimbn]
    constraint = 1.0 - gtoc/gtocmin
    penfac = 1.0*parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2
    
    # Ensure fuel volume makes sense
    Wfmax = parg[igWfmax]
    Wf    = parg[igWfuel]
    constraint = Wf/Wfmax - 1.0
    penfac = 10*parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # Max Fan diameters
    dfanmax = 2.0
    dfan = parg[igdfan]
    constraint = dfan/dfanmax - 1.0
    penfac = parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # Ensure fans will fit within fuselage
    daftfanmax = min(3.0, 0.9*parg[igRfuse])
    daftfan = parg[igdaftfan]
    constraint = daftfan/daftfanmax - 1.0
    penfac = parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2

    # Ensure Fans will fit on wing
    bo = parg[igbo]
    b  = parg[igb]
    lmax = b/2*4/5 - (bo/2+2*parg[igdfan])
    lfans = parg[igneng]/2*parg[igdfan]*1.25
    constraint = lfans/lmax - 1.0
    penfac = parg[igWpay]
    f = f + penfac*max(0.0, constraint)^2
    
    println("X̄ = $x  ⇨  PFEI = $(parm[imPFEI]) f = $f, MTOW = $(parg[igWMTO])")
    push!(farray, f)
    
    return f
end
```
## Running the optimization
```julia
opt_time = @elapsed (optf, optx, ret) = NLopt.optimize(opt, initial)
numevals = opt.numevals # the number of function evaluations

println("got $optf at $optx after $numevals iterations which took $(opt_time/60) min (returned $ret)")

```

## Plotting resulting data
#### Plot aircraft model details:
```julia
savedir = "./example/optimization/"
figname = "Opt_tutorial_ac_details"
global track_fig = TASOPT.plot_details(ac; ax = track_fig)
plt.savefig(savedir*figname*".png")
```
#### Plot optimization outputs over iterations:
```julia
fig, ax = plt.subplots(2,2, figsize = (12,8))
ax[1].plot(PFEIarray)
ax[1].set_xlabel("Iterations")
ax[1].set_ylabel("PFEI (J/Nm)")
ax[2].plot(farray)
ax[2].set_xlabel("Iterations")
ax[2].set_ylabel("Objective f")
ax[3].plot(CDarray)
ax[3].set_xlabel("Iterations")
ax[3].set_ylabel("CD")
ax[4].plot(WMTOarray./1000)
ax[4].set_xlabel("Iterations")
ax[4].set_ylabel("WMTO (1000kg)")
plt.suptitle("Optimization outputs")
fig.savefig("./example/Optimization/Opt_tutorial_iterations.png")
```
