
using Pkg
Pkg.activate("../")
Pkg.instantiate()

# 1) Load TASOPT
include("../tasopt.jl")

# Inputs for testing runs
nmisx = 1
pari = zeros(Int64, iitotal)
parg = zeros(Float64, igtotal)
parm = zeros(Float64, (imtotal, nmisx))
para = zeros(Float64, (iatotal, iptotal, nmisx))
pare = zeros(Float64, (ietotal, iptotal, nmisx))

ft_to_m = 0.3048
in_to_m = 0.0254
nmi_to_m = 1852.0
deg_to_rad = π / 180.0
lbf_to_N = 4.45
kts_to_mps = 0.514
hp_to_W = 745.7

pari[iiaircraftclass] = 737 #specifies aircraft class in order to run appropriate NPSS model
pari[iiopt] = 0 # 0 run sizing loop only; 1 run optimization procedure
pari[iifwing] = 1 # 0 = all fuel stored in tanks; 1 = all fuel stored in wings
pari[iiengtype] = 1  # 0 = Turboelectric engine; 1 = Turbofan engine

#-----------------------------------------
#- parameters for each mission 
parm[imRange, :] .= 3000.0 * nmi_to_m       # m
Npax = 180
Wpax = 215.0 * lbf_to_N
parm[imWpay, :] .= Npax * Wpax # [N]
altTO = 0.0 * ft_to_m
parm[imaltTO, :] .= altTO # altTO    takeoff/landing altitude [m]
T0TO = 288.0
parm[imT0TO, :] .= T0TO # ambient takeoff/landing temp [K]

#-----------------------------------------
# Design cruise-start altitude (with max fuel, max payload)
altCR = 35000.0 * ft_to_m

#-----------------------------------------
#- Takeoff and initial climb parameters
clpmax = 2.25

para[iaalt, ipstatic:ipcutback, :] .= altTO
para[iaclpmax, ipstatic:ipcutback, :] .= clpmax

para[iaalt, ipclimb1, :] .= altTO
para[iaclpmax, ipclimb1, :] .= clpmax

para[iaalt, ipcruise1, :] .= altCR

para[iaalt, ipdescentn, :] .= altTO
para[iaclpmax, ipdescentn, :] .= clpmax

parg[igcdefan] = 0.500
parg[igCDgear] = 0.015
parg[igCDspoil] = 0.10
parg[igmuroll] = 0.025
parg[igmubrake] = 0.35

parg[ighobst] = 35.0 * ft_to_m

parg[iglBFmax] = 8000.0 * ft_to_m
parg[iggtocmin] = 0.015
parg[igdBSLmax] = 90.0
parg[igdBCBmax] = 75.0

parm[imthCB, :] .= 40.0 * π / 180.0
parm[imgamVCB, :] .= 3.0 * π / 180.0
parm[imgamVDE1, :] .= -3.0 * π / 180.0
parm[imgamVDEn, :] .= -3.0 * π / 180.0

#-----------------------------------------
#- sizing-load parameters
parg[igNlift] = 3.0           # Nlift  max vertical load factor for wing bending loads
parg[igNland] = 6.0           # Nland  max vertical load factor for fuse bending loads
parg[igVne] = 280.0 * kts_to_mps # Vne    never-exceed IAS, for tail loads

cabinPressureAlt = 8000.0 # Altitude for cabin pressure [ft]
_, p_cabin, _, _, _ = atmos(cabinPressureAlt * ft_to_m / 1000.0)
parg[igpcabin] = p_cabin

#-----------------------------------------
#- cruise aero parameters
para[iaCL, ipclimb1+1:ipdescentn-1, :] .= 0.57
para[iaMach, ipclimbn:ipdescent1, :] .= 0.80

#----------------------------------------
#- basic wing parameters  
parg[igsweep] = 26.0         # λ wing
parg[igAR] = 10.1        # Aspect ratio
parg[igsweeph] = parg[igsweep]

parg[igbmax] = 117.5 * ft_to_m # Max span for span constraint ICAO Code D/ FAA Group IV

parg[iglambdas] = 0.7
parg[iglambdat] = 0.25

pari[iiwplan] = 1 # wing cantilever with engine 
pari[iifwcen] = 1 # Fuel in center box 

parg[igrWfmax] = 0.90

parg[igfLo] = -0.3   # fLo   fuselage lift carryover loss factor
parg[igfLt] = -0.05  # fLt   tip lift rolloff factor

parg[igzs] = 154.0 * in_to_m

parg[igbo] = 2 * (71.0 * in_to_m) # 2 × wing centerbox halfspan
parg[igetas] = 0.285 # ηs panel break eta location  (strut-attach if iwplan=2)
parg[igrVstrut] = 1.0   # Strut local/free sream velocity

#----------------------------------------
#- tail download parameter at max load case (download looks like added weight)
parg[igCLhNrat] = -0.5 #CLh/ CLmax

#----------------------------------------
#- wing spanwise cl and cm distributions over mission
#- ( rclo = clo/clo = 1.0  by definition, so it's not specified )
#- takeoff, initial climb 
para[iarcls, 1:ipclimb1, :] .= 1.1     #  rcls    break/root cl ratio = cls/clo
para[iarclt, 1:ipclimb1, :] .= 0.6     #  rclt    tip  /root cl ratio = clt/clo
para[iacmpo, 1:ipclimb1, :] .= -0.20   #  cmpo    root  cm
para[iacmps, 1:ipclimb1, :] .= -0.20   #  cmps    break cm
para[iacmpt, 1:ipclimb1, :] .= -0.02   #  cmpt    tip   cm

# Clean climb cruise descent and for wing structure sizing
para[iarcls, ipclimb1+1:ipdescentn-1, :] .= 1.238   #  rcls    
para[iarclt, ipclimb1+1:ipdescentn-1, :] .= 0.90    #  rclt    
para[iacmpo, ipclimb1+1:ipdescentn-1, :] .= -0.06   #  cmpo    
para[iacmps, ipclimb1+1:ipdescentn-1, :] .= -0.06   #  cmps    
para[iacmpt, ipclimb1+1:ipdescentn-1, :] .= -0.06   #  cmpt   

# Landing, forward CG tail sizing case
para[iarcls, ipdescentn, :] .= 1.1     #  rcls  
para[iarclt, ipdescentn, :] .= 0.5     #  rclt  
para[iacmpo, ipdescentn, :] .= -0.35   #  cmpo      
para[iacmps, ipdescentn, :] .= -0.35   #  cmps      
para[iacmpt, ipdescentn, :] .= -0.02   #  cmpt     

#-----------------------------------------
#- wing and tail structural box parameters
parg[igwbox] = 0.50   # wbox    box width/c
parg[ighboxo] = 0.1268 # hboxo   box height/c  (airfoil t/c) at root
parg[ighboxs] = 0.1266 # hboxs   box height/c  (airfoil t/c) at break and tip
parg[igrh] = 0.75   # rh      web-height/hbox ratio
parg[igXaxis] = 0.40   # Xaxis   spar box axis x/c location

parg[ighstrut] = 0.15   # hstrut  strut t/c    (used only if iwplan=2)


#----------------------------------------
#- weight fractions of flight surfaces and secondary wing components, 
parg[igfflap] = 0.20
parg[igfslat] = 0.10
parg[igfaile] = 0.04
parg[igflete] = 0.10
parg[igfribs] = 0.15 #+ 0.01 #Adding this for folding wingtips
parg[igfspoi] = 0.02
parg[igfwatt] = 0.03

#-----------------------------------------
#- horizontal and vertical tail parameters
pari[iiHTsize] = 1 # 1 = set Sh via Vh; 2 = set Sh via igCLhCGfwd at max-forward CG during landing
parg[igVh] = 1.45 # HT volume coeff
parg[igCLhCGfwd] = -0.7

pari[iiVTsize] = 1
parg[igVv] = 0.1
parg[igCLveout] = 0.5

pari[iixwmove] = 2
parg[igCLhspec] = -0.02
parg[igSMmin] = 0.05

parg[igdepsda] = 0.60

parg[igdCLnda] = 3.8

parg[igARh] = 6.0
parg[igARv] = 2.0
parg[iglambdah] = 0.25
parg[iglambdav] = 0.30
# parg[igsweeph] = 25.0
parg[igsweepv] = 25.0

parg[igboh] = 2 * (2.5 * ft_to_m) # 2 × half span
parg[igbov] = 0.0

parg[igfCDhcen] = 0.1

parg[igCLhmax] = 2.0
parg[igCLvmax] = 2.6

parg[igfhadd] = 0.30
parg[igfvadd] = 0.40

parg[igwboxh] = 0.50
parg[igwboxv] = 0.50
parg[ighboxh] = 0.14
parg[ighboxv] = 0.14
parg[igrhh] = 0.75
parg[igrhv] = 0.75

parg[ignvtail] = 1.0

#----------------------------------------
#- cabin and fuselage geometry parameters
parg[igRfuse] = 77 * in_to_m
parg[igdRfuse] = 15.0 * in_to_m
parg[igwfb] = 0.0 * in_to_m
parg[ighfloor] = 5.0 * in_to_m

parg[iganose] = 1.65
parg[igbtail] = 2.0

parg[igxnose] = 0.0 * ft_to_m
parg[igxend] = 124.0 * ft_to_m
parg[igxblend1] = 20.0 * ft_to_m
parg[igxblend2] = 97.0 * ft_to_m
parg[igxshell1] = 17.0 * ft_to_m
parg[igxshell2] = 102.0 * ft_to_m
parg[igxconend] = 117.0 * ft_to_m
parg[igxwbox] = 57.0 * ft_to_m
parg[igxhbox] = 114.5 * ft_to_m
parg[igxvbox] = 110.0 * ft_to_m

parg[igzwing] = -5.5 * ft_to_m
parg[igzhtail] = 0.0 * ft_to_m

pari[iiengloc] = 1 # engines on fuselage

parg[igxeng] = 52.0 * ft_to_m
parg[igyeng] = 16.0 * ft_to_m
parg[igneng] = 2.0

parg[iglambdac] = 0.3 # Tail cone taper ratio

parg[igfstring] = 0.35
parg[igfframe] = 0.25
parg[igffadd] = 0.20

parg[igWfix] = 3000.0 * 4.45  # cockpit, pilots etc converted to [N]
parg[igxfix] = 7.0 * ft_to_m

parg[igWpwindow] = 145.0 * 3.0 #[N/m]
parg[igWppinsul] = 22.0       #[N/m²]
parg[igWppfloor] = 60.0       #[N/m²]

parg[igrMh] = 0.4
parg[igrMv] = 0.7

pari[iifclose] = 0 # 0 = fuse tapers to a point; 1 = fuse end is flat

parg[igCMVf1] = 2390.0 * 0.0283 #60.0  
parg[igCLMf0] = 0.185            # CLMf1  CL where Mfuse = 0

para[iafduo, :, :] .= 0.018# 0.019    # fduo   fuselage velocity overspeed at wing root
para[iafdus, :, :] .= 0.014#0.011    # fdus   fuselage velocity overspeed at wing break
para[iafdut, :, :] .= 0.0045#0.004   # fdut   fuselage velocity overspeed at wing tip

#-----------------------------------------
#- power systems and landing gear locations and weight fractions
parg[igxhpesys] = 62.0 * ft_to_m   #  xhpesys   hyd/pneu/ele system location
parg[igxlgnose] = 14.0 * ft_to_m   #  xlgnose   nose LG location
parg[igdxlgmain] = 1.0 * ft_to_m   # dxlgmain   main LG offset behind wing lift centroid

parg[igfhpesys] = 0.010     # fhpesys    Whpesys/WMTO
parg[igflgnose] = 0.011     # flgnose    Wlgnose/WMTO
parg[igflgmain] = 0.044     # flgmain    Wlgmain/WMTO

#-----------------------------------------
#- other added-weight fractions
parg[igxapu] = 120.0 * ft_to_m # xapu      APU location
parg[igfapu] = 0.035          # fapu   Wapu/Wpay     APU weight fraction

parg[igfseat] = 0.10  # fseat    Wseat/Wpay    seat weight fraction
parg[igfpadd] = 0.35  # fpadd    Wpadd/Wpay    other payload-proportional fraction

parg[igfeadd] = 0.10  # feadd    Weadd/Wbare   engine accessories, fuel system fraction 
parg[igfpylon] = 0.10  # fpylon   Wpylon/We+a+n engine pylon weight fraction   
parg[igfreserve] = 0.20  # freserve Wfreserve/Wburn

#----------------------------------------
#- allowable stresses at sizing cases
parg[igsigfac] = 1.0   #  sigfac   convenient multiplier on all the stress values below                      

parg[igsigskin] = 15000.0 / 0.000145   # sigskin   fuselage pressurization skin stress                      
parg[igsigbend] = 30000.0 / 0.000145   # sigbend   fuselage bending skin+stringer stress                      

parg[igsigcap] = 30000.0 / 0.000145   # sigcap    wing,tail bending caps                      
parg[igtauweb] = 20000.0 / 0.000145   # tauweb    wing,tail shear webs                      
parg[igsigstrut] = 30000.0 / 0.000145   # sigstrut  strut       

#----------------------------------------
#- fuselage shell modulus ratio, for bending material sizing
parg[igrEshell] = 1.0                  # rEshell   Ebend/Eskin  ratio

#----------------------------------------
#- moduli, for strut-induced buckling load estimation
parg[igEcap] = 10.0e6 / 0.000145   # Ecap     wing sparcap
parg[igEstrut] = 10.0e6 / 0.000145   # Estrut   strut

#----------------------------------------
#- structural material densities
parg[igrhoskin] = 2700.0  #  rhoskin     fuselage skin
parg[igrhobend] = 2700.0  #  rhobend     fuselage bending stringers 
parg[igrhocap] = 2700.0  #  rhocap  	wing, tail bending caps	 
parg[igrhoweb] = 2700.0  #  rhoweb  	wing, tail shear webs	 
parg[igrhostrut] = 2700.0  #  rhostrut	strut   


#------------------------------------------------------------------
#- database for wing profile cd in transonic cruise, high climb
para[iacdfw, 1:iptotal, :] .= 0.0085  #  cdfw    wing profile cd for low speed (takeoff, initial climb)
para[iacdpw, 1:iptotal, :] .= 0.0035  #  cdpw    
para[iaRerefw, 1:iptotal, :] .= 20.0e6  #  Rerefw

para[iacdft, 1:iptotal, :] .= 0.0060  #  cdft    tail profile cd
para[iacdpt, 1:iptotal, :] .= 0.0030  #  cdpt    
para[iaRereft, 1:iptotal, :] .= 10.0e6  #  Rereft  

para[iacdfs, 1:iptotal, :] .= 0.0085  #  cdfs    strut profile cd (not used if there's no strut)
para[iacdps, 1:iptotal, :] .= 0.0035  #  cdps    
para[iaRerefs, 1:iptotal, :] .= 1.0e6   #  Rerefs  

para[iaaRexp, 1:iptotal, :] .= -0.15   #  aRexp   exponent for Re-scaling:  CD = cd * (Re/Re_ref)^aRexp

para[iafexcdw, 1:iptotal, :] .= 1.02     #  fexcdw   # wing excrescence drag factor
para[iafexcdt, 1:iptotal, :] .= 1.02     #  fexcdt   # tail excrescence drag factor
para[iafexcdf, 1:iptotal, :] .= 1.030    #  fexcdf   # fuse excrescence drag factor

parg[igfBLIw] = 0.0
parg[igfBLIf] = 0.0
pari[iiBLIc] = 0 # core in clean flow

#-----------------------------------------
#- fuel parameters
pari[iifuel] = 24
parg[igrhofuel] = 817.0
Tfuel = 280.0
pare[ieTfuel, :, :] .= Tfuel


#----------------------------------------
#- engine temperatures
parg[igTmetal] = 1280.0
Tt4TO = 1833.0
parg[igfTt4CL1] = 0.2
parg[igfTt4CLn] = 0.2
Tt4CR = 1587.0

dTstrk = 200.0
Mtexit = 1.0
StA = 0.09
efilm = 0.7
tfilm = 0.30
M4a = 0.9
ruc = 0.15

pare[ieTt4, :, :] .= Tt4CR
pare[ieM4a, :, :] .= M4a
pare[ieruc, :, :] .= ruc
pare[iedTstrk, :, :] .= dTstrk
pare[ieMtexit, :, :] .= Mtexit
pare[ieStA, :, :] .= StA
pare[ieefilm, :, :] .= efilm
pare[ietfilm, :, :] .= tfilm


pare[ieT0, ipstatic, :] .= T0TO
pare[ieT0, iprotate, :] .= T0TO
pare[ieT0, iptakeoff, :] .= T0TO

pare[ieTt4, ipstatic, :] .= Tt4TO
pare[ieTt4, iprotate, :] .= Tt4TO
pare[ieTt4, iptakeoff, :] .= Tt4TO

#------------------------------
#- design pressure ratios, efficiencies, etc.
OPR = 30.0
pilc = 8.0
pif = 1.685
pid = 0.998
pib = 0.94
pifn = 0.98
pitn = 0.989
epolf = 0.8948
epollc = 0.88
epolhc = 0.87
epolht = 0.889
epollt = 0.899

etab = 0.985

pifK = 1.685
epfK = -0.077

BPR = 5.1
Gearf = 1.0
HTRf = 0.30
HTRlc = 0.60
HTRhc = 0.80

M2 = 0.60
M25 = 0.60

#- - - - - - - - - - - - - - - - - - - - -
# power loss fractions
epsl = 0.01
epsh = 0.022

#- - - - - - - - - - - - - - - - - - - - -
# offtakes (total for both engines)
pihc = OPR / pilc
pare[iepid, :, :] .= pid
pare[iepib, :, :] .= pib
pare[iepifn, :, :] .= pifn
pare[iepitn, :, :] .= pitn
pare[iepif, :, :] .= pif
pare[iepilc, :, :] .= pilc
pare[iepihc, :, :] .= pihc
pare[ieepolf, :, :] .= epolf
pare[ieepollc, :, :] .= epollc
pare[ieepolhc, :, :] .= epolhc
pare[ieepolht, :, :] .= epolht
pare[ieepollt, :, :] .= epollt
pare[ieetab, :, :] .= etab
pare[iepifK, :, :] .= pifK
pare[ieepfK, :, :] .= epfK
pare[ieBPR, :, :] .= BPR
pare[ieM2, :, :] .= M2
pare[ieM25, :, :] .= M25
pare[ieepsl, :, :] .= epsl
pare[ieepsh, :, :] .= epsh

#- - - - - - - - - - - - - - - - - - - - -
# offtakes (total for both engines)
mofftpax = 0.00633
mofftmMTO = 0.0 * 0.001

Pofftpax = 200.0
PofftmMTO = 1.8

parg[igmofWpay] = mofftpax / Wpax
parg[igmofWMTO] = mofftmMTO / 9.81
parg[igPofWpay] = Pofftpax / Wpax
parg[igPofWMTO] = PofftmMTO / 9.81

Tt9 = 300.0
pt9 = 30000.0

pare[ieTt9, :, :] .= Tt9
pare[iept9, :, :] .= pt9

#------------------------------
#- fan nozzle area factors relative to cruise design area
A7static = 1.0
A7takeoff = 1.0
A7cutback = 1.0
A7climb1 = 1.0
A7climbn = 1.0
A7descent1 = 1.0
A7descentn = 1.0

#------------------------------
#- core nozzle area factors relative to cruise design area
A5static = 1.0
A5takeoff = 1.0
A5cutback = 1.0
A5climb1 = 1.0
A5climbn = 1.0
A5descent1 = 1.0
A5descentn = 1.0

pare[ieA7fac, ipstatic, :] .= A7static
pare[ieA7fac, iprotate, :] .= A7takeoff
pare[ieA7fac, iptakeoff, :] .= A7takeoff
pare[ieA7fac, ipcutback, :] .= A7cutback

pare[ieA5fac, ipstatic, :] .= A5static
pare[ieA5fac, iprotate, :] .= A5takeoff
pare[ieA5fac, iptakeoff, :] .= A5takeoff
pare[ieA5fac, ipcutback, :] .= A5cutback

for ip = ipclimb1:ipclimbn

    frac = convert(AbstractFloat, ip - ipclimb1) / convert(AbstractFloat, ipclimbn - ipclimb1)

    pare[ieA7fac, ip, :] .= A7climb1 * (1.0 - frac) + A7climbn * frac
    pare[ieA5fac, ip, :] .= A5climb1 * (1.0 - frac) + A5climbn * frac

end

pare[ieA7fac, ipcruise1:ipcruisen, :] .= 1.0
pare[ieA5fac, ipcruise1:ipcruisen, :] .= 1.0

for ip = ipdescent1:ipdescentn

    frac = convert(AbstractFloat, ip - ipdescent1) / convert(AbstractFloat, ipdescentn - ipdescent1)

    pare[ieA7fac, ip, :] .= A7descent1 * (1.0 - frac) + A7descentn * frac
    pare[ieA5fac, ip, :] .= A5descent1 * (1.0 - frac) + A5descentn * frac

end

pare[ieA7fac, iptest, :] .= A7static
pare[ieA5fac, iptest, :] .= A5static

parg[igGearf] = Gearf
parg[igHTRf] = HTRf
parg[igHTRlc] = HTRlc
parg[igHTRhc] = HTRhc

#------------------------------
#- nacelle drag stuff
parg[igrSnace] = 16.0
parg[igrVnace] = 1.02

#------------------------------
#- engine weight model
pari[iiengwgt] = 1

#------------------------------
# Other
parg[ignfweb] = 1.0

# importing the module
using DelimitedFiles


println("pari", pari)
writedlm("pari_jl.txt", pari)
writedlm("parg_jl.txt", parg)
writedlm("parm_jl.txt", parm)
writedlm("para_jl.txt", para)
writedlm("pare_jl.txt", pare)



println(size(pari))
println(size(parg))
println(size(parm))
println(size(para))
println(size(pare))

