"""
# wsize - Main weight sizing section

Calls on various sub-functions to calculate weight of fuselage, wings, tails etc,

and iterates until the MTOW converges to within a specified tolerance.

### Inputs:
- Array of flags that control design choices - fuel types, where to store fuel etc
- Geometrical and structural parameters - dimensions primarily
- Aerodynamic paramters - CL, CD, KE dissipation, etc...
- Mission specific paramters - alt, mach, P, T etc...
- Engine specific parameters 
"""
function wsize()

    
    include("atmos.jl")
    include("fuseW.jl")
    include("fusebl.jl")
    include("wingsc.jl")
    include("surfcm.jl")
    include("surfdx.jl")
    include("wingpo.jl")
    include("tailpo.jl")


## Load parameters
# Set specifed quantites during weight iteration loop
# Store design-mission parameters

# 

## Initial guess section
# Wing panel weights and moments after estimating span
# HT sizing estimate
# Initial weight guesses
# Estimate wing centroid
# Estimate tail centroid
# Any intial guesses for engine cycle - thrust, Win, fan size etc (based on OEI?)
# Estimate the mission fuel fraction from Breguet range equation
# Estimate initial cruize climb angle
# 

## Weight loop
# Fuselage sizing

fusebl!(pari, parg, para[1,ipcruise1])

# Fixed weights and locations -> moments

# Calculate cruise altitude atmospheric conditions
T,p,ρ,a,μ = atmos(para[iaalt, ipcruise1]/1000.0)

# Calc x-offset of wing centroid from wingbox
surfdx()

# Initial fuel fraction estimate from BRE
LoD  = 18.0
TSFC = 1.0/ 7000.0
V    = 240
ffburn = (1.0 - exp(-Rangetot*TSFC/(V*LoD))) # ffburn = Wfuel/Wstart


(tskin, tcone, tfweb, tfloor, xhbend, xvbend,
EIhshell,EIhbend, EIvshell,EIvbend, GJshell ,GJcone,
Wshell, Wcone, Wwindow, Winsul, Wfloor, Whbend, Wvbend,
Wfuse, xWfuse, cabVol) = fuseW(gee, Nland, Wfix, Wpay, Wpadd, Wseat, Wapu, Weng, 
                      fstring, fframe, ffadd, deltap, 
                      Wpwindow, Wppinsul, Wppfloor, 
                      Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax, 
                      bv, lambdav, nvtail, 
                      Rfuse, dRfuse, wfb, nfweb, lambdac, 
                      xnose, xshell1, xshell2, xconend, 
                      xhtail, xvtail, 
                      xwing, xwbox, cbox, 
                      xfix, xapu, xeng, 
                      hfloor, 
                      sigskin, sigbend,  rhoskin, rhobend, 
                      Eskin, Ebend, Gskin)

parg[igtskin ]  = tskin
parg[igtcone ]  = tcone
parg[igtfweb ]  = tfweb
parg[igtfloor]  = tfloor
parg[igxhbend]  = xhbend
parg[igxvbend]  = xvbend

parg[igEIhshell] = EIhshell
parg[igEIhbend ] = EIhbend 
parg[igEIvshell] = EIvshell
parg[igEIvbend ] = EIvbend 
parg[igGJshell ] = GJshell 
parg[igGJcone  ] = GJcone  

parg[igWshell ]  = Wshell
parg[igWcone  ]  = Wcone
parg[igWwindow]  = Wwindow
parg[igWinsul ]  = Winsul
parg[igWfloor ]  = Wfloor

parg[igWhbend]  = Whbend
parg[igWvbend]  = Wvbend

parg[igWfuse ] = Wfuse
parg[igxWfuse] = xWfuse

parg[igcabVol] = cabVol

# Use cabin volume to get actual buoyancy weight
ρcab = max(parg[igpcabin], pambient)/ RSL*TSL
WbuoyCR = (ρcab - ρ0)*gee*cabVol

# Engine weights

# Update weights


#----------------------
## Wing sizing section
#----------------------

# Size wing area and chords at start-of-cruise
ip = ipcruise1
W = WMTO*para[iafracW, ip]
CL = para[iaCL, ip]
ρ0 = #ambient rho
u0 = velocity
qinf = 0.5*ρ0*u0^2
BW = W + WbuoyCR #Weight including buoyancy

# Initial size of the wing area and chords
S, b, bs, co = wingsc(BW, CL, qinf, ηsi, bo, λt, λs)
parg[[igS, igb, igbs, igco]] = [S, b, bs, co]

cbox = co*wbox #Updating wing box chord for fuseW in next iteration

# x-offset of the wing centroid from wingbox
dxwing, macco = surfdx(b, bs, bo, λt, λs, sweep)
xwing = xwbox + dxwing
cma   = macco * co
para[[igxwing, igcma]] = [xwing, cma]

# Calculate wing pitching moment constants
#------------------------------------------
## Takeoff
ip = iptakeoff
cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
γt, γs = parg[iglambdat]*para[iarclt, ip], parg[iglambdas]*para[iarcls, ip]

CMw0, CMw1 = surfcm(b, bs, b0, sweep, Xaxis,
                    λt,λs,γt,γs, 
                    AR,fLo,fLt,cmpo,cmps,cmpt)

para[iaCMw0, ipstatic:ipclimb1] .= CMw0
para[iaCMw1, ipstatic:ipclimb1] .= CMw1

## Cruise
ip = ipcruise1
cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
γt, γs = parg[iglambdat]*para[iarclt, ip], parg[iglambdas]*para[iarcls, ip]

CMw0, CMw1 = surfcm(b, bs, b0, sweep, Xaxis,
                    λt,λs,γt,γs, 
                    AR,fLo,fLt,cmpo,cmps,cmpt)

para[iaCMw0, ipclimb1+1:ipdescentn-1] .= CMw0
para[iaCMw1, ipclimb1+1:ipdescentn-1] .= CMw1

## Descent
ip = ipdescentn
cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
γt, γs = parg[iglambdat]*para[iarclt, ip], parg[iglambdas]*para[iarcls, ip]

CMw0, CMw1 = surfcm(b, bs, b0, sweep, Xaxis,
                λt,λs,γt,γs, 
                AR,fLo,fLt,cmpo,cmps,cmpt)

para[iaCMw0, ipdescentn] .= CMw0
para[iaCMw1, ipdescentn] .= CMw1
#------------------------------------------

# Wing center load po calculation using cruise spanload cl(y)
ip = ipcruise1
γt, γs = parg[iglambdat]*para[iarclt, ip], parg[iglambdas]*para[iarcls, ip]
Lhtail = WMTO * parg[igCLhNrat]*parg[igSh]/parg[igS]

po = wingpo(b,bs,bo,
        λt,λs,γt,γs,
        AR,N,W,Lhtail,fLo,fLt)


results = surfw(gee,po,b,bs,bo,co,zs,
                lambdat,lambdas,gammat,gammas,
                Nload,iwplan,We,
                Winn,Wout,dyWinn,dyWout,
                sweep,wbox,hboxo,hboxs,rh, fLt,
                tauweb,sigcap,sigstrut,Ecap,Eweb,Gcap,Gweb,
                rhoweb,rhocap,rhostrut,rhofuel)

# [TODO] note this assumes wings have some fuel, so need to ensure that is addressed

# -------------------------------
#      Tail sizing section
# -------------------------------

# Set tail CL derivative 
dϵdα   = parg[igdepsda]
sweeph = parg[igsweeph]
tanL   = tan(sweep *π/180.)
tanLh  = tan(sweeph*π/180.)

ip = ipcruise1
Mach = para[iaMach, ip]
β  = sqrt(1.0 - Mach^2) #Prandtl-Glauert factor √(1-M²)
# Calculate the tail lift-curve slope
dCLhdCL = (β + 2.0/AR)/(β + 2.0/ARh) * sqrt(β^2 + tanL^2)/sqrt(β^2 + tanLh^2) * (1.0 - dϵdα)

parg[igdCLhdCL] = dCLhdCL

# Set Nacelle CL derivative fraction
dCLnda  = parg[igdCLnda]
dCLndCL = dCLnda * (β + 2.0/AR) * sqrt(β^2 + tanL^2)/ (2.0*π*(1.0 + 0.5*hboxo))
parg[dCLndCL] = dCLndCL

# Size HT

#if initial iterations or intiial weight =0 then just get tail volume coeff
lhtail = xhtail - xwing
Vh = parg[igVh]
Sh = Vh*S*cma/lhtail
parg[igSh] = Sh

# for subsequent iterations:
htsize(pari, parg, para[1, ipdescentn], para[1, ipcruise1], para[1, ipcruise1])
xwbox, xwing = parg[igxwbox], parg[igxwing]
lhtail = xhtail - xwing
Sh = parg[igSh]
parg[igVh] = Sh*lhtail/(S*cma)
#endif

# set HT max loading magnitude
bh, coh, poh = tailpo(Sh, ARh, λh, qne, CLhmax)
# set VT max loading magnitude, based on singel tail + its bottom image
bv2, cov, pov = tailpo(2.0*Sv/nvtail, 2.0*ARv,lambdav,qne,CLvmax)

# HT weight
# HT centroid x-offset
# HT pitching moment coeff

# VT weight
# VT centroid x-offset


#calculate for start-of-cruise point
ip = ipcruise1

# Pitch trim by adjusting Clh or by moving wing
Wzero = WMTO - parg[igWfuel] #Zero fuel weight
Wf    = para[iafracW, ip]*WMTO - Wzero
rfuel = Wf/parg[igWfuel]
rpay  = 1.0
ξpay  = 0.
itrim = 1
balance(pari,parg,para[1,ip],rfuel,rpay, ξpay, itrim)

# Drag buildup cdsum()
cdsum(pari, parg, para, pare, 1)
LoD = para[iaCL, ip]/para[iaCD, ip]

# Size engine for TOC

# Size PCEC - estimate weights 

# Engine weight section
#  Drela's weight model? Nate Fitszgerald - geared TF weight model


# Fly mission
# Get mission fuel burn (check if fuel capacity is sufficent)

# Recalculate weight wupdate()

# Set previous iteration weights 

# END weight sizing loop

# BFL calculations/ Noise? / Engine perf 

end

