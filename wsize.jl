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
    include("wingsc.jl")
    include("surfcm.jl")
    include("surfdx.jl")
    include("wingpo.jl")
    include("tailpo.jl")


## Load parameters
# Set specifed quantites during weight iteration loop
# Store design-mission parameters
# Fixed weights and locations -> moments
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

# Use cabin volume to get actual buoyancy weight

# Engine weights

# Update weights


#----------------------
## Wing sizing section
#----------------------
# Initial size of the wing area and chords
S, b, bs, co = wingsc(W, CL, qinf, ηsi, bo, λt, λs)

# x-offset of the wing centroid from wingbox
dxwing, macco = surfdx(b, bs, bo, λt, λs, sweep)

# Calculate wing pitching moment constants
CMw0, CMw1 = surfcm(b, bs, b0, sweep, Xaxis,
                 λt,λs,γt,γs, 
                 AR,fLo,fLt,cmpo,cmps,cmpt)

 #repeat above for climb, cruise and descent

# Wing center load po calculation using cruise spanload cl(y)
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

##Tail sizing section

# Size HT
htsize(pari, parg, para[1, ipdescentn], para[1, ipcruise1], para[1, ipcruise1])
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
balance(pari,parg,para[1,ip],rfuel,rpay,xipay, itrim)


# Drag buildup cdsum()
cdsum(pari, parg, para, pare, 1)

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

