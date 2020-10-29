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

# Wing sizing 
results = surfw(gee,po,b,bs,bo,co,zs,
	lambdat,lambdas,gammat,gammas,
	Nload,iwplan,We,
	Winn,Wout,dyWinn,dyWout,
	sweep,wbox,hboxo,hboxs,rh, fLt,
	tauweb,sigcap,sigstrut,Ecap,Eweb,Gcap,Gweb,
	rhoweb,rhocap,rhostrut,rhofuel)

# [TODO] note this assumes wings have some fuel, so need to ensure that is addressed

#Tail sizing section

# Pitch trim by adjusting Clh or by moving wing

# Drag buildup cdsum()

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

