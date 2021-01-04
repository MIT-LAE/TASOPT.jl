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
function wsize(pari, parg, parm, para, pare,
            itermax, wrlx, initwgt, initeng, iairf)

    #include all functions required [TOCHECK] might not include here but only in MAIN program
    include("atmos.jl")
    include("fuseW.jl")
    include("fusebl.jl")
    include("wingsc.jl")
    include("surfcm.jl")
    include("surfdx.jl")
    include("wingpo.jl")
    include("tailpo.jl")

    # Weight convergence tolerance 
    tolerW = 1.0e-10

    # Atmospheric conditions at sea-level
    TSL, pSL, ρSL, aSL, μSL = atmos(0.0)
## Set quantities that are fixed during weight iteration

    # Unpack payload and range for design mission - this is the mission that the structures are sized for
        Rangetot = parm[imRange]
        Wpay     = parm[imWpay ]
    # Store the design mission in the geometry array as well
        parg[igRange] = Rangetot
        parg[igWpay ] = Wpay

    # Fixed weight and location of fixed weight
        Wfix = parg[igWfix]
        xfix = parg[igxfix]

    # Weight fractions
        fapu  = parg[igfapu ]
        fpadd = parg[igfpadd]
        fseat = parg[igfseat]
        feadd = parg[igfeadd]
        fnace = parg[igfnace]
        fhadd = parg[igfhadd]
        fvadd = parg[igfvadd]
        fwadd = parg[igfflap] + parg[igfslat] + parg[igfaile] + parg[igflete] + parg[igfribs] + parg[igfspoi] + parg[igfwatt]

        fstring = parg[igfstring]
        fframe  = parg[igfframe ]
        ffadd   = parg[igffadd  ]
        fpylon  = parg[igfpylon ]

        fhpesys = parg[igfhpesys]
        flgnose = parg[igflgnose]
        flgmain = parg[igflgmain]

        freserve = parg[igfreserve]

    # fuselage lift carryover loss, tip lift loss fractions
        fLo = parg[igfLo]
        fLt = parg[igfLt]
    
    # fuselage dimensions and coordinates
        Rfuse   = parg[igRfuse ] 
        dRfuse  = parg[igdRfuse]
        wfb     = parg[igwfb  ]
        nfweb   = parg[ignfweb]
        hfloor  = parg[ighfloor]
        xnose   = parg[igxnose]
        xend    = parg[igxend ]
        xshell1 = parg[igxshell1]
        xshell2 = parg[igxshell2]
        xconend = parg[igxconend]
        xwbox   = parg[igxwbox]
        xhbox   = parg[igxhbox]
        xvbox   = parg[igxvbox]
        xapu    = parg[igxapu]
        xeng    = parg[igxeng]

    # calculate payload proportional weights from weight fractions
        Wapu  = Wpay * fapu
        Wpadd = Wpay * fpadd
        Wseat = Wpay * fseat

    # window and insulation densities per length and per area
        Wpwindow = parg[igWpwindow]
        Wppinsul = parg[igWppinsul]
        Wppfloor = parg[igWppfloor]

    # fuselage-bending inertial relief factors
        rMh = parg[igrMh]
        rMv = parg[igrMv]

    # tail CL's at structural sizing cases
        CLhmax = parg[igCLhmax]
        CLvmax = parg[igCLvmax]

    # wing break, wing tip taper ratios
        λs = parg[iglambdas]
        λt = parg[iglambdat]

    # tail surface taper ratios (no inner panel, so λs=1)
        λhs = 1.0
        λh  = parg[iglambdah]
        λvs = 1.0
        λv  = parg[iglambdav]

    # tailcone taper ratio
        λc = parg[iglambdac]

    # wing geometry parameters
        sweep = parg[igsweep]
        wbox  = parg[igwbox ]
        hboxo = parg[ighboxo]
        hboxs = parg[ighboxs]
        rh    = parg[igrh]
        AR    = parg[igAR]
        bo    = parg[igbo]
        ηs    = parg[igetas]
        Xaxis = parg[igXaxis]

    # tail geometry parameters
        sweeph = parg[igsweeph]
        wboxh  = parg[igwboxh]
        hboxh  = parg[ighboxh]
        rhh    = parg[igrhh]
        ARh    = parg[igARh]
        boh    = parg[igboh]

        sweepv = parg[igsweepv]
        wboxv  = parg[igwboxv]
        hboxv  = parg[ighboxv]
        rhv    = parg[igrhv]
        ARv    = parg[igARv]
        bov    = parg[igbov]

    # number of vertical tails
        nvtail = parg[ignvtail]

    # strut vertical base height, h/c, strut shell t/h
        zs     = parg[igzs]
        hstrut = parg[ighstrut]
        tohstrut = 0.05

    # assume no struts on tails
        zsh = 0.0
        zsv = 0.0

    # max g load factors for wing, fuselage
        Nlift = parg[igNlift]
        Nland = parg[igNland]

    # never-exceed dynamic pressure for sizing tail structure
        Vne = parg[igVne]
        qne = 0.5*ρSL*Vne^2

    # wingbox stresses and densities [section 3.1.9 taspot.pdf [prash]] #reference
        σcap    = parg[igsigcap ] * parg[igsigfac]
        tauweb  = parg[igtauweb ] * parg[igsigfac]
        rhoweb  = parg[igrhoweb ]
        rhocap  = parg[igrhocap ]

    # fuselage stresses and densities
        σskin = parg[igsigskin] * parg[igsigfac]
        σbend = parg[igsigbend] * parg[igsigfac]
        rhoskin = parg[igrhoskin]
        rhobend = parg[igrhobend]

    # fuselage shell bending/skin modulus ratio
        rEshell = parg[igrEshell]

    # strut stress and density
        σstrut = parg[igsigstrut] * parg[igsigfac]
        rhostrut = parg[igrhostrut]

    # assume tail stresses and densities are same as wing's (keeps it simpler)
        σcaph   = σcap
        tauwebh = tauweb
        rhowebh = rhoweb
        rhocaph = rhocap

        σcapv   = σcap
        tauwebv = tauweb
        rhowebv = rhoweb
        rhocapv = rhocap

    # number of engines, y-position of outermost engine
        neng = parg[igneng]
        yeng = parg[igyeng]

    # fan hub/tip ratio
        HTRf = parg[igHTRf]

    # nacelle wetted area / fan area ratio
        rSnace = parg[igrSnace]
# -------------------------------------------------------    
## Initial guess section [Section 3.2 of TASOPT docs]
# -------------------------------------------------------
    # Allow first iteration
    if(initwgt == 0)
        Whtail = 0.05 * Wpay/parg[igsigfac]
        Wvtail = Whtail
        Wwing  = 0.5  * Wpay/parg[igsigfac]
        Wstrut = 0.0
        Weng   = 0.3 * Wpay
        feng   = 0.08

        dxWhtail = 0.0
        dxWvtail = 0.0

        # Wing panel weights and moments (after estimating span first)
            ip = ipcruise1
            W = 5.0*Wpay
            S = W / (0.5* pare[ierho0,ip] * pare[ieu0,ip]^2 * para(iaCL,ip))
            b  = sqrt(S*parg[igAR])
            bs = b*ηs
            Winn = 0.15 * Wpay / parg[igsigfac]
            Wout = 0.05 * Wpay / parg[igsigfac]
            dyWinn = Winn*0.30*(0.5*(bs-bo))
            dyWout = Wout*0.25*(0.5*(b -bs))

            parg[igWhtail] = Whtail
            parg[igWvtail] = Wvtail
            parg[igWwing ] = Wwing
            parg[igWstrut] = Wstrut
            parg[igWeng  ] = Weng
            parg[igWinn  ] = Winn
            parg[igWout  ] = Wout
            parg[igdxWhtail] = dxWhtail
            parg[igdxWvtail] = dxWvtail
            parg[igdyWinn] = dyWinn
            parg[igdyWout] = dyWout

        # wing centroid x-offset form wingbox
            dxwing, macco = surfdx(b, bs, bo, λt, λs, sweep)
            xwing = xwbox * dxwing
            parg[igxwing] = xwing

        # tail area centroid locations (assume no offset from sweep initially)
            parg[igxhtail], parg[igxvtail] = xhbox, xvbox
        # center wingbox chord extent for fuse weight calcs (small effect)
            cbox = 0.0
        # nacelle, fan duct, core, cowl lengths ℛℯ calcs
            parg[iglnace] = 0.5*S/b

        # nacelle Awet/S
            fSnace = 0.2
            parg[igfSnace] = fSnace

        # Initial fuel fraction estimate from BRE
            LoD  = 18.0
            TSFC = 1.0/ 7000.0
            V    = pare[ieu0, ipcruise1]
            ffburn = (1.0 - exp(-Rangetot*TSFC/(V*LoD))) # ffburn = Wfuel/Wstart
            ffburn = min(ffburn, 0.8/(1.0 + freserve))
        # mission-point fuel fractions
            ffuelb = ffburn*(1.0  + freserve)  # start of climb
            ffuelc = ffburn*(0.90 + freserve)  # start of cruise
            ffueld = ffburn*(0.02 + freserve)  # start of descent
            ffuele = ffburn*(0.0  + freserve)  # end of descent (landing) 

        # max fuel fraction is at start of climb
            ffuel = ffuelb # ffuel is max fuel fraction

        # Set initial climb γ = 0 to force intial guesses
            para[iagamV, :] .= 0.0

        # Put initial-guess weight fractions in mission-point array. 
        
            # These are points before climb starts
            para[iafracW, ipstatic ] = 1.0
            para[iafracW, iprotate ] = 1.0
            para[iafracW, iptakeoff] = 1.0
            para[iafracW, ipcutback] = 1.0

            # Interpolate for the rest
            #--------- Simple linear interpolation ---------
            #     y-y1   y2-y1
            #     ---- = ----- => y = y1*(1-frac) + y2*frac
            #     x-x1   x2-x1
            #-----------------------------------------------
            # Climb
            for ip = ipclimb1:ipclimbn
                frac = float(ip - ipclimb1)/float(ipclimbn - ipclimb1)
                ffp  = ffuelb*(1.0 - frac) + ffuelc*frac
                para[iafracW, ip] = 1.0 - ffuel + ffp
            end
            # Cruise
            for ip = ipcruise1:ipcruisen
                frac = float(ip - ipcruise1)/float(ipcruisen - ipcruise1)
                ffp  = ffuelc*(1.0 - frac) + ffueld*frac
                para[iafracW, ip] = 1.0 - ffuel + ffp
            end
            # Descent
            for ip = ipdescent1:ipdescentn
                frac = float(ip - ipdescent1)/float(ipdescentn - ipdescent1)
                ffp  = ffueld*(1.0 - frac) + ffuele*frac
                para[iafracW, ip] = 1.0 - ffuel + ffp
            end
        
        # Initial tail info for sizing of fuselage bending and torsion added material 
            Sh = (2.0*Wpay) / (qne*CLhmax)
            Sv = (2.0*Wpay) / (qne*CLvmax)
            bv = sqrt(Sv*ARv)
    
            parg[igSh] = Sh
            parg[igSv] = Sv

        # Initial wing and tail pitching moments (including sweep)
            para[iaCMw0, :] .= 0.0
            para[iaCMw1, :] .= 0.0
            para[iaCMh0, :] .= 0.0
            para[iaCMh1, :] .= 0.0
            para[iaCLh , :] .= 0.0
        
        # Initial cruise-climb angle gamVcr needed to estimate end-of-cruise altitude
            gamVcr = 0.0002
            para[iaCD, ipcruise1] = para[iaCL, ipcruise1]/LoD
            para[iagamV, ipcruise1] = gamVcr
        
        # Pressure and altitude at start of cruise
            
        # Guess for OEI [TODO] This needs some thinking about what is "One engine out" mean for a turbo-electric aircraft

        # Guess fan face mach numbers for nacelle CD calcs

    else #Second iteration onwards use previously calcualted values

        # Wing parameters
            S  = parg[igS]
            b  = parg[igb]
            bs = parg[igbs]
            bo = parg[igbo]
            cbox = parg[igco]*parg[igwbox]

            xwing  = parg[igxwing]
            dxwing = parg[igxwing] - parg[igxwbox]

        # Tail parameters
            bh = parg[igbh]
            bv = parg[igbv]

            coh = parg[igcoh]
            cov = parg[igcov]

            Sh  = parg[igSh]
            Sv  = parg[igSv]
            ARh = parg[igARh]
            ARv = parg[igARv]

            xhtail = parg[igxhtail]
            xvtail = parg[igxvtail]
    
            xhbox = parg[igxhbox]
            xvbox = parg[igxvbox]

            dxWhtail = parg[igdxWhtail]
            dxWvtail = parg[igdxWvtail]

        # Weights
            Whtail = parg[igWhtail]
            Wvtail = parg[igWvtail]
            Wwing  = parg[igWwing ]
            Wstrut = parg[igWstrut]
            Weng   = parg[igWeng  ]
            Winn   = parg[igWinn  ]
            Wout   = parg[igWout  ]

        # Others       
            dyWinn = parg[igdyWinn]
            dyWout = parg[igdyWout]    

            WMTO  = parg[igWMTO]
            feng  = parg[igWeng]/WMTO
            ffuel = parg[igWfuel]/WMTO

            fSnace = parg[igfSnace]
    endif
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

    # Calculate fuselage B.L. development at start of cruise: ipcruise1
    fusebl!(pari, parg, para[1,ipcruise1])
    KAfTE   = para[iaKAfTE  , ipcruise1] # Kinetic energy area at T.E.
    DAfsurf = para[iaDAfsurf, ipcruise1] # Surface dissapation area 
    DAfwake = para[iaDAfwake, ipcruise1] # Wake dissapation area
    PAfinf  = para[iaPAfinf , ipcruise1] # Momentum area at ∞

    # Assume K.E., Disspation and momentum areas are const. for all mission points:
    para[iaKAfTE  , :] .= KAfTE
    para[iaDAfsurf, :] .= DAfsurf
    para[iaDAfwake, :] .= DAfwake
    para[iaPAfinf , :] .= PAfinf

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
                      bv, λv, nvtail, 
                      Rfuse, dRfuse, wfb, nfweb, λc, 
                      xnose, xshell1, xshell2, xconend, 
                      xhtail, xvtail, 
                      xwing, xwbox, cbox, 
                      xfix, xapu, xeng, 
                      hfloor, 
                      σskin, σbend,  rhoskin, rhobend, 
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
                λt,λs,gammat,gammas,
                Nload,iwplan,We,
                Winn,Wout,dyWinn,dyWout,
                sweep,wbox,hboxo,hboxs,rh, fLt,
                tauweb,σcap,σstrut,Ecap,Eweb,Gcap,Gweb,
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
bv2, cov, pov = tailpo(2.0*Sv/nvtail, 2.0*ARv,λv,qne,CLvmax)

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

