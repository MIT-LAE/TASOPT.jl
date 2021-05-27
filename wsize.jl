using Printf
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
            itermax, wrlx1, wrlx2, wrlx3,
            initwgt, initeng, iairf, Ldebug, printiter)

time_propsys = 0.0
    # Weight convergence tolerance 
        # tolerW = 1.0e-10
        tolerW = 1.0e-6
        errw   = 1.0
    # Initialze some variables
    fsum = 0.0
    ifirst = true
    NPSS_TS = Base.Process
    NPSS_Fan = Base.Process
    NPSS_AftFan = Base.Process

    # Flags
    ifuel   = pari[iifuel  ] # Fuel type 24 == kerosene TODO need to update this for LH2
    ifwcen  = pari[iifwcen ]
    iwplan  = pari[iiwplan ]
    iengloc = pari[iiengloc]
    iengwgt = pari[iiengwgt]
    iBLIc   = pari[iiBLIc  ]
    ifclose = pari[iifclose]
    iHTsize = pari[iiHTsize]
    iVTsize = pari[iiVTsize]
    ixwmove = pari[iixwmove]

    # Unpack number of powertrain elements
       nfan    = parpt[ipt_nfan]
       ngen    = parpt[ipt_ngen]
       nTshaft = parpt[ipt_nTshaft]

    # Atmospheric conditions at sea-level
    TSL, pSL, ρSL, aSL, μSL = atmos(0.0)

    # Calculate fuselage B.L. development at start of cruise: ipcruise1
    time_fusebl = @elapsed fusebl!(pari, parg, para, ipcruise1)
    # println("Fuse bl time = $time_fusebl")
    KAfTE   = para[iaKAfTE  , ipcruise1] # Kinetic energy area at T.E.
    DAfsurf = para[iaDAfsurf, ipcruise1] # Surface dissapation area 
    DAfwake = para[iaDAfwake, ipcruise1] # Wake dissapation area
    PAfinf  = para[iaPAfinf , ipcruise1] # Momentum area at ∞

    # Assume K.E., Disspation and momentum areas are const. for all mission points:
    para[iaKAfTE  , :] .= KAfTE
    para[iaDAfsurf, :] .= DAfsurf
    para[iaDAfwake, :] .= DAfwake
    para[iaPAfinf , :] .= PAfinf
    
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
        fwadd = parg[igfflap] + parg[igfslat] + 
                parg[igfaile] + parg[igflete] + parg[igfribs] + parg[igfspoi] + parg[igfwatt]

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

    # Set atmos conditions
        ip = ipcruise1
        altkm = para[iaalt, ip]/1000.0
        T0, p0, ρ0, a0, μ0 = atmos(altkm)
        Mach  = para[iaMach, ip]
        pare[iep0  , ip] = p0 
        pare[ieT0  , ip] = T0
        pare[ierho0, ip] = ρ0
        pare[iemu0 , ip] = μ0
        pare[iea0  , ip] = a0

        pare[ieM0, ip] = Mach
        pare[ieu0, ip] = Mach*a0
        para[iaReunit, ip] = Mach*a0 *ρ0/μ0

        #Rotation condition
        ip = iprotate
        altkm = para[iaalt, ip]/1000.0
        T0, p0, ρ0, a0, μ0 = atmos(altkm)
        Mach  = para[iaMach, ip]
        pare[iep0  , ip] = p0 
        pare[ieT0  , ip] = T0
        pare[ierho0, ip] = ρ0
        pare[iemu0 , ip] = μ0
        pare[iea0  , ip] = a0

        pare[ieM0, ip] = Mach
        pare[ieu0, ip] = Mach*a0
        para[iaReunit, ip] = Mach*a0 *ρ0/μ0


# -------------------------------------------------------    
## Initial guess section [Section 3.2 of TASOPT docs]
# -------------------------------------------------------
    # Allow first iteration
    if(initwgt == 0)
        Whtail = 0.05 * Wpay/parg[igsigfac]
        Wvtail = Whtail
        Wwing  = 0.5  * Wpay/parg[igsigfac]
        Wstrut = 0.0
        Weng   = 0.0 * Wpay
        feng   = 0.0

        # Wfan    = Weng*0.2
        # Wmot    = Weng*0.3
        # Winv    = Weng*0.2
        # Wgen    = Weng*0.5
        # Wtshaft = Weng*0.5

        dxWhtail = 0.0
        dxWvtail = 0.0

        # Wing panel weights and moments (after estimating span first)
            ip = ipcruise1
            W = 5.0*Wpay
            S = W / (0.5* pare[ierho0,ip] * pare[ieu0,ip]^2 * para[iaCL,ip])
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

        # Turbo-electric weights
            parg[igWtesys ] = 0.32*Wpay
            parg[igWftank ] = 0.0

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
            ffburn = (1.0 - exp(-Rangetot*TSFC/(V*LoD))) # ffburn = Wfuel/WMTO
            ffburn = min(ffburn, 0.8/(1.0 + freserve))   # 0.8 is the fuel useablilty? 
        # mission-point fuel fractions ⇾ ffuel = Wfuel/WMTO
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
            @inbounds for ip = ipclimb1:ipclimbn
                frac = float(ip - ipclimb1)/float(ipclimbn - ipclimb1)
                ffp  = ffuelb*(1.0 - frac) + ffuelc*frac
                para[iafracW, ip] = 1.0 - ffuel + ffp
            end
            # Cruise
            @inbounds for ip = ipcruise1:ipcruisen
                frac = float(ip - ipcruise1)/float(ipcruisen - ipcruise1)
                ffp  = ffuelc*(1.0 - frac) + ffueld*frac
                para[iafracW, ip] = 1.0 - ffuel + ffp
            end
            # Descent
            @inbounds for  ip = ipdescent1:ipdescentn
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
            Mach = para[iaMach, ipcruise1]
            p0c  = pare[iep0  , ipcruise1]
            altc = para[iaalt , ipcruise1]
            # Guess pressure at end-of-cruise (scales with weight)
            p0d  = p0c * (1.0-ffuel+ffueld)/(1.0-ffuel+ffuelc)
            pare[iep0, ipcruisen] = p0d
            
        # Guess for OEI [TODO] This needs some thinking about what is "One engine out" mean for a turbo-electric aircraft
            # pare[ieFe, iprotate] = 2.0*Wpay/neng
            pare[ieFe, iprotate] = 2.0*Wpay # ieFe now stores total thrust
            pare[ieu0, iprotate] = 70.0
            Afan = 3.0e-5 *Wpay/neng
            parg[igdfan] = sqrt(Afan*4.0/π)

        # Guess fan face mach numbers for nacelle CD calcs
            M2des = pare[ieM2, ipcruise1] = 0.6
            pare[ieM2, ipstatic: ipcruisen   ] .= M2des
            pare[ieM2, ipdescent1: ipdescentn] .= 0.8*M2des

    else #Second iteration onwards use previously calculated values

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

        # Turbo-electric system
            Wtesys = parg[igWtesys]
            Wftank = parg[igWftank]
    
    end

# Initialize previous weight iterations
    WMTO1, WMTO2, WMTO3 = zeros(Float64, 3) #1st-previous to 3rd previous iteration weight for convergence criterion

Lconv = false # no convergence yet

# -------------------------------------------------------    
#                   Weight loop
# -------------------------------------------------------    

    @inbounds for  iterw = 1:itermax
        if iterw == itermax
            println("Reached max iterations in weight sizing loop!")
        end
        if(initwgt == 0)
            #Current weight iteration started from an initial guess so be cautious
            itrlx = 5
        else
            #Current weight started from previously converged solution
            itrlx = 2
        end

        if(iterw <= itrlx)
            # under-relax first n itrlx iterations
            rlx = wrlx1
        elseif(iterw >= 3/4*itermax)
            # under-relax after 3/4th of max iterations
            rlx = wrlx3
        else
            # default is no under-relaxation for weight update
            rlx = wrlx2
        end

        # Fuselage sizing

            # Max tail lifts at maneuver qne
                Lhmax = qne*Sh*CLhmax
                Lvmax = qne*Sv*CLvmax/nvtail
            # Max Δp (fuselage pressure) at end of cruise-climb, assumes p ~ W/S
                wcd = para[iafracW, ipcruisen]/ para[iafracW, ipcruise1]
                Δp = parg[igpcabin] - pare[iep0, ipcruise1]*wcd
                parg[igdeltap] = Δp
            
            # Engine weight mounted on tailcone, if any
                if (iengloc==1)
                    Wengtail = 0.0 
                else
                    Wengtail = (parg[igWtshaft] + parg[igWcat])*nTshaft +
                                parg[igWgen]*ngen + parg[igWftank]
                end
            
            Whtail = parg[igWhtail]
            Wvtail = parg[igWvtail]
            xhtail = parg[igxhtail]
            xvtail = parg[igxvtail]
            xwbox  = parg[igxwbox ]    
            xwing  = parg[igxwing ]

            Wtesys = parg[igWtesys]
            Wftank = parg[igWftank]

            # Call fusew
                Eskin = parg[igEcap]
                Ebend = Eskin * rEshell
                Gskin = Eskin * 0.5/(1.0 +0.3)

                (tskin, tcone, tfweb, tfloor, xhbend, xvbend,
                EIhshell,EIhbend, EIvshell,EIvbend, GJshell ,GJcone,
                Wshell, Wcone, Wwindow, Winsul, Wfloor, Whbend, Wvbend,
                Wfuse, xWfuse, cabVol) = fusew(gee, Nland, Wfix, Wpay, Wpadd, Wseat, Wapu, Wengtail, 
                                                fstring, fframe, ffadd, Δp, 
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
                ρcab = max(parg[igpcabin], pare[iep0, ipcruise1])/ (RSL*TSL)
                WbuoyCR = (ρcab - pare[ierho0, ipcruise1])*gee*cabVol
            # Total max Takeoff weight (MTOW)
                
                # WMTO = Wpay + Wfuse + Wwing + Wstrut + Wtesys + Wftank
                #        Whtail + Wvtail +
                #        Weng + Wfuel + 
                #        Whpesys + Wlgnose + Wlgmain
                if (iterw == 1 && initwgt == 0)
                    feng = 0.0 # Set feng to be zero since we are not using the TFan but a TE system

                    # To allow performing aerodynamic and weight-burn calculations on the first iteration, 
                    # an interim MTOW is computed: 
                    fsum = feng + ffuel + fhpesys + flgnose + flgmain
                    # WMTO = (Wpay + Wfuse + Wwing + Wstrut + Whtail + Wvtail)/(1.0 - fsum)
                    WMTO = (Wpay + Wfuse + Wwing + Wstrut + Whtail + Wvtail +
                            Wtesys + Wftank)/(1.0 - fsum)

                    Weng, Wfuel, Whpesys, Wlgnose, Wlgmain = WMTO .* [feng, ffuel, fhpesys, flgnose, flgmain] 
                    parg[igWMTO] = WMTO
                    parg[igWeng] = Weng
                    parg[igWfuel]= Wfuel
                    println("Wfuel initial = ",(ffuel*WMTO))

                else 
                    # Call a better Wupdate function
                    Wupdate0!(parg, rlx, fsum)
                    if(fsum >= 1.0) 
                        println("Something is wrong!! fsum ≥ 1.0")
                        break
                    end

                    parm[imWTO]   = parg[igWMTO]
                    parm[imWfuel] = parg[igWfuel]

                end
                # this calculated WMTO is the design-mission WTO
                parm[imWTO] = parg[igWMTO]
            # Convergence tests
                WMTO = parg[igWMTO]
                errw1 = (WMTO - WMTO1)/WMTO
                errw2 = (WMTO - WMTO2)/WMTO
                errw3 = (WMTO - WMTO3)/WMTO

                errw = max(abs(errw1), abs(errw2), abs(errw3))

            # Print weight/ convergnce started
            if(printiter && iterw == 1)
            @printf("%5s  %11s  %11s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s \n",
            "iterw", "errW", "errW1", "WMTO", "Wfuel", "Wftank", "Wtesys", "Wgen", "Wtshaft", "Wwing", "span", "area", "HTarea", "xwbox" )
            end
           if printiter
            @printf("%5d  %+9.4e  %+9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e\n",
             iterw, errw, errw1, parm[imWTO], parg[igWfuel], parg[igWftank],
             parg[igWtesys], parg[igWgen], parg[igWtshaft],
             parg[igWwing], parg[igb], parg[igS], 
             parg[igSh], parg[igxwbox])
           end
            if(errw <= tolerW)
                Lconv = true
                break
            end

        #--------------------------------
        ##  Wing sizing section
        #--------------------------------
            WMTO = parg[igWMTO]

            # Size wing area and chords at start-of-cruise
            ip = ipcruise1
            W = WMTO*para[iafracW, ip]
            CL = para[iaCL, ip]
            ρ0 = pare[ierho0, ip]
            u0 = pare[ieu0  , ip]
            qinf = 0.5*ρ0*u0^2
            BW = W + WbuoyCR #Weight including buoyancy

            # Initial size of the wing area and chords
            S, b, bs, co = wingsc(BW, CL, qinf, AR, ηs, bo, λt, λs)
            parg[[igS, igb, igbs, igco]] = [S, b, bs, co]

            #Updating wing box chord for fuseW in next iteration
            cbox = co*wbox 

            # x-offset of the wing centroid from wingbox
            dxwing, macco = surfdx(b, bs, bo, λt, λs, sweep)
            xwing = xwbox + dxwing
            cma   = macco * co
            parg[igxwing] = xwing
            parg[igcma]   = cma

            # Calculate wing pitching moment constants
            #------------------------------------------
            ## Takeoff
            ip = iptakeoff
            cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
            γt = parg[iglambdat]*para[iarclt, ip]
            γs = parg[iglambdas]*para[iarcls, ip]

            CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                                λt,λs,γt,γs, 
                                AR,fLo,fLt,cmpo,cmps,cmpt)

            para[iaCMw0, ipstatic:ipclimb1] .= CMw0
            para[iaCMw1, ipstatic:ipclimb1] .= CMw1

            ## Cruise
            ip = ipcruise1
            cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
            γt = parg[iglambdat]*para[iarclt, ip]
            γs = parg[iglambdas]*para[iarcls, ip]

            CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                                λt,λs,γt,γs, 
                                AR,fLo,fLt,cmpo,cmps,cmpt)

            para[iaCMw0, ipclimb1+1:ipdescentn-1] .= CMw0
            para[iaCMw1, ipclimb1+1:ipdescentn-1] .= CMw1

            ## Descent
            ip = ipdescentn
            cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
            γt = parg[iglambdat]*para[iarclt, ip]
            γs = parg[iglambdas]*para[iarcls, ip]
            
            CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                            λt,λs,γt,γs, 
                            AR,fLo,fLt,cmpo,cmps,cmpt)

            para[iaCMw0, ipdescentn] = CMw0
            para[iaCMw1, ipdescentn] = CMw1
            #------------------------------------------

            # Wing center load po calculation using cruise spanload cl(y)
            # -----------------------------------------
            ip = ipcruise1
            γt, γs = parg[iglambdat]*para[iarclt, ip], parg[iglambdas]*para[iarcls, ip] # Lift "taper ratios"
            Lhtail = WMTO * parg[igCLhNrat]*parg[igSh]/parg[igS]

            po = wingpo(b,bs,bo,
                        λt, λs, γt, γs,
                        AR, Nlift, WMTO, Lhtail, fLo, fLt)
         
            # Wing structure
            # -----------------------------------------
            
            # Fans, mot, inv centre of mass:
            # -----
            dy = 1.0 # space to leave near wing root and tip [m]
            yi = LinRange(bo/2 + dy , b/2 - dy, Int(parg[igneng]/2))
            ηi = yi/(b/2)
            ηo = bo/b
            ci = zero(yi)
            for (i, η)  in enumerate(ηi)
                if η <=ηs
                    ci[i] = co*(1  + (λs -  1)*(η - ηo)/(ηs - ηo))
                else
                    ci[i] = co*(λs + (λt - λs)*(η - ηs)/(1  - ηs))
                end
            end

            tanL = tan(parg[igsweep]*π/180.0)
            parg[igxfan] = mean(tanL * (yi .- bo/2) - 0.4ci) + parg[igxwbox] - 1.0
            parg[igxmot] = parg[igxfan] + 0.5
            # -----
            
            if (iwplan == 1)
                Weng1 = parg[igWfan] + parg[igWmot] + parg[igWinv]
            else
                Weng1 = 0.0
            end

            Winn = parg[igWinn]
            Wout = parg[igWout]
            dyWinn  = parg[igdyWinn]
            dyWout  = parg[igdyWout]
            if (pari[iifwing] == 0) 
                rhofuel = 0.0 # tell surfw that there is no fuel in wings
            else
                rhofuel = parg[igrhofuel]
            end
            
            Ecap = parg[igEcap]
            Eweb = Ecap
            Gcap = Ecap*0.5/(1.0+0.3)
            Gweb = Ecap*0.5/(1.0+0.3)


            Ss,Ms,tbwebs,tbcaps,EIcs,EIns,GJs,
            So,Mo,tbwebo,tbcapo,EIco,EIno,GJo,
            Astrut,lstrutp,cosLs,
            Wscen,Wsinn,Wsout,dxWsinn,dxWsout,dyWsinn,dyWsout,
            Wfcen,Wfinn,Wfout,dxWfinn,dxWfout,dyWfinn,dyWfout,
            Wweb,  Wcap,  Wstrut,
            dxWweb,dxWcap,dxWstrut = surfw(gee,po,b,bs,bo,co,zs,
                                            λt,λs,γt, γs,
                                            Nlift,iwplan,Weng1,
                                            Winn,Wout,dyWinn,dyWout,
                                            sweep,wbox,hboxo,hboxs,rh, fLt,
                                            tauweb,σcap,σstrut,Ecap,Eweb,Gcap,Gweb,
                                            rhoweb,rhocap,rhostrut,rhofuel)
            # println([Ss,Ms,tbwebs,tbcaps,EIcs,EIns,GJs,
            # So,Mo,tbwebo,tbcapo,EIco,EIno,GJo,
            # Astrut,lstrutp,cosLs,
            # Wscen,Wsinn,Wsout,dxWsinn,dxWsout,dyWsinn,dyWsout,
            # Wfcen,Wfinn,Wfout,dxWfinn,dxWfout,dyWfinn,dyWfout,
            # Wweb,  Wcap,  Wstrut,
            # dxWweb,dxWcap,dxWstrut])
            # Wsinn is the in-board skin weight, Wsout is the outboard skin weight. 
            # Multiply by 2 to account for the two wing-halves: 
            Wwing   = 2.0 * (Wscen + Wsinn +   Wsout) * (1.0 + fwadd)
            dxWwing = 2.0 * (      dxWsinn + dxWsout) * (1.0 + fwadd)

            # Note this assumes wings have some fuel, so additional check is performed to see if iifwing is 1
            #Calculate volume limited fuel weight depending on if wing center box has fuel or not
            Wfmax = 0.0
            dxWfmax = 0.0
            rfmax = 0.0
            if (pari[iifwing] == 1) # if fuel is stored in wings only then do this
                if (pari[iifwcen] == 0)
                    Wfmax   = 2.0*(         Wfinn +   Wfout)
                    dxWfmax = 2.0*(       dxWfinn + dxWfout)
                else
                    Wfmax   = 2.0*(Wfcen +  Wfinn +   Wfout)
                    dxWfmax = 2.0*(       dxWfinn + dxWfout)
                end
                rfmax = parg[igWfuel]/Wfmax
            end

            
            # Save wing details into geometry array
            parg[igWwing] = Wwing * rlx + parg[igWwing]*(1.0-rlx)
            parg[igWfmax] = Wfmax
            parg[igdxWwing] = dxWwing
            parg[igdxWfuel] = dxWfmax*rfmax
        
            parg[igtbwebs] = tbwebs
            parg[igtbcaps] = tbcaps
            parg[igtbwebo] = tbwebo
            parg[igtbcapo] = tbcapo
            parg[igAstrut] = Astrut
            parg[igcosLs ] = cosLs 
            parg[igWweb  ] = Wweb  
            parg[igWcap  ] = Wcap  
            parg[igWstrut] = Wstrut
            parg[igSomax ] = So
            parg[igMomax ] = Mo
            parg[igSsmax ] = Ss
            parg[igMsmax ] = Ms
            parg[igEIco] = EIco
            parg[igEIcs] = EIcs
            parg[igEIno] = EIno
            parg[igEIns] = EIns
            parg[igGJo]  = GJo
            parg[igGJs]  = GJs
        
            parg[igWstrut] = Wstrut
            parg[igdxWstrut] = dxWstrut

            # Strut chord (perpendicular to strut)
            cstrut = sqrt(0.5*Astrut/(tohstrut*hstrut))
            Ssturt = 2.0*cstrut*lstrutp
            parg[igcstrut] = cstrut
            parg[igSstrut] = Ssturt

            # Individual panel weights
            Winn = Wsinn*(1.0 + fwadd) + rfmax*Wfinn
            Wout = Wsout*(1.0 + fwadd) + rfmax*Wfout

            dyWinn = dyWsinn*(1.0 + fwadd) + rfmax*dyWfinn
            dyWout = dyWsout*(1.0 + fwadd) + rfmax*dyWfout
      
            parg[igWinn] = Winn
            parg[igWout] = Wout
            parg[igdyWinn] = dyWinn
            parg[igdyWout] = dyWout

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
                parg[igdCLndCL] = dCLndCL

            # Size HT
                if(iterw<=2 && initwgt == 0)
                    #if initial iterations or intiial weight =0 then just get tail volume coeff Vh
                    lhtail = xhtail - xwing
                    Vh = parg[igVh]
                    Sh = Vh*S*cma/lhtail
                    parg[igSh] = Sh
                else
                    # for subsequent iterations:
                    htsize(pari, parg, view(para, :, ipdescentn), view(para,:, ipcruise1), view(para,:, ipcruise1))

                    xwbox, xwing = parg[igxwbox], parg[igxwing]

                    lhtail = xhtail - xwing
                    Sh = parg[igSh]

                    parg[igVh] = Sh*lhtail/(S*cma)
                end
                
            # Vertical tail sizing 
                
                # Estimate thrust at take-off and the resultant moment if OEI - to estimate Vertical tail size
                #   Same logic should hold if outboard electric motors fail
                # [Section 2.13.2 of TASOPT docs]

                    ip = iprotate
                    qstall = 0.5 * pare[ierho0, ip] *(pare[ieu0, ip]/1.2)^2
                    CDAe = parg[igcdefan] * 0.25π *parg[igdfan]^2
                    De = qstall*CDAe
                    Fe = pare[ieFe, ip]
                    Me = (Fe + De)*yeng

                #
                if(iVTsize == 1)
                    lvtail = xvtail - xwing
                    Vv = parg[igVv]
                    Sv = Vv*S*b/lvtail
                    parg[igSv] = Sv
                    parg[igCLveout] = Me/(qstall*Sv*lvtail) # Max lift coeff oc vertical tail with some yaw control when OEI [Eqn. 312 of TASOPT docs]
                else
                    lvtail = xvtail - xwing
                    CLveout = parg[igCLveout]
                    Sv = Me/(qstall * CLveout *lvtail)
                    parg[igSv] = Sv
                    parg[igVv] = Sv*lvtail/(S*b)
                end

            # set HT max loading magnitude
            bh, coh, poh = tailpo(Sh, ARh, λh, qne, CLhmax)
            parg[igbh ] = bh
            parg[igcoh] = coh

            # set VT max loading magnitude, based on single tail + its bottom image
            bv2, cov, pov = tailpo(2.0*Sv/nvtail, 2.0*ARv,λv,qne,CLvmax)
            parg[igbv ] = bv2/2
            parg[igcov] = cov

            # HT weight
            γh  = λh
            γhs = λhs

            ihplan = 0
            Wengh = 0.0
            Ecap = parg[igEcap]
            Eweb = Ecap
            Gcap = Ecap*0.5/(1.0+0.3)
            Gweb = Ecap*0.5/(1.0+0.3)
            
            #[TODO] Use multiple dispatch here to create a new method for surfw
            Ssh,Msh,tbwebsh,tbcapsh,EIcsh,EInsh,GJsh,
            Soh,Moh,tbweboh,tbcapoh,EIcoh,EInoh,GJoh,
            _, _, _,
            Wscenh,Wsinnh,Wsouth,dxWsinnh,dxWsouth,dyWsinnh,dyWsouth,
            Wfcenh,Wfinnh,Wfouth,dxWfinnh,dxWfouth,dyWfinnh,dyWfouth,
            Wwebh,  Wcaph,  Wstruth,
            dxWwebh,dxWcaph,dxWstruth = surfw(gee,poh,bh,boh,boh,coh,zsh,
                                            λh,λhs,γh,γhs,
                                            1, ihplan, Wengh,
                                            0.0, 0.0, 0.0, 0.0,
                                            sweeph,wboxh,hboxh,hboxh,rhh, fLt,
                                            tauwebh,σcaph,σstrut,Ecap,Eweb,Gcap,Gweb,
                                            rhowebh, rhocaph, rhostrut, rhofuel)

              Whtail = 2.0 * (Wscenh + Wsinnh +   Wsouth)*(1.0 + fhadd)
            dxWhtail = 2.0 * (       dxWsinnh + dxWsouth)*(1.0 + fhadd)
            parg[igWhtail] = Whtail
            parg[igdxWhtail] = dxWhtail
            
            parg[igtbwebh] = tbweboh
            parg[igtbcaph] = tbcapoh
            parg[igEIch] = EIcoh
            parg[igEInh] = EInoh
            parg[igGJh]  = GJoh

            # HT centroid x-offset
            dxh, macco = surfdx(bh, boh, boh, λh, λhs, sweeph)
            parg[igxhtail] = xhbox + dxh

            # HT pitching moment coeff
            fLoh = 0.
            fLth = fLt
            cmph = 0.

            CMh0, CMh1 = surfcm(bh, boh, boh, sweeph, Xaxis, λh, 1.0, λh, 1.0,
                                ARh, fLoh, fLth, 0.0, 0.0, 0.0)
            
            para[iaCMh0, :] .= CMh0
            para[iaCMh1, :] .= CMh1

            # VT weight

            γv  = λv
            γvs = λvs
            ivplan = 0
            Wengv = 0.0
            Ecap = parg[igEcap]
            Eweb = Ecap
            Gcap = Ecap*0.5/(1.0+0.3)
            Gweb = Ecap*0.5/(1.0+0.3)

            Ssv,Msv,tbwebsv,tbcapsv,EIcsv,EInsv,GJsv,
            Sov,Mov,tbwebov,tbcapov,EIcov,EInov,GJov,
            _, _, _,
            Wscenv,Wsinnv,Wsoutv,dxWsinnv,dxWsoutv,dyWsinnv,dyWsoutv,
            Wfcenv,Wfinnv,Wfoutv,dxWfinnv,dxWfoutv,dyWfinnv,dyWfoutv,
            Wwebv2,  Wcapv2,  Wstrutv,
            dxWwebv2,dxWcapv2,dxWstrutv = surfw(gee,pov, bv2, bov, bov,cov,zsv,
                                                λv, λvs, γv, γvs,
                                                1.0,ivplan,Wengv,
                                                0.0, 0.0, 0.0, 0.0,
                                                sweepv,wboxv,hboxv,hboxv,rhv, fLt,
                                                tauwebv,σcapv,σstrut,Ecap,Eweb,Gcap,Gweb,
                                                rhowebv,rhocapv,rhostrut,rhofuel)

            Wvtail   = (Wscenv +  Wsinnv +   Wsoutv)*(1.0 + fvadd) * nvtail
            dxWvtail = (        dxWsinnv + dxWsoutv)*(1.0 + fvadd) * nvtail
            parg[igWvtail] = Wvtail
            parg[igdxWvtail] = dxWvtail

            parg[igtbwebv] = tbwebov
            parg[igtbcapv] = tbcapov
            parg[igEIcv] = EIcov
            parg[igEInv] = EInov
            parg[igGJv]  = GJov

            # VT centroid x-offset
            dxv, _ = surfdx(bv2, bov, bov, λv, λvs, sweepv)
            parg[igxvtail] = xvbox + dxv

        # -----------------------------
        # Drag and engine calculations
        # ------------------------------
            # Total Drag

                WMTO = parg[igWMTO]
                #calculate for start-of-cruise point
                ip = ipcruise1

                # Pitch trim by adjusting Clh or by moving wing
                Wzero = WMTO - parg[igWfuel] #Zero fuel weight
                Wf    = para[iafracW, ip]*WMTO - Wzero
                rfuel = Wf/parg[igWfuel]
                rpay  = 1.0
                ξpay  = 0.
                itrim = 1
                balance(pari,parg,view(para, :,ip),rfuel,rpay, ξpay, itrim)
                # Set N.P. at cruise
                parg[igxNP] = para[iaxNP, ip]

                # Drag buildup cdsum()
                cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), 1)

            # L/D and Design point thrust
                LoD = para[iaCL, ip]/para[iaCD, ip]
                gamV = para[iagamV, ip]
                W   = para[iafracW, ip] * WMTO
                BW  = W + WbuoyCR
                Fdes = BW*(1/LoD + gamV)

                pare[ieFe, ip] = Fdes # Let ieFe store total thrust 
                # println("Cruise Ftotal, des = ", Fdes)

            # Size engine for TOC
            ρAir = pare[ierho0, ipcruise1]
            μAir = pare[iemu0 , ipcruise1]

        if (iterw==1)
            NPSS_Fan = startNPSS("NPSS_Turboshaft/", "Fan.bat")
            NPSS_AftFan = startNPSS("NPSS_Turboshaft/", "Fan.bat")
            NPSS_TS  = startNPSS("NPSS_Turboshaft/", "TP.bat" )
        end
            ρ0 = pare[ierho0, ipcruise1]
            u0 = pare[ieu0  , ipcruise1] 
            fBLIf = parg[igfBLIf]

            Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
            Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines

           time_propsys += @elapsed  ηpt, Ppt, Hpt, heatexcess, mpt, SPpt,
            mdotf_tot, BSFC,
            deNOx, EINOx1, EINOx2,
            FAR, Tt3, OPR, Wc3, FanNozArea, Snace1, AftSnace1 =  PowerTrain(NPSS_TS, NPSS_Fan, NPSS_AftFan, para[iaalt, ipcruise1], para[iaMach, ipcruise1], Fdes,
                                        Kinl, Φinl, parg, parpt, parmot, pargen, ifirst)
            ifirst = false

            pare[iedeNOx , ip] = deNOx
            pare[ieEINOx1, ip] = EINOx1
            pare[ieEINOx2, ip] = EINOx2
            
            pare[ieOPR, ip] = OPR
            pare[ieTt3, ip] = Tt3
            pare[ieWc3, ip] = Wc3
            parg[igWc3des]  = Wc3
            pare[ieFAR, ip] = FAR
            pare[iemdotf, ip] = mdotf_tot
            pare[ieemot:ieethermal, ip] .= ηpt[2:end]
            pare[ieHrejmot:ieHrejtot, ip] .= Hpt
            pare[ieHexcess, ip] = heatexcess
            # parg[igWtesys] = Wtesys * rlx + parg[igWtesys]*(1.0 - rlx)
            # Engine weight section
                #  Drela's weight model? Nate Fitszgerald - geared TF weight model

            # Snace = Snace1*neng
            # fSnace = Snace/S
            # parg[igfSnace] = fSnace
            lnace = parg[igdfan]*parg[igrSnace]*0.15
            parg[iglnace] = lnace
            
            #Aft fan
            lnace = parg[igdaftfan]*parg[igrSnace]*0.15
            parg[iglnaceaft] = lnace

        # ----------------------
        #     Fly mission
        # ----------------------
        time_propsys += mission!(pari, parg, parm, para, pare, NPSS_TS, NPSS_Fan, NPSS_AftFan, Ldebug)
        parg[igWfuel] = parm[imWfuel] # This is the design mission fuel

        # ----------------------
        #     LH₂ Tank weight
        # ----------------------
        hconvgas = 0.0
        h_LH2 = 210.0
        Tfuel = 20.0
        Tair  = 288.0 #Heated cabin temp
        h_v = 447000.0
        t_cond = [0.05, 1.524e-5, 0.05, 1.524e-5, 1.57e-2] #assumed from energies
        k = ones(length(t_cond)).*5.0e-3#foam conductivities
        hconvair = 15.0 #from sciencedirect.com https://www.sciencedirect.com/topics/engineering/convection-heat-transfer-coefficient
        time_flight = para[iatime, ipdescent1]
        sigskin = 172.4e6 #AL 2219 Brewer / energies stress for operating conditions (290e6 ultimate operatoin)
        rho_insul = [35.24, 14764, 35.24, 14764, 83] #energies
        rhoskintank =  2825.0 #Al 2219 / energies
        max_boiloff = 0.1
        ARtank = 2.0
        clearance_fuse = 0.10
        rhofuel = parg[igrhofuel] 
        ptank = 2.0 #atm

        Wtank_total, thickness_insul, ltank, mdot_boiloff, Vfuel, Wfuel_tot,
        m_boiloff, tskin, t_head, Rtank, Whead, Wcyl,
        Winsul_sum, Winsul, l_tank, Wtank = tanksize(gee, rhofuel, ptank*101325.0,
                      Rfuse, dRfuse, hconvgas, h_LH2, Tfuel, Tair,
                      h_v, t_cond, k, hconvair, time_flight, fstring,ffadd,
                      wfb, nfweb, sigskin, rho_insul, rhoskintank, 
                      parg[igWfuel], max_boiloff, clearance_fuse, ARtank)
        
        parg[igWfmax] = Vfuel*rhofuel*9.81
        parg[igWftank] = Wtank
        parg[igxWftank] = Wtank * parg[igxftank]
        parg[iglftank] = ltank
        parg[igRftank] = Rtank
        parg[igWinsftank] = Winsul_sum


# Get mission fuel burn (check if fuel capacity is sufficent)

# Recalculate weight wupdate()
    ip = ipcruise1
    Wupdate!(parg, rlx, fsum)

    parm[imWTO] = parg[igWMTO]
    parm[imWfuel] = parg[igWfuel]
    # printstyled("Wfuel = $(parg[igWfuel]) \n", color=:blue)


# Set previous iteration weights 
    WMTO3 = WMTO2
    WMTO2 = WMTO1
    WMTO1 = parg[igWMTO]

# END weight sizing loop

# BFL calculations/ Noise? / Engine perf 

    end
 
    endNPSS(NPSS_TS)
    endNPSS(NPSS_Fan)
    endNPSS(NPSS_AftFan)

    # println("Propsys time = ", time_propsys)
end

"""
Wupdate0 updates the weight of the aircraft
"""
function Wupdate0!(parg, rlx, fsum)
    WMTO = parg[igWMTO]

    ftotadd = sum(parg[[igfhpesys, igflgnose, igflgmain]])
    fsum = 0.0

    Wsum = parg[igWpay] + 
            parg[igWfuse ] +
            parg[igWwing ] +
            parg[igWstrut] +
            parg[igWhtail] +
            parg[igWvtail] +
            parg[igWeng  ] +
            parg[igWfuel ] +
            parg[igWtesys] +
            parg[igWftank]

    WMTO = rlx*Wsum/(1.0 - ftotadd) + (1.0 - rlx)*WMTO
    parg[igWMTO] = WMTO

end


"""
Wupdate
"""
function Wupdate!(parg, rlx, fsum)
    
    WMTO = parg[igWMTO]

    fwing  = parg[igWwing ]/WMTO
    fstrut = parg[igWstrut]/WMTO
    fhtail = parg[igWhtail]/WMTO
    fvtail = parg[igWvtail]/WMTO
    feng   = parg[igWeng  ]/WMTO
    ffuel  = parg[igWfuel ]/WMTO
    fhpesys = parg[igfhpesys]
    flgnose = parg[igflgnose]
    flgmain = parg[igflgmain]

    ftank   = parg[igWftank]/WMTO
    ftesys = parg[igWtesys]/WMTO

    Wtesys  = parg[igWtesys]
    Wftank  = parg[igWftank]
    Wpay    = parg[igWpay]
    Wfuse   = parg[igWfuse]

    fsum = fwing + fstrut + fhtail + fvtail + feng + ffuel + fhpesys +
     flgnose + flgmain + ftank + ftesys

    if(fsum ≥ 1.0)
        println("SOMETHING IS WRONG fsum ≥ 1")
    end

    # WMTO = rlx*(Wpay + Wfuse + Wtesys + Wftank)/(1.0-fsum) + (1.0-rlx)*WMTO
    WMTO = rlx*(Wpay + Wfuse )/(1.0-fsum) + (1.0-rlx)*WMTO

    parg[igWMTO  ]  = WMTO
    parg[igWwing ]  = WMTO*fwing
    parg[igWstrut]  = WMTO*fstrut
    parg[igWhtail]  = WMTO*fhtail
    parg[igWvtail]  = WMTO*fvtail
    parg[igWeng  ]  = WMTO*feng
    parg[igWfuel ]  = WMTO*ffuel
    parg[igWftank]  = WMTO*ftank
    parg[igWtesys]  = WMTO*ftesys


end