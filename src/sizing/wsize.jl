using Printf
"""
    wsize(ac; itermax=35,
    wrlx1=0.5, wrlx2=0.9, wrlx3=0.5, initwgt=false, initeng=0, 
    iairf=1, Ldebug=false, printiter=true, saveODperf=false)

Main weight sizing function. Calls on various sub-functions to calculate weight of fuselage, wings, tails, etc.,
and iterates until the MTOW converges to within a specified tolerance.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - Array of flags that control design choices - fuel types, where to store fuel, etc.
    - Geometric and structural parameters - dimensions primarily
    - Aerodynamic parameters - CL, CD, KE dissipation, etc.
    - Mission-specific parameters - alt, mach, P, T etc.
    - Engine-specific parameters 

    **Outputs:**
    - No explicit outputs. Computed quantities are saved to `par` arrays of `aircraft` model.
"""
function wsize(ac; itermax=35,
    wrlx1=0.5, wrlx2=0.9, wrlx3=0.5, initwgt=false, initeng=0, 
    iairf=1, Ldebug=false, printiter=true, saveODperf=false)

    #Unpack data storage arrays
    pari = ac.pari
    parg = ac.parg
    parm = ac.parmd
    para = ac.parad
    pare = ac.pared      
    
    fuse_tank = ac.fuse_tank #Unpack struct with tank parameters
    
    fuse = ac.fuselage 
    wing = ac.wing

    time_propsys = 0.0

    inite1 = 0

    ichoke5 = zeros(iptotal)
    ichoke7 = zeros(iptotal)
    Tmrow = zeros(ncrowx)
    epsrow = zeros(ncrowx)
    epsrow_Tt3 = zeros(ncrowx)
    epsrow_Tt4 = zeros(ncrowx)
    epsrow_Trr = zeros(ncrowx)

    # Weight convergence tolerance 
    # tolerW = 1.0e-10
    tolerW = 1.0e-8
    # tolerW = 1.0e-6
    errw = 1.0
    # Initialze some variables
    fsum = 0.0
    ifirst = true

    # update_fuse!(parg)

    # Flags
    # Fuel type 24 == kerosene #TODO need to update this for LH2
    ifuel = pari[iifuel]
    # ifwcen = pari[iifwcen]
    # iwplan = pari[iiwplan]
    iengloc = pari[iiengloc]
    iengwgt = pari[iiengwgt]
    # iBLIc = pari[iiBLIc]
    # ifclose = pari[iifclose]
    # iHTsize = pari[iiHTsize]
    iVTsize = pari[iiVTsize]
    # ixwmove = pari[iixwmove]
    ifwing = pari[iifwing]

    # Unpack number of powertrain elements
    # nfan = parpt[ipt_nfan]
    ngen = parpt[ipt_ngen]
    nTshaft = parpt[ipt_nTshaft]

    # Calculate fuselage B.L. development at start of cruise: ipcruise1
    time_fusebl = @elapsed fusebl!(fuse,pari, parg, para, ipcruise1)
    # println("Fuse bl time = $time_fusebl")
    # Kinetic energy area at T.E.
    KAfTE = para[iaKAfTE, ipcruise1]
    # Surface dissapation area 
    DAfsurf = para[iaDAfsurf, ipcruise1]
    # Wake dissapation area
    DAfwake = para[iaDAfwake, ipcruise1]
    # Momentum area at ∞
    PAfinf = para[iaPAfinf, ipcruise1]

    # Assume K.E., Disspation and momentum areas are const. for all mission points:
    para[iaKAfTE, :] .= KAfTE
    para[iaDAfsurf, :] .= DAfsurf
    para[iaDAfwake, :] .= DAfwake
    para[iaPAfinf, :] .= PAfinf

    # Set quantities that are fixed during weight iteration

    # Unpack payload and range for design mission - this is the mission that the structures are sized for
    Rangetot = parm[imRange]
    #Typical payload
    Wpay = parm[imWpay]

    Wpaymax = parg[igWpaymax] # Max payload
    # if Wpay or Wpaymax is unset
    if (Wpaymax == 0)
        println("Max payload weight was not set, setting Wpaymax = Wpay")
        Wpaymax = parg[igWpaymax] = max(Wpay, Wpaymax)
    end

    # Store the design mission in the geometry array as well
    parg[igRange] = Rangetot
    parg[igWpay] = Wpay

    # Fixed weight and location of fixed weight
    Wfix = parg[igWfix]
    xfix = parg[igxfix]

    # Weight fractions
    fapu = parg[igfapu]
    fpadd = parg[igfpadd]
    fseat = parg[igfseat]
    feadd = parg[igfeadd]
    fnace = parg[igfnace]
    fhadd = parg[igfhadd]
    fvadd = parg[igfvadd]
    fwadd = wing_additional_weight(wing)
            
    ffadd = parg[igffadd]

    fpylon = parg[igfpylon]

    fhpesys = parg[igfhpesys]
    flgnose = parg[igflgnose]
    flgmain = parg[igflgmain]

    freserve = parg[igfreserve]

    # fuselage lift carryover loss, tip lift loss fractions
    fLo = parg[igfLo]
    fLt = parg[igfLt]

    xwbox = wing.layout.x_wing_box
    xhbox = parg[igxhbox]
    xvbox = parg[igxvbox]
    xapu = parg[igxapu]
    xeng = parg[igxeng]

    # calculate payload proportional weights from weight fractions
    Wapu = Wpaymax * fapu
    Wpadd = Wpaymax * fpadd
    Wseat = Wpaymax * fseat

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

    # tail surface taper ratios (no inner panel, so λs=1)
    λhs = 1.0
    λh = parg[iglambdah]
    λvs = 1.0
    λv = parg[iglambdav]


    # tail geometry parameters
    sweeph = parg[igsweeph]
    wboxh = parg[igwboxh]
    hboxh = parg[ighboxh]
    rhh = parg[igrhh]
    ARh = parg[igARh]
    boh = parg[igboh]

    sweepv = parg[igsweepv]
    wboxv = parg[igwboxv]
    hboxv = parg[ighboxv]
    rhv = parg[igrhv]
    ARv = parg[igARv]
    bov = parg[igbov]

    # number of vertical tails
    nvtail = parg[ignvtail]

    # strut vertical base height, h/c, strut shell t/h
    # zs = parg[igzs]
    # hstrut = parg[ighstrut]
    tohstrut = 0.05

    # assume no struts on tails
    zsh = 0.0
    zsv = 0.0

    # max g load factors for wing, fuselage
    Nlift = parg[igNlift]
    Nland = parg[igNland]

    # never-exceed dynamic pressure for sizing tail structure
    Vne = parg[igVne]
    qne = 0.5 * ρSL * Vne^2

    # wingbox stresses and densities [section 3.1.9 taspot.pdf [prash]] #reference
    σcap = wing.inboard.caps.σ * parg[igsigfac]
    tauweb =  wing.inboard.webs.material.τmax * parg[igsigfac]
    rhoweb = wing.inboard.webs.ρ
    rhocap = wing.inboard.caps.ρ

    # strut stress and density
    σstrut = ac.wing.strut.material.σmax * parg[igsigfac]
    rhostrut = ac.wing.strut.material.ρ

    # assume tail stresses and densities are same as wing's (keeps it simpler)
    σcaph = σcap
    tauwebh = tauweb
    rhowebh = rhoweb
    rhocaph = rhocap

    σcapv = σcap
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

    # set cruise-altitude atmospheric conditions
    set_ambient_conditions!(ac, ipcruise1)

    # set takeoff-altitude atmospheric conditions
    set_ambient_conditions!(ac, iprotate, 0.25)

    # Set atmos conditions for top of climb
    set_ambient_conditions!(ac, ipclimbn)

    nftanks = pari[iinftanks] #Number of fuel tanks in fuselage

    #Fuselage fuel tank placement and moment
    #Initialize parameters

    xfuel = 0.0
    ltank = 0.0

    if (pari[iifwing] == 1) #If fuel is stored in the wings
        xftank = 0.0
        xftankaft = 0.0
    else
        xftank = fuse.layout.x_start_cylinder + 1.0*ft_to_m
        xftankaft = fuse.layout.x_end_cylinder
    end
        
    parg[igxftank] = xftank
    parg[igxftankaft] = xftankaft
   

    # -------------------------------------------------------    
    ## Initial guess section [Section 3.2 of TASOPT docs]
    # -------------------------------------------------------
    # Allow first iteration
    if (initwgt == 0)

        Whtail = 0.05 * Wpay / parg[igsigfac]
        Wvtail = Whtail
        Wwing = 0.5 * Wpay / parg[igsigfac]
        Wstrut = 0.0
        Wftank = 0.0
        Weng = 0.0 * Wpay
        feng = 0.0

        dxWhtail = 0.0
        dxWvtail = 0.0

        # Wing panel weights and moments (after estimating span first)
        ip = ipcruise1
        W = 5.0 * Wpay
        S = W / (0.5 * pare[ierho0, ip] * pare[ieu0, ip]^2 * para[iaCL, ip])
        b = sqrt(S * wing.layout.AR)
        bs = b * wing.layout.ηs
        # TODO: remove above
        # wing.layout.b = sqrt(S * wing.layout.AR)
        # wing.layout.b_inner = wing.layout.b * wing.layout.ηs
        Winn = 0.15 * Wpay / parg[igsigfac]
        Wout = 0.05 * Wpay / parg[igsigfac]
        dyWinn = Winn * 0.30 * (0.5 * (bs - wing.layout.box_halfspan))
        dyWout = Wout * 0.25 * (0.5 * (b - bs))

        parg[igWhtail] = Whtail
        parg[igWvtail] = Wvtail
        wing.weight = Wwing
        wing.strut.weight = Wstrut
        parg[igWeng] = Weng
        parg[igWinn] = Winn
        parg[igWout] = Wout
        parg[igWftank] = Wftank
        parg[igdxWhtail] = dxWhtail
        parg[igdxWvtail] = dxWvtail
        parg[igdyWinn] = dyWinn
        parg[igdyWout] = dyWout


        # wing centroid x-offset form wingbox
        dxwing, macco = surfdx(b, bs, wing.layout.box_halfspan, wing.layout.λt, wing.layout.λs, wing.layout.sweep)
        xwing = xwbox + dxwing
        wing.layout.x = xwing

        # tail area centroid locations (assume no offset from sweep initially)
        parg[igxhtail], parg[igxvtail] = xhbox, xvbox
        # center wingbox chord extent for fuse weight calcs (small effect)
        cbox = 0.0
        # nacelle, fan duct, core, cowl lengths ℛℯ calcs
        parg[iglnace] = 0.5 * S / b

        # nacelle Awet/S
        fSnace = 0.2
        parg[igfSnace] = fSnace

        # Initial fuel fraction estimate from BRE
        LoD = 18.0
        TSFC = 1.0 / 7000.0
        V = pare[ieu0, ipcruise1]
        ffburn = (1.0 - exp(-Rangetot * TSFC / (V * LoD))) # ffburn = Wfuel/WMTO
        ffburn = min(ffburn, 0.8 / (1.0 + freserve))   # 0.8 is the fuel useability? 

        # mission-point fuel fractions ⇾ ffuel = Wfuel/WMTO
        ffuelb = ffburn * (1.0 + freserve)  # start of climb
        ffuelc = ffburn * (0.90 + freserve)  # start of cruise
        ffueld = ffburn * (0.02 + freserve)  # start of descent
        ffuele = ffburn * (0.0 + freserve)  # end of descent (landing) 

        # max fuel fraction is at start of climb
        ffuel = ffuelb # ffuel is max fuel fraction

        # Set initial climb γ = 0 to force intial guesses
        para[iagamV, :] .= 0.0

        # Put initial-guess weight fractions in mission-point array. 

        # These are points before climb starts
        para[iafracW, ipstatic] = 1.0
        para[iafracW, iprotate] = 1.0
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
            frac = float(ip - ipclimb1) / float(ipclimbn - ipclimb1)
            ffp = ffuelb * (1.0 - frac) + ffuelc * frac
            para[iafracW, ip] = 1.0 - ffuel + ffp
        end

        # Cruise
        @inbounds for ip = ipcruise1:ipcruisen
            frac = float(ip - ipcruise1) / float(ipcruisen - ipcruise1)
            ffp = ffuelc * (1.0 - frac) + ffueld * frac
            para[iafracW, ip] = 1.0 - ffuel + ffp
        end

        # Descent
        @inbounds for ip = ipdescent1:ipdescentn
            frac = float(ip - ipdescent1) / float(ipdescentn - ipdescent1)
            ffp = ffueld * (1.0 - frac) + ffuele * frac
            para[iafracW, ip] = 1.0 - ffuel + ffp
        end

        # Initial tail info for sizing of fuselage bending and torsion added material 
        Sh = (2.0 * Wpaymax) / (qne * CLhmax)
        Sv = (2.0 * Wpaymax) / (qne * CLvmax)
        bv = sqrt(Sv * ARv)

        parg[igSh] = Sh
        parg[igSv] = Sv

        # Initial wing and tail pitching moments (including sweep)
        para[iaCMw0, :] .= 0.0
        para[iaCMw1, :] .= 0.0
        para[iaCMh0, :] .= 0.0
        para[iaCMh1, :] .= 0.0
        para[iaCLh, :] .= 0.0

        # Initial cruise-climb angle gamVcr needed to estimate end-of-cruise altitude
        gamVcr = 0.0002
        para[iaCD, ipcruise1] = para[iaCL, ipcruise1] / LoD
        para[iagamV, ipcruise1] = gamVcr

        # Pressure and altitude at start of cruise
        Mach = para[iaMach, ipcruise1]
        p0c = pare[iep0, ipcruise1]
        altc = para[iaalt, ipcruise1]

        # Guess pressure at end-of-cruise (scales with weight)
        p0d = p0c * (1.0 - ffuel + ffueld) / (1.0 - ffuel + ffuelc)
        pare[iep0, ipcruisen] = p0d

        # Guess for OEI #TODO This needs some thinking about what is "One engine out" mean for a turbo-electric aircraft
        # pare[ieFe, iprotate] = 2.0*Wpay/neng
        pare[ieFe, iprotate] = 2.0 * Wpay # ieFe now stores total thrust
        pare[ieu0, iprotate] = 70.0
        Afan = 3.0e-5 * Wpay / neng
        parg[igdfan] = sqrt(Afan * 4.0 / π)

        # Guess fan face mach numbers for nacelle CD calcs
        M2des = pare[ieM2, ipcruise1] = 0.6
        pare[ieM2, ipstatic:ipcruisen] .= M2des
        pare[ieM2, ipdescent1:ipdescentn] .= 0.8 * M2des

        # calculate initial guesses for cooling mass flow ratios epsrow(.)
        ip = iprotate
        cpc = 1080.0
        cp4 = 1340.0
        Rgc = 288.0
        Rg4 = 288.0
        M0to = pare[ieu0, ip] / pare[iea0, ip]
        T0to = pare[ieT0, ip]
        epolhc = pare[ieepolhc, ip]
        OPRto = pare[iepilc, ipcruise1] * pare[iepihc, ipcruise1]
        Tt4to = pare[ieTt4, ip]
        dTstrk = pare[iedTstrk, ip]
        Mtexit = pare[ieMtexit, ip]
        efilm = pare[ieefilm, ip]
        tfilm = pare[ietfilm, ip]
        StA = pare[ieStA, ip]
        for icrow = 1:ncrowx
            Tmrow[icrow] = parg[igTmetal]
        end

        Tt2to = T0to * (1.0 + 0.5 * (gamSL - 1.0) * M0to^2)
        Tt3to = Tt2to * OPRto^(Rgc / (epolhc * cpc))
        Trrat = 1.0 / (1.0 + 0.5 * Rg4 / (cp4 - Rg4) * Mtexit^2)

        ncrow, epsrow, epsrow_Tt3, epsrow_Tt4, epsrow_Trr = mcool(ncrowx, Tmrow,
            Tt3to, Tt4to, dTstrk, Trrat, efilm, tfilm, StA)

        epstot = 0.0
        for icrow = 1:ncrow
            epstot = epstot + epsrow[icrow]
        end
        fo = pare[iemofft, ip] / pare[iemcore, ip]
        fc = (1.0 - fo) * epstot

        for jp = 1:iptotal
            pare[iefc, jp] = fc
            for icrow = 1:ncrowx
                pare[ieepsc1+icrow-1, jp] = epsrow[icrow]
                pare[ieTmet1+icrow-1, jp] = Tmrow[icrow]
            end
        end
        # end

    else #Second iteration onwards use previously calculated values

        bh = parg[igbh]
        bv = parg[igbv]

        coh = parg[igcoh]
        cov = parg[igcov]

        cbox = wing.layout.chord * wing.layout.box_width_chord

        Whtail = parg[igWhtail]
        Wvtail = parg[igWvtail]
        Wwing = wing.weight
        Wstrut = wing.strut.weight
        Weng = parg[igWeng]
        Winn = parg[igWinn]
        Wout = parg[igWout]
        dxWhtail = parg[igdxWhtail]
        dxWvtail = parg[igdxWvtail]
        dyWinn = parg[igdyWinn]
        dyWout = parg[igdyWout]

        WMTO = parg[igWMTO]
        feng = parg[igWeng] / WMTO
        ffuel = parg[igWfuel] / WMTO

        xwing = wing.layout.x
        dxwing = wing.layout.x - wing.layout.x_wing_box

        xhtail = parg[igxhtail]
        xvtail = parg[igxvtail]

        xhbox = parg[igxhbox]
        xvbox = parg[igxvbox]

        Wftank = parg[igWftank]

        fSnace = parg[igfSnace]

        Sh = parg[igSh]
        Sv = parg[igSv]
        ARh = parg[igARh]
        ARv = parg[igARv]

        # Turbo-electric system
        Wtesys = parg[igWtesys]
        Wftank = parg[igWftank]

    end

    # Initialize previous weight iterations
    WMTO1, WMTO2, WMTO3 = zeros(Float64, 3) #1st-previous to 3rd previous iteration weight for convergence criterion

    # no convergence yet
    Lconv = false

    # set these to zero for first-iteration info printout
    wing.layout.b = 0.0
    wing.layout.S = 0.0

    for ip = 1:iptotal
        ichoke5[ip] = 0
        ichoke7[ip] = 0
    end


    # -------------------------------------------------------    
    #                   Weight loop
    # -------------------------------------------------------    

    @inbounds for iterw = 1:itermax
        if iterw == itermax
            println("Reached max iterations in weight sizing loop!")
        end
        if (initwgt == 0)
            #Current weight iteration started from an initial guess so be cautious
            itrlx = 5
        else
            #Current weight started from previously converged solution
            itrlx = 2
        end

        if (iterw <= itrlx)
            # under-relax first n itrlx iterations
            rlx = wrlx1
        elseif (iterw >= 3 / 4 * itermax)
            # under-relax after 3/4th of max iterations
            rlx = wrlx3
        else
            # default is no under-relaxation for weight update
            rlx = wrlx2
        end
        # Fuselage sizing

        # Max tail lifts at maneuver qne
        Lhmax = qne * Sh * CLhmax
        Lvmax = qne * Sv * CLvmax / nvtail
        # Max Δp (fuselage pressure) at end of cruise-climb, assumes p ~ W/S
        wcd = para[iafracW, ipcruisen] / para[iafracW, ipcruise1]
        Δp = parg[igpcabin] - pare[iep0, ipcruise1] * wcd
        parg[igdeltap] = Δp

        # Engine weight mounted on tailcone, if any
        if (iengloc == 1) # 1: Eng on wing. 2: Eng on aft fuse
            Wengtail = 0.0
            Waftfuel = 0.0
        else
            Wengtail = (parg[igWtshaft] + parg[igWcat]) * nTshaft +
                       parg[igWgen] * ngen
        end

        Whtail = parg[igWhtail]
        Wvtail = parg[igWvtail]
        xhtail = parg[igxhtail]
        xvtail = parg[igxvtail]
        xwbox = wing.layout.x_wing_box
        xwing = wing.layout.x
        xapu = parg[igxapu]
        xeng = parg[igxeng]

        Wtesys = parg[igWtesys]
        nftanks = pari[iinftanks]

        ifwing = pari[iifwing]
        if ifwing == 0 #fuselage fuel store
            tank_placement = fuse_tank.placement
            Wftank_single = parg[igWftank] / nftanks #Weight of a single tank

            #Calculate the weight of the fuel near the tail depending on the tank location
            if tank_placement == "rear"
                Waftfuel = parg[igWfuel]
                xftank_fuse = parg[igxftankaft] #assumed location of tank for fuselage sizing
            elseif tank_placement == "both"
                Waftfuel = parg[igWfuel] / 2.0
                xftank_fuse = parg[igxftankaft]

            elseif tank_placement == "front" #The case when the fuel is at the front is treated specially
                #The code assumes that the fuel is located at the back for the purpose of sizing of the symmetric fuselage
                Waftfuel = parg[igWfuel]
                xftank_fuse = fuse.layout.x_end - parg[igxftank]
            end
        else
            tank_placement = ""
            xftank_fuse = 0.0
            Wftank_single = 0.0
        end
        
        (cabVol) = fusew!(fuse, Nland, Wfix, Wpaymax, Wpadd, Wseat, Wapu, Wengtail, 
            ifwing, nftanks,
            Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
            ffadd, Δp,
            Wpwindow, Wppinsul, Wppfloor,
            Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
            bv, λv, nvtail,
            
            xhtail, xvtail,
            xwing, xwbox, cbox,
            xfix, xapu, xeng, xfuel)

        parg[igcabVol] = cabVol

        # Use cabin volume to get actual buoyancy weight
        ρcab = max(parg[igpcabin], pare[iep0, ipcruise1]) / (RSL * TSL)
        WbuoyCR = (ρcab - pare[ierho0, ipcruise1]) * gee * cabVol

        if (iterw == 1 && initwgt == 0)

            feng = 0.08
            fsum = feng + ffuel + fhpesys + flgnose + flgmain
            WMTO = (Wpay + fuse.weight + Wwing + Wstrut + Whtail + Wvtail) / (1.0 - fsum)

            Weng, Wfuel, Whpesys, Wlgnose, Wlgmain = WMTO .* [feng, ffuel, fhpesys, flgnose, flgmain]
            parg[igWMTO] = WMTO
            parg[igWeng] = Weng
            parg[igWfuel] = Wfuel
            println("Wfuel initial = ", (ffuel * WMTO))

        else
            # Call a better Wupdate function
            Wupdate0!(parg, fuse, wing, rlx, fsum)
            if (fsum >= 1.0)
                println("Something is wrong!! fsum ≥ 1.0")
                break
            end

            parm[imWTO] = parg[igWMTO]
            parm[imWfuel] = parg[igWfuel]

        end
        # this calculated WMTO is the design-mission WTO
        parm[imWTO] = parg[igWMTO]
        # Convergence tests
        WMTO = parg[igWMTO]
        errw1 = (WMTO - WMTO1) / WMTO
        errw2 = (WMTO - WMTO2) / WMTO
        errw3 = (WMTO - WMTO3) / WMTO

        errw = max(abs(errw1), abs(errw2), abs(errw3))

        # Print weight/ convergnce started
        if (printiter && iterw == 1)
            @printf("%5s %15s %15s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s \n",
                "iterw", "errW", "errW1", "WMTO", "Wpay", "Wfuel", "Weng", "Wfuse", "Wwing", "span", "area", "HTarea", "xwbox")
        end
        if printiter
            @printf("%5d %+13.8e %+13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",
                iterw, errw, errw1, parm[imWTO], parg[igWpaymax], parg[igWfuel], parg[igWeng],
                fuse.weight, wing.weight, wing.layout.b, wing.layout.S,
                parg[igSh], wing.layout.x_wing_box)
        end
        if (errw <= tolerW)
            Lconv = true
            break
        end

        #--------------------------------
        ##  Wing sizing section
        #--------------------------------
        WMTO = parg[igWMTO]

        # Size wing area and chords at start-of-cruise
        ip = ipcruise1
        W = WMTO * para[iafracW, ip]
        CL = para[iaCL, ip]
        ρ0 = pare[ierho0, ip]
        u0 = pare[ieu0, ip]
        qinf = 0.5 * ρ0 * u0^2
        BW = W + WbuoyCR #Weight including buoyancy

        # Initial size of the wing area and chords
        wingsc!(BW, CL, qinf, wing)
        # wing.layout.S = S
        # wing.layout.b = b
        # wing.layout.b_inner = bs

        #Updating wing box chord for fuseW in next iteration
        cbox = wing.layout.chord * wing.layout.box_width_chord

        # x-offset of the wing centroid from wingbox
        dxwing, macco = surfdx(wing.layout.b, wing.layout.b_inner, wing.layout.box_halfspan, wing.layout.λt, wing.layout.λs, wing.layout.sweep)
        xwing = xwbox + dxwing
        cma = macco * wing.layout.chord
        wing.layout.x = xwing
        parg[igcma] = cma

        # Calculate wing pitching moment constants
        #------------------------------------------
        ## Takeoff
        ip = iptakeoff
        cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
        γt = wing.layout.λt * para[iarclt, ip]
        γs = wing.layout.λs * para[iarcls, ip]

        CMw0, CMw1 = surfcm(wing.layout.b, wing.layout.b_inner, wing.layout.box_halfspan, wing.layout.sweep, wing.layout.spar_box_x_c,
        wing.layout.λt, wing.layout.λs, γt, γs,
            wing.layout.AR, fLo, fLt, cmpo, cmps, cmpt)

        para[iaCMw0, ipstatic:ipclimb1] .= CMw0
        para[iaCMw1, ipstatic:ipclimb1] .= CMw1

        ## Cruise
        ip = ipcruise1
        cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
        γt = wing.layout.λt * para[iarclt, ip]
        γs = wing.layout.λs * para[iarcls, ip]

        CMw0, CMw1 = surfcm(wing.layout.b, wing.layout.b_inner, wing.layout.box_halfspan, wing.layout.sweep, wing.layout.spar_box_x_c,
        wing.layout.λt, wing.layout.λs, γt, γs,
            wing.layout.AR, fLo, fLt, cmpo, cmps, cmpt)

        para[iaCMw0, ipclimb1+1:ipdescentn-1] .= CMw0
        para[iaCMw1, ipclimb1+1:ipdescentn-1] .= CMw1

        ## Descent
        ip = ipdescentn
        cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
        γt = wing.layout.λt * para[iarclt, ip]
        γs = wing.layout.λs * para[iarcls, ip]

        CMw0, CMw1 = surfcm(wing.layout.b, wing.layout.b_inner, wing.layout.box_halfspan, wing.layout.sweep, wing.layout.spar_box_x_c,
        wing.layout.λt, wing.layout.λs, γt, γs,
            wing.layout.AR, fLo, fLt, cmpo, cmps, cmpt)

        para[iaCMw0, ipdescentn] = CMw0
        para[iaCMw1, ipdescentn] = CMw1
        #------------------------------------------

        # Wing center load po calculation using cruise spanload cl(y)
        # -----------------------------------------
        ip = ipcruise1
        γt, γs = wing.layout.λt * para[iarclt, ip], wing.layout.λs * para[iarcls, ip] # Lift "taper ratios"
        Lhtail = WMTO * parg[igCLhNrat] * parg[igSh] / wing.layout.S

        po = wingpo(wing.layout.b, wing.layout.b_inner, wing.layout.box_halfspan,
        wing.layout.λt, wing.layout.λs, γt, γs,
            wing.layout.AR, Nlift, BW, Lhtail, fLo, fLt)

        if (wing.planform == 1)
            # engines on wing, at ys=wing.layout.b_inner/2
            Weng1 = parg[igWeng] / neng
        else
            # engines not mounted on wing
            Weng1 = 0.0
        end

        # HACK not tested
        nout = 0
        yout = 0.0 # avg. moment arm of outboard engines
        nin = 0
        yinn = 0.0

        if (wing.planform == 1)
            if pari[iiengtype] == 0
                Weng1 = parg[igWfan] + parg[igWmot] + parg[igWfanGB]
            else
                Weng1 = parg[igWeng] / parg[igneng]
            end

        else
            Weng1 = 0.0
        end

        Winn = parg[igWinn]
        Wout = parg[igWout]
        dyWinn = parg[igdyWinn]
        dyWout = parg[igdyWout]
        if (pari[iifwing] == 0)
            rhofuel = 0.0 # tell surfw that there is no fuel in wings
        else
            rhofuel = parg[igrhofuel]
        end
        Ecap = parg[igEcap]
        Eweb = Ecap
        Gcap = Ecap * 0.5 / (1.0 + 0.3)
        Gweb = Ecap * 0.5 / (1.0 + 0.3)


        Ss, Ms, tbwebs, tbcaps, EIcs, EIns, GJs,
        So, Mo, tbwebo, tbcapo, EIco, EIno, GJo,
        Astrut, lstrutp, cosLs,
        Wscen, Wsinn, Wsout, dxWsinn, dxWsout, dyWsinn, dyWsout,
        Wfcen, Wfinn, Wfout, dxWfinn, dxWfout, dyWfinn, dyWfout,
        Wweb, Wcap, Wstrut,
        dxWweb, dxWcap, dxWstrut = surfw(po, wing.layout.b, wing.layout.b_inner, wing.layout.box_halfspan, wing.layout.chord, wing.strut.z,
        wing.layout.λt, wing.layout.λs, γt, γs,
            Nlift, wing.planform, Weng1,
            nout, yout, nin, yinn,
            Winn, Wout, dyWinn, dyWout,
            wing.layout.sweep, wing.layout.box_width_chord, wing.layout.root_chord_thickness, wing.layout.spanbreak_chord_thickness, wing.layout.hweb_to_hbox, fLt,
            tauweb, σcap, σstrut, Ecap, Eweb, Gcap, Gweb,
            rhoweb, rhocap, rhostrut, rhofuel)
        # println([Ss,Ms,tbwebs,tbcaps,EIcs,EIns,GJs,
        # So,Mo,tbwebo,tbcapo,EIco,EIno,GJo,
        # Astrut,lstrutp,cosLs,
        # Wscen,Wsinn,Wsout,dxWsinn,dxWsout,dyWsinn,dyWsout,
        # Wfcen,Wfinn,Wfout,dxWfinn,dxWfout,dyWfinn,dyWfout,
        # Wweb,  Wcap,  Wstrut,
        # dxWweb,dxWcap,dxWstrut])
        # Wsinn is the in-board skin weight, Wsout is the outboard skin weight. 
        # Multiply by 2 to account for the two wing-halves: 
        Wwing = 2.0 * (Wscen + Wsinn + Wsout) * (1.0 + fwadd)
        dxWwing = 2.0 * (dxWsinn + dxWsout) * (1.0 + fwadd)

        # Note this assumes wings have some fuel, so additional check is performed to see if iifwing is 1
        #Calculate volume limited fuel weight depending on if wing center box has fuel or not
        Wfmax = 0.0
        dxWfmax = 0.0
        rfmax = 0.0
        if (pari[iifwing] == 1) # if fuel is stored in wings only then do this
            if (pari[iifwcen] == 0)
                Wfmax = 2.0 * (Wfinn + Wfout)
                dxWfmax = 2.0 * (dxWfinn + dxWfout)
            else
                Wfmax = 2.0 * (Wfcen + Wfinn + Wfout)
                dxWfmax = 2.0 * (dxWfinn + dxWfout)
            end

            #at full payload, the fuel tank cannot be full, so less bending relief from fuel
            Wfuelmp = Wpay - Wpaymax + parg[igWfuel] # Wfuelmp == parg[igWfuel] if Wpaymax = Wpay
            rfmax = Wfuelmp / Wfmax
        end


        # Save wing details into geometry array
        wing.weight = Wwing * rlx + wing.weight * (1.0 - rlx)
        
        parg[igWfmax] = Wfmax
        parg[igdxWwing] = dxWwing
        parg[igdxWfuel] = dxWfmax * rfmax

        parg[igtbwebs] = tbwebs
        parg[igtbcaps] = tbcaps
        parg[igtbwebo] = tbwebo
        parg[igtbcapo] = tbcapo
        parg[igAstrut] = Astrut
        parg[igcosLs] = cosLs
        parg[igWweb] = Wweb
        parg[igWcap] = Wcap
        wing.strut.weight = Wstrut
        parg[igSomax] = So
        parg[igMomax] = Mo
        parg[igSsmax] = Ss
        parg[igMsmax] = Ms
        parg[igEIco] = EIco
        parg[igEIcs] = EIcs
        parg[igEIno] = EIno
        parg[igEIns] = EIns
        parg[igGJo] = GJo
        parg[igGJs] = GJs

        wing.strut.weight = Wstrut
        parg[igdxWstrut] = dxWstrut

        # Strut chord (perpendicular to strut)
        cstrut = sqrt(0.5 * Astrut / (tohstrut * wing.strut.toc))
        Ssturt = 2.0 * cstrut * lstrutp
        parg[igcstrut] = cstrut
        wing.strut.area = Ssturt

        # Individual panel weights
        Winn = Wsinn * (1.0 + fwadd) + rfmax * Wfinn
        Wout = Wsout * (1.0 + fwadd) + rfmax * Wfout

        dyWinn = dyWsinn * (1.0 + fwadd) + rfmax * dyWfinn
        dyWout = dyWsout * (1.0 + fwadd) + rfmax * dyWfout

        parg[igWinn] = Winn
        parg[igWout] = Wout
        parg[igdyWinn] = dyWinn
        parg[igdyWout] = dyWout

        # -------------------------------
        #      Tail sizing section
        # -------------------------------

        # Set tail CL derivative 
        dϵdα = parg[igdepsda]
        sweeph = parg[igsweeph]
        tanL = tan(wing.layout.sweep * π / 180.0)
        tanLh = tan(sweeph * π / 180.0)

        ip = ipcruise1
        Mach = para[iaMach, ip]
        β = sqrt(1.0 - Mach^2) #Prandtl-Glauert factor √(1-M²)
        # Calculate the tail lift-curve slope
        dCLhdCL = (β + 2.0 / wing.layout.AR) / (β + 2.0 / ARh) * sqrt(β^2 + tanL^2) / sqrt(β^2 + tanLh^2) * (1.0 - dϵdα)
        parg[igdCLhdCL] = dCLhdCL

        # Set Nacelle CL derivative fraction
        dCLnda = parg[igdCLnda]
        dCLndCL = dCLnda * (β + 2.0 / wing.layout.AR) * sqrt(β^2 + tanL^2) / (2.0 * π * (1.0 + 0.5 * wing.layout.root_chord_thickness))
        parg[igdCLndCL] = dCLndCL

        # Fuselage pitching moment
        #       Use this with caution - slender body theory is used here to estimate the fuselage 
        #       pitching moment - this ofc isn't true if the aircraft fuselage isn't "slender"
        #       Drela used a 3D panel method to actually calculate the CMVf1 and CMV0  for the aircraft studied in the N+3 work
        #       If sizes are roughly that of the 737/ 777 or D8 perhaps best to use those values and comment out the following bits of code
        #TODO: Add switch to either calculate fuse pitching moment online or use offline specified values
        cosL = cos(wing.layout.sweep * π / 180.0)
        Mperp = Mach * cosL
        βn = sqrt(1 - Mperp^2) # PG correction factor with M⟂ 
        # Estimate finite wing ∂CL/∂α from thin airfoil lift-slope 2π and 
        #  corrections for sweep and compressibility:
        CLα = 2π * cosL / (sqrt(βn^2 + (2 * cosL / wing.layout.AR)^2) + 2 * cosL / wing.layout.AR)
        # Estimate CMVf1 via slender body theory: dM/dα = 𝒱 ⟹ dM/dCL = dM/dα × dα/dCL = 𝒱/(dCL/dα)
        # parg[igCMVf1] = parg[igfuseVol]/CLα

        # Size HT
        if (iterw <= 2 && initwgt == 0)
            #if initial iterations or intiial weight =0 then just get tail volume coeff Vh
            lhtail = xhtail - xwing
            Vh = parg[igVh]
            Sh = Vh * wing.layout.S * cma / lhtail
            parg[igSh] = Sh
        else
            # for subsequent iterations:
            htsize(pari, parg, view(para, :, ipdescentn), view(para, :, ipcruise1), view(para, :, ipcruise1), fuse, wing)

            xwbox, xwing = wing.layout.x_wing_box, wing.layout.x

            lhtail = xhtail - xwing
            Sh = parg[igSh]

            parg[igVh] = Sh * lhtail / (wing.layout.S * cma)
        end

        # Vertical tail sizing 

        # Estimate thrust at take-off and the resultant moment if OEI - to estimate Vertical tail size
        #   Same logic should hold if outboard electric motors fail, however electrical power can be redirected to maintain symmetric Thrust
        #   even if turboshaft fails. Therefore we assume that the worst yaw moment is if the aft engine fails the thrust lost is from the BLI ducted fan 
        #   connected mechanically to the turboshaft engines.
        # [Section 2.13.2 of TASOPT docs]

        ip = iprotate
        qstall = 0.5 * pare[ierho0, ip] * (pare[ieu0, ip] / 1.2)^2
        dfan = parg[igdfan]
        CDAe = parg[igcdefan] * 0.25π * dfan^2
        De = qstall * CDAe
        Fe = pare[ieFe, ip]

        # Calcualte max eng out moment
        Me = (Fe + De) * yeng

        if (iVTsize == 1)
            lvtail = xvtail - xwing
            Vv = parg[igVv]
            Sv = Vv * wing.layout.S * wing.layout.b / lvtail
            parg[igSv] = Sv
            parg[igCLveout] = Me / (qstall * Sv * lvtail) # Max lift coeff oc vertical tail with some yaw control when OEI [Eqn. 312 of TASOPT docs]
        else
            lvtail = xvtail - xwing
            CLveout = parg[igCLveout]
            Sv = Me / (qstall * CLveout * lvtail)
            parg[igSv] = Sv
            parg[igVv] = Sv * lvtail / (wing.layout.S * wing.layout.b)
        end

        # set HT max loading magnitude
        bh, coh, poh = tailpo(Sh, ARh, λh, qne, CLhmax)
        parg[igbh] = bh
        parg[igcoh] = coh

        # set VT max loading magnitude, based on single tail + its bottom image
        bv2, cov, pov = tailpo(2.0 * Sv / nvtail, 2.0 * ARv, λv, qne, CLvmax)
        bv = bv2/2
        parg[igbv] = bv
        parg[igcov] = cov

        # HT weight
        γh = λh
        γhs = λhs
        ihplan = 0
        Wengh = 0.0
        Ecap = parg[igEcap]
        Eweb = Ecap
        Gcap = Ecap * 0.5 / (1.0 + 0.3)
        Gweb = Ecap * 0.5 / (1.0 + 0.3)
        Ssh, Msh, tbwebsh, tbcapsh, EIcsh, EInsh, GJsh,
        Soh, Moh, tbweboh, tbcapoh, EIcoh, EInoh, GJoh,
        _, _, _,
        Wscenh, Wsinnh, Wsouth, dxWsinnh, dxWsouth, dyWsinnh, dyWsouth,
        Wfcenh, Wfinnh, Wfouth, dxWfinnh, dxWfouth, dyWfinnh, dyWfouth,
        Wwebh, Wcaph, Wstruth,
        dxWwebh, dxWcaph, dxWstruth = surfw(poh, bh, boh, boh, coh, zsh,
            λh, λhs, γh, γhs,
            1, ihplan, Wengh, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            sweeph, wboxh, hboxh, hboxh, rhh, fLt,
            tauwebh, σcaph, σstrut, Ecap, Eweb, Gcap, Gweb,
            rhowebh, rhocaph, rhostrut, rhofuel)

        Whtail = 2.0 * (Wscenh + Wsinnh + Wsouth) * (1.0 + fhadd)
        dxWhtail = 2.0 * (dxWsinnh + dxWsouth) * (1.0 + fhadd)
        parg[igWhtail] = Whtail
        parg[igdxWhtail] = dxWhtail

        parg[igtbwebh] = tbweboh
        parg[igtbcaph] = tbcapoh
        parg[igEIch] = EIcoh
        parg[igEInh] = EInoh
        parg[igGJh] = GJoh

        # HT centroid x-offset
        xhbox = parg[igxhbox]
        dxh, macco = surfdx(bh, boh, boh, λh, λhs, sweeph)
        parg[igxhtail] = xhbox + dxh

        # HT pitching moment coeff
        fLoh = 0.0
        fLth = fLt
        cmph = 0.0
        CMh0, CMh1 = surfcm(bh, boh, boh, sweeph, wing.layout.spar_box_x_c, λh, 1.0, λh, 1.0,
            ARh, fLoh, fLth, 0.0, 0.0, 0.0)
        para[iaCMh0, :] .= CMh0
        para[iaCMh1, :] .= CMh1

        # VT weight

        γv = λv
        γvs = λvs
        ivplan = 0
        Wengv = 0.0
        Ecap = parg[igEcap]
        Eweb = Ecap
        Gcap = Ecap * 0.5 / (1.0 + 0.3)
        Gweb = Ecap * 0.5 / (1.0 + 0.3)
        Ssv, Msv, tbwebsv, tbcapsv, EIcsv, EInsv, GJsv,
        Sov, Mov, tbwebov, tbcapov, EIcov, EInov, GJov,
        _, _, _,
        Wscenv, Wsinnv, Wsoutv, dxWsinnv, dxWsoutv, dyWsinnv, dyWsoutv,
        Wfcenv, Wfinnv, Wfoutv, dxWfinnv, dxWfoutv, dyWfinnv, dyWfoutv,
        Wwebv2, Wcapv2, Wstrutv,
        dxWwebv2, dxWcapv2, dxWstrutv = surfw(pov, bv2, bov, bov, cov, zsv,
            λv, λvs, γv, γvs,
            1.0, ivplan, Wengv, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            sweepv, wboxv, hboxv, hboxv, rhv, fLt,
            tauwebv, σcapv, σstrut, Ecap, Eweb, Gcap, Gweb,
            rhowebv, rhocapv, rhostrut, rhofuel)

        Wvtail = (Wscenv + Wsinnv + Wsoutv) * (1.0 + fvadd) * nvtail
        dxWvtail = (dxWsinnv + dxWsoutv) * (1.0 + fvadd) * nvtail
        parg[igWvtail] = Wvtail
        parg[igdxWvtail] = dxWvtail

        parg[igtbwebv] = tbwebov
        parg[igtbcapv] = tbcapov
        parg[igEIcv] = EIcov
        parg[igEInv] = EInov
        parg[igGJv] = GJov

        # VT centroid x-offset
        xvbox = parg[igxvbox]
        dxv, _ = surfdx(bv2, bov, bov, λv, λvs, sweepv)
        parg[igxvtail] = xvbox + dxv

        # ----------------------
        #     Fuselage Fuel Tank weight
        # ----------------------
        if (pari[iifwing] == 0) #If fuel is stored in the fuselage
            #Unpack parameters
            time_flight = para[iatime, ipdescent1]
            tank_placement = fuse_tank.placement
            rhofuel = fuse_tank.rhofuel

            #Convective cooling
            if tank_placement == "rear"
                xftank_heat = parg[igxftankaft]
            else
                xftank_heat = parg[igxftank]
            end
            ifuel = pari[iifuel]
            M_inf = para[iaMach, ipcruise1]
            z_alt = para[iaalt, ipcruise1]
            
            #Fuel tank design
            fuse_tank.Wfuelintank = parg[igWfuel] / nftanks #Each fuel tank carries 1/nftanks of the fuel
            
            mdot_boiloff, Vfuel, Rtank, Winsul_sum, ltank, Wtank = tanksize!(fuse_tank, z_alt, M_inf, xftank_heat,
            time_flight, ifuel)

            parg[igWfmax] = Vfuel * rhofuel * gee * nftanks #If more than one tank, max fuel capacity is nftanks times that of one tank
            parg[igWftank] = nftanks * Wtank #total weight of fuel tanks (including insulation)
            parg[iglftank] = ltank
            parg[igRftank] = Rtank
            parg[igWinsftank] = nftanks * Winsul_sum #total weight of insulation in fuel tanks
            parg[igmdotboiloff] = nftanks * mdot_boiloff #store total fuel boiloff rate

            #Tank placement and weight moment
            lcabin = dx_cabin(fuse)
            if tank_placement == "front"
                flag_front = 1
                flag_aft = 0
                xftank = fuse.layout.x_start_cylinder + 1.0*ft_to_m + ltank/2.0
                xftankaft = 0.0
            elseif tank_placement == "rear"
                flag_front = 0
                flag_aft = 1
                xftank = 0.0
                xftankaft = fuse.layout.x_start_cylinder + lcabin + 1.0*ft_to_m + ltank/2.0
            elseif tank_placement == "both"
                flag_front = 1
                flag_aft = 1
                xftank = fuse.layout.x_start_cylinder + 1.0*ft_to_m + ltank/2.0
                xftankaft = fuse.layout.x_start_cylinder + 1.0*ft_to_m + ltank + 1.0*ft_to_m + lcabin + 1.0*ft_to_m + ltank/2.0
            end
            
            parg[igxftank] = xftank
            parg[igxftankaft] = xftankaft
            parg[igxWftank] = Wtank * (flag_front * xftank + flag_aft * xftankaft) 
            xfuel = (flag_front * xftank + flag_aft * xftankaft) / (flag_front + flag_aft)
            parg[igxWfuel] = parg[igWfuel] * xfuel

            # Update fuselage according to tank requirements
            update_fuse!(pari, parg) #update fuselage length to accommodate tank
            fusebl!(fuse,pari, parg, para, ipcruise1) #Recalculate fuselage bl properties

            #Update fuselage BL properties
            # Kinetic energy area at T.E.
            KAfTE = para[iaKAfTE, ipcruise1]
            # Surface dissapation area 
            DAfsurf = para[iaDAfsurf, ipcruise1]
            # Wake dissapation area
            DAfwake = para[iaDAfwake, ipcruise1]
            # Momentum area at ∞
            PAfinf = para[iaPAfinf, ipcruise1]

            # Assume K.E., Disspation and momentum areas are const. for all mission points:
            para[iaKAfTE, :] .= KAfTE
            para[iaDAfsurf, :] .= DAfsurf
            para[iaDAfwake, :] .= DAfwake
            para[iaPAfinf, :] .= PAfinf

        end

        # -----------------------------
        # Heat exchanger design and operation
        # ------------------------------
        ipdes = ipcruise1 #Design point: start of cruise

        if iterw > 2 #Only include heat exchangers after second iteration
            global HXs = hxdesign!(pare, pari, ipdes, HXs)
            #global HXs_prev = deepcopy(HXs) #Store current heat exchange vector as previous for debugging

        else
            global HXs = []
            #global HXs_prev = deepcopy(HXs) #Store current heat exchange vector as previous
            
        end

        # -----------------------------
        # Drag and engine calculations
        # ------------------------------
        # Total Drag

        WMTO = parg[igWMTO]
        #calculate for start-of-cruise point
        # ip = ipclimbn
        ip = ipcruise1

        #Calculate fuel weight moment for balance
        if (pari[iifwing] == 1) #If fuel is stored in the wings
            xfuel = wing.layout.x_wing_box + parg[igdxWfuel] / parg[igWfuel]
            parg[igxWfuel] = parg[igWfuel] * wing.layout.x_wing_box + parg[igdxWfuel] #Store fuel weight moment
        end

        # Pitch trim by adjusting Clh or by moving wing
        Wzero = WMTO - parg[igWfuel] #Zero fuel weight
        Wf = para[iafracW, ip] * WMTO - Wzero
        rfuel = Wf / parg[igWfuel]
        rpay = 1.0
        ξpay = 0.0
        itrim = 1
        balance(pari, parg, view(para, :, ip), fuse, wing, rfuel, rpay, ξpay, itrim)

        # Set N.P. at cruise
        parg[igxNP] = para[iaxNP, ip]

        para[iaalt, ipclimbn] = para[iaalt, ipcruise1]

        # Drag buildup cdsum()
        cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), wing,  1)

        # L/D and Design point thrust
        # println("CD = ", para[iaCD,ip])
        LoD = para[iaCL, ip] / para[iaCD, ip]
        gamV = para[iagamV, ip]
        W = para[iafracW, ip] * WMTO
        BW = W + WbuoyCR
        # Fdes = BW * (1 / LoD + gamV) * 1.05 #Ad-hoc 5% addition for OEI
        Fdes = BW * (1 / LoD + gamV)

        pare[ieFe, ip] = Fdes / neng

        # Size engine for TOC
        icall = 0
        icool = 1
        if (iterw == 1 || initeng == 0)
            # initialize engine state
            inite1 = 0
        else
            # start with current engine state
            inite1 = 1
        end

        ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), wing, ip, icall, icool, inite1)

        # store engine design-point parameters for all operating points
        parg[igA5] = pare[ieA5, ip] / pare[ieA5fac, ip]
        parg[igA7] = pare[ieA7, ip] / pare[ieA7fac, ip]
        for jp = 1:iptotal
            pare[ieA2, jp] = pare[ieA2, ip]
            pare[ieA25, jp] = pare[ieA25, ip]
            pare[ieA5, jp] = parg[igA5] * pare[ieA5fac, jp]
            pare[ieA7, jp] = parg[igA7] * pare[ieA7fac, jp]

            pare[ieNbfD, jp] = pare[ieNbfD, ip]
            pare[ieNblcD, jp] = pare[ieNblcD, ip]
            pare[ieNbhcD, jp] = pare[ieNbhcD, ip]
            pare[ieNbhtD, jp] = pare[ieNbhtD, ip]
            pare[ieNbltD, jp] = pare[ieNbltD, ip]

            pare[iembfD, jp] = pare[iembfD, ip]
            pare[iemblcD, jp] = pare[iemblcD, ip]
            pare[iembhcD, jp] = pare[iembhcD, ip]
            pare[iembhtD, jp] = pare[iembhtD, ip]
            pare[iembltD, jp] = pare[iembltD, ip]

            pare[iepifD, jp] = pare[iepifD, ip]
            pare[iepilcD, jp] = pare[iepilcD, ip]
            pare[iepihcD, jp] = pare[iepihcD, ip]
            pare[iepihtD, jp] = pare[iepihtD, ip]
            pare[iepiltD, jp] = pare[iepiltD, ip]
        end

        dfan = parg[igdfan]
        dlcomp = parg[igdlcomp]
        dhcomp = parg[igdhcomp]

        Mach = para[iaMach, ip]
        CL = para[iaCL, ip]
        CD = para[iaCD, ip]

        # bare weight for one engine [Newtons]
        mdotc = pare[iemblcD, ip] * sqrt(Tref / TSL) * (pSL / pref)
        BPR = pare[ieBPR, ip]
        OPR = pare[iepilc, ip] * pare[iepihc, ip]

        # weight of engine and related stuff
        Gearf = parg[igGearf]
        Weng, Wnace, Webare, Snace1 = tfweight(iengwgt, Gearf, OPR, BPR, mdotc, dfan, rSnace,
            dlcomp, neng, feadd, fpylon, HXs)

        parg[igWeng] = Weng
        parg[igWebare] = Webare
        parg[igWnace] = Wnace
        parg[igWeng] = Weng

        # set new nacelle area / reference area  fraction fSnace
        Snace = Snace1 * neng
        fSnace = Snace / wing.layout.S
        parg[igfSnace] = fSnace
        lnace = parg[igdfan] * parg[igrSnace] * 0.15
        parg[iglnace] = lnace

        ipc1 = 1
        time_propsys += mission!(pari, parg, parm, para, pare, fuse, wing, Ldebug)

        # this calculated fuel is the design-mission fuel 
        parg[igWfuel] = parm[imWfuel]
        
        # size cooling mass flow at takeoff rotation condition (at Vstall)
        ip = iprotate

        # must define CDwing for this point in case there's wing BLI
        cdfw = para[iacdfw, ip] * para[iafexcdw, ip]
        cdpw = para[iacdpw, ip] * para[iafexcdw, ip]
        cosL = cos(wing.layout.sweep * pi / 180.0)
        para[iaCDwing, ip] = cdfw + cdpw * cosL^3

        icall = 1
        icool = 2
        ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), wing, ip, icall, icool, inite1)

        # Tmetal was specified... set blade row cooling flow ratios for all points
        for jp = 1:iptotal
            for icrow = 1:ncrowx
                pare[ieepsc1+icrow-1, jp] = pare[ieepsc1+icrow-1, ip]
            end
            # also set first estimate of total cooling mass flow fraction
            pare[iefc, jp] = pare[iefc, ip]
        end

        # Recalculate weight wupdate()
        ip = ipcruise1
        Wupdate!(parg, fuse, wing,rlx, fsum)

        parm[imWTO] = parg[igWMTO]
        parm[imWfuel] = parg[igWfuel]

        # Set previous iteration weights 
        WMTO3 = WMTO2
        WMTO2 = WMTO1
        WMTO1 = parg[igWMTO]

        ifirst = false

        # Get mission fuel burn (check if fuel capacity is sufficent)

        # Recalculate weight wupdate()
        ip = ipcruise1
        Wupdate!(parg, fuse, wing, rlx, fsum)

        parm[imWTO] = parg[igWMTO]
        parm[imWfuel] = parg[igWfuel]
        # printstyled("Wfuel = $(parg[igWfuel]) \n", color=:blue)

        # Set previous iteration weights 
        WMTO3 = WMTO2
        WMTO2 = WMTO1
        WMTO1 = parg[igWMTO]

        # BFL calculations/ Noise? / Engine perf 

    end

    # normal takeoff and balanced-field takeoff calculations
    # set static thrust for takeoff routine
    ip = ipstatic
    icall = 1
    icool = 1

    ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), wing, ip, icall, icool, inite1)

    # set rotation thrust for takeoff routine
    # (already available from cooling calculations)
    ip = iprotate
    icall = 1
    icool = 1
    ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), wing, ip, icall, icool, inite1)

    # calculate takeoff and balanced-field lengths
    takeoff!(pari, parg, parm, para, pare, wing, initeng, ichoke5, ichoke7)

    # calculate CG limits from worst-case payload fractions and packings
    rfuel0, rfuel1, rpay0, rpay1, xCG0, xCG1 = cglpay(pari, parg,fuse, wing)
    parg[igxCGfwd] = xCG0
    parg[igxCGaft] = xCG1
    parg[igrpayfwd] = rpay0
    parg[igrpayaft] = rpay1

    # set neutral point at cruise
    ip = ipcruise1
    Wzero = WMTO - parg[igWfuel]
    Wf = para[iafracW, ip] * WMTO - Wzero
    rfuel = Wf / parg[igWfuel]
    rpay = 1.0
    ξpay = 0.0
    itrim = 0
    balance(pari, parg, view(para, :, ip), fuse, wing, rfuel, rpay, ξpay, itrim)
    
end

"""
Wupdate0 updates the weight of the aircraft
"""
function Wupdate0!(parg, fuse, wing, rlx, fsum)
    WMTO = parg[igWMTO]
    

    ftotadd = sum(parg[[igfhpesys, igflgnose, igflgmain]])
    fsum = 0.0

    Wsum = parg[igWpay] +
           fuse.weight +
           wing.weight +
           wing.strut.weight +
           parg[igWhtail] +
           parg[igWvtail] +
           parg[igWeng] +
           parg[igWfuel] +
           parg[igWtesys] +
           parg[igWftank]

    WMTO = rlx * Wsum / (1.0 - ftotadd) + (1.0 - rlx) * WMTO
    parg[igWMTO] = WMTO

end


"""
Wupdate
"""
function Wupdate!(parg, fuse, wing, rlx, fsum)

    WMTO = parg[igWMTO]

    fwing = wing.weight / WMTO
    fstrut = wing.strut.weight / WMTO
    fhtail = parg[igWhtail] / WMTO
    fvtail = parg[igWvtail] / WMTO
    feng = parg[igWeng] / WMTO
    ffuel = parg[igWfuel] / WMTO
    fhpesys = parg[igfhpesys]
    flgnose = parg[igflgnose]
    flgmain = parg[igflgmain]

    ftesys = parg[igWtesys] / WMTO

    Wtesys = parg[igWtesys]
    Wftank = parg[igWftank]
    Wpay = parg[igWpay]
    Wfuse = fuse.weight

    ftank = parg[igWftank] / WMTO

    fsum = fwing + fstrut + fhtail + fvtail + feng + ffuel + fhpesys +
           flgnose + flgmain + ftank + ftesys

    if (fsum ≥ 1.0)
        println("SOMETHING IS WRONG fsum ≥ 1")
    end

    # WMTO = rlx*(Wpay + Wfuse + Wtesys + Wftank)/(1.0-fsum) + (1.0-rlx)*WMTO
    WMTO = rlx * (Wpay + fuse.weight) / (1.0 - fsum) + (1.0 - rlx) * WMTO

    parg[igWMTO] = WMTO
    wing.weight = WMTO * fwing
    wing.strut.weight = WMTO * fstrut
    parg[igWhtail] = WMTO * fhtail
    parg[igWvtail] = WMTO * fvtail
    parg[igWeng] = WMTO * feng
    parg[igWfuel] = WMTO * ffuel
   
    
    parg[igWftank] = WMTO * ftank 


    parg[igWtesys] = WMTO * ftesys


end

"""
    set_ambient_conditions!(ac, mis_point, Mach=NaN)

Sets ambient condition at the given mission point `mis_point`.
"""
function set_ambient_conditions!(ac, mis_point, Mach=NaN)
    mis_point = mis_point
    altkm = ac.parad[iaalt, mis_point]/1000.0
    T0, p0, ρ0, a0, μ0 = atmos(altkm)
    if Mach === NaN
        Mach = ac.parad[iaMach, mis_point]
    end
    ac.pared[iep0, mis_point] = p0
    ac.pared[ieT0, mis_point] = T0
    ac.pared[iea0, mis_point] = a0
    ac.pared[ierho0, mis_point] = ρ0
    ac.pared[iemu0, mis_point] = μ0
    ac.pared[ieM0, mis_point] = Mach
    ac.pared[ieu0, mis_point] = Mach * a0
    ac.parad[iaReunit, mis_point] = Mach * a0 * ρ0 / μ0

end  # function set_ambient_conditions
