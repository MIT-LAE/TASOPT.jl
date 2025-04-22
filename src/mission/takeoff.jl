#TODO: takeoff doc page needed; docstrings need updating
"""
    takeoff!(ac)

Calculates takeoff parameters and balanced field length.
The aircraft must be defined in parg array. The ipstatic and iprotate points are assumed to exist.

"""
function takeoff!(ac; printTO = true)
    parg  = ac.parg
    parm  = ac.parmd
    para  = ac.parad
    pare  = ac.pared  
    wing  = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    imission = 1

    #---- Newton convergence tolerance
    toler = 1.0e-7

    #---- zero-fuel weight for this mission
    Wzero = parg[igWMTO] -
            parg[igWfuel] -
            parg[igWpay] +
            parm[imWpay]

    #---- mission TO fuel weight
    WfTO = parm[imWTO] - Wzero

    #---- unpack parameters passed in via global data arrays parg,pare
    W = parm[imWTO]    # total takeoff weight
    S = wing.layout.S   # reference (wing) area
    sweep = wing.layout.sweep # sweep angle, degrees
    dfan = parg[igdfan]   # fan diameter , for engine-out CD_eng estimate
    HTRf = parg[igHTRf]   # hub/tip ratio, for engine-out CD_eng estimate
    neng = parg[igneng]   # number of engines
    muroll = parg[igmuroll]   # rolling resistance coefficient
    mubrake = parg[igmubrake]  # max braking resistance coefficient
    hobst = parg[ighobst]    # obstacle height (FAA spec, ~10m I think)

    Vstall = pare[ieu0, iprotate]
    V2 = pare[ieu0, iptakeoff]

    cosL = cosd(sweep)
    Afan = 0.25 * pi * dfan^2 * (1.0 - HTRf^2)
    CDgear = parg[igCDgear]
    CDeng = parg[igcdefan] * (0.25 * pi * dfan^2) / S
    CDspoiler = parg[igCDspoil]
    CDivert = 0.002

    Fmax = pare[ieFe, ipstatic]
    pare[ieFe, iptakeoff] = Fmax
    Fref = pare[ieFe, iprotate]

    #---- single-engine thrust-curve constants for takeoff roll calculations
    F01 = (Fmax + Fref) / 2.0
    KV1 = (Fmax - Fref) / Vstall^2

    #---- set max takeoff thrust for output
    FTO = Fmax * neng

    #---- normal takeoff
    ip = iprotate
    rho0 = pare[ierho0, ip]
    CLroll = 0.0
    para[iaCL, ip] = CLroll
    para[iaCLh, ip] = 0.0

    #cc      write(*,*) '^ 3a', Fmax, Fref

    #---- total CD during roll
    computes_wing_direct = false
    # iairf = 1
    aircraft_drag!(ac, imission, ip, computes_wing_direct)
    CDroll = para[iaCD, ip] + parg[igCDgear]

    #---- thrust constants for all engines operating
    KV = KV1 * neng
    F0 = F01 * neng

    #---- calculate takeoff length and takeoff time
    kA = (KV + rho0 * S * CDroll) / (W / gee)
    VAlimsq = 2.0 * (F0 - W * muroll) / (KV + rho0 * S * CDroll)

    #      write(*,*) 'F0 W muroll KV', F0, W, muroll, KV
    #      write(*,*) 'VAlimsq', VAlimsq

    VAlim = sqrt(VAlimsq)
    if (VAlim <= V2)
        println("Normal takeoff impossible. VAlim < V2:", VAlim, V2)
        lTO = 10.0 / kA
        tTO = 2.0 * lTO / V2
    else
        Vrat = V2 / VAlim
        lTO = -log(1.0 - Vrat^2) / kA
        tTO = log((1.0 + Vrat) / (1.0 - Vrat)) / (kA * VAlim)
    end

    #---- initial normal climbout sin[angle], gear up
    CDclimb = para[iaCD, ipclimb1]
    F2 = F0 - KV * 0.5 * V2^2
    singTO = (F2 - 0.5 * V2^2 * rho0 * S * CDclimb) / W
    singTO = max(0.01, min(0.99, singTO))

    #---- balanced field length

    #---- total CD during ground roll with one engine out
    CDb = CDroll + CDeng + CDivert

    #---- total CD during braking
    CDc = CDroll + CDeng * neng + CDspoiler

    #---- thrust constants with one engine out
    KV = KV1 * (neng - 1.0)
    F0 = F01 * (neng - 1.0)

    #---- velocity segment constants
    kB = (KV + rho0 * S * CDb) / (W / gee)
    kC = (rho0 * S * CDc) / (W / gee)

    VBlimsq = 2.0 * (F0 - W * muroll) / (KV + rho0 * S * CDb)
    VClimsq = 2.0 * (-W * mubrake) / (rho0 * S * CDc)

    VBlim = sqrt(VBlimsq)

    if (VBlim <= V2)
        println("Engine-out takeoff impossible. VBlim < V2:", VBlim, V2)
        l1 = 0.7 * lTO
        lBF = 10.0 * lTO
        V1 = Vstall
        V2 = Vstall
        singBF = 0.0

        #---- calculate fuel burn during takeoff run
        ip1 = ipstatic
        ip2 = iprotate
        mdotf1 = pare[iemcore, ip1] * pare[ieff, ip1] * neng
        mdotf2 = pare[iemcore, ip2] * pare[ieff, ip2] * neng
        WfTO = 0.5 * (mdotf1 + mdotf2) * tTO * gee
        para[iafracW, ipstatic] = para[iafracW, iptakeoff] +
                                  WfTO / parg[igWMTO]
        para[iatime, ipstatic] = para[iatime, iptakeoff] - tTO
        para[iaRange, ipstatic] = para[iaRange, iptakeoff] - lTO

        #---- store output quantities for returning
        parm[imV1] = V1
        parm[imV2] = V2
        parm[imlTO] = lTO
        parm[iml1] = l1
        parm[imlBF] = lBF
        parm[imtTO] = tTO
        parm[imFTO] = FTO
        parm[imgamVTO] = asin[singTO]
        parm[imgamVBF] = asin[singBF]

        return
    end

    #---- initial guesses for Newton iteration
    l1 = 0.8 * lTO
    lBF = 1.3 * lTO
    V2sq = V2^2

    if printTO
        @printf("\nTakeoff:\n%2s %10s %10s %10s %10s\n", 
        "#", "lTO", "l1", "lBF", "dmax")
    end
    
    #---- Newton iteration loop
    for iter = 1:15
        exA = exp(kA * (-l1))
        exB = exp(kB * (lBF - l1))
        exC = exp(kC * (lBF - l1))
        exA_l1 = -kA * exA
        exB_l1 = -kB * exB
        exC_l1 = -kC * exC
        exA_lBF = 0.0
        exB_lBF = kB * exB
        exC_lBF = kC * exC

        r1 = VAlimsq * (1.0 - exA) + (VBlimsq - V2sq) * exB - VBlimsq
        a11 = VAlimsq * (-exA_l1) + (VBlimsq - V2sq) * exB_l1
        a12 = VAlimsq * (-exA_lBF) + (VBlimsq - V2sq) * exB_lBF

        r2 = VAlimsq * (1.0 - exA) - VClimsq * (1.0 - exC)
        a21 = VAlimsq * (-exA_l1) - VClimsq * (-exC_l1)
        a22 = VAlimsq * (-exA_lBF) - VClimsq * (-exC_lBF)

        det = a11 * a22 - a12 * a21
        dl1 = -(r1 * a22 - a12 * r2) / det
        dlBF = -(a11 * r2 - r1 * a21) / det

        dmax = max(abs(dl1), abs(dlBF))

        #  print convergence history for debugging
        if printTO
            @printf("%2d %10.3f %10.3f %10.3f %10.3f\n", 
            iter, lTO * 3.28, l1 * 3.28, lBF * 3.28, dmax * 3.28)
        end

        l1 = l1 + dl1
        lBF = lBF + dlBF

        if (dmax < toler * lTO)

            V1 = VAlim * sqrt(1.0 - exp(-kA * l1))

            #---- engine-out climbout, gear down
            F2 = F0 - KV * 0.5 * V2^2
            CD = CDclimb + CDgear + CDeng
            #cc   singBF = (F2 - 0.5*V2^2 * rho0*S*CDb) / W
            singBF = (F2 - 0.5 * V2^2 * rho0 * S * CD) / W
            singBF = max(0.01, min(0.99, singBF))

            #---- calculate fuel burn during takeoff run
            ip1 = ipstatic
            ip2 = iprotate
            mdotf1 = pare[iemcore, ip1] * pare[ieff, ip1] * neng
            mdotf2 = pare[iemcore, ip2] * pare[ieff, ip2] * neng
            WfTO = 0.5 * (mdotf1 + mdotf2) * tTO * gee
            para[iafracW, ipstatic] = para[iafracW, iptakeoff] +
                                      WfTO / parg[igWMTO]
            para[iatime, ipstatic] = para[iatime, iptakeoff] - tTO
            para[iaRange, ipstatic] = para[iaRange, iptakeoff] - lTO

            #---- store output quantities for returning
            parm[imV1] = V1
            parm[imV2] = V2
            parm[imlTO] = lTO
            parm[iml1] = l1
            parm[imlBF] = lBF
            parm[imtTO] = tTO
            parm[imFTO] = FTO
            parm[imgamVTO] = asin(singTO)
            parm[imgamVBF] = asin(singBF)

            return
        end
    end
    println("TAKEOFF: BFL convergence failed.  dl1,dlBF=", dl1, dlBF)
    println("         Continuing anyway...")

end

