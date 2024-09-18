function ductedfancalc!(pari, parg, para, pare, ip,
    icall, icool, initeng)

    Lprint = false

    if (Lprint)
            println("entering TFCALC", icall, icool, initeng)
    end

    iBLIc = pari[iiBLIc]

    neng = parg[igneng]
    S = parg[igS]

    Fsp = pare[ieFsp]
    pif = pare[iepif]
    pid = pare[iepid]
    pifn = pare[iepifn]
    
    epolf = pare[ieepolf]
    
    pifK = pare[iepifK]
    epfK = pare[ieepfK]
    M2 = pare[ieM2]
    M0 = pare[ieM0]
    Tt0 = pare[ieTt0]
    ht0 = pare[ieht0]
    pt0 = pare[iept0]
    cpt0 = pare[iecpt0]
    Rt0 = pare[ieRt0]
    p0 = pare[iep0]
    a0 = pare[iea0]
    rho0 = pare[ierho0]
    mu0 = pare[iemu0]
    T0 = pare[ieT0]
    u0 = pare[ieu0]
    
    #Heat exchanger variables
    Δh_PreC = pare[iePreCDeltah]
    Δh_InterC = pare[ieInterCDeltah]
    Δh_Regen = pare[ieRegenDeltah]
    Δh_TurbC = pare[ieTurbCDeltah]
    Δp_PreC = pare[iePreCDeltap]
    Δp_InterC = pare[ieInterCDeltap]
    Δp_Regen = pare[ieRegenDeltap]

    #- - - - - - - - - - - - - - - - - - - - - - - 
    #---- set BL ingestion parameters
    if (M0 == 0.0)
            #----- no ingestion for static case
            Phiinl = 0.0
            Kinl = 0.0

    else
            #----- assume engine is at TE of fuselage
            DAfsurf = para[iaDAfsurf]
            KAfTE = para[iaKAfTE]
            fBLIf = parg[igfBLIf]

            #----- assume 85% of wing dissipation is on surface
            fDwake = 0.15
            CDAwing = para[iaCDwing] * S
            DAwsurf = CDAwing * (1.0 - fDwake)
            KAwTE = DAwsurf

            fBLIw = parg[igfBLIw]

            #----- set ingested  PKinl-PVinl = Phiinl  for one engine
            Phiinl = 0.5 * rho0 * u0^3 * (DAfsurf * fBLIf + DAwsurf * fBLIw) / neng
            Kinl = 0.5 * rho0 * u0^3 * (KAfTE * fBLIf + KAwTE * fBLIw) / neng
    end
    #- - - - - - - - - - - - - - - - - - - - - - - 

    # #--------------------------------------------------------------------------
    if (icall == 0)
        #----- engine sizing case

        Fe = pare[ieFe]

        Fsp, Pfan,
        Tt0, ht0, pt0, cpt0, Rt0,
        Tt18, ht18, pt18, cpt18, Rt18,
        Tt2, ht2, pt2, cpt2, Rt2,
        Tt21, ht21, pt21, cpt21, Rt21,
        Tt7, ht7, pt7, cpt7, Rt7,
        u0,
        T2, u2, p2, cp2, R2, A2,
        T7, u7, p7, cp7, R7, A7,
        T8, u8, p8, cp8, R8, A8,
        epf,
        etaf, Lconv = ductedfansize!(gee, M0, T0, p0, a0, M2,
                        Fe, Phiinl, Kinl, iBLIc,
                        pif,
                        pid, pifn, 
                        epolf,
                        pifK, epfK
                        )
    else
        #----- fixed parameters
        A2 = pare[ieA2]
        A7 = pare[ieA7]

        mbfD = pare[iembfD]

        pifD = pare[iepifD]

        if (initeng == 0)
                #------ force TFOPER to initialize these state variables
                mbf = 0.0
                pif = 0.0
                M2 = 1.0
        else
                #------ use existing state variables as initial guesses
                mbf = pare[iembf]
                pif = pare[iepif]
                M2 = pare[ieM2]
        end
        if (icall == 1) #Power is specified, thrust will be computed
            Feng = 0.0
            Peng = 1e6
            iPspec = true

            TSEC, Fsp, Feng, Peng, mfan, 
            pif, mbf, 
            Tt0, ht0, pt0, cpt0, Rt0,
            Tt18, ht18, pt18, cpt18, Rt18,
            Tt2, ht2, pt2, cpt2, Rt2,
            Tt21, ht21, pt21, cpt21, Rt21,
            Tt7, ht7, pt7, cpt7, Rt7,
            u0, T2, u2, p2, cp2, R2, M2,
            T7, u7, p7, cp7, R7, M7,
            T8, u8, p8, cp8, R8, M8, A8,
            epf, etaf = ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                            Phiinl, Kinl, iBLIc,
                            pid, pifn, 
                            pifD, 
                            mbfD, 
                            A2, A7,
                            epolf,
                            pifK, epfK,
                            Feng, Peng,
                            M2, pif, mbf, 
                            iPspec)
        elseif (icall == 2) #Thrust is specified, power to be computed
            Feng = 1e6
            Peng = 0.0
            iPspec = false

            TSEC, Fsp, Feng, Peng, mfan, 
            pif, mbf, 
            Tt0, ht0, pt0, cpt0, Rt0,
            Tt18, ht18, pt18, cpt18, Rt18,
            Tt2, ht2, pt2, cpt2, Rt2,
            Tt21, ht21, pt21, cpt21, Rt21,
            Tt7, ht7, pt7, cpt7, Rt7,
            u0, T2, u2, p2, cp2, R2, M2,
            T7, u7, p7, cp7, R7, M7,
            T8, u8, p8, cp8, R8, M8, A8,
            epf, etaf = ductedfanoper!(M0, T0, p0, a0, Tref, pref,
                            Phiinl, Kinl, iBLIc,
                            pid, pifn, 
                            pifD, 
                            mbfD, 
                            A2, A7,
                            epolf,
                            pifK, epfK,
                            Feng, Peng,
                            M2, pif, mbf, 
                            iPspec)
        end
    end
end