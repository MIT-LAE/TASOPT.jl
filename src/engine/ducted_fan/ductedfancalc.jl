""" 
    ductedfancalc!(ac, case::String, imission::Int64, ip::Int64, initializes_engine::Bool, iterw::Int64 = 0)

Calls function ductedfansize! or ductedfanoper! for one operating point.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object
    - `case::String`: case identifier, e.g. "sizing" or "off_design"
    - `imission::Int64`: mission index
    - `ip::Int64`: mission point index
    - `initializes_engine::Bool`: flag to initialize engine
      - `true`: initialize variables for iteration in engine
      - `false`: use current variables as initial guesses in engine
    - `iterw::Int64`: sizing loop iteration

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine parameters.
"""
function ductedfancalc!(ac, case::String, imission::Int64, ip::Int64, initializes_engine::Bool, iterw::Int64 = 0)
    #Unpack data storage arrays
    parg, _, para, pare, options, _, _, wing, _, _, eng, _ = unpack_ac(ac, imission, ip=ip)
    iBLIc = eng.model.has_BLI_cores

    neng = parg[igneng]
    S = wing.layout.S

    Fsp = pare[ieFsp]
    pif = pare[iepif]
    pid = pare[iepid]
    pifn = pare[iepifn]
    
    epolf = pare[ieepolf]
    
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
    Δh_radiator = pare[ieRadiatorDeltah]
    Δp_radiator = pare[ieRadiatorDeltap]

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
    if (case == "design")
        #----- engine sizing case

        Fe = pare[ieFe] #ducted fan sized for a given thrust

        TSEC, Fsp, Pfan, mfan,
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
                        Δh_radiator, Δp_radiator
                        )

        mbf = mfan * sqrt(Tt2 / Tref) / (pt2 / pref)
        M7 = u7 / sqrt(T7 * R7 * cp7 / (cp7 - R7))

        pifD = pif
        mbfD = mbf

        Nf = 1.0 #Arbitrarily set to 1 as only ratios matter
        Nbf = Nf / sqrt(Tt2 / Tref)
        NbfD = Nbf

        #----- store design-point parameters
        pare[ieA2] = A2
        pare[ieA7] = A7

        pare[iembfD] = mbfD
        pare[iepifD] = pifD
        pare[ieNbfD] = Nbf

    else
        #----- fixed parameters
        A2 = pare[ieA2]
        A7 = pare[ieA7]

        mbfD = pare[iembfD]
        pifD = pare[iepifD]
        NbfD = pare[ieNbfD]
        Δh_radiator = pare[ieRadiatorDeltah]
        Δp_radiator = pare[ieRadiatorDeltap]

        if initializes_engine
                #------ initialize these state variables
                mbf = 0.0
                pif = 0.0
                M2 = 1.0
        else
                #------ use existing state variables as initial guesses
                mbf = pare[iembf]
                pif = pare[iepif]
                M2 = pare[ieM2]
        end
        if ip in range(ipstatic, iprotate) 
            #Power is specified, thrust will be computed -- in takeoff
            Feng = 0.0
            Peng = pare[iePfanmax]
            iPspec = true

            TSEC, Fsp, Feng, Pfan, mfan, 
            pif, mbf, Nbf,
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
                            mbfD, NbfD,
                            A2, A7,
                            epolf,
                            Feng, Peng,
                            M2, pif, mbf, 
                            Δh_radiator, Δp_radiator,
                            iPspec)
            pare[ieFe] = Feng

        else #Thrust is specified, power to be computed
            Feng = pare[ieFe]
            Peng = 0.0
            iPspec = false

            TSEC, Fsp, Feng, Pfan, mfan, 
            pif, mbf, Nbf,
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
                            mbfD, NbfD,
                            A2, A7,
                            epolf,
                            Feng, Peng,
                            M2, pif, mbf, 
                            Δh_radiator, Δp_radiator,
                            iPspec)
        end
    end
    HTRf = parg[igHTRf]
    dfan = sqrt(4.0 * A2 / (pi * (1.0 - HTRf^2)))
    parg[igdfan] = dfan

    pare[iePfan] = Pfan
    pare[ieTSEC] = TSEC
    pare[ieTSFC] = TSEC / 120e6  #TODO change when fuel is not burnt
    pare[ieFsp] = Fsp
    pare[iemfan] = mfan
    pare[iePhiinl] = Phiinl
    pare[ieKinl] = Kinl

    pare[iembf] = mbf
    pare[iepif] = pif
    pare[ieNbf] = Nbf
    pare[ieNf] = Nbf * sqrt(Tt2 / Tref)

    pare[ieTt0] = Tt0
    pare[ieht0] = ht0
    pare[iept0] = pt0
    pare[iecpt0] = cpt0
    pare[ieRt0] = Rt0

    pare[ieTt18] = Tt18
    pare[ieht18] = ht18
    pare[iept18] = pt18
    pare[iecpt18] = cpt18
    pare[ieRt18] = Rt18

    pare[ieTt2] = Tt2
    pare[ieht2] = ht2
    pare[iept2] = pt2
    pare[iecpt2] = cpt2
    pare[ieRt2] = Rt2

    pare[ieTt21] = Tt21
    pare[ieht21] = ht21
    pare[iept21] = pt21
    pare[iecpt21] = cpt21
    pare[ieRt21] = Rt21

    pare[ieTt7] = Tt7
    pare[ieht7] = ht7
    pare[iept7] = pt7
    pare[iecpt7] = cpt7
    pare[ieRt7] = Rt7

    pare[ieu0] = u0

    pare[iep2] = p2
    pare[ieT2] = T2
    pare[ieR2] = R2
    pare[iecp2] = cp2
    pare[ieu2] = u2

    pare[iep7] = p7
    pare[ieT7] = T7
    pare[ieR7] = R7
    pare[iecp7] = cp7
    pare[ieu7] = u7

    pare[iep8] = p8
    pare[ieT8] = T8
    pare[ieR8] = R8
    pare[iecp8] = cp8
    pare[ieu8] = u8

    pare[ieA8] = A8

    pare[ieepf] = epf

    pare[ieetaf] = etaf


    if (M7 <= 0.999999)
            ichoke7 = 0
    else
            ichoke7 = 1
    end

end