""" 
    tfcalc!(wing, engine, parg, para, pare, ip, ifuel, opt_calc_call, opt_cooling, initializes_engine)

Calls on-design sizing function [`tfsize!`](@ref) or off-design analysis function
[`tfoper!`](@ref) for one operating point `ip`.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `opt_calc_call`:
      - "sizing": call on-design sizing routine `tfsize!`
      - `"oper_fixedTt4"`: call off-design analysis routine `tfoper!` with specified Tt4
      - `"oper_fixedFe"`: call off-design analysis routine `tfoper!` with specified
        net thrust (`Fe`)

    - `opt_cooling`: turbine cooling flag
      - "none": no cooling mass flow
      - `"fixed_coolingflowratio"`: use specified cooling flow ratios `epsrow(.)`;
        calculate `Tmrow(.)`
      - `"fixed_Tmetal"`: use specified metal temperatures `Tmrow(.)`; calculate
        `epsrow(.)`
    - `initializes_engine`:
      - `true`: initialize variables for iteration in `tfoper!`
      - `false`: use current variables as initial guesses in `tfoper!`
"""
function tfcalc!(wing, engine, parg::Vector{Float64}, para, pare, ip::Int64, ifuel::Int64, 
        opt_calc_call::String, opt_cooling::String, initializes_engine::Bool)

        Lprint = false

        if (Lprint)
                println("entering TFCALC", opt_calc_call, opt_cooling, initializes_engine)
        end

        eng_has_BLI_cores = engine.model.has_BLI_cores

        Gearf = parg[igGearf]
        Tmetal = parg[igTmetal]
        neng = parg[igneng]

        mofWpay = parg[igmofWpay]
        mofWMTO = parg[igmofWMTO]
        PofWpay = parg[igPofWpay]
        PofWMTO = parg[igPofWMTO]

        Wpay = parg[igWpay]
        WMTO = parg[igWMTO]

        TSFC = pare[ieTSFC]
        Fsp = pare[ieFsp]
        hfuel = pare[iehfuel]
        Tfuel = pare[ieTfuel]
        Tt4 = pare[ieTt4]
        BPR = pare[ieBPR]
        pif = pare[iepif]
        pilc = pare[iepilc]
        pihc = pare[iepihc]
        pid = pare[iepid]
        pib = pare[iepib]
        pifn = pare[iepifn]
        pitn = pare[iepitn]
        epolf = pare[ieepolf]
        epollc = pare[ieepollc]
        epolhc = pare[ieepolhc]
        epolht = pare[ieepolht]
        epollt = pare[ieepollt]
        etab = pare[ieetab]
        M2 = pare[ieM2]
        M25 = pare[ieM25]
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
        mcore = pare[iemcore]
        dTstrk = pare[iedTstrk]
        StA = pare[ieStA]
        Mtexit = pare[ieMtexit]
        M4a = pare[ieM4a]
        ruc = pare[ieruc]
        efilm = pare[ieefilm]
        tfilm = pare[ietfilm]
        epsl = pare[ieepsl]
        epsh = pare[ieepsh]
        epsrow = zeros(ncrowx)
        Tmrow = zeros(ncrowx)
        hvap = pare[iehvapcombustor] #Enthalpy of vaporization of fuel

        #Effect of cooling on HPT efficiency
        epht_fc = pare[iedehtdfc]
        fc0 = fc0 = pare[iefc0]

        #Heat exchanger variables
        Δh_PreC = pare[iePreCDeltah]
        Δh_InterC = pare[ieInterCDeltah]
        Δh_Regen = pare[ieRegenDeltah]
        Δh_TurbC = pare[ieTurbCDeltah]
        Δp_PreC = pare[iePreCDeltap]
        Δp_InterC = pare[ieInterCDeltap]
        Δp_Regen = pare[ieRegenDeltap]

        if compare_strings(opt_cooling, "fixed_coolingflowratio")
                ncrow = ncrowx
                for icrow = 1:ncrowx
                        epsrow[icrow] = pare[ieepsc1+icrow-1]
                end
        elseif compare_strings(opt_cooling, "fixed_Tmetal")
                ncrow = ncrowx
                for icrow = 1:ncrowx
                        #cc      Tmrow[icrow]  = pare(ieTmet1+icrow-1)
                        Tmrow[icrow] = parg[igTmetal]
                end
        end

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
                CDAwing = para[iaCDwing] * wing.layout.S
                DAwsurf = CDAwing * (1.0 - fDwake)
                KAwTE = DAwsurf

                fBLIw = parg[igfBLIw]

                #----- set ingested  PKinl-PVinl = Phiinl  for one engine
                Phiinl = 0.5 * rho0 * u0^3 * (DAfsurf * fBLIf + DAwsurf * fBLIw) / neng
                Kinl = 0.5 * rho0 * u0^3 * (KAfTE * fBLIf + KAwTE * fBLIw) / neng
        end
        #- - - - - - - - - - - - - - - - - - - - - - - 

        #---- mass and power offtakes
        Pofft_HX = pare[ieHXrecircP] #power offtakes to drive heat exchanger recirculation per engine
        mofft = (mofWpay * Wpay + mofWMTO * WMTO) / neng
        Pofft = (PofWpay * Wpay + PofWMTO * WMTO) / neng + Pofft_HX

        Tt9 = pare[ieTt9]
        pt9 = pare[iept9]

        #--------------------------------------------------------------------------
        #Engine model convergence
        pare[ieConvFail] = 0.0 #Converged by default

        # #--------------------------------------------------------------------------
        if compare_strings(opt_calc_call, "sizing")
                #----- engine sizing case

                Fe = pare[ieFe]

                if (Lprint)
                        println("TFSIZE  M0 p0 =", M0, p0)
                        println("     pif pilc =", pif, pilc)
                        println("        Tt4 F =", Tt4, Fe)
                        println("ncrow =", ncrow)
                        println("        altkm =", para[iaalt] / 1000.0)
                end

                epsrow, Tmrow,
                TSFC, Fsp, hfuel, ff, mcore,
                Tt0, ht0, pt0, cpt0, Rt0,
                Tt18, ht18, pt18, cpt18, Rt18,
                Tt19, ht19, pt19, cpt19, Rt19,
                Tt19c, ht19c, pt19c, cpt19c, Rt19c,
                Tt2, ht2, pt2, cpt2, Rt2,
                Tt21, ht21, pt21, cpt21, Rt21,
                Tt25, ht25, pt25, cpt25, Rt25,
                Tt25c, ht25c, pt25c, cpt25c, Rt25c,
                Tt3, ht3, pt3, cpt3, Rt3,
                ht4, pt4, cpt4, Rt4,
                Tt41, ht41, pt41, cpt41, Rt41,
                Tt45, ht45, pt45, cpt45, Rt45,
                Tt49, ht49, pt49, cpt49, Rt49,
                Tt5, ht5, pt5, cpt5, Rt5,
                Tt7, ht7, pt7, cpt7, Rt7,
                u0,
                T2, u2, p2, cp2, R2, A2,
                T25, u25, p25, cp25, R25, A25,
                T5, u5, p5, cp5, R5, A5,
                T6, u6, p6, cp6, R6, A6,
                T7, u7, p7, cp7, R7, A7,
                T8, u8, p8, cp8, R8, A8,
                u9, A9,
                epf, eplc, ephc, epht, eplt,
                etaf, etalc, etahc, etaht, etalt,
                Lconv = tfsize!(gee, M0, T0, p0, a0, M2, M25,
                        Fe, Phiinl, Kinl, eng_has_BLI_cores,
                        BPR, pif, pilc, pihc,
                        pid, pib, pifn, pitn,
                        Tfuel, ifuel, hvap, etab,
                        epolf, epollc, epolhc, epolht, epollt,
                        mofft, Pofft,
                        Tt9, pt9, Tt4,
                        epsl, epsh,
                        opt_cooling,
                        Mtexit, dTstrk, StA, efilm, tfilm,
                        fc0, epht_fc,
                        M4a, ruc,
                        ncrowx, ncrow,
                        epsrow, Tmrow,
                        Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                        Δp_PreC, Δp_InterC, Δp_Regen)

                #        tadd(time0,t_tfsize)

                M5 = u5 / sqrt(T5 * R5 * cp5 / (cp5 - R5))
                M6 = u6 / sqrt(T6 * R6 * cp6 / (cp6 - R6))
                M7 = u7 / sqrt(T7 * R7 * cp7 / (cp7 - R7))
                M8 = u8 / sqrt(T8 * R8 * cp8 / (cp8 - R8))

                #      pare(ieM5 ) = M5   
                #      pare(ieM6 ) = M6   
                #      pare(ieM7 ) = M7   
                #      pare(ieM8 ) = M8   

                fo = mofft / mcore

                if (Lprint)
                        println("exited TFSIZE", fo)
                end

                #----- corrected mass flows
                mbf = mcore * sqrt(Tt2 / Tref) / (pt2 / pref) * BPR
                mblc = mcore * sqrt(Tt19c / Tref) / (pt19c / pref)
                mbhc = mcore * sqrt(Tt25c / Tref) / (pt25c / pref) * (1.0 - fo)
                mbht = mcore * sqrt(Tt41 / Tref) / (pt41 / pref) * (1.0 - fo + ff)
                mblt = mcore * sqrt(Tt45 / Tref) / (pt45 / pref) * (1.0 - fo + ff)

                #----- spool speed fractions for design case are unity by definition
                Nf = 1.0 / Gearf
                N1 = 1.0
                N2 = 1.0

                Nbf = Nf / sqrt(Tt2 / Tref)
                Nblc = N1 / sqrt(Tt19c / Tref)
                Nbhc = N2 / sqrt(Tt25c / Tref)
                Nbht = N2 / sqrt(Tt41 / Tref)
                Nblt = N1 / sqrt(Tt45 / Tref)

                #----- set quantities fixed by this design case for all operating points
                mbfD = mbf
                mblcD = mblc
                mbhcD = mbhc
                mbhtD = mbht
                mbltD = mblt

                pifD = pif
                pilcD = pilc
                pihcD = pihc
                pihtD = pt41 / pt45
                piltD = pt45 / pt49

                NbfD = Nbf
                NblcD = Nblc
                NbhcD = Nbhc
                NbhtD = Nbht
                NbltD = Nblt

                #----- but recalculate turbine pressure ratios using slightly approximate form,
                #-      to be fuly consistent with TFOPER's turbine efficiency function
                Trh = Tt41 / (Tt41 + (ht45 - ht41) / cpt41)
                Trl = Tt45 / (Tt45 + (ht49 - ht45) / cpt45)
                gexh = cpt41 / (Rt41 * epolht)
                gexl = cpt45 / (Rt45 * epollt)
                pihtD = Trh^gexh
                piltD = Trl^gexl

                #----- store design-point parameters
                pare[ieA2] = A2
                pare[ieA25] = A25
                pare[ieA5] = A5
                pare[ieA7] = A7

                pare[ieNbfD] = NbfD
                pare[ieNblcD] = NblcD
                pare[ieNbhcD] = NbhcD
                pare[ieNbhtD] = NbhtD
                pare[ieNbltD] = NbltD

                pare[iembfD] = mbfD
                pare[iemblcD] = mblcD
                pare[iembhcD] = mbhcD
                pare[iembhtD] = mbhtD
                pare[iembltD] = mbltD

                pare[iepifD] = pifD
                pare[iepilcD] = pilcD
                pare[iepihcD] = pihcD
                pare[iepihtD] = pihtD
                pare[iepiltD] = piltD

                #Fuel mass flow rate
                pare[iemfuel] = ff * mcore * neng

                HTRf = parg[igHTRf]
                HTRlc = parg[igHTRlc]
                HTRhc = parg[igHTRhc]
                Alc = A2 / (1.0 + BPR)
                dfan = sqrt(4.0 * A2 / (pi * (1.0 - HTRf^2)))
                dlcomp = sqrt(4.0 * Alc / (pi * (1.0 - HTRlc^2)))
                dhcomp = sqrt(4.0 * A25 / (pi * (1.0 - HTRhc^2)))
                parg[igdfan] = dfan
                parg[igdlcomp] = dlcomp
                parg[igdhcomp] = dhcomp

                #--------------------------------------------------------------------------
        else
                #----- off-design operation case

                #----- fixed parameters
                A2 = pare[ieA2]
                A25 = pare[ieA25]
                A5 = pare[ieA5]
                A7 = pare[ieA7]

                NbfD = pare[ieNbfD]
                NblcD = pare[ieNblcD]
                NbhcD = pare[ieNbhcD]
                NbhtD = pare[ieNbhtD]
                NbltD = pare[ieNbltD]

                mbfD = pare[iembfD]
                mblcD = pare[iemblcD]
                mbhcD = pare[iembhcD]
                mbhtD = pare[iembhtD]
                mbltD = pare[iembltD]

                pifD = pare[iepifD]
                pilcD = pare[iepilcD]
                pihcD = pare[iepihcD]
                pihtD = pare[iepihtD]
                piltD = pare[iepiltD]

                if (initializes_engine)
                        #------ force TFOPER to initialize these state variables
                        mbf = 0.0
                        mblc = 0.0
                        mbhc = 0.0
                        pif = 0.0
                        pilc = 0.0
                        pihc = 0.0
                        pt5 = 0.0
                        M2 = 1.0
                        M25 = 1.0
                else
                        #------ use existing state variables as initial guesses
                        mbf = pare[iembf]
                        mblc = pare[iemblc]
                        mbhc = pare[iembhc]
                        pif = max(pare[iepif], 1.1)
                        pilc = max(pare[iepilc], 1.1)
                        pihc = max(pare[iepihc], 1.1)
                        pt5 = pare[iept5]
                        M2 = pare[ieM2]
                        M25 = pare[ieM25]
                end

                Fe = 0.0
                Tt4 = pare[ieTt4]

                if compare_strings(opt_calc_call, "oper_fixedTt4")
                        #------ specified Tt4 -- Fe will be computed
                        nothing; #nothing special is done
                elseif compare_strings(opt_calc_call, "oper_fixedFe")
                        #------ specified Fe -- Tt4 will be computed (set initial guess here)
                        Fe = pare[ieFe]
                end

                if (Lprint)
                        println(cplab[ip], opt_calc_call, Tt4, Fe)
                        println("Calling TFOPER...", opt_calc_call, opt_cooling, ip)
                        println(DAwsurf, para[iagamV])
                        println(rho0, u0, parg[igWMTO])
                        println(mcore, M2, M25)
                        println("Phiinl, Kinl", Phiinl, Kinl)

                end

                TSFC, Fsp, hfuel, ff,
                Fe, mcore,
                pif, pilc, pihc,
                mbf, mblc, mbhc,
                Nbf, Nblc, Nbhc,
                Tt0, ht0, pt0, cpt0, Rt0,
                Tt18, ht18, pt18, cpt18, Rt18,
                Tt19, ht19, pt19, cpt19, Rt19,
                Tt19c, ht19c, pt19c, cpt19c, Rt19c,
                Tt2, ht2, pt2, cpt2, Rt2,
                Tt21, ht21, pt21, cpt21, Rt21,
                Tt25, ht25, pt25, cpt25, Rt25,
                Tt25c, ht25c, pt25c, cpt25c, Rt25c,
                Tt3, ht3, pt3, cpt3, Rt3,
                Tt4, ht4, pt4, cpt4, Rt4,
                Tt41, ht41, pt41, cpt41, Rt41,
                Tt45, ht45, pt45, cpt45, Rt45,
                Tt49, ht49, pt49, cpt49, Rt49,
                Tt5, ht5, pt5, cpt5, Rt5,
                Tt7, ht7, pt7, cpt7, Rt7,
                u0,
                T2, u2, p2, cp2, R2, M2,
                T25, u25, p25, cp25, R25, M25,
                T5, u5, p5, cp5, R5, M5,
                T6, u6, p6, cp6, R6, M6, A6,
                T7, u7, p7, cp7, R7, M7,
                T8, u8, p8, cp8, R8, M8, A8,
                u9, A9,
                epf, eplc, ephc, epht, eplt,
                etaf, etalc, etahc, etaht, etalt,
                Lconv = tfoper!(gee, M0, T0, p0, a0, Tref, pref,
                        Phiinl, Kinl, eng_has_BLI_cores,
                        pid, pib, pifn, pitn,
                        Gearf,
                        pifD, pilcD, pihcD, pihtD, piltD,
                        mbfD, mblcD, mbhcD, mbhtD, mbltD,
                        NbfD, NblcD, NbhcD, NbhtD, NbltD,
                        A2, A25, A5, A7,
                        opt_calc_call,
                        Tfuel, ifuel, hvap, etab,
                        epolf, epollc, epolhc, epolht, epollt,
                        mofft, Pofft,
                        Tt9, pt9,
                        epsl, epsh,
                        opt_cooling,
                        Mtexit, dTstrk, StA, efilm, tfilm,
                        fc0, epht_fc,
                        M4a, ruc,
                        ncrowx, ncrow,
                        epsrow, Tmrow,
                        Fe,
                        M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
                        Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                        Δp_PreC, Δp_InterC, Δp_Regen)


                if (Lprint)
                        println("exited TFOPER", Lconv)
                end

                if (!Lconv)
                        #@warn "Convergence failed on operating point: $ip"
                        pare[ieConvFail] = 1.0 #Store convergence failure
                end

                fo = mofft / mcore

                Nf = Nbf * sqrt(Tt2 / Tref)
                N1 = Nblc * sqrt(Tt19c / Tref)
                N2 = Nbhc * sqrt(Tt25c / Tref)

                pare[ieM2] = M2
                pare[ieM25] = M25

                if compare_strings(opt_calc_call, "oper_fixedTt4")
                        pare[ieFe] = Fe
                elseif compare_strings(opt_calc_call, "oper_fixedFe")
                        pare[ieTt4] = Tt4
                end

                pare[ieBPR] = mbf / mblc * sqrt(Tt19c / Tt2) * pt2 / pt19c

                # println("exited TFOPER call")

        end
        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        if compare_strings(opt_cooling, "fixed_coolingflowratio")
                #------ cooling flow was specified... set metal temperatures
                for icrow = 1:ncrowx
                        pare[ieTmet1+icrow-1] = Tmrow[icrow]
                end

        elseif compare_strings(opt_cooling, "fixed_Tmetal")
                #------ Tmetal was specified... set cooling flow ratios
                epstot = 0.0
                for icrow = 1:ncrowx
                        epstot = epstot + epsrow[icrow]
                        pare[ieepsc1+icrow-1] = epsrow[icrow]
                end
                pare[iefc] = (1.0 - fo) * epstot

        end

        pare[ieTSFC] = TSFC
        pare[ieFsp] = Fsp
        pare[iehfuel] = hfuel
        pare[ieff] = ff
        pare[iemcore] = mcore
        pare[iemofft] = mofft
        pare[iePofft] = Pofft
        pare[iePhiinl] = Phiinl
        pare[ieKinl] = Kinl

        pare[ieNf] = Nbf * sqrt(Tt2 / Tref)
        pare[ieN1] = Nblc * sqrt(Tt19c / Tref)
        pare[ieN2] = Nbhc * sqrt(Tt25c / Tref)
        pare[ieNbf] = Nbf
        pare[ieNblc] = Nblc
        pare[ieNbhc] = Nbhc
        pare[iembf] = mbf
        pare[iemblc] = mblc
        pare[iembhc] = mbhc
        pare[iepif] = pif
        pare[iepilc] = pilc
        pare[iepihc] = pihc

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

        pare[ieTt19] = Tt19
        pare[ieht19] = ht19
        pare[iept19] = pt19
        pare[iecpt19] = cpt19
        pare[ieRt19] = Rt19

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

        pare[ieTt25] = Tt25
        pare[ieht25] = ht25
        pare[iept25] = pt25
        pare[iecpt25] = cpt25
        pare[ieRt25] = Rt25

        pare[ieTt3] = Tt3
        pare[ieht3] = ht3
        pare[iept3] = pt3
        pare[iecpt3] = cpt3
        pare[ieRt3] = Rt3

        #cc   pare[ieTt4 ] = Tt4
        pare[ieht4] = ht4
        pare[iept4] = pt4
        pare[iecpt4] = cpt4
        pare[ieRt4] = Rt4

        pare[ieTt41] = Tt41
        pare[ieht41] = ht41
        pare[iept41] = pt41
        pare[iecpt41] = cpt41
        pare[ieRt41] = Rt41

        pare[ieTt45] = Tt45
        pare[ieht45] = ht45
        pare[iept45] = pt45
        pare[iecpt45] = cpt45
        pare[ieRt45] = Rt45

        pare[ieTt49] = Tt49
        pare[ieht49] = ht49
        pare[iept49] = pt49
        pare[iecpt49] = cpt49
        pare[ieRt49] = Rt49

        pare[ieTt5] = Tt5
        pare[ieht5] = ht5
        pare[iept5] = pt5
        pare[iecpt5] = cpt5
        pare[ieRt5] = Rt5

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

        pare[iep25] = p25
        pare[ieT25] = T25
        pare[ieR25] = R25
        pare[iecp25] = cp25
        pare[ieu25] = u25

        pare[iep5] = p5
        pare[ieT5] = T5
        pare[ieR5] = R5
        pare[iecp5] = cp5
        pare[ieu5] = u5

        pare[iep6] = p6
        pare[ieT6] = T6
        pare[ieR6] = R6
        pare[iecp6] = cp6
        pare[ieu6] = u6

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

        pare[ieu9] = u9

        pare[ieA6] = A6
        pare[ieA8] = A8
        pare[ieA9] = A9

        pare[ieepf] = epf
        pare[ieeplc] = eplc
        pare[ieephc] = ephc
        pare[ieepht] = epht
        pare[ieeplt] = eplt

        pare[ieetaf] = etaf
        pare[ieetalc] = etalc
        pare[ieetahc] = etahc
        pare[ieetaht] = etaht
        pare[ieetalt] = etalt

        #Fuel mass flow rate
        pare[iemfuel] = ff * mcore * neng

        if (M5 <= 0.999999)
                ichoke5 = 0
        else
                ichoke5 = 1
        end

        if (M7 <= 0.999999)
                ichoke7 = 0
        else
                ichoke7 = 1
        end

        if (Lprint)
                println(" exiting TFCALC")
                println("Tt3 Tt4 u0 u6 u8 fo fc", Tt3, Tt4, u0, u6, u8, fo, pare[iefc])
        end
        
        return ichoke5, ichoke7
end # tfcalc

function check_engine_convergence_failure(pare)
        if sum(pare[ieConvFail, :]) > 0.0 #If any operating point failed to converge
                return true
        else
                return false #All operating points converged
        end
end