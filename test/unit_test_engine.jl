using Zygote
using DelimitedFiles
using ForwardDiff

isGradient = false

@testset "Engine models" begin
    @testset "gasfun.jl" begin

        # =========================
        # gas_N2
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_N2(t1)

        @test s1 == 2358.21669383134
        @test h1 == 2.43487432826e6
        @test cp1 == 1301.485606111
        @test r1 == 296.94

        if isGradient
            ds_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[1], t1)[1]
            dh_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[2], t1)[1]
            dcp_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[3], t1)[1]
            dr_dt = gradient(t1 -> TASOPT.engine.gas_N2(t1)[4], t1)[1]

            @test ds_dt == 0.5578891260784105
            @test dh_dt == 1301.4526599999983
            @test dcp_dt == 0.04276110100000096
            @test dr_dt == nothing
        end

        # =========================
        # gas_N2
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_Ar(t1)

        @test s1 == 1070.0677477147947
        @test h1 == 1.0582e6
        @test cp1 == 520.0
        @test r1 == 208.0

        # =========================
        # gas_CO2
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_CO2(t1)

        @test s1 == 2382.871389635858
        @test h1 == -6.405596226739999e6
        @test cp1 == 1389.688099952
        @test r1 == 188.96

        # =========================
        # gas_H2O
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_H2O(t1)

        @test s1 == 4656.047912649002
        @test h1 == -8.43504044348e6
        @test cp1 == 2943.545113578
        @test r1 == 461.91


        # =========================
        # gas_CH4
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_CH4(t1)

        @test s1 == 8467.631968582504
        @test h1 == 4.574569283279987e6
        @test cp1 == 5288.091482799951
        @test r1 == 519.65


        # =========================
        # gas_C2H6
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C2H6(t1)

        @test s1 == 7408.861495898641
        @test h1 == 5.381058419279992e6
        @test cp1 == 4498.240840800101
        @test r1 == 277.15


        # =========================
        # gas_C3H8
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C3H8(t1)

        @test s1 == 7209.141132729834
        @test h1 == 6.1354732452e6
        @test cp1 == 4411.240840800101
        @test r1 == 188.5


        # =========================
        # gas_C4H10
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C4H10(t1)

        @test s1 == 7230.2763377713845
        @test h1 == 5.908308419279992e6
        @test cp1 == 4448.240840800101
        @test r1 == 143.3

        # =========================
        # gas_C8H18
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C8H18(t1)

        @test s1 == 7262.8799183106685
        @test h1 == 6.234133419279992e6
        @test cp1 == 4443.240840800101
        @test r1 == 72.9

        # =========================
        # gas_C14H30
        # =========================

        t1 = 2333.00e0
        s1, h1, cp1, r1 = TASOPT.engine.gas_C14H30(t1)

        @test s1 == 7253.868041333698
        @test h1 == 6.369750419279992e6
        @test cp1 == 4432.240840800101
        @test r1 == 167.0
        #    


        # =========================
        # gasfun
        # =========================

        igas = 14
        t = 2333.00e0
        s, s_t, h, h_t, cp, r = TASOPT.engine.gasfun(igas, t)

        @test s == 7230.2763377713845
        @test s_t == 1.9066613119588947
        @test h == 5.908308419279992e6
        @test h_t == 4448.240840800101
        @test cp == 4448.240840800101
        @test r == 143.3

        if isGradient
            igas = 14
            t1 = 2333.00e0
            ds_dt = gradient(t1 -> TASOPT.engine.gasfun(igas, t1)[1], t1)[1]
            dh_dt = gradient(t1 -> TASOPT.engine.gasfun(igas, t1)[3], t1)[1]

            epsilon = 1e-6
            s_perturb, s_t_perturb, h_perturb, h_t_perturb, cp_perturb, r_perturb = TASOPT.engine.gasfun(igas, t + epsilon)
            ds_dt_FD = (s_perturb - s) / epsilon
            dh_dt_FD = (h_perturb - h) / epsilon
        end

        # println([s_t, h_t]) # HACK: Not matching
        # println([ds_dt, dh_dt])
        # println([ds_dt_FD, dh_dt_FD])

        # =========================
        # gaschem
        # =========================

        igas = 13
        nchon = TASOPT.engine.gaschem(igas)

        @test nchon == [3, 8, 0, 0]

    end


    @testset "gascalc.jl" begin

        # =========================
        # gassum, gassumd
        # =========================

        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        n = 6
        n_air = n - 1
        t = 215.0
        s, s_t, h, h_t, cp, r = TASOPT.engine.gassum(alpha, n_air, t)


        @test s == -329.0381537005463
        @test s_t == 4.684344004111629
        @test h == -87237.76845599999
        @test h_t == 1007.1339608840001
        @test cp == 1007.1339608840001
        @test r == 288.29530400000004

        s, s_t, h, h_t, cp, cp_t, r = TASOPT.engine.gassumd(alpha, n_air, t)

        @test s == -329.0381537005463
        @test s_t == 4.684344004111629
        @test h == -87237.76845599999
        @test h_t == 1007.1339608840001
        @test cp == 1007.1339608840001
        @test cp_t == 0.009382799699317275
        @test r == 288.29530400000004

        # =========================
        # gas_tset, gas_tsetd
        # =========================

        tguess = 200.0
        t = TASOPT.engine.gas_tset(alpha, n_air, h, tguess)

        @test t == 214.99999914397367

        t, t_hspec, t_al = TASOPT.engine.gas_tsetd(alpha, n_air, h, tguess)
        @test t == 214.99999914397367
        @test t_hspec == 0.000992916572022094
        @test t_al == [85.64080275342843, 75.32364484439914, 8945.435851261927, 13490.80321923344, 42.8542796904542]

        # =========================
        # gas_prat, gas_pratd
        # =========================
        # TODO []: gas_prat, gas_pratd not verified
        # TODO:
        # Fix fortran compiling issue


        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        n = 5
        p = 0.22e5
        t = 215.0

        # =========================
        # gas_delh, gas_delhd
        # =========================
        s, s_t, h, h_t, cp, r = TASOPT.engine.gassum(alpha, n, t)

        delh = 10000.0
        epol = 0.99

        p, t, h, s, cp, r = TASOPT.engine.gas_delh(alpha, n, p, t, h, s, cp, r, delh, epol)

        @test p == 25717.186436176951
        @test t == 224.92715529368544
        @test h == -77237.768461931657
        @test s == -283.57571753601650
        @test cp == 1007.2429752832797
        @test r == 288.29530400000004

        # =========================
        # gas_burn, gas_burnd
        # =========================

        n = 6
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        gamma = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0]
        ifuel = 13
        to = 450.0
        tf = 300.0
        t = 1400.0
        hvap = 0.0

        f, lambda = TASOPT.engine.gas_burn(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

        @test f == -0.45350655959892078
        @test lambda[1] == 1.4291113895654686

        f, lambda, f_to, f_tf, f_t, l_to, l_tf, l_t = TASOPT.engine.gas_burnd(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)
        @test f_to == 4.3531820788987320E-004
        @test f_tf == -3.2351348482849527E-004
        @test f_t == -5.1186081193040599E-004
        @test l_to[1] == -1.1383818413704377E-003
        @test l_tf[1] == 8.4600613962923770E-004
        @test l_t[1] == 1.3385450988489609E-003

        # =========================
        # gas_mach, gas_machd
        # =========================
        # TODO [] gas_machd not verified
        n = 6
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        to = 450.0
        po = 200000.0
        so, s_to, ho, h_to, cpo, ro = TASOPT.engine.gassum(alpha, n - 1, to)
        mo = 1.1
        m = 1.3
        epol = 0.99

        p, t, h, s, cp, r = TASOPT.engine.gas_mach(alpha, n - 1, po, to, ho, so, cpo, ro, mo, m, epol)

        @test p ≈ 154338.88081676030 rtol = 1e-10
        @test t ≈ 417.96446969196779 rtol = 1e-10
        @test h ≈ 117997.67846481415 rtol = 1e-10
        @test s ≈ 342.75772292744381 rtol = 1e-10
        @test cp ≈ 1019.5701545538500 rtol = 1e-10
        @test r ≈ 288.29530400000004 rtol = 1e-10

        # =========================
        # gas_mass
        # =========================

        n = 6
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]
        to = 450.0
        po = 200000.0
        so, s_to, ho, h_to, cpo, ro = TASOPT.engine.gassum(alpha, n - 1, to)

        Mguess = 1.15
        u = sqrt(1.4 * ro * to) * 1.1
        rho = (po / (ro * to))
        mflux = rho * u

        p, t, h, s, cp, r = TASOPT.engine.gas_mass(alpha, n - 1, po, to, ho + rho * u^2 / 2, so, cpo, ro, mflux, Mguess)

        @test p ≈ 118662.19814638462 rtol = 1e-10
        @test t ≈ 388.25628795926565 rtol = 1e-10
        @test h ≈ 87768.12641472857 rtol = 1e-10
        @test s ≈ 267.72827064815743 rtol = 1e-10
        @test cp ≈ 1015.8338429023191 rtol = 1e-10
        @test r ≈ 288.29530400000004 rtol = 1e-10

        n = 6
        alpha = [0.781, 0.209, 0.0004, 0.0, 0.00965, 0.0]

        # =========================
        # gasfuel
        # =========================

        ifuel = 13
        gamma = TASOPT.engine.gasfuel(ifuel, n)
        @test gamma[2] == -3.6283226981894474

    end

    @testset "tfmap.jl" begin

        # TODO []
        # ecmap1, Ncmap1, etmap not verified

        # =========================
        # Ncmap
        # =========================

        pratio = 10.0
        mb = 100.0
        piD = 15.0
        mbD = 120.0
        NbD = 10.0

        Cmap = zeros(9)
        Cmap[1] = 3.50
        Cmap[2] = 0.80
        Cmap[3] = 0.03
        Cmap[4] = 0.75
        Cmap[5] = -0.5
        Cmap[6] = 3.0
        Cmap[7] = 6.0
        Cmap[8] = 2.5
        Cmap[9] = 15.0

        Nb, Nb_pi, Nb_mb = TASOPT.engine.Ncmap(pratio, mb, piD, mbD, NbD, Cmap)

        @test Nb == 8.3597298887214055
        @test Nb_pi == 0.26149475347651785
        @test Nb_mb == 2.4236392203562176E-002

        # =========================
        # ecmap
        # =========================

        pratio = 10.0
        mb = 100.0
        piD = 15.0
        mbD = 120.0
        NbD = 10.0

        Cmap = zeros(9)
        Cmap[1] = 3.50
        Cmap[2] = 0.80
        Cmap[3] = 0.03
        Cmap[4] = 0.75
        Cmap[5] = -0.5
        Cmap[6] = 3.0
        Cmap[7] = 6.0
        Cmap[8] = 2.5
        Cmap[9] = 15.0

        effo = 0.95
        piK = 0.9 * pratio
        effK = 0.93

        eff, eff_pi, eff_mb = TASOPT.engine.ecmap(pratio, mb, piD, mbD, Cmap, effo, piK, effK)

        @test eff == 1.8781007331361292
        @test eff_pi == 0.92374562099125368
        @test eff_mb == 1.6164204096982837E-003

        # =========================
        # Pimap
        # =========================

        pratio = 10.0
        mb = 100.0
        Nb = 10.0
        piD = 15.0
        mbD = 120.0
        NbD = 10.0

        Cmap = zeros(9)
        Cmap[1] = 3.50
        Cmap[2] = 0.80
        Cmap[3] = 0.03
        Cmap[4] = 0.75
        Cmap[5] = -0.5
        Cmap[6] = 3.0
        Cmap[7] = 6.0
        Cmap[8] = 2.5
        Cmap[9] = 15.0

        eff, eff_pi, eff_mb = TASOPT.engine.Pimap(mb, Nb, piD, mbD, NbD, Cmap)

        @test eff == 16.579462807918382
        @test eff_pi == -3.5593220338983059E-002
        @test eff_mb == 4.419641196046076

    end

    @testset "tfcool.jl" begin

        ncrowx = 2
        ncrow = 2
        Tt3 = 300.0
        Tt4 = 1300.0
        dTstreak = 300.0
        Trrat = 0.9
        efilm = 0.5
        tfilm = 0.3
        StA = 0.2

        epsrow = zeros(ncrowx)
        epsrow[1] = 0.1
        epsrow[2] = 0.1

        Tmrow = zeros(ncrowx)
        Tmrow[1] = 1100.0
        Tmrow[2] = 1100.0

        ncrow, epsrow, epsrow_Tt3, epsrow_Tt4, epsrow_Trr = TASOPT.engine.mcool(ncrowx,
            Tmrow, Tt3, Tt4, dTstreak, Trrat,
            efilm, tfilm, StA)

        ncrow == 1
        epsrow[1] == 0.10313901345291480
        epsrow_Tt3[1] == 1.7595366888535866E-004
        epsrow_Tt4[1] == 2.8152587021657390E-004
        epsrow_Trr[1] == 0.0000000000000000


        ncrowx = 2
        ncrow = 2
        Tt3 = 300.0
        Tt4 = 1300.0
        dTstreak = 300.0
        Trrat = 0.9
        efilm = 0.5
        tfilm = 0.3
        StA = 0.2

        epsrow = zeros(ncrowx)
        epsrow[1] = 0.1
        epsrow[2] = 0.1

        Tmrow = TASOPT.engine.Tmcalc(ncrowx, ncrow,
            Tt3, Tt4, dTstreak, Trrat,
            efilm, tfilm, StA, epsrow)

        @test Tmrow[1] == 1106.8965517241379
        @test Tmrow[2] == 840.00000000000000

        if isGradient
            # AD
            d_Tmrow_d_Tt3 = gradient(Tt3 -> TASOPT.engine.Tmcalc(ncrowx, ncrow,
                    Tt3, Tt4, dTstreak, Trrat,
                    efilm, tfilm, StA, epsrow)[1], Tt3)[1]

            # FD
            epsilon = 1e-6
            Tmrow_d = TASOPT.engine.Tmcalc(ncrowx, ncrow,
                Tt3 + epsilon, Tt4, dTstreak, Trrat,
                efilm, tfilm, StA, epsrow)

            d_Tmrow_d_Tt3_FD = (Tmrow_d[1] - Tmrow[1]) / epsilon

            # FD AD consistency test
            @test d_Tmrow_d_Tt3_FD ≈ d_Tmrow_d_Tt3 rtol = 1e-5
        end

        #Test cooled HPT efficiency
        epht0 = 0.9
        epht_fc = -0.36
        fc0 = 0.16
        fc = 0.2
        epht = TASOPT.engine.find_cooled_hpt_efficiency(epht0, epht_fc, fc0, fc)
        @test epht ≈ 0.8856
    end

    @testset "tfsize.jl" begin

        # TODO []: icool change the function signature. Only works for icool = 1.
        # note: icool::Int -> opt_cooling::String
        # TODO []: the output shall cover most of the functions. Ideally, test all other outputs.

        gee = 9.8100000000000005
        M0 = 0.80000000000000004
        T0 = 219.43067572699252
        p0 = 23922.608843328788
        a0 = 296.85578884697560
        M2 = 0.59999999999999998
        M25 = 0.59999999999999998
        Feng = 22182.101361240744
        Phiinl = 0.0000000000000000
        Kinl = 0.0000000000000000
        eng_has_BLI_cores = false
        BPR = 5.0999999999999996
        pif = 1.6850000000000001
        pilc = 8.0000000000000000
        pihc = 3.7500000000000000
        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 89407.867373646965
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = "fixed_coolingflowratio"
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow = zeros(ncrowx)
        epsrow[1] = 0.12811308404512714
        epsrow[2] = 5.4411331284501797E-002
        epsrow[3] = 1.5791188045605239E-002
        epsrow[4] = 0.0000000000000000
        Tt4 = 1587.0000000000000
        Tmrow = zeros(ncrowx)
        fc0 = 0.0
        epht_fc = 0.0
        hvap = 0.0
        Δh_PreC = 0.0
        Δh_InterC = 0.0
        Δh_Regen = 0.0
        Δh_TurbC = 0.0
        Δp_PreC = 0.0
        Δp_InterC = 0.0
        Δp_Regen = 0.0

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
        T25c, u25c, p25c, cp25c, R25c, A25,
        T5, u5, p5, cp5, R5, A5,
        T6, u6, p6, cp6, R6, A6,
        T7, u7, p7, cp7, R7, A7,
        T8, u8, p8, cp8, R8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfsize!(gee, M0, T0, p0, a0, M2, M25,
            Feng, Phiinl, Kinl, eng_has_BLI_cores,
            BPR, pif, pilc, pihc,
            pid, pib, pifn, pitn,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
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



        @test etaf == 0.8539899545024271
        @test etalc ==  0.8298895182207285
        @test etahc == 0.8384795893726779
        @test etaht ==  0.8978789812518558
        @test etalt == 0.9191845671925591
        @test Tmrow[1] ≈ 1119.1584455133655  rtol = 1e-10

        if isGradient
            # AD
            d_etaf_dM0 = gradient(M0 -> TASOPT.engine.tfsize(gee, M0, T0, p0, a0, M2, M25,
                    Feng, Phiinl, Kinl, eng_has_BLI_cores,
                    BPR, pif, pilc, pihc,
                    pid, pib, pifn, pitn,
                    Ttf, ifuel, hvap, etab,
                    epf0, eplc0, ephc0, epht0, eplt0,
                    mofft, Pofft,
                    Tt9, pt9, Tt4,
                    epsl, epsh,
                    opt_cooling,
                    Mtexit, dTstrk, StA, efilm, tfilm,
                    M4a, ruc,
                    ncrowx, ncrow,
                    epsrow)[115], M0, Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                    Δp_PreC, Δp_InterC, Δp_Regen)[1]

            # FD
            epsilon = 1e-6
            etaf_d = TASOPT.engine.tfsize(gee, M0 + epsilon, T0, p0, a0, M2, M25,
                Feng, Phiinl, Kinl, eng_has_BLI_cores,
                BPR, pif, pilc, pihc,
                pid, pib, pifn, pitn,
                Ttf, ifuel, hvap, etab,
                epf0, eplc0, ephc0, epht0, eplt0,
                mofft, Pofft,
                Tt9, pt9, Tt4,
                epsl, epsh,
                opt_cooling,
                Mtexit, dTstrk, StA, efilm, tfilm,
                fc0, epht_fc,
                M4a, ruc,
                ncrowx, ncrow,
                epsrow, Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
                Δp_PreC, Δp_InterC, Δp_Regen)[115]

            d_etaf_dM0_FD = (etaf_d - etaf) / epsilon
            # FD AD consistency test

            @test d_etaf_dM0_FD ≈ d_etaf_dM0 rtol = 1e-4
        end

        gee = 9.8100000000000005
        M0 = 0.80000000000000004
        T0 = 219.43067572699252
        p0 = 23922.608843328788
        a0 = 296.85578884697560
        M2 = 0.59999999999999998
        M25 = 0.59999999999999998
        Feng = 22182.101361240744
        Phiinl = 0.0000000000000000
        Kinl = 0.0000000000000000
        eng_has_BLI_cores = false
        BPR = 5.0999999999999996
        pif = 1.6850000000000001
        pilc = 8.0000000000000000
        pihc = 3.7500000000000000
        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 89407.867373646965
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = "fixed_Tmetal"
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow = zeros(ncrowx)
        Tt4 = 1587.0000000000000
        Tmrow = [1114.7762513333382, 1102.6358896742904, 1098.9326716275868, 1021.0788483883676]


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
        T25c, u25c, p25c, cp25c, R25c, A25,
        T5, u5, p5, cp5, R5, A5,
        T6, u6, p6, cp6, R6, A6,
        T7, u7, p7, cp7, R7, A7,
        T8, u8, p8, cp8, R8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfsize!(gee, M0, T0, p0, a0, M2, M25,
            Feng, Phiinl, Kinl, eng_has_BLI_cores,
            BPR, pif, pilc, pihc,
            pid, pib, pifn, pitn,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
            mofft, Pofft,
            Tt9, pt9, Tt4,
            epsl, epsh,
            opt_cooling,
            Mtexit, dTstrk, StA, efilm, tfilm,
            fc0, epht_fc,
            M4a, ruc,
            ncrowx, ncrow,
            epsrow, Tmrow, Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
            Δp_PreC, Δp_InterC, Δp_Regen)

        @test epsrow[1] ≈ 0.1303162243242347 rtol = 1e-10

    end

    @testset "gaussn.jl" begin
        nsiz = 3
        nn = 3
        nrhs = 1

        z = zeros((nsiz, nsiz))
        z[1, 1] = 0.52762351
        z[1, 2] = 0.85729295
        z[1, 3] = 0.32073675
        z[2, 1] = 0.70046643
        z[2, 2] = 0.43823310
        z[2, 3] = 0.60231747
        z[3, 1] = 0.86206149
        z[3, 2] = 0.35227634
        z[3, 3] = 0.61320909

        r = zeros((nsiz, nrhs))
        r[1, 1] = 0.52873024
        r[2, 1] = 0.78246347
        r[3, 1] = 0.10212863

        r = TASOPT.engine.gaussn(nsiz, nn, z, r, nrhs)

        @test r[1, 1] == -3.9352200465675033
        @test r[2, 1] == 1.1548320796061180
        @test r[3, 1] == 5.0353302305150418

    end

    @testset "tfoper.jl mode 1" begin

        gee = 9.8100000000000005
        M0 = 0.26302467815397762
        T0 = 288.00000000000000
        p0 = 101320.00000000000
        a0 = 340.08940001123227
        Tref = 288.19999999999999
        pref = 101320.00000000000

        Phiinl = 0.0
        Kinl = 0.0
        eng_has_BLI_cores = false

        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Gearf = 1.0000000000000000
        pifD = 1.6850000000000001
        pilcD = 8.0000000000000000
        pihcD = 3.7500000000000000
        pihtD = 2.1601257635200488
        piltD = 6.2886975330083716
        mbfD = 235.16225770724063
        mblcD = 46.110246609262873
        mbhcD = 7.8056539219349039
        mbhtD = 4.3594697284253883
        mbltD = 8.7016090343744406
        NbfD = 1.0790738309310697
        NblcD = 1.0790738309310697
        NbhcD = 0.77137973563891493
        NbhtD = 0.44698693289691338
        NbltD = 0.48396724306758404
        A2 = 1.3863121762890294
        A25 = 3.8585338087708761E-002
        A5 = 0.19210855588408102
        A7 = 0.64211443204484309
        opt_calc_call = "oper_fixedTt4"
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 77800.595231538944
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = "fixed_coolingflowratio"
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow3 = [0.12061791584226822, 5.1292591721870069E-002, 1.5478853228971187E-002, 0.0000000000000000]
        Tmrow3 = [1000.0, 1000.0, 1000.0, 1000.0]
        fc0 = 0.0
        epht_fc = 0.0

        M2 = 1.0
        pif = 0.0000000000000000
        pilc = 0.0000000000000000
        pihc = 0.0000000000000000
        mbf = 0.0000000000000000
        mblc = 0.0000000000000000
        mbhc = 0.0000000000000000
        Tt4 = 1783.8000000000002
        pt5 = 0.0000000000000000
        mcore = 0.0
        Feng = 0.0
        M25 = 0.0
        hvap = 0.0
        Δh_PreC = 0.0
        Δh_InterC = 0.0
        Δh_Regen = 0.0
        Δh_TurbC = 0.0
        Δp_PreC = 0.0
        Δp_InterC = 0.0
        Δp_Regen = 0.0

        # Fix mass flow find temperature
        TSFC, Fsp, hfuel, ff,
        Feng, mcore,
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
        T25c, u25c, p25c, cp25c, R25c, M25c,
        T5, u5, p5, cp5, R5, M5,
        T6, u6, p6, cp6, R6, M6, A6,
        T7, u7, p7, cp7, R7, M7,
        T8, u8, p8, cp8, R8, M8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfoper!(gee, M0, T0, p0, a0, Tref, pref,
            Phiinl, Kinl, eng_has_BLI_cores,
            pid, pib, pifn, pitn,
            Gearf,
            pifD, pilcD, pihcD, pihtD, piltD,
            mbfD, mblcD, mbhcD, mbhtD, mbltD,
            NbfD, NblcD, NbhcD, NbhtD, NbltD,
            A2, A25, A5, A7,
            opt_calc_call,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
            mofft, Pofft,
            Tt9, pt9,
            epsl, epsh,
            opt_cooling,
            Mtexit, dTstrk, StA, efilm, tfilm,
            fc0, epht_fc,
            M4a, ruc,
            ncrowx, ncrow,
            epsrow3, Tmrow3,
            Feng,
            M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
            Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
            Δp_PreC, Δp_InterC, Δp_Regen)


        @test etaf ≈ 0.8741690445868673 atol = 1e-8
        @test etalc ≈ 0.8367067552587594 atol = 1e-8
        @test etahc ≈ 0.8393327258267146 atol = 1e-8
        @test etaht ≈ 0.8974130931530157  atol = 1e-8
        @test etalt ≈ 0.9176389341085217 atol = 1e-8

        Tt4_ref = 1783.8000000000002
        ht4_ref = 470172.6605599733
        pt4_ref = 2.7406297732608365e6
        cpt4_ref = 1312.0669933280028
        Rt4_ref =  288.12027917359467
        @test Tt4 ≈ Tt4_ref atol = 1e-8 * Tt4_ref
        @test ht4 ≈ ht4_ref atol = 1e-8 * ht4_ref
        @test pt4 ≈ pt4_ref atol = 1e-8 * pt4_ref
        @test cpt4 ≈ cpt4_ref atol = 1e-8 * cpt4_ref
        @test Rt4 ≈ Rt4_ref atol = 1e-8 * Rt4_ref

    end

    @testset "tfoper.jl mode 2" begin

        gee = 9.8100000000000005
        M0 = 0.26302467815397762
        T0 = 288.00000000000000
        p0 = 101320.00000000000
        a0 = 340.08940001123227
        Tref = 288.19999999999999
        pref = 101320.00000000000

        Phiinl = 0.0
        Kinl = 0.0
        eng_has_BLI_cores = false

        pid = 0.99800000000000000
        pib = 0.93999999999999995
        pifn = 0.97999999999999998
        pitn = 0.98899999999999999
        Gearf = 1.0000000000000000
        pifD = 1.6850000000000001
        pilcD = 8.0000000000000000
        pihcD = 3.7500000000000000
        pihtD = 2.1601257635200488
        piltD = 6.2886975330083716
        mbfD = 235.16225770724063
        mblcD = 46.110246609262873
        mbhcD = 7.8056539219349039
        mbhtD = 4.3594697284253883
        mbltD = 8.7016090343744406
        NbfD = 1.0790738309310697
        NblcD = 1.0790738309310697
        NbhcD = 0.77137973563891493
        NbhtD = 0.44698693289691338
        NbltD = 0.48396724306758404
        A2 = 1.3863121762890294
        A25 = 3.8585338087708761E-002
        A5 = 0.19210855588408102
        A7 = 0.64211443204484309
        opt_calc_call = "oper_fixedTt4"
        Ttf = 280.00000000000000
        ifuel = 24
        etab = 0.98499999999999999
        epf0 = 0.89480000000000004
        eplc0 = 0.88000000000000000
        ephc0 = 0.87000000000000000
        epht0 = 0.88900000000000001
        eplt0 = 0.89900000000000002
        mofft = 0.56969999999999998
        Pofft = 77800.595231538944
        Tt9 = 300.00000000000000
        pt9 = 30000.000000000000
        epsl = 1.0000000000000000E-002
        epsh = 2.1999999999999999E-002
        opt_cooling = "fixed_coolingflowratio"
        Mtexit = 1.0000000000000000
        dTstrk = 200.00000000000000
        StA = 8.9999999999999997E-002
        efilm = 0.69999999999999996
        tfilm = 0.29999999999999999
        M4a = 0.90000000000000002
        ruc = 0.14999999999999999
        ncrowx = 4
        ncrow = 4
        epsrow3 = [0.12061791584226822, 5.1292591721870069E-002, 1.5478853228971187E-002, 0.0000000000000000]
        Tmrow3 = [1000.0, 1000.0, 1000.0, 1000.0]
        fc0 = 0.0
        epht_fc = 0.0

        M2 = 1.0
        pif = 0.0000000000000000
        pilc = 0.0000000000000000
        pihc = 0.0000000000000000
        mbf = 0.0000000000000000
        mblc = 0.0000000000000000
        mbhc = 0.0000000000000000
        Tt4 = 1783.8000000000002
        pt5 = 0.0000000000000000
        mcore = 0.0
        M25 = 0.0
        hvap = 0.0
        Feng = 0.0
        Δh_PreC = 0.0
        Δh_InterC = 0.0
        Δh_Regen = 0.0
        Δh_TurbC = 0.0
        Δp_PreC = 0.0
        Δp_InterC = 0.0
        Δp_Regen = 0.0

        # Fix temperature find mass flow
        opt_cooling = "fixed_Tmetal"
        TSFC, Fsp, hfuel, ff,
        Feng, mcore,
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
        T25c, u25c, p25c, cp25c, R25c, M25c,
        T5, u5, p5, cp5, R5, M5,
        T6, u6, p6, cp6, R6, M6, A6,
        T7, u7, p7, cp7, R7, M7,
        T8, u8, p8, cp8, R8, M8, A8,
        u9, A9,
        epf, eplc, ephc, epht, eplt,
        etaf, etalc, etahc, etaht, etalt,
        Lconv = TASOPT.engine.tfoper!(gee, M0, T0, p0, a0, Tref, pref,
            Phiinl, Kinl, eng_has_BLI_cores,
            pid, pib, pifn, pitn,
            Gearf,
            pifD, pilcD, pihcD, pihtD, piltD,
            mbfD, mblcD, mbhcD, mbhtD, mbltD,
            NbfD, NblcD, NbhcD, NbhtD, NbltD,
            A2, A25, A5, A7,
            opt_calc_call,
            Ttf, ifuel, hvap, etab,
            epf0, eplc0, ephc0, epht0, eplt0,
            mofft, Pofft,
            Tt9, pt9,
            epsl, epsh,
            opt_cooling,
            Mtexit, dTstrk, StA, efilm, tfilm,
            fc0, epht_fc,
            M4a, ruc,
            ncrowx, ncrow,
            epsrow3, Tmrow3,
            Feng,
            M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
            Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
            Δp_PreC, Δp_InterC, Δp_Regen)

        @test etaf ≈ 0.8862602081467577 rtol = 1e-6
        @test etalc ≈ 0.6539900828671803 rtol = 1e-6
        @test etahc ≈ 0.8390003200177845 rtol = 1e-6
        @test etaht ≈ 0.8977614695292404 rtol = 1e-6
        @test etalt ≈ 0.8649073213143893 rtol = 1e-6

        @test Tt4 ≈ 1783.8000000000002 rtol = 1e-6
        @test ht4 ≈ 292168.48621470056  rtol = 1e-6
        @test pt4 ≈ 1.0428916387810945e6  rtol = 1e-6
        @test cpt4 ≈ 1323.0363088269542 rtol = 1e-6
        @test Rt4 ≈ 288.2159811164396  rtol = 1e-6
    end

    @testset "tfweight.jl" begin

        ac = TASOPT.load_default_model()
        ip = ipcruise1

        ac.parg[igGearf] = 1.0
        ac.pared[iemblcD, ip] = 46.110246609262873
        ac.pared[ieBPR, ip] = 5.0999999999999996
        ac.pared[iepilc, ip] = 3.0
        ac.pared[iepihc, ip] = 10.000000000000000 
        ac.parg[igdfan] = 1.3927234305722356
        ac.parg[igdlcomp] = 0.67240459668963337
        ac.parg[igrSnace] = 16.000000000000000
        ac.parg[igneng] = 2.0
        ac.parg[igfeadd] = 0.10000000000000001
        ac.parg[igfpylon] = 0.10000000000000001

        Weng, Wnac, Webare, W_HXs, Snace1 = TASOPT.engine.tfweight(ac)

        @test Weng ≈ 46847.51154286845 rtol = 1e-10
        @test Wnac ≈ 9411.055803345604 rtol = 1e-10
        @test Webare ≈ 30161.446412552305 rtol = 1e-10
        @test Snace1 ≈ 24.374719583103083 rtol = 1e-10

    end

    @testset "map_functions.jl" begin
        #Check reverse interpolation function
        map = TASOPT.engine.FanMap
        Wc_target = 559.314
        PR_target = 1.31
        outps_NR = TASOPT.engine.find_NR_inverse_with_derivatives(map.itp_Wc, map.itp_PR, 
            Wc_target, PR_target)
        outp_NR_check = (0.7000000000000002, 1.799999999999995, 0.0002760733345430258, 0.6100706271038377, 0.00430742991535455, -2.9037427622059644)
        
        for (i,outp_NR) in enumerate(outps_NR)
            @test outp_NR ≈ outps_NR[i]
        end

        #Check function to compute speed and efficiency
        pratio = 1.6
        mb = 0.8
        piD = 1.7
        mbD = 1.0
        NbD = 1.0
        ep0 = 0.9
        outps_Ne = TASOPT.engine.calculate_compressor_speed_and_efficiency(map, pratio, mb, piD, mbD, NbD, ep0)
        outps_Ne_check = (0.8826430459793184, 0.8251883153075635, 0.49351112427246263, -0.003692654225532323, -0.39245022849891215, 0.929263500483822, 0.8738166155195252, 1.4836439672220545)
        
        for (i,outp_Ne) in enumerate(outps_Ne)
            @test outp_Ne ≈ outps_Ne_check[i]
        end

        #Use finite difference to validate derivatives of the maps
        eps = 1e-6
        outps_pi = TASOPT.engine.calculate_compressor_speed_and_efficiency(map, pratio + eps, mb, piD, mbD, NbD, ep0)
        outps_mb = TASOPT.engine.calculate_compressor_speed_and_efficiency(map, pratio, mb + eps, piD, mbD, NbD, ep0)
        
        dN_dpi = (outps_pi[1] - outps_Ne[1])/eps
        depol_dpi = (outps_pi[2] - outps_Ne[2])/eps

        dN_dmb = (outps_mb[1] - outps_Ne[1])/eps
        depol_dmb = (outps_mb[2] - outps_Ne[2])/eps

        @test dN_dpi ≈ outps_Ne[3] rtol = 1e-4
        @test dN_dmb ≈ outps_Ne[4] rtol = 1e-4
        @test depol_dpi ≈ outps_Ne[5] rtol = 1e-4
        @test depol_dmb ≈ outps_Ne[6] rtol = 1e-4
    end 
end