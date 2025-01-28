"""
function check_AD(gee, M0, T0, p0, a0, Tref, pref, Phiinl, Kinl, pid, pib, pifn, pitn, Gearf,
    pifD, pilcD, pihcD, pihtD, piltD, mbfD, mblcD, mbhcD, mbhtD, mbltD, NbfD, NblcD, NbhcD,
    NbhtD, NbltD, A2, A25, A5, A7, Ttf, ifuel, etab, epf0, eplc0, ephc0, epht0, eplt0, pifK,
    epfK, mofft, Pofft, Tt9, pt9, epsl, epsh, Mtexit, dTstrk, StA, efilm, tfilm, M4a, ruc, M25,
    mcore, pt5, Tt4, mbhc, mblc, mbf, pihc, pilc, pif, M2, epsrow)

Newton matrix type checker

    Calculates the most demanding data type required to capture turbofan operations. See tfoper
    inputs/outputs for more detail on inputs and outputs

"""
function check_AD(gee, M0, T0, p0, a0, Tref, pref, Phiinl, Kinl, pid, pib, pifn, pitn, Gearf,
    pifD, pilcD, pihcD, pihtD, piltD, mbfD, mblcD, mbhcD, mbhtD, mbltD, NbfD, NblcD, NbhcD,
    NbhtD, NbltD, A2, A25, A5, A7, Ttf, ifuel, etab, epf0, eplc0, ephc0, epht0, eplt0, pifK,
    epfK, mofft, Pofft, Tt9, pt9, epsl, epsh, Mtexit, dTstrk, StA, efilm, tfilm, M4a, ruc, M25,
    mcore, pt5, Tt4, mbhc, mblc, mbf, pihc, pilc, pif, M2, epsrow)

    prod = gee * M0 * T0 * p0 * a0 * Tref * pref * Phiinl * Kinl * pid * pib * pifn * pitn * Gearf
    prod *= pifD * pilcD * pihcD * pihtD * piltD
    prod *= mbfD * mblcD * mbhcD * mbhtD * mbltD
    prod *= NbfD * NblcD * NbhcD * NbhtD * NbltD
    prod *= A2 * A25 * A5 * A7
    prod *= Ttf * ifuel * etab
    prod *= epf0 * eplc0 * ephc0 * epht0 * eplt0
    prod *= pifK * epfK * mofft * Pofft * Tt9 * pt9 * epsl * epsh
    prod *= Mtexit * dTstrk * StA * efilm * tfilm
    prod *= M4a * ruc
    prod *= epsrow[1] * M2 * pif * pilc * pihc * mbf * mblc * mbhc * Tt4 * pt5 * mcore * M25

    T = typeof(prod)
    return T

end # check_AD


"""
function create_newton_arrays(T)

Create Newton solver arrays for the turbofan sizing
"""
function create_newton_arrays(T)
    res = zeros(T, 9, 1)
    resPlusOne = zeros(T, 9, 1)
    a = zeros(T, 9, 9)
    rrel = zeros(T, 9)
    rsav = zeros(T, 9)
    asav = zeros(T, 9, 10)

    res_dlls = zeros(T, 9)
    a_dlls = zeros(T, 9, 9)

    return res, resPlusOne, a, rrel, rsav, asav, res_dlls, a_dlls
end # create_newton_arrays'


function create_massfrac_arrays(T, n)
    alpha = zeros(T, n)    # air
    beta = zeros(T, n)     # fuel
    gamma = zeros(T, n)    # combustion-caused change in air
    lambda = zeros(T, n)   # combustion product gas
    lambdap = zeros(T, n)  # combustion product gas with cooling flow added

    lam_Tt3 = zeros(T, n)
    lam_Ttf = zeros(T, n)
    lam_pl = zeros(T, n)
    lam_ph = zeros(T, n)
    lam_ml = zeros(T, n)
    lam_mh = zeros(T, n)
    lam_Tb = zeros(T, n)
    lamp_pl = zeros(T, n)
    lamp_ph = zeros(T, n)
    lamp_mf = zeros(T, n)
    lamp_ml = zeros(T, n)
    lamp_mh = zeros(T, n)
    lamp_Tb = zeros(T, n)
    lamp_Mi = zeros(T, n)

    T_al = zeros(T, n)
    p_al = zeros(T, n)
    h_al = zeros(T, n)
    s_al = zeros(T, n)
    R_al = zeros(T, n)
    cp_al = zeros(T, n)

    return alpha, beta, gamma, lambda, lambdap, lam_Tt3, lam_Ttf, 
        lam_pl, lam_ph, lam_ml, lam_mh, lam_Tb, lamp_pl, lamp_ph,
        lamp_mf, lamp_ml, lamp_mh, lamp_Tb, lamp_Mi, T_al, p_al, 
        h_al, s_al, R_al, cp_al

end # create_massfrac_arrays


"""
function freestream_gasprop(alpha, nair, T0, p0, M0, a0)

Calculate the static/total properties of freestream air
"""
function freestream_gasprop(alpha, nair, t, M0, a0)
    #---- freestream static quantities
    s0, dsdt0, h0, dhdt0, cp0, R0 = gassum(alpha, nair, t)
    gam0 = cp0 / (cp0 - R0)
    u0 = M0 * a0
    #
    # ===============================================================
    #---- freestream total quantities
    hspec = h0 + 0.5 * u0^2
    Tguess = T0 * (1.0 + 0.5 * (gam0 - 1.0) * M0^2)
    Tt0 = gas_tset(alpha, nair, hspec, Tguess)
    st0, dsdt, ht0, dhdt, cpt0, Rt0 = gassum(alpha, nair, Tt0)
    pt0 = p0 * exp((st0 - s0) / Rt0)
    gamt0 = cpt0 / (cpt0 - Rt0)
    at0 = sqrt(Tt0 * Rt0 * gamt)

    fs = flowStation(
        Mi = M0,

        ps = p0,
        Ts = T0,
        ss = s0,
        ss_Ts = dsdt0,
        hs = h0,
        hs_Ts = dhdt0,
        cps = cp0,
        gams = gam0,
        Rs = R0,
        as = a0,

        pt = pt0,
        Tt = Tt0,
        st = st0,
        st_Tt = dsdt,
        ht = ht0,
        ht_Tt = dhdt,
        cpt = cpt0,
        gamt = gamt0,
        Rt = Rt0,
        at = at0
    )

    return fs

end # freestream_gasprop


"""
function offtake_plume(fs0, fs9)

    Calculates key offtake plume flowstation properties

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fs0::flowStation` : freestream flowStation
    - `fs9::flowStation` : offtake plume flowStation

    **Output**
    ------
    - `fs9` : offtake plume flowStation (updated) 
"""
function offtake_plume(fs0::flowStation, pt9, Tt9)
    Trat = (fs0.ps / pt9)^(fs0.Rt / fs0.cpt)
    if (Trat < 1.0)
          u9 = sqrt(2.0 * fs0.cpt * Tt9 * (1.0 - Trat))
          rhos9 = fs0.ps / (fs0.Rt * fs0.Tt * Trat)
    else
          u9 = 0.0
          rhos9 = fs0.ps / (fs0.Rt * fs0.Tt)
    end

    fs9 = flowStation(
        u = u9,
        rhos = rhos9
    )

    return fs9
end # offtake_plume


"""
function cooling_mass_flow_ratio(icool, mcore, mofft, icore, ncrow, epsrow)

Calculates cooling mass flow ratios
"""
function cooling_mass_flow_ratio(icool, mcore, mofft, icore, ncrow, epsrow)
    if (icool == 1)
        if (mcore == 0.0)
              fo = 0.0
        else
              fo = mofft / mcore
        end
        for icrow = 1:ncrow
              fc = fc + (1.0 - fo) * epsrow[icrow]
        end
    end

    return fo, fc
end # cooling_mass_flow_ratio


"""
function boundary_layer_calculations(at0, gam0, Mi, iBLIc, mf, ml, Tref, pref, Tt0, pt0, Kinl)

Boundary layer calcualtor
    
    Calculates pressure corrections due to boundary layers in the fan and core streams

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `at0`:   freestream total speed of sound [m/s]
    - `gam0`:  freestream specific heat ratio
    - `Mi`:    inlet mach number
    - `iBLIc`: 0=core in clear flow, 1=core sees Phiinl
    - `mf`:    fan corrected mass flow [kg/s]
    - `ml`:    lpc corrected mass flow [kg/s]
    - `Tref`:  reference temperature [K]
    - `pref`:  reference pressure [Pa]
    - `Tt0`:   freestream total temperature [K]
    - `pt0`:   freestream total pressure [Pa]

    **Output**
    ------
    - `sbfan`:     fan stream pressure correction
    - `sbfan_?`:   fan correction derivatives
    - `sbcore`:    core stream pressure correction
    - `sbcore_?`:  core correction derivatives
"""
function boundary_layer_calculations(fs2, fs19, fs18, fs0, Mi, iBLIc, mf, ml, Tref, pref, Kinl)
    if (fs0.u == 0.0)
        #----- static case... no BL defect
        sbfan = 0.0
        sbfan_mf = 0.0
        sbfan_ml = 0.0
        sbfan_Mi = 0.0

        sbcore = 0.0
        sbcore_mf = 0.0
        sbcore_ml = 0.0
        sbcore_Mi = 0.0

    else
            #----- account for inlet BLI defect via mass-averaged entropy
            a2sq = fs0.at^2 / (1.0 + 0.5 * (fs0.gams - 1.0) * Mi^2)
            a2sq_Mi = -a2sq / (1.0 + 0.5 * (fs0.gams - 1.0) * Mi^2) *
                    (fs0.gams - 1.0) * Mi

            if (iBLIc == 0)
                #------ BL mixes with fan flow only
                #c      mmix    = mf*sqrt(Tref/Tt2) * pt2   /pref
                #c      mmix_mf =    sqrt(Tref/Tt2) * pt2   /pref
                #c      mmix_Mi = mf*sqrt(Tref/Tt2) * pt2_Mi/pref

                mmix = mf * sqrt(Tref / fs0.Tt) * fs0.pt / pref
                mmix_mf = sqrt(Tref / fs0.Tt) * fs0.pt / pref
                mmix_Mi = 0.0

                sbfan = Kinl * fs0.gams / (mmix * a2sq)
                sbfan_mf = (-sbfan / mmix) * mmix_mf
                sbfan_ml = 0.0
                sbfan_Mi = (-sbfan / mmix) * mmix_Mi +
                            (-sbfan / a2sq) * a2sq_Mi

                sbcore = 0.0
                sbcore_mf = 0.0
                sbcore_ml = 0.0
                sbcore_Mi = 0.0

            else
                #------ BL mixes with fan + core flow
                #c      mmix    = mf*sqrt(Tref/Tt2 ) * pt2    /pref
                #c             + ml*sqrt(Tref/Tt19) * pt19   /pref
                #c      mmix_mf =    sqrt(Tref/Tt2 ) * pt2    /pref
                #c      mmix_ml =    sqrt(Tref/Tt19) * pt19   /pref
                #c      mmix_Mi = mf*sqrt(Tref/Tt2 ) * pt2_Mi /pref
                #c             + ml*sqrt(Tref/Tt19) * pt19_Mi/pref

                mmix = mf * sqrt(Tref / fs0.Tt) * fs0.pt / pref +
                        ml * sqrt(Tref / fs0.Tt) * fs0.pt / pref
                mmix_mf = sqrt(Tref / fs0.Tt) * fs0.pt / pref
                mmix_ml = sqrt(Tref / fs0.Tt) * fs0.pt / pref
                mmix_Mi = 0.0

                sbfan = Kinl * fs0.gams / (mmix * a2sq)
                sbfan_mf = (-sbfan / mmix) * mmix_mf
                sbfan_ml = (-sbfan / mmix) * mmix_ml
                sbfan_Mi = (-sbfan / mmix) * mmix_Mi +
                            (-sbfan / a2sq) * a2sq_Mi

                sbcore = sbfan
                sbcore_mf = sbfan_mf
                sbcore_ml = sbfan_ml
                sbcore_Mi = sbfan_Mi
            end

    end

    #---- note: BL is assumed adiabatic, 
    #-     so Tt2,ht2,st2,cpt2,Rt2  will not change due to BL ingestion

    fs2.Tt = fs18.Tt
    fs2.ht = fs18.ht
    fs2.st = fs18.st
    fs2.cpt = fs18.cpt
    fs2.Rt = fs18.Rt
    fs2.pt = fs18.pt * exp(-sbfan)
    fs2.pt_mf = fs2.pt * (-sbfan_mf)
    fs2.pt_ml = fs2.pt * (-sbfan_ml)
    fs2.pt_M  = fs2.pt * (-sbfan_Mi)

    fs19.Tt = fs18.Tt
    fs19.ht = fs18.ht
    fs19.st = fs18.st
    fs19.cpt = fs18.cpt
    fs19.Tt_ht = 1 / fs19.cpt
    fs19.Rt = fs18.Rt
    fs19.pt = fs18.pt * exp(-sbcore)
    fs19.pt_mf = fs19.pt * (-sbcore_mf)
    fs19.pt_ml = fs19.pt * (-sbcore_ml)
    fs19.pt_M  = fs19.pt * (-sbcore_Mi)

    return fs2, fs19
end # boundary_layer_calculations


"""
function comp_eff(comp::compressor, piK, epfK; useTbl=true)

    Wrapper for ecmap and ecTblMap for use with the compressor struct
"""
function comp_eff(pr, mb, comp::compressor, pifK, epfK; useTbl=true)
    ep = 0.0
    ep_mb = 0.0
    ep_pr = 0.0
    
    if useTbl
        ep, ep_pr, ep_mb = ecTblMap(pr, mb, comp.prD, comp.mbD, comp.map, comp.oob_map, comp.epD; pifK=pifK, epfK=epfK)
    else
        ep, ep_pr, ep_mb = ecmap(pr, mb, comp.prD, comp.mbD, comp.oob_map, comp.epD, pifK, epfK)
    end
    return ep, ep_pr, ep_mb
end # comp_eff


"""
function comp_spd(comp::compressor, piK, epfK; useTbl=true)

    Wrapper for Ncmap and NcTblMap for use with the compressor struct
"""
function comp_spd(pr, mb, comp::compressor; useTbl=true)
    Nb = 0.0
    Nb_mb = 0.0
    Nb_pr = 0.0

    if useTbl
        Nb, Nb_pr, Nb_mb = NcTblMap(pr, mb, comp.prD, comp.mbD, comp.NbD, comp.map, comp.oob_map)
    else
        Nb, Nb_pr, Nb_mb = Ncmap(pr, mb, comp.prD, comp.mbD, comp.NbD, comp.oob_map)
    end
    return Nb, Nb_pr, Nb_mb
end # comp_eff


"""
function check_in_bounds_reggrid(m, p, compTbl)

    Checks whether a given (mb, pr) point in within the associated map defined as a
    compressorTbl. Uses properties of regular grids to prevent array allocation
"""
function check_in_bounds_reggrid(m, p, comp::compressor)
    # ---- Calculate map scaling factors
    s_pr = (comp.prD - 1) / (comp.map.pr_des - 1)
    s_mb = comp.mbD / comp.map.mb_des

    p_descl = (p - 1) / s_pr + 1
    m_descl = m / s_mb

    miLow, miHigh, mclamp = regGridBoundingIndices(m_descl, comp.map.dm, comp.map.Nm, comp.map.mbGrid[1])
    piLow, piHigh, pclamp = regGridBoundingIndices(p_descl, comp.map.dp, comp.map.Np, comp.map.prGrid[1])

    if (mclamp == 1 || mclamp == -1) || (pclamp == 1 || pclamp == -1)
        return false
    end

    p1 = comp.map.mask[miLow, piLow]
    p2 = comp.map.mask[miLow, piHigh]
    p3 = comp.map.mask[miHigh, piHigh]
    p4 = comp.map.mask[miHigh, piLow]

    return (p1||p2)||(p3||p4)

    # mlow = comp.map.mbGrid[miLow]
    # plow = comp.map.prGrid[piLow]

    # frac_m = (m - mlow) / comp.map.dm
    # frac_p = (p - plow) / comp.map.dp

    # miNear = 0
    # piNear = 0

    # if frac_m < 0.5
    #     miNear = miLow
    # else
    #     miNear = miHigh
    # end

    # if frac_p < 0.5
    #     piNear = piLow
    # else
    #     piNear = piHigh
    # end

    # return comp.map.mask[miNear, piNear]
end # check_in_bounds_reggrid()


"""
function tfoper_residual(res, pf, pl, ph, mf, ml, mh, Tb, Pc, Mi,
    Phiinl, Kinl, iBLIc,
    pid, pib, pifn, pitn,
    Gearf,
    A2, A25, A5, A7,
    iTFspec,
    Ttf, ifuel, hvap, etab,
    pifK, epfK,
    mofft, Pofft,
    Tt9, pt9,
    epsl, epsh,
    icool,
    Mtexit, dTstrk, StA, efilm, tfilm,
    M4a, ruc,
    ncrowx, ncrow,
    epsrow, Tmrow,
    Feng,
    M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
    Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
    Δp_PreC, Δp_InterC, Δp_Regen,
    Tref, pref,
    alpha, beta, gamma, lambda, lambdap,
    fan::compressor, lpc::compressor, hpc::compressor, hpt::turbine, lpt::turbine,
    pt0, Tt0, ht0, st0, cpt0, Rt0, at0, p0, T0, h0, s0, cp0, R0, gam0, u0
    pt18, Tt18, ht18, st18, cpt18, Rt18,
    pt19, Tt19, ht19, st19, cpt19, Rt19)

Used to run the engine cycle without the jacobian calculation
"""
function tfoper_residual(res, pf, pl, ph, mf, ml, mh, Tb, Pc, Mi,
    Phiinl, Kinl, iBLIc,
    pid, pib, pifn, pitn,
    Gearf,
    A2, A25, A5, A7,
    iTFspec,
    Ttf, ifuel, hvap, etab,
    pifK, epfK,
    mofft, Pofft,
    Tt9, pt9,
    epsl, epsh,
    icool,
    Mtexit, dTstrk, StA, efilm, tfilm,
    M4a, ruc,
    ncrowx, ncrow,
    epsrow, Tmrow,
    Feng,
    M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4set, pt5, mcore, M25, 
    Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
    Δp_PreC, Δp_InterC, Δp_Regen,
    Tref, pref,
    alpha, beta, gamma, lambda, lambdap,
    fan::compressor, lpc::compressor, hpc::compressor, hpt::turbine, lpt::turbine,
    pt0, Tt0, ht0, st0, cpt0, Rt0, at0, p0, T0, h0, s0, cp0, R0, gam0, u0,
    pt18, Tt18, ht18, st18, cpt18, Rt18,
    pt19, Tt19, ht19, st19, cpt19, Rt19,
    useTbl, epfmin, n, nair, ncrowy)

    if (iTFspec == 1)
        Tt4spec = Tt4set
    else
        Fspec = Feng
    end

    epf = fan.epD
    eplc = lpc.epD
    ephc = hpc.epD
    epht = hpt.epD
    eplt = lpt.epD

    # ===============================================================
    #---- set fan inlet total pressure pt21 corrected for BLI
    #-    (Tt2,pt2 approximated with Tt0,pt0 here to avoid circular definition)
    if (u0 == 0.0)
            #----- static case... no BL defect
            sbfan = 0.0
            sbcore = 0.0
    else
            #----- account for inlet BLI defect via mass-averaged entropy
            a2sq = at0^2 / (1.0 + 0.5 * (gam0 - 1.0) * Mi^2)

            if (iBLIc == 0)
                #------ BL mixes with fan flow only
                mmix = mf * sqrt(Tref / Tt0) * pt0 / pref

                sbfan = Kinl * gam0 / (mmix * a2sq)
                sbcore = 0.0
            else
                #------ BL mixes with fan + core flow
                mmix = mf * sqrt(Tref / Tt0) * pt0 / pref +
                        ml * sqrt(Tref / Tt0) * pt0 / pref

                sbfan = Kinl * gam0 / (mmix * a2sq)
                sbcore = sbfan
            end

    end

    #---- note: BL is assumed adiabatic, 
    #-     so Tt2,ht2,st2,cpt2,Rt2  will not change due to BL ingestion
    Tt2 = Tt18
    ht2 = ht18
    st2 = st18
    cpt2 = cpt18
    Rt2 = Rt18
    pt2 = pt18 * exp(-sbfan)

    Tt19 = Tt18
    ht19 = ht18
    st19 = st18
    cpt19 = cpt18
    Rt19 = Rt18
    pt19 = pt18 * exp(-sbcore)

    p2, T2, h2, s2, cp2, R2,
    p2_st2,
    p2_pt2,
    dum,
    p2_Tt2, T2_Tt2, h2_Tt2, s2_Tt2,
    p2_ht2, T2_ht2, h2_ht2, s2_ht2,
    p2_Mi, T2_Mi, h2_Mi, s2_Mi,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_machd(alpha, nair,
            pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, Mi, 1.0)
    u2 = sqrt(2.0 * (ht2 - h2))

    rho2 = p2 / (R2 * T2)


    p19, T19, h19, s19, cp19, R19,
    p19_st19,
    p19_pt19,
    dum,
    p19_Tt19, T19_Tt19, h19_Tt19, s19_Tt19,
    p19_ht19, T19_ht19, h19_ht19, s19_ht19,
    p19_Mi, T19_Mi, h19_Mi, s19_Mi,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_machd(alpha, nair,
            pt19, Tt19, ht19, st19, cpt19, Rt19, 0.0, Mi, 1.0)
    u19 = sqrt(2.0 * (ht19 - h19))

    rho19 = p19 / (R19 * T19)

    # ===============================================================
    #---- fan flow 2-7
    epf, epf_pf, epf_mf = comp_eff(pf, mf, fan, pifK, epfK; useTbl=useTbl)
    # epf, epf_pf, epf_mf = ecTblMap(pf, mf, pifD, mbfD, E3fan, epf0)

    if (epf < epfmin)
            epf = epfmin
    end
    if (pf < 1.0)
            epf = 1.0 / epf
    end

    pt21, Tt21, ht21, st21, cpt21, Rt21,
    pt21_pt2,
    pt21_st2, Tt21_st2, ht21_st2, st21_st2,
    pt21_pf, Tt21_pf, ht21_pf, st21_pf,
    pt21_epf, Tt21_epf, ht21_epf, st21_epf,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_pratd(alpha, nair,
            pt2, Tt2, ht2, st2, cpt2, Rt2, pf, epf)

    #---- fan duct nozzle total quantities
    pt7 = pt21 * pifn
    Tt7 = Tt21
    ht7 = ht21
    st7 = st21
    cpt7 = cpt21
    Rt7 = Rt21

    # ===============================================================
    #---- Compressor precooler 19-19c
    pt19c = pt19 - Δp_PreC
    ht19c = ht19 + Δh_PreC
    Tt19c, Tt19c_ht19c, _ = gas_tsetd(alpha, nair, ht19c, Tt19)
    st19c, st19c_Tt19c, ht19c, ht29c_Tt29c, cpt19c, cpt19c_Tt19c, Rt19c = gassumd(alpha, nair, Tt19c)

    p19c, T19c, h19c, s19c, cp19c, R19c,
    p19c_st19c,
    p19c_pt19c,
    dum,
    p19c_Tt19c, T19c_Tt19c, h19c_Tt19c, s19c_Tt19c,
    p19c_ht19c, T19c_ht19c, h19c_ht19c, s19c_ht19c,
    p19c_Mi, T19c_Mi, h19c_Mi, s19c_Mi,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_machd(alpha, nair,
            pt19c, Tt19c, ht19c, st19c, cpt19c, Rt19c, 0.0, Mi, 1.0)
    u19c = sqrt(2.0 * (ht19c - h19c))

    rho19c = p19c / (R19c * T19c)

    #--------------------------------------------------------------
    #---- offtake mass ratio
    #cc   fo = mofft / mcore
    fo = mofft / ml * sqrt(Tt19c / Tref) * pref / pt19c

    #---- normalized power offtake Pofft / mcore
    Pom = Pofft / ml * sqrt(Tt19c / Tref) * pref / pt19c

    # ===============================================================
    #---- LP compressor flow 2-25
    eplc, eplc_pl, eplc_ml = comp_eff(pl, ml, lpc, 0.0, 0.0; useTbl=useTbl)
    # eplc, eplc_pl, eplc_ml = ecTblMap(pl, ml, pilcD, mblcD, E3lpc, eplc0)

    if (eplc < 0.70)
            eplc = 0.70
    end

    pt25, Tt25, ht25, st25, cpt25, Rt25,
    pt25_pt19c,
    pt25_st19c, Tt25_st19c, ht25_st19c, st25_st19c,
    pt25_pl, Tt25_pl, ht25_pl, st25_pl,
    pt25_eplc, Tt25_eplc, ht25_eplc, st25_eplc,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_pratd(alpha, nair,
            pt19c, Tt19c, ht19c, st19c, cpt19c, Rt19c, pl, eplc)

    # ===============================================================
    #---- Compressor intercooler 25-25c
    pt25c = pt25 - Δp_InterC
    ht25c = ht25 + Δh_InterC
    Tt25c, Tt25c_ht25c, _ = gas_tsetd(alpha, nair, ht25c, Tt25)
    st25c, st25c_Tt25c, ht25c, ht25c_Tt25c, cpt25c, cpt25c_Tt25c, Rt25c = gassumd(alpha, nair, Tt25c)

    # ===============================================================
    #---- HP compressor flow 25-3
    ephc, ephc_ph, ephc_mh = comp_eff(ph, mh, hpc, 0.0, 0.0; useTbl=useTbl)
    # ephc, ephc_ph, ephc_mh = ecTblMap(ph, mh, pihcD, mbhcD, E3hpc, ephc0)
    if (ephc < 0.70)
            ephc = 0.70
    end

    pt3, Tt3, ht3, st3, cpt3, Rt3,
    pt3_pt25c,
    pt3_st25c, Tt3_st25c, ht3_st25c, st3_st25c,
    pt3_ph, Tt3_ph, ht3_ph, st3_ph,
    pt3_ephc, Tt3_ephc, ht3_ephc, st3_ephc,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_pratd(alpha, nair,
            pt25c, Tt25c, ht25c, st25c, cpt25c, Rt25c, ph, ephc)

    # ===============================================================
    #---- burner fuel flow from Tt3-Tb difference   (ffb = m_fuel/m_burner)
    ffb, lambda,
    ffb_Tt3, ffb_Ttf, ffb_Tb,
    lam_Tt3, lam_Ttf, lam_Tb = gas_burnd(alpha, beta, gamma, n, ifuel, Tt3, Ttf, Tb, hvap)

    #---- all station 4 quantities
    Tt4 = Tb
    Tt4_Tb = 1.0
    st4, st4_Tt4,
    ht4, ht4_Tt4,
    cpt4, cpt4_Tt4, Rt4 = gassumd(lambda, nair, Tt4)

    pt4 = pib * pt3

    # HSC: SEEMS TO BE FINE
    # ===============================================================
    if (icool == 0)
            #----- no cooling air present... station 41 is same as 4
            pt41 = pt4
            Tt41 = Tt4
            ht41 = ht4
            st41 = st4
            cpt41 = cpt4
            Rt41 = Rt4
            for i = 1:nair
                lambdap[i] = lambda[i]
            end

            #
            fc = 0.0

            #----- set ff = m_fuel/m_core = ffb * m_burner/m_core
            ff = (1.0 - fo) * ffb

            #----------------------------------------------------------------
    else
            #----- cooling air is present... 

            #----- hot-section temperature ratio for each blade row (for cooling model)
            gmi4 = Rt4 / (cpt4 - Rt4)
            Trrat = 1.0 / (1.0 + 0.5 * gmi4 * Mtexit^2)

            # Heat exchanger to cool turbine cooling air
            ht_tc = ht3 + Δh_TurbC #Specific enthalpy of turbine cooling air
            Tt_tc, Tttc_httc, _ = gas_tsetd(alpha, nair, ht_tc, Tt3) #Temperature of turbine cooling air

            if (icool == 1)
                #------ epsrow(.) is assumed to be passed in.. calculate Tmrow(.)
                Tmrow_copy = Tmcalc(ncrowx, ncrow,
                        Tt_tc, Tb, dTstrk, Trrat,
                        efilm, tfilm, StA, epsrow)
                Tmrow[:] = Tmrow_copy[:]

                #------ total cooling flow fraction
                fc = 0.0
                for icrow = 1:ncrow
                        fc = fc + (1.0 - fo) * epsrow[icrow]
                end

            else
                #------ calculate cooling mass flow ratios epsrow(.) to get specified Tmrow(.)
                ncrow, epsrow_copy, epsrow_Tttc, epsrow_Tb, epsrow_Trr = mcool(ncrowx,
                        Tmrow, Tt_tc, Tb, dTstrk, Trrat,
                        efilm, tfilm, StA)
                epsrow[:] = epsrow_copy[:]

                #------ total cooling flow fraction
                fc = 0.0
                for icrow = 1:ncrow
                        fc = fc + (1.0 - fo) * epsrow[icrow]
                end

            end

            #----- set ff = m_fuel/m_core = ffb * m_burner/m_core
            ff = (1.0 - fo - fc) * ffb

            #----- calculate station 41
            pt4a = pt4
            Tt4a = Tt4
            ht4a = ht4
            st4a = st4
            cpt4a = cpt4
            rt4a = Rt4

            #----- speed at start-of-mixing station 4a
            p4a, t4a, h4a, s4a, cp4a, r4a,
            p4a_st4a,
            p4a_pt4a,
            dum,
            p4a_Tt4a, T4a_Tt4a, h4a_Tt4a, s4a_Tt4a,
            p4a_ht4a, T4a_ht4a, h4a_ht4a, s4a_ht4a,
            p4a_M4a, T4a_M4a, h4a_M4a, s4a_M4a,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_machd(lambda, nair,
                pt4a, Tt4a, ht4a, st4a, cpt4a, rt4a, 0.0, M4a, 1.0)

            if (ht4a > h4a)
                u4a = sqrt(2.0 * (ht4a - h4a))
            else
                u4a = 0.0
            end

            #----- exit speed of cooling air at station 4a
            uc = ruc * u4a

            #----- IGV exit mixing
            frac4 = (1.0 - fo - fc + ff) / (1.0 - fo + ff)
            fracm = fc / (1.0 - fo + ff)

            #----- mixed constituent fraction vector from mass equation
            for i = 1:nair
                lambdap[i] = frac4 * lambda[i] + fracm * alpha[i]
            end

            #----- mixed total enthalpy from enthalpy equation
            ht41 = frac4 * ht4 + fracm * ht_tc

            #----- total temperature from total enthalpy
            Tguess = frac4 * Tt4 + fracm * Tt_tc
            Tt41, Tt41_ht41, T_al = gas_tsetd(lambdap, nair, ht41, Tguess)

            #----- will also need st41,cpt41,Rt41 derivatives
            st41, st41_Tt41,
            ht41, ht41_Tt41,
            cpt41, cpt41_Tt41, Rt41 = gassumd(lambdap, nair, Tt41)

            #----- mixed velocity from momentum equation, assuming constant static pressure
            p41 = p4a
            u41 = frac4 * u4a + fracm * uc

            #----- static temperature from static enthalpy
            h41 = ht41 - 0.5 * u41^2

            Tguess = t4a + (h41 - h4a) / cp4a
            T41, T41_h41, T_al = gas_tsetd(lambdap, nair, h41, Tguess)
            s41, s41_T41, h41, dum, cp41, R41 = gassum(lambdap, nair, T41)

            #----- all stagnation quantities, from total-static enthalpy difference
            dhb = ht41 - h41
            epi = 1.0
            pt41, Tt41, ht41, st41, cpt41, Rt41,
            pt41_s41,
            pt41_p41,
            dum,
            pt41_h41, Tt41_h41, ht41_h41, st41_h41,
            pt41_dhb, Tt41_dhb, ht41_dhb, st41_dhb,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_delhd(lambdap, nair,
                p41, T41, h41, s41, cp41, R41, dhb, epi)
    end

    # ===============================================================
    #---- fan corrected speed
    Nf, Nf_pf, Nf_mf = comp_spd(pf, mf, fan; useTbl=useTbl)

    #---- LPC corrected speed
    Nl, Nl_pl, Nl_ml = comp_spd(pl, ml, lpc; useTbl=useTbl)

    #---- HPC corrected speed
    Nh, Nh_ph, Nh_mh = comp_spd(ph, mh, hpc; useTbl=useTbl)

    # ===============================================================
    #---- HPT and LPT work

    #---- bypass ratio
    BPR = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c

    #---- HPT work
    dhfac = -(1.0 - fo) / (1.0 - fo + ff) / (1.0 - epsh)
    dhht = (ht3 - ht25c) * dhfac

    #---- LPT work
    dlfac = -1.0 / (1.0 - fo + ff) / (1.0 - epsl)
    dhlt = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac

    #---- HPT corrected mass flow, using LPC corrected mass flow and fuel/air ratio
    mbht = ml * (1.0 - fo + ff) * sqrt(Tt41 / Tt19c) * pt19c / pt41

    #---- HPT corrected speed
    Nbht = Nh * sqrt(Tt25c / Tt41)
    Nbht_Nh = sqrt(Tt25c / Tt41)

    #---- HPT efficiency
    # HACK: HSC
    epht, epht_dhht, epht_mbht, epht_Nbht, epht_Tt41,
    epht_cpt41, epht_Rt41 = etmap(dhht, mbht, Nbht, hpt.prD, hpt.mbD, hpt.NbD, hpt.epD, hpt.map,
            Tt41, cpt41, Rt41)

    if (epht < 0.80)
            epht = 0.80
    end

    #---- HPT work to determine station 45
    epi = 1.0 / epht
    pt45, Tt45, ht45, st45, cpt45, Rt45,
    pt45_st41,
    pt45_pt41,
    pt45_epi,
    pt45_ht41, Tt45_ht41, ht45_ht41, st45_ht41,
    pt45_dhht, Tt45_dhht, ht45_dhht, st45_dhht,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_delhd(lambdap, nair,
            pt41, Tt41, ht41, st41, cpt41, Rt41, dhht, epi)

    #---- will also need cpt45,Rt45 derivatives
    st45, st45_Tt45,
    ht45, ht45_Tt45,
    cpt45, cpt45_Tt45, Rt45 = gassumd(lambdap, nair, Tt45)

    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #---- LPT corrected mass flow (in terms of Newton variables)
    mblt = ml * (1.0 - fo + ff) * sqrt(Tt45 / Tt19c) * pt19c / pt45

    #---- LPT corrected speed
    Nblt = Nl * sqrt(Tt19c / Tt45)
    Nblt_Nl = sqrt(Tt19c / Tt45)

    #---- LPT efficiency
    eplt,
    eplt_dhlt, eplt_mblt, eplt_Nblt,
    eplt_Tt45, eplt_cpt45, eplt_Rt45 = etmap(dhlt, mblt, Nblt, lpt.prD, lpt.mbD, lpt.NbD, lpt.epD, lpt.map,
            Tt45, cpt45, Rt45)

    if (eplt < 0.80)
        eplt = 0.80
    end

    if (Pc == 0.0)
            #----- initial guess for pt5
            epi = 1.0 / eplt
            pt49, Tt49, ht49, st49, cpt49, Rt49 = gas_delh(lambdap, nair,
                pt45, Tt45, ht45, st45, cpt45, Rt45, dhlt, epi)
            Pc = pt49 * pitn
            Pc = max(Pc, p0 * (1.0 + 0.2 * M0^2)^3.5)
    end

    Pc = max(Pc, 1.000001 * p0)

    pilt = Pc / pt45 / pitn
    pilt_Pc = 1.0 / pt45 / pitn
    pilt_pt45 = -pilt / pt45

    epi = 1.0 / eplt
    epi_eplt = -epi / eplt
    pt49, Tt49, ht49, st49, cpt49, Rt49,
    pt49_pt45,
    pt49_st45, Tt49_st45, ht49_st45, st49_st45,
    pt49_pilt, Tt49_pilt, ht49_pilt, st49_pilt,
    pt49_epi, Tt49_epi, ht49_epi, st49_epi,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_pratd(lambdap, nair,
            pt45, Tt45, ht45, st45, cpt45, Rt45, pilt, epi)

    # ===============================================================
    #---- Regenerative cooling heat exchanger 49-49c
    pt49c = pt49 - Δp_Regen
    ht49c = ht49 + Δh_Regen
    Tt49c, Tt49c_ht49c, _ = gas_tsetd(lambdap, nair, ht49c, Tt49)
    st49c, st49c_Tt49c, ht49c, ht49c_Tt49, cpt49c, cpt49c_Tt49c, Rt49c = gassumd(lambdap, nair, Tt49c)

    #---- apply core nozzle loss
    pt5 = pt49c * pitn
    Tt5 = Tt49c
    ht5 = ht49c
    st5 = st49c
    cpt5 = cpt49c
    Rt5 = Rt49c

    # ===============================================================
    #---- fan nozzle flow 7-8, use alpha mass fraction (air)
    pfn = p0 / pt7
    pfn_pt7 = -pfn / pt7
    p7, T7, h7, s7, cp7, R7,
    p7_pt7,
    p7_st7, T7_st7, h7_st7, s7_st7,
    p7_pfn, T7_pfn, h7_pfn, s7_pfn,
    dum, dum, dum, dum,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_pratd(alpha, nair,
            pt7, Tt7, ht7, st7, cpt7, Rt7, pfn, 1.0)
    u7 = sqrt(2.0 * max(ht7 - h7, 0.0))
    M7 = u7 / sqrt(T7 * cp7 * R7 / (cp7 - R7))

    if (M7 >= 1.0)
        #----- fan nozzle is choked... specify M7 = 1 instead
        M7 = 1.0
        p7, T7, h7, s7, cp7, R7,
        p7_st7,
        p7_pt7,
        dum,
        p7_Tt7, T7_Tt7, h7_Tt7, s7_Tt7,
        p7_ht7, T7_ht7, h7_ht7, s7_ht7,
        p7_M7, T7_M7, h7_M7, s7_M7,
        p_al, T_al, h_al, s_al,
        cp_al, R_al = gas_machd(alpha, nair,
            pt7, Tt7, ht7, st7, cpt7, Rt7, 0.0, M7, 1.0)
    end

    if (ht7 > h7)
            u7 = sqrt(2.0 * (ht7 - h7))
    else
            u7 = 0.0
    end

    rho7 = p7 / (R7 * T7)

    # ===============================================================
    #---- core nozzle flow 5-6, use lambda mass fraction (combustion products)
    pcn = p0 / pt5
    pcn_pt5 = -pcn / pt5

    p5, T5, h5, s5, cp5, R5,
    p5_pt5,
    p5_st5, T5_st5, h5_st5, s5_st5,
    p5_pcn, T5_pcn, h5_pcn, s5_pcn,
    dum, dum, dum, dum,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_pratd(lambdap, nair,
            pt5, Tt5, ht5, st5, cpt5, Rt5, pcn, 1.0)

    u5 = sqrt(2.0 * max(ht5 - h5, 0.0))
    M5 = u5 / sqrt(T5 * cp5 * R5 / (cp5 - R5))

    if (M5 >= 1.0)
        #----- core nozzle is choked... specify M5 = 1 instead

        M5 = 1.0
        p5, T5, h5, s5, cp5, R5,
        p5_st5,
        p5_pt5,
        dum,
        p5_Tt5, T5_Tt5, h5_Tt5, s5_Tt5,
        p5_ht5, T5_ht5, h5_ht5, s5_ht5,
        p5_M5, T5_M5, h5_M5, s5_M5,
        p_al, T_al, h_al, s_al,
        cp_al, R_al = gas_machd(lambdap, nair,
            pt5, Tt5, ht5, st5, cpt5, Rt5, 0.0, M5, 1.0)
    end

    if (ht5 > h5)
            u5 = sqrt(2.0 * (ht5 - h5))
    else
            u5 = 0.0
    end

    rho5 = p5 / (R5 * T5)

    #---- Calculate Residuals

    #---- fan/LPC speed constraint
    trf = sqrt(Tt2 / Tref)
    trl = sqrt(Tt19c / Tref)
    res[1, 1] = Gearf * trf * Nf - trl * Nl

    #---- fixed corrected mass flow at LPT IGV (vertical-line LPT map)
    res[2, 1] = mblt - lpt.mbD

    #---- fixed corrected mass flow at HPT IGV (vertical-line HPT map)
    res[3, 1] = mbht - hpt.mbD

    #---- fan nozzle mass flow, choked or unchoked
    mdotf = mf * sqrt(Tref / Tt2) * pt2 / pref
    res[4, 1] = mdotf - rho7 * u7 * A7

    #---- core nozzle mass flow, choked or unchoked
    mdotc = (1.0 - fo + ff) * ml * sqrt(Tref / Tt19c) * pt19c / pref
    res[5, 1] = mdotc - rho5 * u5 * A5

    #---- LPC mass flow = HPC mass flow + mass offtake
    mdotl = ml * sqrt(Tref / Tt19c) * pt19c / pref
    mdoth = mh * sqrt(Tref / Tt25c) * pt25c / pref
    res[6, 1] = mdoth - mdotl + mofft

    if (iTFspec == 1)
        #----- specified Tt4 constraint
        res[7, 1] = Tt4 - Tt4spec
    else
        #----- specified thrust constraint

        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        #----- fan nozzle flow 7-8, use alpha mass fraction (air)
        pfn = p0 / pt7
        epi = 1.0
        p8, T8, h8, s8, cp8, R8,
        p8_pt7,
        p8_st7, T8_st7, h8_st7, s8_st7,
        p8_pfn, T8_pfn, h8_pfn, s8_pfn,
        p8_epi, T8_epi, h8_epi, s8_epi,
        p_al, T_al, h_al, s_al,
        cp_al, R_al = gas_pratd(alpha, nair,
              pt7, Tt7, ht7, st7, cpt7, Rt7, pfn, epi)

        if (ht7 > h8)
              u8 = sqrt(2.0 * (ht7 - h8))
        else
              u8 = 0.0
        end
        rho8 = p8 / (R8 * T8)

        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        #----- core nozzle flow 5-6, use lambda mass fraction (combustion products)
        pcn = p0 / pt5
        epi = 1.0
        p6, T6, h6, s6, cp6, R6,
        p6_pt5,
        p6_st5, T6_st5, h6_st5, s6_st5,
        p6_pcn, T6_pcn, h6_pcn, s6_pcn,
        p6_epi, T6_epi, h6_epi, s6_epi,
        p_al, T_al, h_al, s_al,
        cp_al, R_al = gas_pratd(lambdap, nair,
              pt5, Tt5, ht5, st5, cpt5, Rt5, pcn, epi)

        if (ht5 > h6)
              u6 = sqrt(2.0 * (ht5 - h6))
        else
              u6 = 0.0
        end

        BPR = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c
        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        #----- overall thrust
        if (u0 == 0.0)
              Finl = 0.0
        else
              Finl = Phiinl / u0
        end
        F = ((1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9) * mdotl + Finl

        #----- specified-F constraint
        res[7, 1] = F - Fspec
    end

    #---- LPT power matched to fan+LPC power
    epi = 1.0 / eplt
    pt5h, Tt5h, ht5h, st5h, cpt5h, Rt5h,
    pt5h_st45,
    pt5h_pt45,
    pt5h_epi,
    pt5h_ht45, Tt5h_ht45, ht5h_ht45, st5h_ht45,
    pt5h_dhlt, Tt5h_dhlt, ht5h_dhlt, st5h_dhlt,
    p_al, T_al, h_al, s_al,
    cp_al, R_al = gas_delhd(lambdap, nair,
          pt45, Tt45, ht45, st45, cpt45, Rt45, dhlt, epi)

    res[8, 1] = pt49 - pt5h

    #---- inlet Mach - area constraint
    mfA = mf * sqrt(Tref / Tt2) * pt2 / pref / (rho2 * u2)
    mlA = ml * sqrt(Tref / Tt19c) * pt19c / pref / (rho19c * u19c)
    res[9, 1] = mfA + mlA - A2


    return res

end # tfoper_residual