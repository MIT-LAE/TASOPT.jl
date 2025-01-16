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
    a = zeros(T, 9, 9)
    rrel = zeros(T, 9)
    rsav = zeros(T, 9)
    asav = zeros(T, 9, 10)

    res_dlls = zeros(T, 9)
    a_dlls = zeros(T, 9, 9)

    return res, a, rrel, rsav, asav, res_dlls, a_dlls
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
function freestream_gasprop(alpha, nair, T0, p0, M0, a0)
    #---- freestream static quantities
    s0, dsdt, h0, dhdt, cp0, R0 = gassum(alpha, nair, T0)
    gam0 = cp0 / (cp0 - R0)
    #      a0 = sqrt(gam0*R0*T0)
    u0 = M0 * a0

    # ===============================================================
    #---- freestream total quantities
    hspec = h0 + 0.5 * u0^2
    Tguess = T0 * (1.0 + 0.5 * (gam0 - 1.0) * M0^2)
    Tt0 = gas_tset(alpha, nair, hspec, Tguess)
    st0, dsdt, ht0, dhdt, cpt0, Rt0 = gassum(alpha, nair, Tt0)
    pt0 = p0 * exp((st0 - s0) / Rt0)
    at0 = sqrt(Tt0 * Rt0 * cpt0 / (cpt0 - Rt0))

    return s0,  dsdt, h0,  dhdt, cp0,  R0,
           st0, dsdt, ht0, dhdt, cpt0, Rt0,
           gam0, u0, hspec, Tguess, Tt0, pt0, at0

end # freestream_gasprop


"""
function design_initialization(pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, M2,
    pifD, pilcD, pihcD, mbfD, mblcD, mbhcD)

Generate guesses for primary Newton variables
"""
function design_initialization(pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, M2,
    pifD, pilcD, pihcD, mbfD, mblcD, mbhcD)

    pf = pif
    pl = pilc
    ph = pihc
    mf = mbf
    ml = mblc
    mh = mbhc
    Tb = Tt4
    Pc = pt5
    Mi = M2

    if (pf == 0.0)
        pf = pifD
    end
    if (pl == 0.0)
        pl = pilcD
    end
    if (ph == 0.0)
        ph = pihcD
    end
    if (mf == 0.0)
        mf = mbfD
    end
    if (ml == 0.0)
        ml = mblcD
    end
    if (mh == 0.0)
        mh = mbhcD
    end
    if (Mi == 0.0)
        Mi = 0.6
    end
    if (Mi == 1.0)
        Mi = 0.6
    end

    return pf, pl, ph, mf, ml, mh, Tb, Pc, Mi

end # design_initialization


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
function boundary_layer_calculations(at0, gam0, Mi, iBLIc, mf, ml, Tref, pref, Tt0, pt0, Kinl)
    if (u0 == 0.0)
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
            a2sq = at0^2 / (1.0 + 0.5 * (gam0 - 1.0) * Mi^2)
            a2sq_Mi = -a2sq / (1.0 + 0.5 * (gam0 - 1.0) * Mi^2) *
                    (gam0 - 1.0) * Mi

            if (iBLIc == 0)
                #------ BL mixes with fan flow only
                #c      mmix    = mf*sqrt(Tref/Tt2) * pt2   /pref
                #c      mmix_mf =    sqrt(Tref/Tt2) * pt2   /pref
                #c      mmix_Mi = mf*sqrt(Tref/Tt2) * pt2_Mi/pref

                mmix = mf * sqrt(Tref / Tt0) * pt0 / pref
                mmix_mf = sqrt(Tref / Tt0) * pt0 / pref
                mmix_Mi = 0.0

                sbfan = Kinl * gam0 / (mmix * a2sq)
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

                mmix = mf * sqrt(Tref / Tt0) * pt0 / pref +
                        ml * sqrt(Tref / Tt0) * pt0 / pref
                mmix_mf = sqrt(Tref / Tt0) * pt0 / pref
                mmix_ml = sqrt(Tref / Tt0) * pt0 / pref
                mmix_Mi = 0.0

                sbfan = Kinl * gam0 / (mmix * a2sq)
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

    Tt2 = Tt18
    ht2 = ht18
    st2 = st18
    cpt2 = cpt18
    Rt2 = Rt18
    pt2 = pt18 * exp(-sbfan)
    pt2_mf = pt2 * (-sbfan_mf)
    pt2_ml = pt2 * (-sbfan_ml)
    pt2_Mi = pt2 * (-sbfan_Mi)

    Tt19 = Tt18
    ht19 = ht18
    st19 = st18
    cpt19 = cpt18
    Tt19_ht19 = 1 / cpt19
    Rt19 = Rt18
    pt19 = pt18 * exp(-sbcore)
    pt19_mf = pt19 * (-sbcore_mf)
    pt19_ml = pt19 * (-sbcore_ml)
    pt19_Mi = pt19 * (-sbcore_Mi)

    return Tt2,  ht2,  st2,  cpt2,             Rt2,  pt2,  pt2_mf,  pt2_ml,  pt2_Mi,
           Tt19, ht19, st19, cpt19, Tt19_ht19, Rt19, pt19, pt19_mf, pt19_ml, pt19_Mi

end # boundary_layer_calculations


function compressor(pf, mf, pifD, mbfD, Cmap)
    return 0
end # compressor