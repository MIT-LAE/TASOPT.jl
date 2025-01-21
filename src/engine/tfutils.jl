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
function freestream_gasprop(alpha, nair, fs::flowStation)
    #---- freestream static quantities
    fs = gassum_fs(alpha, nair, fs; type="static")
    fs.u = fs.M * fs.as

    # ===============================================================
    #---- freestream total quantities
    hspec = fs.hs + 0.5 * fs.u^2
    Tguess = fs.Ts * (1.0 + 0.5 * (fs.gams - 1.0) * fs.M^2)
    fs.Tt = gas_tset(alpha, nair, hspec, Tguess)
    fs = gassum_fs(alpha, nair, fs; type="total")
    fs.pt = fs.ps * exp((fs.st - fs.ss) / fs.Rt)
    fs.at = sqrt(fs.Tt * fs.Rt * fs.gamt)

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
function offtake_plume(fs0::flowStation, fs9::flowStation)
    Trat = (fs0.ps / fs9.pt)^(fs0.Rt / fs0.cpt)
    if (Trat < 1.0)
          fs9.u = sqrt(2.0 * fs0.cpt * fs9.Tt * (1.0 - Trat))
          fs9.rhos = fs0.ps / (fs0.Rt * fs0.Tt * Trat)
    else
          fs9.u = 0.0
          fs9.rhos = p0 / (Rt0 * Tt0)
    end

    return fs9
end # offtake_plume


"""
function design_initialization(fan, lpc, hpc, Tt4, pt5, M2)

Generate guesses for primary Newton variables
"""
function design_initialization(Tt4, pt5, M2)
    Tb = Tt4
    Pc = pt5
    Mi = M2

    if (Mi == 0.0)
        Mi = 0.6
    end
    if (Mi == 1.0)
        Mi = 0.6
    end

    return Tb, Pc, Mi

end # design_initialization


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
function comp_eff(comp::compressor, piK, epfK; useTbl=true)
    if useTbl
        comp.ep, comp.ep_pr, comp.ep_mb = ecTblMap(comp.pr, comp.mb, comp.prD, comp.mbD, comp.map, comp.epD)
    else
        comp.ep, comp.ep_pr, comp.ep_mb = ecmap(comp.pr, comp.mb, comp.prD, comp.mbD, comp.oob_map, comp.epD, pifK, epfK)
    end
end # comp_eff
