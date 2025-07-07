"""
    Ncmap(pratio, mb, piD, mbD, NbD, Cmap)

Calculates compressor or fan corrected speed as a function of pressure ratio and corrected mass flow

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pratio`: pressure ratio 
    - `mb`:       corrected mass flow
    - `piD`:      design pressure ratio
    - `mbD`:      design corrected mass flow
    - `NbD`:      design corrected speed
    - `Cmap(.)`:  map constants

    **Outputs:**
    - `Nb`:     wheel speed
    - `Nb_?`:   derivatives

"""
function Ncmap(pratio, mb, piD, mbD, NbD, Cmap)

    eps = 1.0e-11

    #----              a     b     k     mo     da    c    d     C    D
    #     data Cmapf / 3.50, 0.80, 0.03, 0.75, -0.50, 3.0, 6.0,  2.5, 15.0 ]
    a = Cmap[1]
    b = Cmap[2]
    k = Cmap[3]
    mrato = Cmap[4]
    da = Cmap[5]
    c = Cmap[6]
    d = Cmap[7]
    CK = Cmap[8]
    DK = Cmap[9]

    #---- corrected mass flow / design corrected mass flow
    m = mb / mbD
    m_mb = 1.0 / mbD

    #---- scaled pressure ratio
    p = (pratio - 1.0) / (piD - 1.0)
    p_pi = 1.0 / (piD - 1.0)

    #---- value of p on spine at specified m
    psm = m^a

    #---- corrected speed ratio initial guess from spine
    N = m^(1.0 / b)

    #---- converge to actual corrected speed by inverting compressor map function
    for iter = 1:20

        ms = N^b
        ms_N = b * ms / N

        ps = N^(a * b)
        ps_N = a * b * ps / N

        if (p >= psm)
            #------- above spine... specify p
            #           write(*,*) 'above'
            res = ps + 2.0 * N * k * log(1.0 - (m - ms) / k) - p
            res_N = ps_N + 2.0 * k * log(1.0 - (m - ms) / k) +
                    2.0 * N * k / (1.0 - (m - ms) / k) * ms_N / k
            res_m = -2.0 * N * k / (1.0 - (m - ms) / k) * 1.0 / k
            res_p = -1.0

        else
            #------- below spine... specify m
            #           write(*,*) 'below'
            res = ms + k * (1.0 - exp((p - ps) / (2.0 * N * k))) - m
            res_N = ms_N + k * (-exp((p - ps) / (2.0 * N * k))) *
                           (-ps_N / (2.0 * N * k) - (p - ps) / (2.0 * N * k) / N)
            res_m = -1.0
            res_p = k * (-exp((p - ps) / (2.0 * N * k))) / (2.0 * N * k)

        end

        dN = -res / res_N

        rlx = 1.0
        dNlim = 0.05
        if (rlx * dN > dNlim)
            rlx = dNlim / dN
        end
        if (rlx * dN < -dNlim)
            rlx = -dNlim / dN
        end

        if (abs(dN) < eps)
            #---- N( m[mb] , p[pratio] )
            N_m = -res_m / res_N
            N_p = -res_p / res_N

            #---- N( mb, pratio )
            N_mb = N_m * m_mb
            N_pi = N_p * p_pi

            #---- Nb( N(mb,pratio) )
            Nb = NbD * N
            Nb_mb = NbD * N_mb
            Nb_pi = NbD * N_pi

            return Nb, Nb_pi, Nb_mb

        end

        N = N + rlx * dN
    end
    println("Ncmap: conv.failed.  N, dN =", N, dN)


end # Ncmap

"""
    ecmap(pratio, mb, piD, mbD, Cmap, effo, piK, effK)

Calculates compressor or fan efficiency as a function of pressure ratio and corrected mass flow

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pratio`:        pressure ratio
    - `mb`:        corrected mass flow
    - `piD`:       design pressure ratio
    - `mbD`:       design corrected mass flow
    - `Cmap(.)`:   map constants
    - `effo`:      maximum efficiency
    - `piK`:       pi-dependence offset  eff = effo + effK*(pi-piK)
    - `effK`:      pi-dependence slope

    **Outputs:**
    - `eff`:        efficiency
    - `eff_?`:      derivatives

"""
function ecmap(pratio, mb, piD, mbD, Cmap, effo, piK, effK)

    eps = 1.0e-8

    #----              a     b     k     mo     da    c    d     C    D
    #     data Cmapf / 3.50, 0.80, 0.03, 0.75, -0.50, 3.0, 6.0,  2.5, 15.0 ]
    a = Cmap[1]
    b = Cmap[2]
    k = Cmap[3]
    mo = Cmap[4]
    da = Cmap[5]
    c = Cmap[6]
    d = Cmap[7]
    CK = Cmap[8]
    DK = Cmap[9]

    adm = a + da - 1.0

    #---- scaled pressure ratio
    p = (pratio - 1.0) / (piD - 1.0)
    p_pi = 1.0 / (piD - 1.0)

    #---- corrected mass flow / design correct mass flow
    m = mb / mbD
    m_mb = 1.0 / mbD

    # psgn = sign(1.0, p / m^adm - m)
    if (p / m^adm - m >= 0.0)
        psgn = 1.0
    else
        psgn = -1.0
    end

    # msgn = sign(1.0, m / mo - 1.0)
    if (m / mo - 1.0 >= 0.0)
        msgn = 1.0
    else
        msgn = -1.0
    end


    #c    eo  = 1.0 - DK*(  abs(1.0 /mo   - 1.0))^d
    eo = 1.0

    #---- eff( m(mdot,Tr,pr) , p[pratio] , pi )
    e1 = 1.0 - CK * (psgn * (p / m^adm - m))^c -
         DK * (msgn * (m / mo - 1.0))^d

    e1_p = -CK * (psgn * (p / m^adm - m))^(c - 1.0) *
           c * psgn / m^adm

    e1_m = -CK * (psgn * (p / m^adm - m))^(c - 1.0) *
           c * psgn * (-adm * p / m^adm / m - 1.0) -
           DK * (msgn * (m / mo - 1.0))^(d - 1.0) *
           d * msgn / mo
    eff = effo * e1 / eo + effK * (pratio - piK)
    eff_p = effo * e1_p / eo
    eff_m = effo * e1_m / eo
    eff_pi = effK

    #---- eff( pi, nb )
    eff_pi = eff_p * p_pi + eff_pi
    eff_mb = eff_m * m_mb

    return eff, eff_pi, eff_mb
end # ecmap

"""
    Ncmap1(pratio, m, piD, mbD, NbD, ABCDm, iabcd, Tr, pr)

Calculates compressor or fan efficiency as a function of pressure ratio and corrected mass flow

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pratio`:        pressure ratio
    - `mb`:        corrected mass flow
    - `piD`:       design pressure ratio
    - `mbD`:       design corrected mass flow
    - `NbD`:      design corrected speed
    - `ABCDm`:   map constants
    - `iabcd`:   map exponents
    - `Tr`:      T/Tref
    - `pr`:      p/pref

    **Outputs:**
    - `N`:        wheel speed
    - `N_?`:      derivatives

"""
function Ncmap1(pratio, m, piD, mbD, NbD, ABCDm, iabcd, Tr, pr)

    A = ABCDm[1]
    B = ABCDm[2]

    ia = iabcd[1]
    ib = iabcd[2]

    binv = 1.0 / float(ib)

    #---- corrected mass flow
    mb = m * sqrt(Tr) / pr
    mb_m = sqrt(Tr) / pr
    mb_Tr = 0.5 * mb / Tr
    mb_pr = -mb / pr

    #---- calculate inverse map
    prat = (pratio - 1.0) / (piD - 1.0)
    mrat = mb / mbD

    fac = (prat - 1.0 + A * (mrat^ia - 1.0)) / B + 1.0
    fac_prat = 1.0 / B
    fac_mrat = A * mrat^(ia - 1) * float(ia) / B

    Nb = NbD * fac^binv
    Nb_prat = binv * Nb / fac * fac_prat
    Nb_mrat = binv * Nb / fac * fac_mrat

    Nb_pi = Nb_prat / (piD - 1.0)
    Nb_mb = Nb_mrat / mbD

    #---- N( mb(m,Tr,pr) , pi, Tr )
    N = Nb * sqrt(Tr)
    N_mb = Nb_mb * sqrt(Tr)
    N_pi = Nb_pi * sqrt(Tr)
    N_Tr = Nb * 0.5 / sqrt(Tr)

    #---- N( m, pratio, Tr, pr )
    N_m = N_mb * mb_m
    #cc   N_pi =              N_pi
    N_Tr = N_mb * mb_Tr + N_Tr
    N_pr = N_mb * mb_pr

    return N, N_pi, N_m, N_Tr, N_pr
end # Ncmap1

"""
    ecmap1(pratio, m, piD, mbD, ABCDm, iabcd, effo, Tr, pr)

Calculates compressor or fan efficiency as a function of pressure ratio and mass flow

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pratio`:        pressure ratio
    - `mb`:        corrected mass flow
    - `piD`:       design pressure ratio
    - `mbD`:       design corrected mass flow
    - `ABCDm`:   map constants
    - `iabcd`:   map exponents
    - `effo`:     maximum efficiency
    - `Tr`:      T/Tref
    - `pr`:      p/pref

    **Outputs:**
    - `eff`:        efficiency
    - `eff_?`:      derivatives

"""
function ecmap1(pratio, m, piD, mbD, ABCDm, iabcd, effo, Tr, pr)

    eps = 1.0e-8

    C = ABCDm[3]
    D = ABCDm[4]
    mrato = ABCDm[5]

    ic = iabcd[3]
    id = iabcd[4]

    #---- corrected mass flow
    mb = m * sqrt(Tr) / pr
    mb_m = sqrt(Tr) / pr
    mb_Tr = 0.5 * mb / Tr
    mb_pr = -mb / pr

    prat = (pratio - 1.0) / (piD - 1.0)
    mrat = mb / mbD

    psgn = sign(1.0, prat / m^2 - 1.0)
    msgn = sign(1.0, mrat / mrato - 1.0)

    #c    eo  = 1.0 - D*(  abs(1.0 /mrato   - 1.0))^id
    eo = 1.0

    e1 = 1.0 - C * (psgn * (prat / mrat^2 - 1.0))^ic -
         D * (msgn * (mrat / mrato - 1.0))^id
    e1_prat = -C * (psgn * (prat / mrat^2 - 1.0))^(ic - 1) *
              psgn / mrat^2 * float(ic)
    e1_mrat = -D * (msgn * (mrat / mrato - 1.0))^(id - 1) *
              msgn / mrato * float(id) +
              C * (psgn * (prat / mrat^2 - 1.0))^(ic - 1) *
              2.0 * psgn * prat / mrat^3 * float(ic)

    eff = effo * e1 / eo
    eff_prat = effo * e1_prat / eo
    eff_mrat = effo * e1_mrat / eo

    #---- eff( mb(m,Tr,pr) , pi )
    eff_pi = eff_prat / (piD - 1.0)
    eff_mb = eff_mrat / mbD

    #---- eff( m, pi, Tr, pr )
    #cc   eff_pi = eff_pi
    eff_m = eff_mb * mb_m
    eff_Tr = eff_mb * mb_Tr
    eff_pr = eff_mb * mb_pr

    return eff, eff_pi, eff_m, eff_Tr, eff_pr
end # ecmap1

"""
    etmap(dh, mb, Nb, piD, mbD, NbD, ept0, Tmap, Tt, cpt, Rt)

Calculates turbine efficiency as a function of work and corrected mass flow

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `dh`:        enthalpy change
    - `mb`:        corrected mass flow
    - `Nb`:        corrected speed
    - `piD`:      design pressure ratio
    - `mbD`:      design corrected mass flow
    - `NbD`:      design corrected speed
    - `ept0`:      turbine polytropic efficiency estimate
    - `Tmap(.)`:   map constants

    **Outputs:**
    - `eff`:        efficiency
    - `eff_?`:      derivatives

"""
function etmap(dh, mb, Nb, piD, mbD, NbD, ept0, Tmap,
    Tt, cpt, Rt)

    eps = 1.0e-8

    pcon = Tmap[1]
    Ncon = Tmap[2]

    Trat = Tt / (Tt + dh / cpt)
    Trat_Tt = (1.0 - Trat) / (Tt + dh / cpt)
    Trat_dh = -Trat / (Tt + dh / cpt) / cpt
    Trat_cpt = Trat / (Tt + dh / cpt) * dh / cpt^2

    gex = cpt / (Rt * ept0)
    gex_cpt = 1.0 / (Rt * ept0)
    gex_Rt = -gex / Rt

    prat = Trat^gex
    prat_Trat = gex * prat / Trat
    prat_gex = prat * log(Trat)

    prat_Tt = prat_Trat * Trat_Tt
    prat_dh = prat_Trat * Trat_dh
    prat_cpt = prat_Trat * Trat_cpt + prat_gex * gex_cpt
    prat_Rt = prat_gex * gex_Rt

    Nmb = Nb * mb
    Nmb_Nb = mb
    Nmb_mb = Nb

    NmbD = NbD * mbD

    ept = ept0 * (1.0 - pcon * (1.0 - prat / piD)^2 -
                  Ncon * (1.0 - Nmb / NmbD)^2)
    ept_prat = ept0 * (pcon * (1.0 - prat / piD) * 2.0 / piD)
    ept_Nmb = ept0 * (Ncon * (1.0 - Nmb / NmbD) * 2.0 / NmbD)

    ept_dh = ept_prat * prat_dh
    ept_mb = ept_Nmb * Nmb_mb
    ept_Nb = ept_Nmb * Nmb_Nb
    ept_Tt = ept_prat * prat_Tt
    ept_cpt = ept_prat * prat_cpt
    ept_Rt = ept_prat * prat_Rt

    return ept,
    ept_dh, ept_mb, ept_Nb,
    ept_Tt, ept_cpt, ept_Rt
end # etmap

"""
    Pimap(mb, Nb, piD, mbD, NbD, Cmap)

Calculates compressor or fan pressure ratio as a function of pressure ratio and corrected mass flow

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mb`:        corrected mass flow
    - `Nb`:        corrected speed
    - `piD`:      design pressure ratio
    - `mbD`:      design corrected mass flow
    - `NbD`:      design corrected speed
    - `Cmap(.)`:  map constants

    **Outputs:**
    - `pratio`:      pressure ratio
    - `pi_?`:    derivatives

"""
function Pimap(mb, Nb, piD, mbD, NbD, Cmap)

    #----              a     b     k     mo     da    c    d     C    D
    #     data Cmapf / 3.50, 0.80, 0.03, 0.75, -0.50, 3.0, 6.0,  2.5, 15.0 /
    a = Cmap[1]
    b = Cmap[2]
    k = Cmap[3]
    #      mrato = Cmap[4]
    #      da = Cmap[5]
    #      c  = Cmap[6]
    #      d  = Cmap[7]
    #      CK = Cmap[8]
    #      DK = Cmap[9]

    #---- corrected mass flow / design corrected mass flow
    m = mb / mbD
    m_mb = 1.0 / mbD

    #---- corrected speed / design corrected speed
    N = Nb / NbD
    N_Nb = 1.0 / NbD

    Nxb = N^b
    Nxb_N = Nxb * b / N

    Nxab = N^(a * b)
    Nxab_N = Nxab * a * b / N

    g = log(1.0 - (m - Nxb) / k)
    g_N = 1.0 / (1.0 - (m - Nxb) / k) * (Nxb_N / k)
    g_m = 1.0 / (1.0 - (m - Nxb) / k) * (-1.0 / k)


    p = Nxab + 2.0 * N * k * g
    p_N = Nxab_N + 2.0 * k * g + 2.0 * N * k * g_N
    p_m = 2.0 * N * k * g_m

    pratio = 1.0 + (piD - 1.0) * p
    pi_N = (piD - 1.0) * p_N
    pi_m = (piD - 1.0) * p_m

    pi_Nb = pi_N * N_Nb
    pi_mb = pi_m * m_mb

    return pratio, pi_mb, pi_Nb
end # Pimap

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

eff, eff_pi, eff_mb = Pimap(mb, Nb, piD, mbD, NbD, Cmap)

"""
    find_cooled_hpt_efficiency(epht0, epht_fc, fc0, fc)

Calculates the cooled efficiency of an HPT based on a linear approximation.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `epht0::Float64`: uncooled polytropic efficiency
    - `epht_fc::Float64`: derivative of polytropic efficiency with respect to cooling fraction
    - `fc0::Float64`: baseline cooling fraction
    - `fc::Float64`: actual cooling fraction

    **Outputs:**
    - `epht::Float64`: cooled polytropic efficiency
"""
function find_cooled_hpt_efficiency(epht0, epht_fc, fc0, fc)
    epht = epht0 + epht_fc * (fc - fc0)
    return epht
end