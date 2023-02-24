include("gasfun.jl")
include("gascalc.jl")

function tfan(gee, M0, T0, p0, Mfan, Afan,
    BPR, pif, pic, pid, pib,
    Tt4, Ttf, ifuel,
    epolf, epolc, epolt,
    icool,
    Mtexit, Tmetal, dTstrk, Stc,
    M4a, ruc)
    #----------------------------------------------------------------
    #     Turbofan performance routine.
    #
    #     Calculation procedure follows that of Kerrebrock,
    #     but the usual gas property formulas are replaced
    #     by function calls, which can therefore implement
    #     more general gas models.  
    #     In addition, a turbine cooling model is added.
    #
    #     The gas routines reside in the following source files:
    #      gascalc.f  Routines for various processes 
    #                 (compressor, turbine, combustor, etc)
    #      gasfun.f   Routines for computing cp[T], h[t], sigma[T], R,
    #                 called by the routines in gascalc.f
    #
    #     Input
    #     -----
    #     gee    gravity acceleration
    #     M0     freestream Mach
    #     T0     freestream temperature  [K]
    #     p0     freestream pressure  [Pa]
    #     Mfan   fan-face Mach number, for computing engine mass flow
    #     Afan   fan-face area [m^2] , for computing engine mass flow
    #     BPR    bypass ratio  = mdot_fan/mdot_core
    #     pif    fan      pressure ratio  ( = pt7/pt2)
    #     pic    overall  pressure ratio  ( = pt3/pt2)
    #     pid    diffuser pressure ratio  ( = pt2/pt0)
    #     pib    burner   pressure ratio  ( = pt4/pt3)
    #     Tt4    turbine-inlet total temperature [K]
    #     Ttf    fuel temperature entering combustor
    #     ifuel  fuel index, see function gasfun (in gasfun.f)
    #     epolf  fan        polytropic efficiency
    #     epolc  compressor polytropic efficiency
    #     epolt  turbine    polytropic efficiency
    #
    #     icool   turbine cooling flag
    #              0 = no cooling, ignore all cooling parameters below
    #              1 = usual cooling, using passed-in BPRc
    #              2 = usual cooling, but set (and return) BPRc from Tmetal
    #     Mtexit  turbine blade-row exit Mach, for setting temperature drops
    #     Tmetal  specified metal temperature  [K], used only if icool=2
    #     dTstrk  hot-streak temperature delta {K}, used only if icool=2
    #     Stc     area-weighted Stanton number    , used only if icool=2
    #     M4a     effective Mach at cooling-flow outlet (start of mixing)
    #     ruc     cooling-flow outlet velocity ratio, u/ue
    #     BPRc    cooling-flow bypass ratio, mdot_cool/mdot_core, input if icool=1
    #
    #     Output
    #     ------
    #     BPRc    cooling-flow bypass ratio, mdot_cool/mdot_core, output if icool=2
    #     TSFC   thrust specific fuel consumption = mdot_fuel g / F   [1/s]
    #     Fsp    specific thrust  = F / (mdot a) = F / ((1+BPR) mdot_core a)
    #     hfuel  fuel heating value   [J / kg K]
    #     ff     fuel mass flow fraction  =  mdot_fuel / mdot_core
    #     mdot   core mass flow = mdot_core  [kg/s]
    #     Tt?    total temperature
    #     ht?    total complete enthalpy (includes heat of formation)
    #     pt?    total pressure
    #     cpt?   specific heat at stagnation temperature  (= dh/dT)
    #     Rt?    gas constant  at stagnation conditions
    #     T?     static temperature
    #     u?     velocity
    #     etaf   fan        overall efficiency
    #     etac   compressor overall efficiency
    #     etat   turbine    overall efficiency
    #
    #     The "?" symbol denotes the station index:
    #       0  freestream
    #       2  fan face
    #       3  compressor exit
    #       4  turbine inlet
    #       5  turbine exit
    #       6  turbine flow downstream
    #       7  fan exit
    #       8  fan flow downstream
    #----------------------------------------------------------------
    n = 6
    # c---- air fractions  N2      O2      CO2    H2O      Ar       fuel
    alpha = [0.781, 0.209, 0.0004, 0.0, 0.0096, 0.0]
    # c---- fuel fractions
    beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    # c

    #---- fan and turbine nozzle flow total-pressure ratios 
    pifnoz = 1.0
    pitnoz = 1.0

    #---- number of air constitutents (all but fuel)
    nair = n - 1

    # =========================
    #---- set combustion-change mass fractions gamma[i] for specified fuel
    gamma = gasfuel(ifuel, n)
    #
    # =========================
    #---- freestream static quantities
    s0, dsdt, h0, dhdt, cp0, R0 = gassum(alpha, nair, T0,)
    gam0 = cp0 / (cp0 - R0)
    a0 = sqrt(gam0 * R0 * T0)
    u0 = M0 * a0
    #
    # =========================
    #---- freestream total quantities
    hspec = h0 + 0.5 * u0^2
    Tguess = T0 * (1.0 + 0.5 * (gam0 - 1.0) * M0^2)
    Tt0 = gas_tset(alpha, nair, hspec, Tguess)
    # 
    st0, dsdt, ht0, dhdt, cpt0, Rt0 = gassum(alpha, nair, Tt0)
    pt0 = p0 * exp((st0 - s0) / Rt0)
    #
    # =========================
    #---- diffuser flow 0..2
    Tt2 = Tt0
    st2 = st0
    ht2 = ht0
    cpt2 = cpt0
    Rt2 = Rt0
    pt2 = pt0 * pid
    #
    # =========================
    #---- fan flow 2..7
    pt7, Tt7, ht7, st7, cpt7, Rt7 = gas_prat(alpha, nair, pt2, Tt2, ht2, st2, cpt2, Rt2, pif, epolf)

    #---- calculate fan efficiency (informative only -- not needed here)
    etaf = 0.0
    pt7i, Tt7i, ht7i, st7i, cpt7i, Rt7i = gas_prat(alpha, nair, pt2, Tt2, ht2, st2, cpt2, Rt2, pif, 1.0,
    )
    etaf = (ht7i - ht2) / (ht7 - ht2)
    #      write(*,*) 'eta_f', etaf, epolf
    #
    # =========================
    #---- compressor flow 2..3
    pt3, Tt3, ht3, st3, cpt3, Rt3 = gas_prat(alpha, nair, pt2, Tt2, ht2, st2, cpt2, Rt2, pic, epolc)
    #c    write(*,*) 'epolc Tt3', epolc, Tt3

    #---- calculate compressor efficiency (informative only -- not needed here)
    etac = 0.0
    pt3i, Tt3i, ht3i, st3i, cpt3i, Rt3i = gas_prat(alpha, nair, pt2, Tt2, ht2, st2, cpt2, Rt2, pic, 1.0)
    etac = (ht3i - ht2) / (ht3 - ht2)
    #c      write(*,*) 'eta_c', etac, epolc

    # =========================
    #---- combustor flow 3..4   (ffb = mdot_fuel/mdot_burner)
    ffb, lambda = gas_burn(alpha, beta, gamma, n, ifuel, Tt3, Ttf, Tt4)
    st4, dsdt, ht4, dhdt, cpt4, Rt4 = gassum(lambda, nair, Tt4)
    pt4 = pt3 * pib
    gam4 = cpt4 / (cpt4 - Rt4)

    #---- set ff = mdot_fuel/mdot_core = ffb * mdot_burner/mdot_core
    BPRc = 0.0
    ff = ffb * (1.0 - BPRc)

    #      write(*,*) 'ff    Tt4', ff, Tt4
    #      write(*,12) 'pt0 pt2 pt3 pt4', p0, pt0, pt2, pt3, pt4


    #     write(*,'(1x,8f9.5)') alpha,
    #     alpha[1]+alpha[2]+alpha[3]+alpha[4]+alpha[5]+alpha[6]
    #     write(*,'(1x,8f9.5)') gamma,
    #     gamma[1]+gamma[2]+gamma[3]+gamma[4]+gamma[5]+gamma[6]
    #     write(*,'(1x,8f9.5)') lambda,
    #     lambda[1]+lambda[2]+lambda[3]+lambda[4]+lambda[5]+lambda[6]

    # =========================
    if (icool == 2)
        #----- set cooling mass flow
        efilm = 0.7
        tfilm = 0.4
        #
        gmi4 = Rt4 / (cpt4 - Rt4)
        Trrat = 1.0 / (1.0 + 0.5 * gmi4 * Mtexit^2)

        #----- calculate cooling mass flow ratio BPRc
        throw, epsrow, BPRc, nrow, nrowx = mcool(Tmetal, Tt3, Tt4, dTstrk, Trrat,
            efilm, tfilm, Stc)

        # else
        #----- BPRc is assumed to be passed in
    end

    lambdap = zeros(nair)
    #----------------------------------------------------------------
    if (icool == 0)
        #----- no cooling air present
        pt41 = pt4
        Tt41 = Tt4
        ht41 = ht4
        st41 = st4
        cpt41 = cpt4
        Rt41 = Rt4
        for i = 1:nair
            lambdap[i] = lambda[i]
        end
        #----------------------------------------------------------------
    else
        pt4a = pt4
        Tt4a = Tt4
        ht4a = ht4
        st4a = st4
        cpt4a = cpt4
        Rt4a = Rt4

        #----- speed at start-of-mixing station 4a
        p4a, T4a, h4a, s4a, cp4a, r4a = gas_mach(lambda, nair, pt4a, Tt4a, ht4a, st4a, cpt4a, rt4a, 0.0, M4a, 1.0)
        u4sq = max(2.0 * (ht4a - h4a), 0.0)
        u4a = sqrt[u4sq]

        uc = ruc * u4a
        #---------------------------------------------------------
        #----- IGV exit mixing
        frac4 = (1.0 - BPRc + ff) / (1.0 + ff)
        fracm = BPRc / (1.0 + ff)

        #----- mixed constituent fraction vector from mass equation
        for i = 1:nair
            lambdap[i] = frac4 * lambda[i] + fracm * alpha[i]
        end

        #----- mixed total enthalpy from enthalpy equation
        ht41 = frac4 * ht4 + fracm * ht3

        #----- total temperature from total enthalpy
        Tguess = frac4 * Tt4 + fracm * Tt3
        Tt41 = gas_tset(lambdap, nair, ht41, Tguess)

        #----- all total quantities, from total temperature
        st41, dsdt, ht41, dhdt, cpt41, Rt41 = gassum(lambdap, nair, Tt41)
        pt41 = pt4

        #----- mixed velocity from momentum equation
        p41 = p4a
        u41 = frac4 * u4a + fracm * uc
        #
        #----- static temperature from static enthalpy
        h41 = ht41 - 0.5 * u41^2
        Tguess = t4a + (h41 - h4a) / cp4a
        T41 = gas_tset(lambdap, nair, h41, Tguess)
        #
        #----- all static quantities, from static temperature
        s41, dsdt, h41, dhdt, cp41, R41 = gassum(lambdap, nair, T41)

        #----- all stagnation quantities, from total-static enthalpy difference
        delh = ht41 - h41
        epi = 1.0
        pt41, Tt41, ht41, st41, cpt41, Rt41 = gas_delh(lambdap, nair, p41, T41, h41, s41, cp41, R41, delh, epi)
    end
    #----------------------------------------------------------------

    # =========================
    #---- turbine flow 4..5
    delh = -(ht3 - ht2 + BPR * (ht7 - ht2)) / (1.0 + ff)
    epi = 1.0 / epolt
    pt5, Tt5, ht5, st5, cpt5, Rt5 = gas_delh(lambdap, nair,
        pt41, Tt41, ht41, st41, cpt41, Rt41, delh, epi)

    #c    write(*,*) 'delh     ', delh
    #c    write(*,*) 'epolt Tt5', epolt, Tt5

    #---- calculate turbine efficiency (informative only -- not needed here)
    etat = 0.0
    pit = pt5 / pt41
    pt5i, Tt5i, ht5i, st5i, cpt5i, Rt5i = gas_prat(lambdap, nair, pt41, Tt41, ht41, st41, cpt41, Rt41, pit, 1.0)
    etat = (ht5 - ht41) / (ht5i - ht41)
    #c      write(*,*) 'eta_t', etat, epolt

    #
    # =========================
    #---- fan nozzle flow 7..8, use alpha mass fraction (air)
    pt8 = pt7 * pifnoz
    ht8 = ht7
    Tt8 = Tt7
    st8 = st7
    cpt8 = cpt7
    Rt8 = Rt7
    pratfn = p0 / pt8
    p8, T8, h8, s8, cp8, R8 = gas_prat(alpha, nair, pt8, Tt8, ht8, st8, cpt8, Rt8, pratfn, 1.0)
    u8 = sqrt(2.0 * (ht8 - h8))

    #
    #---- fan specific thrust = F_fan / mdot_core a
    F8sp = BPR * (u8 - u0) / a0
    #
    # =========================
    #---- turbine nozzle flow 5..6, use lambda mass fraction (combustion products)
    pt6 = pt5 * pitnoz
    ht6 = ht5
    Tt6 = Tt5
    st6 = st5
    cpt6 = cpt5
    Rt6 = Rt5
    prattn = p0 / pt6
    p6, T6, h6, s6, cp6, R6 = gas_prat(lambdap, nair, pt6, Tt6, ht6, st6, cpt6, Rt6, prattn, 1.0)
    h6 = min(h6, ht6)
    u6 = sqrt(2.0 * (ht6 - h6))

    #c    write(*,*) 'Pt6   u6 ', pt6, u6
    #
    #---- core specific thrust = F_core / mdot_core a
    F6sp = ((1.0 + ff) * u6 - u0) / a0
    #
    # =========================
    #---- overall specific thrust  Fsp  =  F / (mdot_fan + mdot_core) a
    Fsp = (F6sp + F8sp) / (1.0 + BPR)

    #---- overall specific impulse and fuel consumption
    Isp = Fsp * a0 * (1.0 + BPR) / (gee * ff)
    TSFC = 1.0 / Isp
    #
    cpa = 0.5 * (cpt3 + cpt41)
    hfuel = ((Tt41 - Tt3) + ff * (Tt41 - Ttf)) * cpa / ff
    #
    # =========================
    #---- actual core mass flow for specified Afan = A2
    if (Afan > 0.0)
        M2 = Mfan
        p2, t2, h2, s2, cp2, r2 = gas_mach(alpha, nair, pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, M2, 1.0)
        u2 = sqrt(2.0 * (ht2 - h2))
        rho2 = p2 / (r2 * t2)
        mdot = rho2 * u2 * Afan / (1.0 + BPR)
    end

    return BPRc,
    TSFC, Fsp, hfuel, ff, mdot,
    Tt0, ht0, pt0, cpt0, Rt0,
    Tt2, ht2, pt2, cpt2, Rt2,
    Tt3, ht3, pt3, cpt3, Rt3,
    ht4, pt4, cpt4, Rt4,
    Tt41, ht41, pt41, cpt41, Rt41,
    Tt5, ht5, pt5, cpt5, Rt5,
    Tt7, ht7, pt7, cpt7, Rt7,
    u0,
    T6, u6,
    T8, u8,
    etaf, etac, etat
end # tfan



function mcool(Tmetal, Tt3, Tt4, dTstreak, Trrat,
    efilm, tfilm, Stc)

    # =========================
    #     Calculates cooling mass flow requirement.
    #     All temperatures are in Kelvin.
    #
    #     Input
    #     -----
    #     Tmetal    design metal temperature
    #     Tt3       cooling flow temperature
    #     Tt4       hot gas temperature from burner
    #     dTstreak  added temperature seen by first IGV
    #     Trrat     static temperature ratio across each blade row, T4.1 / T4
    #     
    #     Output
    #     ------
    #     throw(.)  cooling effectiveness for each blade row 1,2,3...nrow
    #     epsrow(.) cooling mass flow ratio for each blade row, m_c_row/m_core
    #     BPRc      total cooling mass flow for all rows with cooling
    #     nrow      number of blade rows which need cooling
    #     nrowx     dimension of arrays
    #
    # =========================
    #
    for irow = 1:nrowx
        if (irow == 1)
            Tg = Tt4 + dTstreak
        else
            Tg = Tt4 * Trrat^(irow - 1)
        end
        #
        throw[irow] = (Tg - Tmetal) / (Tg - Tt3)
        eps0 = Stc * (throw[irow] * (1.0 - efilm * tfilm) - tfilm * (1.0 - efilm)) / (efilm * (1.0 - tfilm))
        if (eps0 < 0.0)
            BPRc = 0.0
            for irow = 1:nrow
                BPRc = BPRc + epsrow[irow]
            end

            return throw, epsrow, BPRc, nrow, nrowx
        end
        nrow = irow
        epsrow[irow] = eps0 / (1.0 + eps0)
    end


end # mcool