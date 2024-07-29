"""
function tfoper!(gee, M0, T0, p0, a0, Tref, pref,
      Phiinl, Kinl, iBLIc,
      pid, pib, pifn, pitn,
      Gearf,
      pifD, pilcD, pihcD, pihtD, piltD,
      mbfD, mblcD, mbhcD, mbhtD, mbltD,
      NbfD, NblcD, NbhcD, NbhtD, NbltD,
      A2, A25, A5, A7,
      iTFspec,
      Ttf, ifuel, etab,
      epf0, eplc0, ephc0, epht0, eplt0,
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
      M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25)

Turbofan operation routine

    Calculation procedure follows that of Kerrebrock,
    but the usual gas property formulas are replaced
    by function calls, which can therefore implement
    more general gas models.  
    In addition, a turbine cooling model is added.

    The gas routines reside in the following source files:
     gascalc.f  Routines for various processes 
                (compressor, turbine, combustor, etc)
     gasfun.f   Routines for computing cp[T], h[t], sigma[T], R,
                called by the routines in gascalc.f

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gee`:     gravity acceleration
    - `M0`:      freestream Mach
    - `T0`:      freestream temperature  [K]
    - `p0`:      freestream pressure  [Pa]
    - `Tref`:    reference temperature for corrected mass flow and speed
    - `pref`:    reference pressure for corrected mass flow
    - `Phiinl`:  inlet ingested dissipation  Phi_inl
    - `iBLIc`:   0=core in clear flow, 1=core sees Phiinl
    - `pid`:     diffuser pressure ratio  ( = pt2/pt0)
    - `pib`:     burner   pressure ratio  ( = pt4/pt3)
    - `pifn`:    fan     nozzle pressure ratio  ( = pt7/pt6.9)
    - `pitn`:    turbine nozzle pressure ratio  ( = pt5/pt4.9)
    - `Gearf`:   fan gear ratio  ( = Nl/Nf )
    - `pifD`:    design fan pressure ratio  ( = pt21/pt2 )
    - `pilcD`:   design LPC pressure ratio  ( = pt25/pt19)
    - `pihcD`:   design HPC pressure ratio  ( = pt3 /pt25)
    - `pihtD`:   design HPT pressure ratio  ( = pt45/pt41)
    - `piltD`:   design LPT pressure ratio  ( = pt49/pt45)
    - `mbfD`:    design corrected fan mass flow ( = mf*sqrt(Tt2 /Tref)/(pt2 /pref) )
    - `mblcD`:   design corrected LPC mass flow ( = mc*sqrt(Tt19/Tref)/(pt19/pref) )
    - `mbhcD`:   design corrected HLC mass flow ( = mc*sqrt(Tt25/Tref)/(pt25/pref) )
    - `mbhtD`:   design corrected HPT mass flow ( = mt*sqrt(Tt41/Tref)/(pt41/pref) )
    - `mbltD`:   design corrected LPT mass flow ( = mt*sqrt(Tt45/Tref)/(pt45/pref) )
    - `NbfD`:    design corrected fan speed ( = Nf/sqrt(Tt2 /Tref) )
    - `NblcD`:   design corrected LPC speed ( = Nl/sqrt(Tt19/Tref) )
    - `NbhcD`:   design corrected HPC speed ( = Nh/sqrt(Tt25/Tref) )
    - `NbhtD`:   design corrected HPT speed ( = Nh/sqrt(Tt41/Tref) )
    - `NbltD`:   design corrected LPT speed ( = Nl/sqrt(Tt45/Tref) )
    - `A2`:      fan-face area [m^2]                mf = mc*BPR, mt = mc*(1+ff)
    - `A25`:     HPC-face area [m^2]
    - `A5`:      core nozzle area [m^2]
    - `A7`:      fan  nozzle area [m^2]
    - `iTFspec`:   = 1 Tt4  is specified
                   = 2 Feng is specified
    - `Tt4`:     turbine-inlet total temperature [K]
    - `Ttf`:     fuel temperature entering combustor
    - `ifuel`:   fuel index, see function gasfun (in gasfun.f)
    - `etab`:    combustor efficiency (fraction of fuel burned)
    - `epf0`:    max fan polytropic efficiency
    - `eplc0`:   LPC max polytropic efficiency
    - `ephc0`:   HPC max polytropic efficiency
    - `epht0`:   HPT max polytropic efficiency
    - `eplt0`:   LPT max polytropic efficiency
    - `pifK`:    fan efficiency FPR offset:    epolf = epf0 + epfK*(pif-pifK)
    - `epfK`:    fan efficiency pif derivative

    - `mofft`:    mass flow offtake at LPC discharge station 2.5
    - `Pofft`:    low spool power offtake
    - `Tt9`:     offtake air discharge total temperature
    - `pt9`:     offtake air discharge total pressure
    - `epsl`:    low  spool power loss fraction
    - `epsh`:    high spool power loss fraction

    - `icool`:    turbine cooling flag
                  0 = no cooling, ignore all cooling parameters below
                  1 = usual cooling, using passed-in fc
                  2 = usual cooling, but set (and return) fc from Tmetal
    - `Mtexit`:   turbine blade-row exit Mach, for setting temperature drops
    - `Tmetal`:   specified metal temperature  [K], used only if icool=2
    - `dTstrk`:   hot-streak temperature delta {K}, used only if icool=2
    - `StA`:      area-weighted Stanton number    , used only if icool=2
    - `M4a`:      effective Mach at cooling-flow outlet (start of mixing)
    - `ruc`:      cooling-flow outlet velocity ratio, u/ue
    - `ncrowx`:      dimension of epsrow array
    - `ncrow`:       number of blade rows requiring cooling
    - `epsrow(.)`:   input specified  cooling-flow bypass ratio if icool=1
                     output resulting cooling-flow bypass ratio if icool=2
    - `Tmrow(.)`:    input specified  metal temperature  [K]    if icool=2
                     output resulting metal temperature  [K]    if icool=1

      **Output:**
    ------
    - `epsrow(.)`:   see above
    - `Tmrow(.)`:    see above
    - `TSFC`:    thrust specific fuel consumption = mdot_fuel g / F   [1/s]
    - `Fsp`:     specific thrust  = F / (mdot u0) = F / ((1+BPR) mdot_core u0)
    - `hfuel`:   fuel heating value   [J / kg K]
    - `ff`:      fuel mass flow fraction  =  mdot_fuel / mdot_core
    - `Feng`:    net effective thrust  = (PK_inl+PK_out-Phi_jet)/u0  =  sum( mdot u)
    - `mcore`:   core mass flow = mdot_core  [kg/s]
    - `BPR`:     bypass ratio   = mdot_fan/mdot_core
    - `Tt?`:     total temperature
    - `ht?`:     total complete enthalpy (includes heat of formation)
    - `pt?`:     total pressure
    - `cpt?`:    specific heat at stagnation temperature  (= dh/dT)
    - `Rt?`:     gas constant  at stagnation conditions
    - `T?`:     static temperature
    - `u?`:      velocity
    - `etaf`:    fan          overall efficiency
    - `etac`:    compressor   overall efficiency
    - `etatf`:   fan-turbine  overall efficiency
    - `etatc`:   comp-turbine overall efficiency
    - `Lconv`:   T if convergence was successful, F otherwise

    The "?" symbol denotes the station index:
      0   freestream
      18  fan face outside of casing BLs
      19  fan face over LPC portion
      2   fan face over fan portion
      21  fan exit, precooler inlet
      19c precooler outlet, LPC inlet
      25  LPC exit, intercooler inlet 
      25c intercooler exit, HPC inlet
      3   compressor exit
      4   combustor exit before cooling air addition
      41  turbine  inlet after  cooling air addition
      45  HPT exit, LPT inlet
      49  LPT exit, regenerative cooler inlet
      49c regenerative cooler outlet
      5   core nozzle
      6   core flow downstream
      7   fan nozzle
      8   fan flow downstream
"""
function ductedfanoper!(gee, M0, T0, p0, a0, Tref, pref,
      Phiinl, Kinl, iBLIc,
      pid, pifn, 
      pifD, 
      mbfD, 
      NbfD,
      A2, A7,
      epf0,
      pifK,
      Feng,
      M2, pif, mbf)

      sol = nlsolve(residual, guess, ftol = 1e-7) 


        return TSFC, Fsp, hfuel, ff,
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
        Lconv

end




end # tfoper

function res_df(x)
    # ===============================================================
    #---- freestream static quantities
    s0, dsdt, h0, dhdt, cp0, R0 = gassum(alpha, nair, T0)
    gam0 = cp0 / (cp0 - R0)
    #      a0 = sqrt(gam0*R0*T0)
    u0 = M0 * a0
    #
    # ===============================================================
    #---- freestream total quantities
    hspec = h0 + 0.5 * u0^2
    Tguess = T0 * (1.0 + 0.5 * (gam0 - 1.0) * M0^2)
    Tt0 = gas_tset(alpha, nair, hspec, Tguess)
    st0, dsdt, ht0, dhdt, cpt0, Rt0 = gassum(alpha, nair, Tt0)
    pt0 = p0 * exp((st0 - s0) / Rt0)
    at0 = sqrt(Tt0 * Rt0 * cpt0 / (cpt0 - Rt0))
    if (u0 == 0.0)
        #----- static case... no BL defect
        sbfan = 0.0

    else
        #----- account for inlet BLI defect via mass-averaged entropy
        a2sq = at0^2 / (1.0 + 0.5 * (gam0 - 1.0) * Mi^2)

        if (iBLIc == 0)
            #------ BL mixes with fan flow only
            #c      mmix    = mf*sqrt(Tref/Tt2) * pt2   /pref
            #c      mmix_mf =    sqrt(Tref/Tt2) * pt2   /pref
            #c      mmix_Mi = mf*sqrt(Tref/Tt2) * pt2_Mi/pref

            mmix = mf * sqrt(Tref / Tt0) * pt0 / pref
            
            sbfan = Kinl * gam0 / (mmix * a2sq)

        else
            
            mmix = mf * sqrt(Tref / Tt0) * pt0 / pref +
                    ml * sqrt(Tref / Tt0) * pt0 / pref
            
            sbfan = Kinl * gam0 / (mmix * a2sq)
        end
    end

    Tt2 = Tt18
    ht2 = ht18
    st2 = st18
    cpt2 = cpt18
    Rt2 = Rt18
    pt2 = pt18 * exp(-sbfan)

    p2, T2, h2, s2, cp2, R2 = gas_mach(alpha, nair,
            pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, Mi, 1.0)

    u2 = sqrt(2.0 * (ht2 - h2))

    rho2 = p2 / (R2 * T2)

    epf, epf_pf, epf_mf = ecmap(pf, mf, pifD, mbfD, Cmapf, epf0, pifK, epfK)

    if (epf < epfmin)
            epf = epfmin
            epf_pf = 0.0
            epf_mf = 0.0
    end
    if (pf < 1.0)
            epf_pf = (-1.0 / epf^2) * epf_pf
            epf_mf = (-1.0 / epf^2) * epf_mf
            epf = 1.0 / epf
    end

    pt21, Tt21, ht21, st21, cpt21, Rt21 = gas_prat(alpha, nair,
            pt2, Tt2, ht2, st2, cpt2, Rt2, pf, epf)

    #---- fan duct nozzle total quantities
    pt7 = pt21 * pifn
    Tt7 = Tt21
    ht7 = ht21
    st7 = st21
    cpt7 = cpt21
    Rt7 = Rt21

 # ===============================================================
    #---- fan nozzle flow 7-8, use alpha mass fraction (air)
    pfn = p0 / pt7

    p7, T7, h7, s7, cp7, R7 = gas_prat(alpha, nair,
    pt7, Tt7, ht7, st7, cpt7, Rt7, pfn, 1.0)

    u7 = sqrt(2.0 * max(ht7 - h7, 0.0))
    M7 = u7 / sqrt(T7 * cp7 * R7 / (cp7 - R7))

    if (M7 > 1.0)
        #----- fan nozzle is choked... specify M7 = 1 instead
        M7 = 1.0

        p7, T7, h7, s7, cp7, R7 = gas_mach(alpha, nair,
                    pt7, Tt7, ht7, st7, cpt7, Rt7, 0.0, M7, 1.0)

    end

    if (ht7 > h7)
            u7 = sqrt(2.0 * (ht7 - h7))
            
    else
            u7 = 0.0

    end

    rho7 = p7 / (R7 * T7)

# ===============================================================
    #---- #Set up residuals
    res = zeros(3)

    mfan = mf * sqrt(Tref / Tt2) * pt2 / pref
    #Fan nozzle mass flow, choked or unchoked
    res[1] = mfan - rho7 * A7 * u7

    #Inlet Mach area constraint
    mfA = mf * sqrt(Tref / Tt2) * pt2 / pref / (rho2 * u2)
    res[2] = mfA - A2

    if (iPspec == 0) #Specified power constraint
        P = mfan * (ht21 - ht2)
        res[3] = P - Pspec

    else #Specified thrust constraint
        epi = 1.0
        p8, T8, h8, s8, cp8, R8 = gas_prat(alpha, nair,
                        pt7, Tt7, ht7, st7, cpt7, Rt7, pfn, epi)

        if (ht7 > h8)
            u8 = sqrt(2.0 * (ht7 - h8))
        else
            u8 = 0.0
        end

        #----- overall thrust
        if (u0 == 0.0)
            Finl = 0.0
        else
            Finl = Phiinl / u0
        end
        F = mfan * (u8 - u0)  + Finl

        #residual
        res[3] = F - Fspec
    end
    return res
end
