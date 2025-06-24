"""
    function tfoper!(gee, M0, T0, p0, a0, Tref, pref,
      Phiinl, Kinl, eng_has_BLI_cores,
      pid, pib, pifn, pitn,
      Gearf,
      pifD, pilcD, pihcD, pihtD, piltD,
      mbfD, mblcD, mbhcD, mbhtD, mbltD,
      NbfD, NblcD, NbhcD, NbhtD, NbltD,
      A2, A25, A5, A7,
      opt_calc_call,
      Ttf, ifuel, etab,
      epf0, eplc0, ephc0, epht0, eplt0,
      mofft, Pofft,
      Tt9, pt9,
      epsl, epsh,
      opt_cooling,
      Mtexit, dTstrk, StA, efilm, tfilm,
      M4a, ruc,
      ncrowx, ncrow,
      epsrow, Tmrow,
      Feng,
      M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25)

Turbofan operation routine. Calculation procedure follows that of Kerrebrock, but
the usual gas property formulas are replaced by function calls, which can therefore
implement more general gas models. In addition, a turbine cooling model is added.

The gas routines are described in [Gas Calculations](@ref)

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gee`:     gravity acceleration
    - `M0`:      freestream Mach
    - `T0`:      freestream temperature  [K]
    - `p0`:      freestream pressure  [Pa]
    - `Tref`:    reference temperature for corrected mass flow and speed
    - `pref`:    reference pressure for corrected mass flow
    - `Phiinl`:  inlet ingested dissipation  `Phi_inl`
    - `eng_has_BLI_cores`:
      - `false`: core in clear flow
      - `true`: core sees `Phiinl`
    - `pid`:     diffuser pressure ratio  ( = `pt2/pt0`)
    - `pib`:     burner   pressure ratio  ( = `pt4/pt3`)
    - `pifn`:    fan     nozzle pressure ratio  ( = `pt7/pt6.9`)
    - `pitn`:    turbine nozzle pressure ratio  ( = `pt5/pt4.9`)
    - `Gearf`:   fan gear ratio  ( = Nl/Nf )
    - `pifD`:    design fan pressure ratio  ( = `pt21/pt2`)
    - `pilcD`:   design LPC pressure ratio  ( = `pt25/pt19`)
    - `pihcD`:   design HPC pressure ratio  ( = `pt3/pt25`)
    - `pihtD`:   design HPT pressure ratio  ( = `pt45/pt41`)
    - `piltD`:   design LPT pressure ratio  ( = `pt49/pt45`)
    - `mbfD`:    design corrected fan mass flow ( = `mf*sqrt(Tt2/Tref)/(pt2/pref)` )
      where `mf = mc*BPR`
    - `mblcD`:   design corrected LPC mass flow ( = `mc*sqrt(Tt19/Tref)/(pt19/pref)` )
    - `mbhcD`:   design corrected HLC mass flow ( = `mc*sqrt(Tt25/Tref)/(pt25/pref)` )
    - `mbhtD`:   design corrected HPT mass flow ( = `mt*sqrt(Tt41/Tref)/(pt41/pref)` )
      where `mt = mc*(1+ff)`
    - `mbltD`:   design corrected LPT mass flow ( = `mt*sqrt(Tt45/Tref)/(pt45/pref)` )
    - `NbfD`:    design corrected fan speed ( = `Nf/sqrt(Tt2/Tref)` )
    - `NblcD`:   design corrected LPC speed ( = `Nl/sqrt(Tt19/Tref)` )
    - `NbhcD`:   design corrected HPC speed ( = `Nh/sqrt(Tt25/Tref)` )
    - `NbhtD`:   design corrected HPT speed ( = `Nh/sqrt(Tt41/Tref)` )
    - `NbltD`:   design corrected LPT speed ( = `Nl/sqrt(Tt45/Tref)` )
    - `A2`:      fan-face area [m^2]
    - `A25`:     HPC-face area [m^2]
    - `A5`:      core nozzle area [m^2]
    - `A7`:      fan  nozzle area [m^2]
    - `opt_calc_call`:
      - `"oper_fixedTt4"`: `Tt4` is specified
      - `"oper_fixedFe"`: `Feng` is specified
    - `Tt4`:     turbine-inlet total temperature [K]
    - `Ttf`:     fuel temperature entering combustor
    - `ifuel`:   fuel index, see function [`gasfun`](@ref)
    - `hvap`:    fuel enthalpy of vaporization (J/kg)
    - `etab`:    combustor efficiency (fraction of fuel burned)
    - `epf0`:    max fan polytropic efficiency
    - `eplc0`:   LPC max polytropic efficiency
    - `ephc0`:   HPC max polytropic efficiency
    - `epht0`:   HPT max polytropic efficiency
    - `eplt0`:   LPT max polytropic efficiency

    - `mofft`:    mass flow offtake at LPC discharge station 2.5
    - `Pofft`:    low spool power offtake
    - `Tt9`:     offtake air discharge total temperature
    - `pt9`:     offtake air discharge total pressure
    - `epsl`:    low  spool power loss fraction
    - `epsh`:    high spool power loss fraction

    - `opt_cooling`:   turbine cooling flag
      - `"none"`: no cooling, ignore all cooling parameters below
      - `"fixed_coolingflowratio"`: usual cooling, using passed-in `fcool`
      - `"fixed_Tmetal"`: usual cooling, but set (and return) `fcool` from `Tmetal`
    - `Mtexit`:   turbine blade-row exit Mach, for setting temperature drops
    - `Tmetal`:   specified metal temperature  [K], used only if `opt_cooling="fixed_Tmetal"`
    - `dTstrk`:   hot-streak temperature delta {K}, used only if `opt_cooling="fixed_Tmetal"`
    - `StA`:      area-weighted Stanton number    , used only if `opt_cooling="fixed_Tmetal"`
    - `M4a`:      effective Mach at cooling-flow outlet (start of mixing)
    - `ruc`:      cooling-flow outlet velocity ratio, `u/ue`
    - `ncrowx`:      dimension of `epsrow` array
    - `ncrow`:       number of blade rows requiring cooling
    - `epsrow(.)`: input specified  cooling-flow bypass ratio if
      `opt_cooling="fixed_coolingflowratio"`; output resulting cooling-flow bypass ratio
      if `opt_cooling="fixed_Tmetal"`.
    - `Tmrow(.)`: input specified metal temperature [K] if `opt_cooling="fixed_Tmetal"`;
      output resulting metal temperature [K] if `opt_cooling="fixed_coolingflowratio"`

    **Output:**
    - `epsrow(.)`:   see above
    - `Tmrow(.)`:    see above
    - `TSFC`:    thrust specific fuel consumption = `mdot_fuel g / F`   [1/s]
    - `Fsp`:     specific thrust  = `F / (mdot u0) = F / ((1+BPR) mdot_core u0)`
    - `hfuel`:   fuel heating value   [J / kg K]
    - `ff`:      fuel mass flow fraction  =  `mdot_fuel / mdot_core`
    - `Feng`:    net effective thrust  = `(PK_inl+PK_out-Phi_jet)/u0  =  sum(mdot u)`
    - `mcore`:   core mass flow = `mdot_core`  [kg/s]
    - `BPR`:     bypass ratio   = `mdot_fan/mdot_core`
    - `Tt?`:     total temperature
    - `ht?`:     total complete enthalpy (includes heat of formation)
    - `pt?`:     total pressure
    - `cpt?`:    specific heat at stagnation temperature  (= `dh/dT`)
    - `Rt?`:     gas constant  at stagnation conditions
    - `T?`:      static temperature
    - `u?`:      velocity
    - `etaf`:    fan          overall efficiency
    - `etac`:    compressor   overall efficiency
    - `etatf`:   fan-turbine  overall efficiency
    - `etatc`:   comp-turbine overall efficiency
    - `Lconv`:   `true` if convergence was successful, `false` otherwise

    The "?" symbol denotes the station index:
    - 0: freestream
    - 18: fan face outside of casing BLs
    - 19: fan face over LPC portion
    - 2: fan face over fan portion
    - 21: fan exit, precooler inlet
    - 19c: precooler outlet, LPC inlet
    - 25: LPC exit, intercooler inlet
    - 25c: intercooler exit, HPC inlet
    - 3: compressor exit
    - 4: combustor exit before cooling air addition
    - 41: turbine inlet after cooling air addition
    - 45: HPT exit, LPT inlet
    - 49: LPT exit, regenerative cooler inlet
    - 49c: regenerative cooler outlet
    - 5: core nozzle
    - 6: core flow downstream
    - 7: fan nozzle
    - 8: fan flow downstream
"""
function tfoper!(gee, M0, T0, p0, a0, Tref, pref,
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
      epsrow, Tmrow,
      Feng,
      M2, pif, pilc, pihc, mbf, mblc, mbhc, Tt4, pt5, mcore, M25, 
      Δh_PreC, Δh_InterC, Δh_Regen, Δh_TurbC,
      Δp_PreC, Δp_InterC, Δp_Regen)

      #---- ncrowy must be at least as big as ncrowx defined in index.inc
      ncrowy = 8

      # Determine whether we are in AD...
      prod = gee * M0 * T0 * p0 * a0 * Tref * pref * Phiinl * Kinl * pid * pib * pifn * pitn * Gearf
      prod *= pifD * pilcD * pihcD * pihtD * piltD
      prod *= mbfD * mblcD * mbhcD * mbhtD * mbltD
      prod *= NbfD * NblcD * NbhcD * NbhtD * NbltD
      prod *= A2 * A25 * A5 * A7
      prod *= Ttf * ifuel * etab
      prod *= epf0 * eplc0 * ephc0 * epht0 * eplt0
      prod *= mofft * Pofft * Tt9 * pt9 * epsl * epsh
      prod *= Mtexit * dTstrk * StA * efilm * tfilm
      prod *= M4a * ruc
      prod *= epsrow[1] * M2 * pif * pilc * pihc * mbf * mblc * mbhc * Tt4 * pt5 * mcore * M25
      T = typeof(prod)


      #---- Newton system arrays
      res = zeros(T, 9, 1)
      a = zeros(T, 9, 9)
      rrel = zeros(T, 9)
      rsav = zeros(T, 9)
      asav = zeros(T, 9, 10)

      res_dlls = zeros(T, 9)
      a_dlls = zeros(T, 9, 9)

      Random.seed!(1234) #Seed for RNG in relaxation

      #---- number of gas constituents
      n = 6

      #---- mass fractions
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

      #---- minimum allowable fan efficiency
      epfmin = 0.60

      # from "airfrac.inc"
      # air fractions  
      #        N2      O2      CO2    H2O      Ar       fuel
      alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127, 0.0]
      beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

      #---- convergence tolerance
      #      data toler  1.0e-7 
      toler = 1.0e-9

      itmax = 50

      #---- max fan-face Mach number, above which it will be artificially limited
      Mimax = 0.98

      Lprint = false

      if compare_strings(opt_calc_call, "oper_fixedTt4")
            Tt4spec = Tt4
      elseif compare_strings(opt_calc_call, "oper_fixedFe")
            Fspec = Feng
      end

      #---- number of air constitutents (all but fuel)
      nair = n - 1

      # ===============================================================
      #---- set combustion-change mass fractions gamma[i] for specified fuel
      gamma = gasfuel(ifuel, n)

      # Convert type for ForwardDiff
      if (typeof(etab) <: ForwardDiff.Dual)
            gamma = convert(Array{typeof(etab)}, gamma)
      end
      #---- apply combustor efficiency
      for i = 1:nair
            gamma[i] = etab * gamma[i]
      end
      gamma[n] = 1.0 - etab

      #Starting guesses for compressor map non-linear solvers
      Nfg = 0.5
      Rfg = 2.0
      Nlcg = 0.5
      Rlcg = 2.0
      Nhcg = 0.5
      Rhcg = 2.0
      #
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

      #
      # ===============================================================
      #---- Offtake plume flow 9
      Trat = (p0 / pt9)^(Rt0 / cpt0)
      if (Trat < 1.0)
            u9 = sqrt(2.0 * cpt0 * Tt9 * (1.0 - Trat))
            rho9 = p0 / (Rt0 * Tt0 * Trat)
      else
            u9 = 0.0
            rho9 = p0 / (Rt0 * Tt0)
      end

      # ===============================================================
      #---- diffuser flow 0-2
      Tt18 = Tt0
      st18 = st0
      ht18 = ht0
      cpt18 = cpt0
      Rt18 = Rt0
      pt18 = pt0 * pid

      #---- initial guesses for primary Newton variables
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

      #---- total cooling mass flow ratio
      fc = 0.0

      if compare_strings(opt_cooling, "fixed_coolingflowratio")
            if (mcore == 0.0)
                  fo = 0.0
            else
                  fo = mofft / mcore
            end
            for icrow = 1:ncrow
                  fc = fc + (1.0 - fo) * epsrow[icrow]
            end
      end

      fc_pl = 0.0
      fc_ph = 0.0
      fc_ml = 0.0
      fc_mh = 0.0
      fc_Tb = 0.0

      for iter = 1:itmax
      
            if (iter == -1)
                  eps1 = 2.0e-7
                  eps = eps1

                  pf = pf + eps
                  j = 1

            end

            epf = epf0
            eplc = eplc0
            ephc = ephc0
            epht = epht0
            eplt = eplt0

            # println("epf, eplc, ephc, epht, eplt")
            # println(epf, eplc, ephc, epht, eplt)
            # exit()

            # HSC: SEEMS TO BE FINE

            # ===============================================================
            #---- set fan inlet total pressure pt21 corrected for BLI
            #-    (Tt2,pt2 approximated with Tt0,pt0 here to avoid circular definition)
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

                  if eng_has_BLI_cores
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

                  else #clean flow
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
            u2_Mi = (-1.0 / u2) * h2_Mi

            rho2 = p2 / (R2 * T2)
            rho2_p2 = 1.0 / (R2 * T2)
            rho2_T2 = -rho2 / T2

            rho2_Mi = rho2_T2 * T2_Mi +
                      rho2_p2 * (p2_pt2 * pt2_Mi + p2_Mi)
            rho2_mf = rho2_p2 * p2_pt2 * pt2_mf
            rho2_ml = rho2_p2 * p2_pt2 * pt2_ml


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
            u19_Mi = (-1.0 / u19) * h19_Mi

            rho19 = p19 / (R19 * T19)
            rho19_p19 = 1.0 / (R19 * T19)
            rho19_T19 = -rho19 / T19

            rho19_Mi = rho19_T19 * T19_Mi +
                       rho19_p19 * (p19_pt19 * pt19_Mi + p19_Mi)
            rho19_mf = rho19_p19 * p19_pt19 * pt19_mf
            rho19_ml = rho19_p19 * p19_pt19 * pt19_ml

            # ===============================================================
            #---- fan flow 2-7
            Nf, epf, Nf_pf, Nf_mf, epf_pf, epf_mf, Nfg, Rfg = 
                  calculate_compressor_speed_and_efficiency(FanMap, pf, mf, pifD, mbfD, NbfD, epf0, Ng = Nfg, Rg = Rfg)
            
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

            pt21, Tt21, ht21, st21, cpt21, Rt21,
            pt21_pt2,
            pt21_st2, Tt21_st2, ht21_st2, st21_st2,
            pt21_pf, Tt21_pf, ht21_pf, st21_pf,
            pt21_epf, Tt21_epf, ht21_epf, st21_epf,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_pratd(alpha, nair,
                  pt2, Tt2, ht2, st2, cpt2, Rt2, pf, epf)

            pt21_pf = pt21_epf * epf_pf + pt21_pf
            Tt21_pf = Tt21_epf * epf_pf + Tt21_pf
            ht21_pf = ht21_epf * epf_pf + ht21_pf
            st21_pf = st21_epf * epf_pf + st21_pf

            pt21_mf = pt21_epf * epf_mf + pt21_pt2 * pt2_mf
            Tt21_mf = Tt21_epf * epf_mf
            ht21_mf = ht21_epf * epf_mf
            st21_mf = st21_epf * epf_mf

            pt21_ml = pt21_pt2 * pt2_ml
            pt21_Mi = pt21_pt2 * pt2_Mi

            #---- fan duct nozzle total quantities
            pt7 = pt21 * pifn
            Tt7 = Tt21
            ht7 = ht21
            st7 = st21
            cpt7 = cpt21
            Rt7 = Rt21

            pt7_pf = pt21_pf * pifn
            Tt7_pf = Tt21_pf
            ht7_pf = ht21_pf
            st7_pf = st21_pf

            pt7_mf = pt21_mf * pifn
            Tt7_mf = Tt21_mf
            ht7_mf = ht21_mf
            st7_mf = st21_mf

            pt7_ml = pt21_ml * pifn
            ht7_ml = 0.0

            pt7_Mi = pt21_Mi * pifn
            ht7_Mi = 0.0
            # ===============================================================
            #---- Compressor precooler 19-19c
            pt19c = pt19 - Δp_PreC
            ht19c = ht19 + Δh_PreC
            Tt19c, Tt19c_ht19c, _ = gas_tsetd(alpha, nair, ht19c, Tt19)
            st19c, st19c_Tt19c, ht19c, ht29c_Tt29c, cpt19c, cpt19c_Tt19c, Rt19c = gassumd(alpha, nair, Tt19c)

             #Derivatives with respect to pressure are unchanged because pt19c_pt19 = 1
            pt19c_pt19 = 1.0
            ht19c_ht19 = 1.0

            Tt19c_Tt19 = Tt19c_ht19c * ht19c_ht19 / Tt19_ht19
            pt19c_mf = pt19_mf * pt19c_pt19
            pt19c_Mi = pt19_Mi * pt19c_pt19
            pt19c_ml = pt19_ml * pt19c_pt19

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
            u19c_Mi = (-1.0 / u19c) * h19c_Mi

            rho19c = p19c / (R19c * T19c)
            rho19c_p19c = 1.0 / (R19c * T19c)
            rho19c_T19c = -rho19c / T19c

            rho19c_Mi = rho19c_T19c * T19c_Mi +
                       rho19c_p19c * (p19c_pt19c * pt19c_Mi + p19c_Mi)
            rho19c_mf = rho19c_p19c * p19c_pt19c * pt19c_mf
            rho19c_ml = rho19c_p19c * p19c_pt19c * pt19c_ml


            #--------------------------------------------------------------
            #---- offtake mass ratio
            #cc   fo = mofft / mcore
            fo = mofft / ml * sqrt(Tt19c / Tref) * pref / pt19c
            fo_ml = -fo / ml - (fo / pt19c) * pt19c_ml
            fo_mf = -(fo / pt19c) * pt19c_mf
            fo_Mi = -(fo / pt19c) * pt19c_Mi

            #---- normalized power offtake Pofft / mcore
            Pom = Pofft / ml * sqrt(Tt19c / Tref) * pref / pt19c
            Pom_ml = -Pom / ml - (Pom / pt19c) * pt19c_ml
            Pom_mf = -(Pom / pt19c) * pt19c_mf
            Pom_Mi = -(Pom / pt19c) * pt19c_Mi
      
            # ===============================================================
            #---- LP compressor flow 2-25
            Nl, eplc, Nl_pl, Nl_ml, eplc_pl, eplc_ml, Nlcg, Rlcg = 
                  calculate_compressor_speed_and_efficiency(LPCMap, pl, ml, pilcD, mblcD, NblcD, eplc0, Ng = Nlcg, Rg = Rlcg)
                  
            if (eplc < 0.70)
                  eplc = 0.70
                  eplc_pl = 0.0
                  eplc_ml = 0.0
            end

            pt25, Tt25, ht25, st25, cpt25, Rt25,
            pt25_pt19c,
            pt25_st19c, Tt25_st19c, ht25_st19c, st25_st19c,
            pt25_pl, Tt25_pl, ht25_pl, st25_pl,
            pt25_eplc, Tt25_eplc, ht25_eplc, st25_eplc,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_pratd(alpha, nair,
                  pt19c, Tt19c, ht19c, st19c, cpt19c, Rt19c, pl, eplc)

            pt25_pl = pt25_eplc * eplc_pl + pt25_pl
            Tt25_pl = Tt25_eplc * eplc_pl + Tt25_pl
            ht25_pl = ht25_eplc * eplc_pl + ht25_pl
            st25_pl = st25_eplc * eplc_pl + st25_pl

            pt25_ml = pt25_eplc * eplc_ml + pt25_pt19c * pt19c_ml
            Tt25_ml = Tt25_eplc * eplc_ml
            ht25_ml = ht25_eplc * eplc_ml
            st25_ml = st25_eplc * eplc_ml

            pt25_mf = pt25_pt19c * pt19c_mf
            pt25_Mi = pt25_pt19c * pt19c_Mi

            Tt25_ht25 = 1 / cpt25
            st25_Tt25 = cpt25 / Tt25

            # ===============================================================
            #---- Compressor intercooler 25-25c
            pt25c = pt25 - Δp_InterC
            ht25c = ht25 + Δh_InterC
            Tt25c, Tt25c_ht25c, _ = gas_tsetd(alpha, nair, ht25c, Tt25)
            st25c, st25c_Tt25c, ht25c, ht25c_Tt25c, cpt25c, cpt25c_Tt25c, Rt25c = gassumd(alpha, nair, Tt25c)

            #Derivatives with respect to pressure are unchanged
            pt25c_pt25 = 1.0
            ht25c_ht25 = 1.0

            Tt25c_Tt25 = Tt25c_ht25c * ht25c_ht25 / Tt25_ht25
            st25c_st25 = st25c_Tt25c * Tt25c_Tt25 / st25_Tt25

            pt25c_pl = pt25_pl * pt25c_pt25
            Tt25c_pl = Tt25_pl * Tt25c_Tt25
            ht25c_pl = ht25_pl * ht25c_ht25
            st25c_pl = st25_pl * st25c_st25

            pt25c_ml = pt25_ml * pt25c_pt25
            Tt25c_ml = Tt25_ml * Tt25c_Tt25
            ht25c_ml = ht25_ml * ht25c_ht25
            st25c_ml = st25_ml * st25c_st25

            pt25c_mf = pt25_mf * pt25c_pt25
            pt25c_Mi = pt25_Mi * pt25c_pt25
            pt25c_ml = pt25_ml * pt25c_pt25
      
            # ===============================================================
            #---- HP compressor flow 25-3
            Nh, ephc, Nh_ph, Nh_mh, ephc_ph, ephc_mh, Nhcg, Rhcg = 
                  calculate_compressor_speed_and_efficiency(HPCMap, ph, mh, pihcD, mbhcD, NbhcD, ephc0, Ng = Nhcg, Rg = Rhcg)
            
            if (ephc < 0.70)
                  ephc = 0.70
                  ephc_ph = 0.0
                  ephc_mh = 0.0
            end

            pt3, Tt3, ht3, st3, cpt3, Rt3,
            pt3_pt25c,
            pt3_st25c, Tt3_st25c, ht3_st25c, st3_st25c,
            pt3_ph, Tt3_ph, ht3_ph, st3_ph,
            pt3_ephc, Tt3_ephc, ht3_ephc, st3_ephc,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_pratd(alpha, nair,
                  pt25c, Tt25c, ht25c, st25c, cpt25c, Rt25c, ph, ephc)

            pt3_pl = pt3_pt25c * pt25c_pl +
                     pt3_st25c * st25c_pl
            Tt3_pl = Tt3_st25c * st25c_pl
            ht3_pl = ht3_st25c * st25c_pl
            st3_pl = st3_st25c * st25c_pl

            pt3_ph = pt3_ephc * ephc_ph + pt3_ph
            Tt3_ph = Tt3_ephc * ephc_ph + Tt3_ph
            ht3_ph = ht3_ephc * ephc_ph + ht3_ph
            st3_ph = st3_ephc * ephc_ph + st3_ph

            pt3_mf = pt3_pt25c * pt25c_mf

            pt3_ml = pt3_pt25c * pt25c_ml +
                     pt3_st25c * st25c_ml
            Tt3_ml = Tt3_st25c * st25c_ml
            ht3_ml = ht3_st25c * st25c_ml
            st3_ml = st3_st25c * st25c_ml

            pt3_mh = pt3_ephc * ephc_mh
            Tt3_mh = Tt3_ephc * ephc_mh
            ht3_mh = ht3_ephc * ephc_mh
            st3_mh = st3_ephc * ephc_mh

            pt3_Mi = pt3_pt25c * pt25c_Mi

            Tt3_ht3 = 1 / cpt3

            # ===============================================================
            #---- burner fuel flow from Tt3-Tb difference   (ffb = m_fuel/m_burner)
            ffb, lambda,
            ffb_Tt3, ffb_Ttf, ffb_Tb,
            lam_Tt3, lam_Ttf, lam_Tb = gas_burnd(alpha, beta, gamma, n, ifuel, Tt3, Ttf, Tb, hvap)
            ffb_pl = ffb_Tt3 * Tt3_pl
            ffb_ph = ffb_Tt3 * Tt3_ph
            ffb_ml = ffb_Tt3 * Tt3_ml
            ffb_mh = ffb_Tt3 * Tt3_mh
            for i = 1:n
                  lam_pl[i] = lam_Tt3[i] * Tt3_pl
                  lam_ph[i] = lam_Tt3[i] * Tt3_ph
                  lam_ml[i] = lam_Tt3[i] * Tt3_ml
                  lam_mh[i] = lam_Tt3[i] * Tt3_mh
            end

            #---- all station 4 quantities
            Tt4 = Tb
            Tt4_Tb = 1.0
            st4, st4_Tt4,
            ht4, ht4_Tt4,
            cpt4, cpt4_Tt4, Rt4 = gassumd(lambda, nair, Tt4)

            #c     gassum(lambda,nair,Tt4,st4   ,dum,ht4   ,dum,cpt4   ,Rt4   )
            st4_pl, dum, ht4_pl, dum, cpt4_pl, Rt4_pl = gassum(lam_pl, nair, Tt4)
            st4_ph, dum, ht4_ph, dum, cpt4_ph, Rt4_ph = gassum(lam_ph, nair, Tt4)
            st4_ml, dum, ht4_ml, dum, cpt4_ml, Rt4_ml = gassum(lam_ml, nair, Tt4)
            st4_mh, dum, ht4_mh, dum, cpt4_mh, Rt4_mh = gassum(lam_mh, nair, Tt4)
            st4_Tb, dum, ht4_Tb, dum, cpt4_Tb, Rt4_Tb = gassum(lam_Tb, nair, Tt4)

            pt4 = pib * pt3
            pt4_pl = pib * pt3_pl
            pt4_ph = pib * pt3_ph
            pt4_ml = pib * pt3_ml
            pt4_mh = pib * pt3_mh
            pt4_Tb = 0.0
            pt4_Mi = pib * pt3_Mi

            st4_Tb = st4_Tb + st4_Tt4 * Tt4_Tb
            ht4_Tb = ht4_Tb + ht4_Tt4 * Tt4_Tb
            cpt4_Tb = cpt4_Tb + cpt4_Tt4 * Tt4_Tb

            # HSC: SEEMS TO BE FINE
            # ===============================================================
            if compare_strings(opt_cooling, "none")
                  #----- no cooling air present... station 41 is same as 4
                  pt41 = pt4
                  Tt41 = Tt4
                  ht41 = ht4
                  st41 = st4
                  cpt41 = cpt4
                  Rt41 = Rt4
                  for i = 1:nair
                        lambdap[i] = lambda[i]
                        lamp_pl[i] = lam_pl[i]
                        lamp_ph[i] = lam_ph[i]
                        lamp_mf[i] = 0.0
                        lamp_ml[i] = lam_ml[i]
                        lamp_mh[i] = lam_mh[i]
                        lamp_Tb[i] = lam_Tb[i]
                        lamp_Mi[i] = 0.0
                  end

                  pt41_pl = pt4_pl
                  pt41_ph = pt4_ph
                  pt41_ml = pt4_ml
                  pt41_mh = pt4_mh
                  pt41_Tb = pt4_Tb
                  pt41_Mi = pt4_Mi

                  Tt41_pl = 0.0
                  Tt41_ph = 0.0
                  Tt41_mf = 0.0
                  Tt41_ml = 0.0
                  Tt41_mh = 0.0
                  Tt41_Tb = 1.0
                  Tt41_Mi = 0.0

                  ht41_pl = ht4_pl
                  ht41_ph = ht4_ph
                  ht41_mf = 0.0
                  ht41_ml = ht4_ml
                  ht41_mh = ht4_mh
                  ht41_Tb = ht4_Tb
                  ht41_Mi = 0.0

                  st41_pl = st4_pl
                  st41_ph = st4_ph
                  st41_ml = st4_ml
                  st41_mh = st4_mh
                  st41_Tb = st4_Tb

                  cpt41_pl = cpt4_pl
                  cpt41_ph = cpt4_ph
                  cpt41_ml = cpt4_ml
                  cpt41_mh = cpt4_mh
                  cpt41_Tb = cpt4_Tb

                  Rt41_pl = Rt4_pl
                  Rt41_ph = Rt4_ph
                  Rt41_ml = Rt4_ml
                  Rt41_mh = Rt4_mh
                  Rt41_Tb = Rt4_Tb
                  #
                  fc = 0.0
                  fc_pl = 0.0
                  fc_ph = 0.0
                  fc_mf = 0.0
                  fc_ml = 0.0
                  fc_mh = 0.0
                  fc_Tb = 0.0
                  fc_Mi = 0.0

                  #----- set ff = m_fuel/m_core = ffb * m_burner/m_core
                  ff = (1.0 - fo) * ffb
                  ff_pl = (1.0 - fo) * ffb_pl
                  ff_ph = (1.0 - fo) * ffb_ph
                  ff_mf = -fo_mf * ffb
                  ff_ml = (1.0 - fo) * ffb_ml - fo_ml * ffb
                  ff_mh = (1.0 - fo) * ffb_mh
                  ff_Tb = (1.0 - fo) * ffb_Tb
                  ff_Mi = -fo_Mi * ffb

                  #----------------------------------------------------------------
            else
                  #----- cooling air is present... 

                  #----- hot-section temperature ratio for each blade row (for cooling model)
                  gmi4 = Rt4 / (cpt4 - Rt4)
                  gmi4_Rt4 = (1.0 + gmi4) / (cpt4 - Rt4)
                  gmi4_cpt4 = -gmi4 / (cpt4 - Rt4)
                  Trrat = 1.0 / (1.0 + 0.5 * gmi4 * Mtexit^2)
                  Trr_gmi4 = -Trrat / (1.0 + 0.5 * gmi4 * Mtexit^2) * 0.5 * Mtexit^2
                  Trr_pl = Trr_gmi4 * (gmi4_Rt4 * Rt4_pl + gmi4_cpt4 * cpt4_pl)
                  Trr_ph = Trr_gmi4 * (gmi4_Rt4 * Rt4_ph + gmi4_cpt4 * cpt4_ph)
                  Trr_ml = Trr_gmi4 * (gmi4_Rt4 * Rt4_ml + gmi4_cpt4 * cpt4_ml)
                  Trr_mh = Trr_gmi4 * (gmi4_Rt4 * Rt4_mh + gmi4_cpt4 * cpt4_mh)
                  Trr_Tb = Trr_gmi4 * (gmi4_Rt4 * Rt4_Tb + gmi4_cpt4 * cpt4_Tb)

                  # Heat exchanger to cool turbine cooling air
                  ht_tc = ht3 + Δh_TurbC #Specific enthalpy of turbine cooling air
                  Tt_tc, Tttc_httc, _ = gas_tsetd(alpha, nair, ht_tc, Tt3) #Temperature of turbine cooling air

                  httc_ht3 = 1.0
                  Tttc_Tt3 = Tttc_httc * httc_ht3 / Tt3_ht3

                  if compare_strings(opt_cooling, "fixed_coolingflowratio")
                        #------ epsrow(.) is assumed to be passed in.. calculate Tmrow(.)
                        Tmrow_copy = Tmcalc(ncrowx, ncrow,
                              Tt_tc, Tb, dTstrk, Trrat,
                              efilm, tfilm, StA, epsrow)
                        Tmrow[:] = Tmrow_copy[:]

                        #------ total cooling flow fraction
                        fc = 0.0
                        fc_fo = 0.0
                        for icrow = 1:ncrow
                              fc = fc + (1.0 - fo) * epsrow[icrow]
                              fc_fo = fc_fo - epsrow[icrow]
                        end
                        fc_pl = 0.0
                        fc_ph = 0.0
                        fc_mf = fc_fo * fo_mf
                        fc_ml = fc_fo * fo_ml
                        fc_mh = 0.0
                        fc_Tb = 0.0
                        fc_Mi = fc_fo * fo_Mi

                  else
                        #------ calculate cooling mass flow ratios epsrow(.) to get specified Tmrow(.)
                        ncrow, epsrow_copy, epsrow_Tttc, epsrow_Tb, epsrow_Trr = mcool(ncrowx,
                              Tmrow, Tt_tc, Tb, dTstrk, Trrat,
                              efilm, tfilm, StA)
                        epsrow[:] = epsrow_copy[:]

                        #------ total cooling flow fraction
                        fc = 0.0
                        fc_fo = 0.0
                        fc_Tttc = 0.0
                        fc_Tb = 0.0
                        fc_Trr = 0.0
                        for icrow = 1:ncrow
                              fc = fc + (1.0 - fo) * epsrow[icrow]
                              fc_fo = fc_fo - epsrow[icrow]
                              fc_Tttc = fc_Tttc + (1.0 - fo) * epsrow_Tttc[icrow]
                              fc_Tb = fc_Tb + (1.0 - fo) * epsrow_Tb[icrow]
                              fc_Trr = fc_Trr + (1.0 - fo) * epsrow_Trr[icrow]
                        end
                        fc_Tt3 = fc_Tttc * Tttc_Tt3

                        fc_pl = fc_Tt3 * Tt3_pl + fc_Trr * Trr_pl
                        fc_ph = fc_Tt3 * Tt3_ph + fc_Trr * Trr_ph
                        fc_mf = fc_fo * fo_mf
                        fc_ml = fc_Tt3 * Tt3_ml + fc_Trr * Trr_ml + fc_fo * fo_ml
                        fc_mh = fc_Tt3 * Tt3_mh + fc_Trr * Trr_mh
                        fc_Tb = fc_Tb + fc_Trr * Trr_Tb
                        fc_Mi = fc_fo * fo_Mi

                  end

                  #----- set ff = m_fuel/m_core = ffb * m_burner/m_core
                  ff = (1.0 - fo - fc) * ffb
                  ff_pl = (1.0 - fo - fc) * ffb_pl - fc_pl * ffb
                  ff_ph = (1.0 - fo - fc) * ffb_ph - fc_ph * ffb
                  ff_mf = -fc_mf * ffb - fo_mf * ffb
                  ff_ml = (1.0 - fo - fc) * ffb_ml - fc_ml * ffb - fo_ml * ffb
                  ff_mh = (1.0 - fo - fc) * ffb_mh - fc_mh * ffb
                  ff_Tb = (1.0 - fo - fc) * ffb_Tb - fc_Tb * ffb
                  ff_Mi = -fc_Mi * ffb - fo_Mi * ffb

                  #----- calculate station 41
                  pt4a = pt4
                  Tt4a = Tt4
                  ht4a = ht4
                  st4a = st4
                  cpt4a = cpt4
                  rt4a = Rt4

                  Tt4a_Tb = Tt4_Tb

                  pt4a_pl = pt4_pl
                  pt4a_ph = pt4_ph
                  pt4a_ml = pt4_ml
                  pt4a_mh = pt4_mh
                  pt4a_Tb = pt4_Tb
                  pt4a_Mi = pt4_Mi

                  ht4a_pl = ht4_pl
                  ht4a_ph = ht4_ph
                  ht4a_ml = ht4_ml
                  ht4a_mh = ht4_mh
                  ht4a_Tb = ht4_Tb

                  st4a_pl = st4_pl
                  st4a_ph = st4_ph
                  st4a_ml = st4_ml
                  st4a_mh = st4_mh
                  st4a_Tb = st4_Tb

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

                  p4a_pl = p4a_st4a * st4a_pl +
                           p4a_pt4a * pt4a_pl +
                           p4a_ht4a * ht4a_pl
                  T4a_pl = T4a_ht4a * ht4a_pl
                  h4a_pl = h4a_ht4a * ht4a_pl
                  s4a_pl = s4a_ht4a * ht4a_pl

                  p4a_ph = p4a_st4a * st4a_ph +
                           p4a_pt4a * pt4a_ph +
                           p4a_ht4a * ht4a_ph
                  T4a_ph = T4a_ht4a * ht4a_ph
                  h4a_ph = h4a_ht4a * ht4a_ph
                  s4a_ph = s4a_ht4a * ht4a_ph


                  p4a_ml = p4a_st4a * st4a_ml +
                           p4a_pt4a * pt4a_ml +
                           p4a_ht4a * ht4a_ml
                  T4a_ml = T4a_ht4a * ht4a_ml
                  h4a_ml = h4a_ht4a * ht4a_ml
                  s4a_ml = s4a_ht4a * ht4a_ml

                  p4a_mh = p4a_st4a * st4a_mh +
                           p4a_pt4a * pt4a_mh +
                           p4a_ht4a * ht4a_mh
                  T4a_mh = T4a_ht4a * ht4a_mh
                  h4a_mh = h4a_ht4a * ht4a_mh
                  s4a_mh = s4a_ht4a * ht4a_mh

                  p4a_Tb = p4a_st4a * st4a_Tb +
                           p4a_pt4a * pt4a_Tb +
                           p4a_ht4a * ht4a_Tb + p4a_Tt4a * Tt4a_Tb
                  T4a_Tb = T4a_ht4a * ht4a_Tb + T4a_Tt4a * Tt4a_Tb
                  h4a_Tb = h4a_ht4a * ht4a_Tb + h4a_Tt4a * Tt4a_Tb
                  s4a_Tb = s4a_ht4a * ht4a_Tb + s4a_Tt4a * Tt4a_Tb

                  p4a_Mi = p4a_pt4a * pt4a_Mi

                  for i = 1:nair
                        p4a_pl = p4a_pl + p_al[i] * lam_pl[i]
                        p4a_ph = p4a_ph + p_al[i] * lam_ph[i]
                        p4a_ml = p4a_ml + p_al[i] * lam_ml[i]
                        p4a_mh = p4a_mh + p_al[i] * lam_mh[i]
                        p4a_Tb = p4a_Tb + p_al[i] * lam_Tb[i]

                        T4a_pl = T4a_pl + T_al[i] * lam_pl[i]
                        T4a_ph = T4a_ph + T_al[i] * lam_ph[i]
                        T4a_ml = T4a_ml + T_al[i] * lam_ml[i]
                        T4a_mh = T4a_mh + T_al[i] * lam_mh[i]
                        T4a_Tb = T4a_Tb + T_al[i] * lam_Tb[i]

                        h4a_pl = h4a_pl + h_al[i] * lam_pl[i]
                        h4a_ph = h4a_ph + h_al[i] * lam_ph[i]
                        h4a_ml = h4a_ml + h_al[i] * lam_ml[i]
                        h4a_mh = h4a_mh + h_al[i] * lam_mh[i]
                        h4a_Tb = h4a_Tb + h_al[i] * lam_Tb[i]

                        s4a_pl = s4a_pl + s_al[i] * lam_pl[i]
                        s4a_ph = s4a_ph + s_al[i] * lam_ph[i]
                        s4a_ml = s4a_ml + s_al[i] * lam_ml[i]
                        s4a_mh = s4a_mh + s_al[i] * lam_mh[i]
                        s4a_Tb = s4a_Tb + s_al[i] * lam_Tb[i]
                  end

                  if (ht4a > h4a)
                        u4a = sqrt(2.0 * (ht4a - h4a))
                        u4a_pl = (ht4a_pl - h4a_pl) / u4a
                        u4a_ph = (ht4a_ph - h4a_ph) / u4a
                        u4a_ml = (ht4a_ml - h4a_ml) / u4a
                        u4a_mh = (ht4a_mh - h4a_mh) / u4a
                        u4a_Tb = (ht4a_Tb - h4a_Tb) / u4a
                  else
                        u4a = 0.0
                        u4a_pl = 0.0
                        u4a_ph = 0.0
                        u4a_ml = 0.0
                        u4a_mh = 0.0
                        u4a_Tb = 0.0
                  end

                  #----- exit speed of cooling air at station 4a
                  uc = ruc * u4a
                  uc_pl = ruc * u4a_pl
                  uc_ph = ruc * u4a_ph
                  uc_ml = ruc * u4a_ml
                  uc_mh = ruc * u4a_mh
                  uc_Tb = ruc * u4a_Tb

                  #----- IGV exit mixing
                  frac4 = (1.0 - fo - fc + ff) / (1.0 - fo + ff)
                  frac4_fo = -(1.0 - frac4) / (1.0 - fo + ff)
                  frac4_ff = (1.0 - frac4) / (1.0 - fo + ff)
                  frac4_fc = -1.0 / (1.0 - fo + ff)

                  fracm = fc / (1.0 - fo + ff)
                  fracm_fo = fracm / (1.0 - fo + ff)
                  fracm_ff = -fracm / (1.0 - fo + ff)
                  fracm_fc = 1.0 / (1.0 - fo + ff)

                  frac4_pl = frac4_fc * fc_pl + frac4_ff * ff_pl
                  frac4_ph = frac4_fc * fc_ph + frac4_ff * ff_ph
                  frac4_mf = frac4_fc * fc_mf + frac4_ff * ff_mf + frac4_fo * fo_mf
                  frac4_ml = frac4_fc * fc_ml + frac4_ff * ff_ml + frac4_fo * fo_ml
                  frac4_mh = frac4_fc * fc_mh + frac4_ff * ff_mh
                  frac4_Tb = frac4_fc * fc_Tb + frac4_ff * ff_Tb
                  frac4_Mi = frac4_fc * fc_Mi + frac4_ff * ff_Mi + frac4_fo * fo_Mi

                  fracm_pl = fracm_fc * fc_pl + fracm_ff * ff_pl
                  fracm_ph = fracm_fc * fc_ph + fracm_ff * ff_ph
                  fracm_mf = fracm_fc * fc_mf + fracm_ff * ff_mf + fracm_fo * fo_mf
                  fracm_ml = fracm_fc * fc_ml + fracm_ff * ff_ml + fracm_fo * fo_ml
                  fracm_mh = fracm_fc * fc_mh + fracm_ff * ff_mh
                  fracm_Tb = fracm_fc * fc_Tb + fracm_ff * ff_Tb
                  fracm_Mi = fracm_fc * fc_Mi + fracm_ff * ff_Mi + fracm_fo * fo_Mi

                  #----- mixed constituent fraction vector from mass equation
                  for i = 1:nair
                        lambdap[i] = frac4 * lambda[i] + fracm * alpha[i]
                        lamp_pl[i] = frac4_pl * lambda[i] + fracm_pl * alpha[i] +
                                     frac4 * lam_pl[i]
                        lamp_ph[i] = frac4_ph * lambda[i] + fracm_ph * alpha[i] +
                                     frac4 * lam_ph[i]
                        lamp_mf[i] = frac4_mf * lambda[i] + fracm_mf * alpha[i]
                        lamp_ml[i] = frac4_ml * lambda[i] + fracm_ml * alpha[i] +
                                     frac4 * lam_ml[i]
                        lamp_mh[i] = frac4_mh * lambda[i] + fracm_mh * alpha[i] +
                                     frac4 * lam_mh[i]
                        lamp_Tb[i] = frac4_Tb * lambda[i] + fracm_Tb * alpha[i] +
                                     frac4 * lam_Tb[i]
                        lamp_Mi[i] = frac4_Mi * lambda[i] + fracm_Mi * alpha[i]
                  end

                  #derivatives for turbine cooling air when there is a HX
                  httc_pl = ht3_pl * httc_ht3
                  httc_ph = ht3_ph * httc_ht3
                  httc_ml = ht3_ml * httc_ht3
                  httc_mh = ht3_mh * httc_ht3

                  #----- mixed total enthalpy from enthalpy equation
                  ht41 = frac4 * ht4 + fracm * ht_tc
                  ht41_pl = frac4_pl * ht4 + frac4 * ht4_pl + fracm_pl * ht_tc + fracm * httc_pl
                  ht41_ph = frac4_ph * ht4 + frac4 * ht4_ph + fracm_ph * ht_tc + fracm * httc_ph
                  ht41_mf = frac4_mf * ht4 + fracm_mf * ht_tc
                  ht41_ml = frac4_ml * ht4 + frac4 * ht4_ml + fracm_ml * ht_tc + fracm * httc_ml
                  ht41_mh = frac4_mh * ht4 + frac4 * ht4_mh + fracm_mh * ht_tc + fracm * httc_mh
                  ht41_Tb = frac4_Tb * ht4 + frac4 * ht4_Tb + fracm_Tb * ht_tc
                  ht41_Mi = frac4_mf * ht4 + fracm_Mi * ht_tc

                  #----- total temperature from total enthalpy
                  Tguess = frac4 * Tt4 + fracm * Tt_tc
                  Tt41, Tt41_ht41, T_al = gas_tsetd(lambdap, nair, ht41, Tguess)
                  Tt41_pl = Tt41_ht41 * ht41_pl
                  Tt41_ph = Tt41_ht41 * ht41_ph
                  Tt41_mf = Tt41_ht41 * ht41_mf
                  Tt41_ml = Tt41_ht41 * ht41_ml
                  Tt41_mh = Tt41_ht41 * ht41_mh
                  Tt41_Tb = Tt41_ht41 * ht41_Tb
                  Tt41_Mi = Tt41_ht41 * ht41_Mi
                  for i = 1:nair
                        Tt41_pl = Tt41_pl + T_al[i] * lamp_pl[i]
                        Tt41_ph = Tt41_ph + T_al[i] * lamp_ph[i]
                        Tt41_mf = Tt41_mf + T_al[i] * lamp_mf[i]
                        Tt41_ml = Tt41_ml + T_al[i] * lamp_ml[i]
                        Tt41_mh = Tt41_mh + T_al[i] * lamp_mh[i]
                        Tt41_Tb = Tt41_Tb + T_al[i] * lamp_Tb[i]
                        Tt41_Mi = Tt41_Mi + T_al[i] * lamp_Mi[i]
                  end

                  #----- will also need st41,cpt41,Rt41 derivatives
                  st41, st41_Tt41,
                  ht41, ht41_Tt41,
                  cpt41, cpt41_Tt41, Rt41 = gassumd(lambdap, nair, Tt41)
                  st41_pl = st41_Tt41 * Tt41_pl
                  st41_ph = st41_Tt41 * Tt41_ph
                  st41_mf = st41_Tt41 * Tt41_mf
                  st41_ml = st41_Tt41 * Tt41_ml
                  st41_mh = st41_Tt41 * Tt41_mh
                  st41_Tb = st41_Tt41 * Tt41_Tb
                  st41_Mi = st41_Tt41 * Tt41_Mi

                  cpt41_pl = cpt41_Tt41 * Tt41_pl
                  cpt41_ph = cpt41_Tt41 * Tt41_ph
                  cpt41_mf = cpt41_Tt41 * Tt41_mf
                  cpt41_ml = cpt41_Tt41 * Tt41_ml
                  cpt41_mh = cpt41_Tt41 * Tt41_mh
                  cpt41_Tb = cpt41_Tt41 * Tt41_Tb
                  cpt41_Mi = cpt41_Tt41 * Tt41_Mi

                  Rt41_pl = 0.0
                  Rt41_ph = 0.0
                  Rt41_mf = 0.0
                  Rt41_ml = 0.0
                  Rt41_mh = 0.0
                  Rt41_Tb = 0.0
                  Rt41_Mi = 0.0
                  for i = 1:nair
                        si, s_ti, hi, h_ti, cpi, Ri = gasfun(i, Tt41)
                        st41_pl = st41_pl + si * lamp_pl[i]
                        st41_ph = st41_ph + si * lamp_ph[i]
                        st41_mf = st41_mf + si * lamp_mf[i]
                        st41_ml = st41_ml + si * lamp_ml[i]
                        st41_mh = st41_mh + si * lamp_mh[i]
                        st41_Tb = st41_Tb + si * lamp_Tb[i]
                        st41_Mi = st41_Mi + si * lamp_Mi[i]

                        cpt41_pl = cpt41_pl + cpi * lamp_pl[i]
                        cpt41_ph = cpt41_ph + cpi * lamp_ph[i]
                        cpt41_mf = cpt41_mf + cpi * lamp_mf[i]
                        cpt41_ml = cpt41_ml + cpi * lamp_ml[i]
                        cpt41_mh = cpt41_mh + cpi * lamp_mh[i]
                        cpt41_Tb = cpt41_Tb + cpi * lamp_Tb[i]
                        cpt41_Mi = cpt41_Mi + cpi * lamp_Mi[i]

                        Rt41_pl = Rt41_pl + Ri * lamp_pl[i]
                        Rt41_ph = Rt41_ph + Ri * lamp_ph[i]
                        Rt41_mf = Rt41_mf + Ri * lamp_mf[i]
                        Rt41_ml = Rt41_ml + Ri * lamp_ml[i]
                        Rt41_mh = Rt41_mh + Ri * lamp_mh[i]
                        Rt41_Tb = Rt41_Tb + Ri * lamp_Tb[i]
                        Rt41_Mi = Rt41_Mi + Ri * lamp_Mi[i]
                  end

                  #----- mixed velocity from momentum equation, assuming constant static pressure
                  p41 = p4a
                  p41_pl = p4a_pl
                  p41_ph = p4a_ph
                  p41_ml = p4a_ml
                  p41_mh = p4a_mh
                  p41_Tb = p4a_Tb
                  p41_Mi = p4a_Mi

                  u41 = frac4 * u4a + fracm * uc
                  u41_pl = frac4_pl * u4a + frac4 * u4a_pl + fracm_pl * uc + fracm * uc_pl
                  u41_ph = frac4_ph * u4a + frac4 * u4a_ph + fracm_ph * uc + fracm * uc_ph
                  u41_mf = frac4_mf * u4a + fracm_mf * uc
                  u41_ml = frac4_ml * u4a + frac4 * u4a_ml + fracm_ml * uc + fracm * uc_ml
                  u41_mh = frac4_mh * u4a + frac4 * u4a_mh + fracm_mh * uc + fracm * uc_mh
                  u41_Tb = frac4_Tb * u4a + frac4 * u4a_Tb + fracm_Tb * uc + fracm * uc_Tb
                  u41_Mi = frac4_Mi * u4a + fracm_Mi * uc

                  #----- static temperature from static enthalpy
                  h41 = ht41 - 0.5 * u41^2
                  h41_pl = ht41_pl - u41 * u41_pl
                  h41_ph = ht41_ph - u41 * u41_ph
                  h41_mf = ht41_mf - u41 * u41_mf
                  h41_ml = ht41_ml - u41 * u41_ml
                  h41_mh = ht41_mh - u41 * u41_mh
                  h41_Tb = ht41_Tb - u41 * u41_Tb
                  h41_Mi = ht41_Mi - u41 * u41_Mi

                  Tguess = t4a + (h41 - h4a) / cp4a
                  T41, T41_h41, T_al = gas_tsetd(lambdap, nair, h41, Tguess)
                  T41_pl = T41_h41 * h41_pl
                  T41_ph = T41_h41 * h41_ph
                  T41_mf = T41_h41 * h41_mf
                  T41_ml = T41_h41 * h41_ml
                  T41_mh = T41_h41 * h41_mh
                  T41_Tb = T41_h41 * h41_Tb
                  T41_Mi = T41_h41 * h41_Mi
                  for i = 1:nair
                        T41_pl = T41_pl + T_al[i] * lamp_pl[i]
                        T41_ph = T41_ph + T_al[i] * lamp_ph[i]
                        T41_mf = T41_mf + T_al[i] * lamp_mf[i]
                        T41_ml = T41_ml + T_al[i] * lamp_ml[i]
                        T41_mh = T41_mh + T_al[i] * lamp_mh[i]
                        T41_Tb = T41_Tb + T_al[i] * lamp_Tb[i]
                        T41_Mi = T41_Mi + T_al[i] * lamp_Mi[i]
                  end
                  s41, s41_T41, h41, dum, cp41, R41 = gassum(lambdap, nair, T41)
                  s41_pl = s41_T41 * T41_pl
                  s41_ph = s41_T41 * T41_ph
                  s41_mf = s41_T41 * T41_mf
                  s41_ml = s41_T41 * T41_ml
                  s41_mh = s41_T41 * T41_mh
                  s41_Tb = s41_T41 * T41_Tb
                  s41_Mi = s41_T41 * T41_Mi
                  for i = 1:nair
                        si, s_ti, hi, h_ti, cpi, Ri = gasfun(i, T41)
                        s41_pl = s41_pl + si * lamp_pl[i]
                        s41_ph = s41_ph + si * lamp_ph[i]
                        s41_mf = s41_mf + si * lamp_mf[i]
                        s41_ml = s41_ml + si * lamp_ml[i]
                        s41_mh = s41_mh + si * lamp_mh[i]
                        s41_Tb = s41_Tb + si * lamp_Tb[i]
                        s41_Mi = s41_Mi + si * lamp_Mi[i]
                  end

                  #----- all stagnation quantities, from total-static enthalpy difference
                  dhb = ht41 - h41
                  dhb_pl = ht41_pl - h41_pl
                  dhb_ph = ht41_ph - h41_ph
                  dhb_mf = ht41_mf - h41_mf
                  dhb_ml = ht41_ml - h41_ml
                  dhb_mh = ht41_mh - h41_mh
                  dhb_Tb = ht41_Tb - h41_Tb
                  dhb_Mi = ht41_Mi - h41_Mi
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

                  pt41_pl = pt41_s41 * s41_pl +
                            pt41_p41 * p41_pl +
                            pt41_h41 * h41_pl + pt41_dhb * dhb_pl

                  pt41_ph = pt41_s41 * s41_ph +
                            pt41_p41 * p41_ph +
                            pt41_h41 * h41_ph + pt41_dhb * dhb_ph

                  pt41_mf = pt41_s41 * s41_mf +
                            pt41_h41 * h41_mf + pt41_dhb * dhb_mf

                  pt41_ml = pt41_s41 * s41_ml +
                            pt41_p41 * p41_ml +
                            pt41_h41 * h41_ml + pt41_dhb * dhb_ml

                  pt41_mh = pt41_s41 * s41_mh +
                            pt41_p41 * p41_mh +
                            pt41_h41 * h41_mh + pt41_dhb * dhb_mh

                  pt41_Tb = pt41_s41 * s41_Tb +
                            pt41_p41 * p41_Tb +
                            pt41_h41 * h41_Tb + pt41_dhb * dhb_Tb

                  pt41_Mi = pt41_s41 * s41_Mi +
                            pt41_p41 * p41_Mi +
                            pt41_h41 * h41_Mi + pt41_dhb * dhb_Mi

                  for i = 1:nair
                        pt41_pl = pt41_pl + p_al[i] * lamp_pl[i]
                        pt41_ph = pt41_ph + p_al[i] * lamp_ph[i]
                        pt41_mf = pt41_mf + p_al[i] * lamp_mf[i]
                        pt41_ml = pt41_ml + p_al[i] * lamp_ml[i]
                        pt41_mh = pt41_mh + p_al[i] * lamp_mh[i]
                        pt41_Tb = pt41_Tb + p_al[i] * lamp_Tb[i]
                        pt41_Mi = pt41_Mi + p_al[i] * lamp_Mi[i]

                  end

            end

            # ===============================================================
            #---- HPT and LPT work

            #---- bypass ratio
            BPR = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c
            BPR_mf = 1.0 / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c +
                     BPR / pt2 * pt2_mf -
                     BPR / pt19c * pt19c_mf
            BPR_ml = -BPR / ml +
                     BPR / pt2 * pt2_ml -
                     BPR / pt19c * pt19c_ml
            BPR_Mi = BPR / pt2 * pt2_Mi -
                     BPR / pt19c * pt19c_Mi

            #---- HPT work
            dhfac = -(1.0 - fo) / (1.0 - fo + ff) / (1.0 - epsh)
            dhfac_fo = dhfac / (1.0 - fo + ff) +
                       1.0 / (1.0 - fo + ff) / (1.0 - epsh)
            dhfac_ff = -dhfac / (1.0 - fo + ff)

            dhfac_pl = dhfac_ff * ff_pl
            dhfac_ph = dhfac_ff * ff_ph
            dhfac_mf = dhfac_ff * ff_mf + dhfac_fo * fo_mf
            dhfac_ml = dhfac_ff * ff_ml + dhfac_fo * fo_ml
            dhfac_mh = dhfac_ff * ff_mh
            dhfac_Tb = dhfac_ff * ff_Tb
            dhfac_Mi = dhfac_ff * ff_Mi + dhfac_fo * fo_Mi

            dhht = (ht3 - ht25c) * dhfac
            dhht_pl = (ht3 - ht25c) * dhfac_pl + (ht3_pl - ht25c_pl) * dhfac
            dhht_ph = (ht3 - ht25c) * dhfac_ph + (ht3_ph) * dhfac
            dhht_mf = (ht3 - ht25c) * dhfac_mf
            dhht_ml = (ht3 - ht25c) * dhfac_ml + (ht3_ml - ht25c_ml) * dhfac
            dhht_mh = (ht3 - ht25c) * dhfac_mh + (ht3_mh) * dhfac
            dhht_Tb = (ht3 - ht25c) * dhfac_Tb
            dhht_Mi = (ht3 - ht25c) * dhfac_Mi

            #---- LPT work
            dlfac = -1.0 / (1.0 - fo + ff) / (1.0 - epsl)
            dlfac_fo = dlfac / (1.0 - fo + ff)
            dlfac_ff = -dlfac / (1.0 - fo + ff)

            dlfac_pl = dlfac_ff * ff_pl
            dlfac_ph = dlfac_ff * ff_ph
            dlfac_mf = dlfac_ff * ff_mf + dlfac_fo * fo_mf
            dlfac_ml = dlfac_ff * ff_ml + dlfac_fo * fo_ml
            dlfac_mh = dlfac_ff * ff_mh
            dlfac_Tb = dlfac_ff * ff_Tb
            dlfac_Mi = dlfac_ff * ff_Mi + dlfac_fo * fo_Mi

            dhlt = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac
            dhlt_pf = (BPR * ht21_pf) * dlfac
            dhlt_pl = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac_pl +
                      (ht25_pl) * dlfac
            dhlt_ph = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac_ph
            dhlt_mf = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac_mf +
                      (BPR * ht21_mf + Pom_mf +
                       BPR_mf * (ht21 - ht2)) * dlfac
            dhlt_ml = (ht25 - ht19c + BPR * (ht21 - ht2)) * dlfac_ml +
                      (ht25_ml + BPR_ml * (ht21 - ht2) + Pom_ml) * dlfac
            dhlt_mh = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac_mh
            dhlt_Tb = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac_Tb
            dhlt_Mi = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac_Mi +
                      (BPR_Mi * (ht21 - ht2) + Pom_Mi) * dlfac

            #---- HPT corrected mass flow, using LPC corrected mass flow and fuel/air ratio
            mbht = ml * (1.0 - fo + ff) * sqrt(Tt41 / Tt19c) * pt19c / pt41
            mbht_ml = (1.0 - fo + ff) * sqrt(Tt41 / Tt19c) * pt19c / pt41
            mbht_ff = ml * sqrt(Tt41 / Tt19c) * pt19c / pt41
            mbht_fo = -ml * sqrt(Tt41 / Tt19c) * pt19c / pt41
            mbht_Tt41 = 0.5 * mbht / Tt41
            mbht_pt41 = -mbht / pt41
            mbht_pt19c = ml * (1.0 - fo + ff) * sqrt(Tt41 / Tt19c) / pt41

            mbht_pl = mbht_ff * ff_pl + mbht_Tt41 * Tt41_pl + mbht_pt41 * pt41_pl
            mbht_ph = mbht_ff * ff_ph + mbht_Tt41 * Tt41_ph + mbht_pt41 * pt41_ph
            mbht_mf = mbht_ff * ff_mf + mbht_Tt41 * Tt41_mf + mbht_pt41 * pt41_mf +
                      mbht_fo * fo_mf + mbht_pt19c * pt19c_mf
            mbht_ml = mbht_ff * ff_ml + mbht_Tt41 * Tt41_ml + mbht_pt41 * pt41_ml +
                      mbht_fo * fo_ml + mbht_pt19c * pt19c_ml + mbht_ml
            mbht_mh = mbht_ff * ff_mh + mbht_Tt41 * Tt41_mh + mbht_pt41 * pt41_mh
            mbht_Tb = mbht_ff * ff_Tb + mbht_Tt41 * Tt41_Tb + mbht_pt41 * pt41_Tb
            mbht_Mi = mbht_ff * ff_Mi + mbht_Tt41 * Tt41_Mi + mbht_pt41 * pt41_Mi +
                      mbht_fo * fo_Mi + mbht_pt19c * pt19c_Mi

            #---- HPT corrected speed
            Nbht = Nh * sqrt(Tt25c / Tt41)
            Nbht_Nh = sqrt(Tt25c / Tt41)
            Nbht_Tt25c = 0.5 * Nbht / Tt25c
            Nbht_Tt41 = -0.5 * Nbht / Tt41

            Nbht_pl = Nbht_Tt25c * Tt25c_pl + Nbht_Tt41 * Tt41_pl
            Nbht_ph = Nbht_Nh * Nh_ph + Nbht_Tt41 * Tt41_ph
            Nbht_mf = Nbht_Tt41 * Tt41_mf
            Nbht_ml = Nbht_Tt25c * Tt25c_ml + Nbht_Tt41 * Tt41_ml
            Nbht_mh = Nbht_Nh * Nh_mh + Nbht_Tt41 * Tt41_mh
            Nbht_Tb = Nbht_Tt41 * Tt41_Tb
            Nbht_Mi = Nbht_Tt41 * Tt41_Mi

            #---- HPT efficiency
            # HACK: HSC
            #First find uncooled HPT efficiency and derivatives
            epht1, epht1_dhht, epht1_mbht, epht1_Nbht, epht1_Tt41,
            epht1_cpt41, epht1_Rt41 = etmap(dhht, mbht, Nbht, pihtD, mbhtD, NbhtD, epht0, Tmaph,
                  Tt41, cpt41, Rt41)

            #Find cooled HPT efficiency epht
            epht = find_cooled_hpt_efficiency(epht1, epht_fc, fc0, fc)

            if (epht < 0.80)
                  epht = 0.80
                  epht1_dhht = 0.0
                  epht1_mbht = 0.0
                  epht1_Nbht = 0.0
                  epht1_Tt41 = 0.0
                  epht1_cpt41 = 0.0
                  epht1_Rt41 = 0.0
            end

            epht1_pl = epht1_dhht * dhht_pl + epht1_mbht * mbht_pl + epht1_Nbht * Nbht_pl
            epht1_ph = epht1_dhht * dhht_ph + epht1_mbht * mbht_ph + epht1_Nbht * Nbht_ph
            epht1_mf = epht1_dhht * dhht_mf + epht1_mbht * mbht_mf + epht1_Nbht * Nbht_mf
            epht1_ml = epht1_dhht * dhht_ml + epht1_mbht * mbht_ml + epht1_Nbht * Nbht_ml
            epht1_mh = epht1_dhht * dhht_mh + epht1_mbht * mbht_mh + epht1_Nbht * Nbht_mh
            epht1_Tb = epht1_dhht * dhht_Tb + epht1_mbht * mbht_Tb + epht1_Nbht * Nbht_Tb
            epht1_Mi = epht1_dhht * dhht_Mi + epht1_mbht * mbht_Mi + epht1_Nbht * Nbht_Mi

            epht1_pl = epht1_Tt41 * Tt41_pl + epht1_cpt41 * cpt41_pl +
                      epht1_Rt41 * Rt41_pl + epht1_pl
            epht1_ph = epht1_Tt41 * Tt41_ph + epht1_cpt41 * cpt41_ph +
                      epht1_Rt41 * Rt41_ph + epht1_ph
            epht1_mf = epht1_Tt41 * Tt41_mf + epht1_cpt41 * cpt41_mf +
                      epht1_Rt41 * Rt41_mf + epht1_mf
            epht1_ml = epht1_Tt41 * Tt41_ml + epht1_cpt41 * cpt41_ml +
                      epht1_Rt41 * Rt41_ml + epht1_ml
            epht1_mh = epht1_Tt41 * Tt41_mh + epht1_cpt41 * cpt41_mh +
                      epht1_Rt41 * Rt41_mh + epht1_mh
            epht1_Tb = epht1_Tt41 * Tt41_Tb + epht1_cpt41 * cpt41_Tb +
                      epht1_Rt41 * Rt41_Tb + epht1_Tb
            epht1_Mi = epht1_Tt41 * Tt41_Mi + epht1_cpt41 * cpt41_Mi +
                      epht1_Rt41 * Rt41_Mi + epht1_Mi

            #Chain rule for HPT efficiency derivatives
            epht_pl = epht1_pl + epht_fc * fc_pl
            epht_ph = epht1_ph + epht_fc * fc_ph
            epht_mf = epht1_mf + epht_fc * fc_mf
            epht_ml = epht1_ml + epht_fc * fc_ml
            epht_mh = epht1_mh + epht_fc * fc_mh
            epht_Tb = epht1_Tb + epht_fc * fc_Tb
            epht_Mi = epht1_Mi + epht_fc * fc_Mi
           
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

            pt45_epht = pt45_epi * (-epi / epht)

            pt45_pl = pt45_st41 * st41_pl +
                      pt45_pt41 * pt41_pl +
                      pt45_epht * epht_pl +
                      pt45_ht41 * ht41_pl + pt45_dhht * dhht_pl
            Tt45_pl = Tt45_ht41 * ht41_pl + Tt45_dhht * dhht_pl
            ht45_pl = ht45_ht41 * ht41_pl + ht45_dhht * dhht_pl
            st45_pl = st45_ht41 * ht41_pl + st45_dhht * dhht_pl

            pt45_ph = pt45_st41 * st41_ph +
                      pt45_pt41 * pt41_ph +
                      pt45_epht * epht_ph +
                      pt45_ht41 * ht41_ph + pt45_dhht * dhht_ph
            Tt45_ph = Tt45_ht41 * ht41_ph + Tt45_dhht * dhht_ph
            ht45_ph = ht45_ht41 * ht41_ph + ht45_dhht * dhht_ph
            st45_ph = st45_ht41 * ht41_ph + st45_dhht * dhht_ph

            pt45_mf = pt45_st41 * st41_mf +
                      pt45_pt41 * pt41_mf +
                      pt45_epht * epht_mf +
                      pt45_ht41 * ht41_mf + pt45_dhht * dhht_mf
            Tt45_mf = Tt45_ht41 * ht41_mf + Tt45_dhht * dhht_mf
            ht45_mf = ht45_ht41 * ht41_mf + ht45_dhht * dhht_mf
            st45_mf = st45_ht41 * ht41_mf + st45_dhht * dhht_mf

            pt45_ml = pt45_st41 * st41_ml +
                      pt45_pt41 * pt41_ml +
                      pt45_epht * epht_ml +
                      pt45_ht41 * ht41_ml + pt45_dhht * dhht_ml
            Tt45_ml = Tt45_ht41 * ht41_ml + Tt45_dhht * dhht_ml
            ht45_ml = ht45_ht41 * ht41_ml + ht45_dhht * dhht_ml
            st45_ml = st45_ht41 * ht41_ml + st45_dhht * dhht_ml

            pt45_mh = pt45_st41 * st41_mh +
                      pt45_pt41 * pt41_mh +
                      pt45_epht * epht_mh +
                      pt45_ht41 * ht41_mh + pt45_dhht * dhht_mh
            Tt45_mh = Tt45_ht41 * ht41_mh + Tt45_dhht * dhht_mh
            ht45_mh = ht45_ht41 * ht41_mh + ht45_dhht * dhht_mh
            st45_mh = st45_ht41 * ht41_mh + st45_dhht * dhht_mh

            pt45_Tb = pt45_st41 * st41_Tb +
                      pt45_pt41 * pt41_Tb +
                      pt45_epht * epht_Tb +
                      pt45_ht41 * ht41_Tb + pt45_dhht * dhht_Tb
            Tt45_Tb = Tt45_ht41 * ht41_Tb + Tt45_dhht * dhht_Tb
            ht45_Tb = ht45_ht41 * ht41_Tb + ht45_dhht * dhht_Tb
            st45_Tb = st45_ht41 * ht41_Tb + st45_dhht * dhht_Tb

            pt45_Mi = pt45_st41 * st41_Mi +
                      pt45_pt41 * pt41_Mi +
                      pt45_epht * epht_Mi +
                      pt45_ht41 * ht41_Mi + pt45_dhht * dhht_Mi
            Tt45_Mi = Tt45_ht41 * ht41_Mi + Tt45_dhht * dhht_Mi
            ht45_Mi = ht45_ht41 * ht41_Mi + ht45_dhht * dhht_Mi
            st45_Mi = st45_ht41 * ht41_Mi + st45_dhht * dhht_Mi


            for i = 1:nair
                  pt45_pl = pt45_pl + p_al[i] * lamp_pl[i]
                  Tt45_pl = Tt45_pl + T_al[i] * lamp_pl[i]
                  ht45_pl = ht45_pl + h_al[i] * lamp_pl[i]
                  st45_pl = st45_pl + s_al[i] * lamp_pl[i]

                  pt45_ph = pt45_ph + p_al[i] * lamp_ph[i]
                  Tt45_ph = Tt45_ph + T_al[i] * lamp_ph[i]
                  ht45_ph = ht45_ph + h_al[i] * lamp_ph[i]
                  st45_ph = st45_ph + s_al[i] * lamp_ph[i]

                  pt45_mf = pt45_mf + p_al[i] * lamp_mf[i]
                  Tt45_mf = Tt45_mf + T_al[i] * lamp_mf[i]
                  ht45_mf = ht45_mf + h_al[i] * lamp_mf[i]
                  st45_mf = st45_mf + s_al[i] * lamp_mf[i]

                  pt45_ml = pt45_ml + p_al[i] * lamp_ml[i]
                  Tt45_ml = Tt45_ml + T_al[i] * lamp_ml[i]
                  ht45_ml = ht45_ml + h_al[i] * lamp_ml[i]
                  st45_ml = st45_ml + s_al[i] * lamp_ml[i]

                  pt45_mh = pt45_mh + p_al[i] * lamp_mh[i]
                  Tt45_mh = Tt45_mh + T_al[i] * lamp_mh[i]
                  ht45_mh = ht45_mh + h_al[i] * lamp_mh[i]
                  st45_mh = st45_mh + s_al[i] * lamp_mh[i]

                  pt45_Tb = pt45_Tb + p_al[i] * lamp_Tb[i]
                  Tt45_Tb = Tt45_Tb + T_al[i] * lamp_Tb[i]
                  ht45_Tb = ht45_Tb + h_al[i] * lamp_Tb[i]
                  st45_Tb = st45_Tb + s_al[i] * lamp_Tb[i]

                  pt45_Mi = pt45_Mi + p_al[i] * lamp_Mi[i]
                  Tt45_Mi = Tt45_Mi + T_al[i] * lamp_Mi[i]
                  ht45_Mi = ht45_Mi + h_al[i] * lamp_Mi[i]
                  st45_Mi = st45_Mi + s_al[i] * lamp_Mi[i]
            end

            #---- will also need cpt45,Rt45 derivatives
            st45, st45_Tt45,
            ht45, ht45_Tt45,
            cpt45, cpt45_Tt45, Rt45 = gassumd(lambdap, nair, Tt45)
            cpt45_pl = cpt45_Tt45 * Tt45_pl
            cpt45_ph = cpt45_Tt45 * Tt45_ph
            cpt45_mf = cpt45_Tt45 * Tt45_mf
            cpt45_ml = cpt45_Tt45 * Tt45_ml
            cpt45_mh = cpt45_Tt45 * Tt45_mh
            cpt45_Tb = cpt45_Tt45 * Tt45_Tb
            cpt45_Mi = cpt45_Tt45 * Tt45_Mi

            Rt45_pl = 0.0
            Rt45_ph = 0.0
            Rt45_mf = 0.0
            Rt45_ml = 0.0
            Rt45_mh = 0.0
            Rt45_Tb = 0.0
            Rt45_Mi = 0.0
            for i = 1:nair
                  si, s_ti, hi, h_ti, cpi, Ri = gasfun(i, Tt45)
                  cpt45_pl = cpt45_pl + cpi * lamp_pl[i]
                  cpt45_ph = cpt45_ph + cpi * lamp_ph[i]
                  cpt45_mf = cpt45_mf + cpi * lamp_mf[i]
                  cpt45_ml = cpt45_ml + cpi * lamp_ml[i]
                  cpt45_mh = cpt45_mh + cpi * lamp_mh[i]
                  cpt45_Tb = cpt45_Tb + cpi * lamp_Tb[i]
                  cpt45_Mi = cpt45_Mi + cpi * lamp_Mi[i]

                  Rt45_pl = Rt45_pl + Ri * lamp_pl[i]
                  Rt45_ph = Rt45_ph + Ri * lamp_ph[i]
                  Rt45_mf = Rt45_mf + Ri * lamp_mf[i]
                  Rt45_ml = Rt45_ml + Ri * lamp_ml[i]
                  Rt45_mh = Rt45_mh + Ri * lamp_mh[i]
                  Rt45_Tb = Rt45_Tb + Ri * lamp_Tb[i]
                  Rt45_Mi = Rt45_Mi + Ri * lamp_Mi[i]
            end


            #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            #---- LPT corrected mass flow (in terms of Newton variables)
            mblt = ml * (1.0 - fo + ff) * sqrt(Tt45 / Tt19c) * pt19c / pt45
            mblt_ml = (1.0 - fo + ff) * sqrt(Tt45 / Tt19c) * pt19c / pt45
            mblt_ff = ml * sqrt(Tt45 / Tt19c) * pt19c / pt45
            mblt_fo = -ml * sqrt(Tt45 / Tt19c) * pt19c / pt45
            mblt_Tt45 = 0.5 * mblt / Tt45
            mblt_pt45 = -mblt / pt45
            mblt_pt19c = ml * (1.0 - fo + ff) * sqrt(Tt45 / Tt19c) / pt45

            mblt_pl = mblt_ff * ff_pl + mblt_Tt45 * Tt45_pl + mblt_pt45 * pt45_pl
            mblt_ph = mblt_ff * ff_ph + mblt_Tt45 * Tt45_ph + mblt_pt45 * pt45_ph
            mblt_mf = mblt_ff * ff_mf + mblt_Tt45 * Tt45_mf + mblt_pt45 * pt45_mf +
                      mblt_fo * fo_mf + mblt_pt19c * pt19c_mf
            mblt_ml = mblt_ff * ff_ml + mblt_Tt45 * Tt45_ml + mblt_pt45 * pt45_ml +
                      mblt_fo * fo_ml + mblt_pt19c * pt19c_ml + mblt_ml
            mblt_mh = mblt_ff * ff_mh + mblt_Tt45 * Tt45_mh + mblt_pt45 * pt45_mh
            mblt_Tb = mblt_ff * ff_Tb + mblt_Tt45 * Tt45_Tb + mblt_pt45 * pt45_Tb
            mblt_Mi = mblt_ff * ff_Mi + mblt_Tt45 * Tt45_Mi + mblt_pt45 * pt45_Mi +
                      mblt_fo * fo_Mi + mblt_pt19c * pt19c_Mi

            #---- LPT corrected speed
            Nblt = Nl * sqrt(Tt19c / Tt45)
            Nblt_Nl = sqrt(Tt19c / Tt45)
            Nblt_Tt45 = -0.5 * Nblt / Tt45

            Nblt_pl = Nblt_Nl * Nl_pl + Nblt_Tt45 * Tt45_pl
            Nblt_ph = Nblt_Tt45 * Tt45_ph
            Nblt_mf = Nblt_Tt45 * Tt45_mf
            Nblt_ml = Nblt_Nl * Nl_ml + Nblt_Tt45 * Tt45_ml
            Nblt_mh = Nblt_Tt45 * Tt45_mh
            Nblt_Tb = Nblt_Tt45 * Tt45_Tb
            Nblt_Mi = Nblt_Tt45 * Tt45_Mi

            #---- LPT efficiency
            eplt,
            eplt_dhlt, eplt_mblt, eplt_Nblt,
            eplt_Tt45, eplt_cpt45, eplt_Rt45 = etmap(dhlt, mblt, Nblt, piltD, mbltD, NbltD, eplt0, Tmapl,
                  Tt45, cpt45, Rt45)

            if (eplt < 0.80)
                  eplt = 0.80
                  eplt_dhlt = 0.0
                  eplt_mblt = 0.0
                  eplt_Nblt = 0.0
                  eplt_Tt45 = 0.0
                  eplt_cpt45 = 0.0
                  eplt_Rt45 = 0.0
            end

            eplt_pf = eplt_dhlt * dhlt_pf
            eplt_pl = eplt_dhlt * dhlt_pl + eplt_mblt * mblt_pl + eplt_Nblt * Nblt_pl
            eplt_ph = eplt_dhlt * dhlt_ph + eplt_mblt * mblt_ph + eplt_Nblt * Nblt_ph
            eplt_mf = eplt_dhlt * dhlt_mf + eplt_mblt * mblt_mf + eplt_Nblt * Nblt_mf
            eplt_ml = eplt_dhlt * dhlt_ml + eplt_mblt * mblt_ml + eplt_Nblt * Nblt_ml
            eplt_mh = eplt_dhlt * dhlt_mh + eplt_mblt * mblt_mh + eplt_Nblt * Nblt_mh
            eplt_Tb = eplt_dhlt * dhlt_Tb + eplt_mblt * mblt_Tb + eplt_Nblt * Nblt_Tb
            eplt_Mi = eplt_dhlt * dhlt_Mi + eplt_mblt * mblt_Mi + eplt_Nblt * Nblt_Mi

            eplt_pl = eplt_Tt45 * Tt45_pl + eplt_cpt45 * cpt45_pl +
                      eplt_Rt45 * Rt45_pl + eplt_pl
            eplt_ph = eplt_Tt45 * Tt45_ph + eplt_cpt45 * cpt45_ph +
                      eplt_Rt45 * Rt45_ph + eplt_ph
            eplt_mf = eplt_Tt45 * Tt45_mf + eplt_cpt45 * cpt45_mf +
                      eplt_Rt45 * Rt45_mf + eplt_mf
            eplt_ml = eplt_Tt45 * Tt45_ml + eplt_cpt45 * cpt45_ml +
                      eplt_Rt45 * Rt45_ml + eplt_ml
            eplt_mh = eplt_Tt45 * Tt45_mh + eplt_cpt45 * cpt45_mh +
                      eplt_Rt45 * Rt45_mh + eplt_mh
            eplt_Tb = eplt_Tt45 * Tt45_Tb + eplt_cpt45 * cpt45_Tb +
                      eplt_Rt45 * Rt45_Tb + eplt_Tb
            eplt_Mi = eplt_Tt45 * Tt45_Mi + eplt_cpt45 * cpt45_Mi +
                      eplt_Rt45 * Rt45_Mi + eplt_Mi


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

            pt49_Pc = pt49_pilt * pilt_Pc
            Tt49_Pc = Tt49_pilt * pilt_Pc
            ht49_Pc = ht49_pilt * pilt_Pc
            st49_Pc = st49_pilt * pilt_Pc

            pt49_pt45 = pt49_pilt * pilt_pt45 + pt49_pt45
            Tt49_pt45 = Tt49_pilt * pilt_pt45
            ht49_pt45 = ht49_pilt * pilt_pt45
            st49_pt45 = st49_pilt * pilt_pt45

            pt49_eplt = pt49_epi * epi_eplt
            Tt49_eplt = Tt49_epi * epi_eplt
            ht49_eplt = ht49_epi * epi_eplt
            st49_eplt = st49_epi * epi_eplt


            pt49_pf = pt49_eplt * eplt_pf
            Tt49_pf = Tt49_eplt * eplt_pf
            ht49_pf = ht49_eplt * eplt_pf
            st49_pf = st49_eplt * eplt_pf

            pt49_pl = pt49_pt45 * pt45_pl + pt49_st45 * st45_pl + pt49_eplt * eplt_pl
            Tt49_pl = Tt49_pt45 * pt45_pl + Tt49_st45 * st45_pl + Tt49_eplt * eplt_pl
            ht49_pl = ht49_pt45 * pt45_pl + ht49_st45 * st45_pl + ht49_eplt * eplt_pl
            st49_pl = st49_pt45 * pt45_pl + st49_st45 * st45_pl + st49_eplt * eplt_pl

            pt49_ph = pt49_pt45 * pt45_ph + pt49_st45 * st45_ph + pt49_eplt * eplt_ph
            Tt49_ph = Tt49_pt45 * pt45_ph + Tt49_st45 * st45_ph + Tt49_eplt * eplt_ph
            ht49_ph = ht49_pt45 * pt45_ph + ht49_st45 * st45_ph + ht49_eplt * eplt_ph
            st49_ph = st49_pt45 * pt45_ph + st49_st45 * st45_ph + st49_eplt * eplt_ph

            pt49_mf = pt49_pt45 * pt45_mf + pt49_st45 * st45_mf + pt49_eplt * eplt_mf
            Tt49_mf = Tt49_pt45 * pt45_mf + Tt49_st45 * st45_mf + Tt49_eplt * eplt_mf
            ht49_mf = ht49_pt45 * pt45_mf + ht49_st45 * st45_mf + ht49_eplt * eplt_mf
            st49_mf = st49_pt45 * pt45_mf + st49_st45 * st45_mf + st49_eplt * eplt_mf

            pt49_ml = pt49_pt45 * pt45_ml + pt49_st45 * st45_ml + pt49_eplt * eplt_ml
            Tt49_ml = Tt49_pt45 * pt45_ml + Tt49_st45 * st45_ml + Tt49_eplt * eplt_ml
            ht49_ml = ht49_pt45 * pt45_ml + ht49_st45 * st45_ml + ht49_eplt * eplt_ml
            st49_ml = st49_pt45 * pt45_ml + st49_st45 * st45_ml + st49_eplt * eplt_ml

            pt49_mh = pt49_pt45 * pt45_mh + pt49_st45 * st45_mh + pt49_eplt * eplt_mh
            Tt49_mh = Tt49_pt45 * pt45_mh + Tt49_st45 * st45_mh + Tt49_eplt * eplt_mh
            ht49_mh = ht49_pt45 * pt45_mh + ht49_st45 * st45_mh + ht49_eplt * eplt_mh
            st49_mh = st49_pt45 * pt45_mh + st49_st45 * st45_mh + st49_eplt * eplt_mh


            pt49_Tb = pt49_pt45 * pt45_Tb + pt49_st45 * st45_Tb + pt49_eplt * eplt_Tb
            Tt49_Tb = Tt49_pt45 * pt45_Tb + Tt49_st45 * st45_Tb + Tt49_eplt * eplt_Tb
            ht49_Tb = ht49_pt45 * pt45_Tb + ht49_st45 * st45_Tb + ht49_eplt * eplt_Tb
            st49_Tb = st49_pt45 * pt45_Tb + st49_st45 * st45_Tb + st49_eplt * eplt_Tb

            pt49_Mi = pt49_pt45 * pt45_Mi + pt49_st45 * st45_Mi + pt49_eplt * eplt_Mi
            Tt49_Mi = Tt49_pt45 * pt45_Mi + Tt49_st45 * st45_Mi + Tt49_eplt * eplt_Mi
            ht49_Mi = ht49_pt45 * pt45_Mi + ht49_st45 * st45_Mi + ht49_eplt * eplt_Mi
            st49_Mi = st49_pt45 * pt45_Mi + st49_st45 * st45_Mi + st49_eplt * eplt_Mi

            for i = 1:nair
                  pt49_pl = pt49_pl + p_al[i] * lamp_pl[i]
                  Tt49_pl = Tt49_pl + T_al[i] * lamp_pl[i]
                  ht49_pl = ht49_pl + h_al[i] * lamp_pl[i]
                  st49_pl = st49_pl + s_al[i] * lamp_pl[i]

                  pt49_ph = pt49_ph + p_al[i] * lamp_ph[i]
                  Tt49_ph = Tt49_ph + T_al[i] * lamp_ph[i]
                  ht49_ph = ht49_ph + h_al[i] * lamp_ph[i]
                  st49_ph = st49_ph + s_al[i] * lamp_ph[i]

                  pt49_mf = pt49_mf + p_al[i] * lamp_mf[i]
                  Tt49_mf = Tt49_mf + T_al[i] * lamp_mf[i]
                  ht49_mf = ht49_mf + h_al[i] * lamp_mf[i]
                  st49_mf = st49_mf + s_al[i] * lamp_mf[i]

                  pt49_ml = pt49_ml + p_al[i] * lamp_ml[i]
                  Tt49_ml = Tt49_ml + T_al[i] * lamp_ml[i]
                  ht49_ml = ht49_ml + h_al[i] * lamp_ml[i]
                  st49_ml = st49_ml + s_al[i] * lamp_ml[i]

                  pt49_mh = pt49_mh + p_al[i] * lamp_mh[i]
                  Tt49_mh = Tt49_mh + T_al[i] * lamp_mh[i]
                  ht49_mh = ht49_mh + h_al[i] * lamp_mh[i]
                  st49_mh = st49_mh + s_al[i] * lamp_mh[i]

                  pt49_Tb = pt49_Tb + p_al[i] * lamp_Tb[i]
                  Tt49_Tb = Tt49_Tb + T_al[i] * lamp_Tb[i]
                  ht49_Tb = ht49_Tb + h_al[i] * lamp_Tb[i]
                  st49_Tb = st49_Tb + s_al[i] * lamp_Tb[i]

                  pt49_Mi = pt49_Mi + p_al[i] * lamp_Mi[i]
                  Tt49_Mi = Tt49_Mi + T_al[i] * lamp_Mi[i]
                  ht49_Mi = ht49_Mi + h_al[i] * lamp_Mi[i]
                  st49_Mi = st49_Mi + s_al[i] * lamp_Mi[i]
            end
            Tt49_ht49 = 1 / cpt49
            st49_Tt49 = cpt49 / Tt49

            # ===============================================================
            #---- Regenerative cooling heat exchanger 49-49c
            pt49c = pt49 - Δp_Regen
            ht49c = ht49 + Δh_Regen
            Tt49c, Tt49c_ht49c, _ = gas_tsetd(lambdap, nair, ht49c, Tt49)
            st49c, st49c_Tt49c, ht49c, ht49c_Tt49, cpt49c, cpt49c_Tt49c, Rt49c = gassumd(lambdap, nair, Tt49c)

            #Derivatives with respect to pressure are unchanged
            pt49c_pt49 = 1.0
            ht49c_ht49 = 1.0

            Tt49c_Tt49 = Tt49c_ht49c * ht49c_ht49 / Tt49_ht49
            st49c_st49 = st49c_Tt49c * Tt49c_Tt49 / st49_Tt49

            #Apply chain rule
            pt49c_pf = pt49_pf * pt49c_pt49
            Tt49c_pf = Tt49_pf * pt49c_pt49
            ht49c_pf = ht49_pf * pt49c_pt49
            st49c_pf = st49_pf * pt49c_pt49

            pt49c_pl = pt49_pl * pt49c_pt49 
            Tt49c_pl = Tt49_pl * Tt49c_Tt49
            ht49c_pl = ht49_pl * ht49c_ht49
            st49c_pl = st49_pl * st49c_st49

            pt49c_ph = pt49_ph * pt49c_pt49  
            Tt49c_ph = Tt49_ph * Tt49c_Tt49
            ht49c_ph = ht49_ph * ht49c_ht49
            st49c_ph = st49_ph * st49c_st49

            pt49c_mf = pt49_mf * pt49c_pt49  
            Tt49c_mf = Tt49_mf * Tt49c_Tt49
            ht49c_mf = ht49_mf * ht49c_ht49
            st49c_mf = st49_mf * st49c_st49

            pt49c_ml = pt49_ml * pt49c_pt49 
            Tt49c_ml = Tt49_ml * Tt49c_Tt49
            ht49c_ml = ht49_ml * ht49c_ht49
            st49c_ml = st49_ml * st49c_st49

            pt49c_mh = pt49_mh * pt49c_pt49  
            Tt49c_mh = Tt49_mh * Tt49c_Tt49
            ht49c_mh = ht49_mh * ht49c_ht49
            st49c_mh = st49_mh * st49c_st49

            pt49c_Tb = pt49_Tb * pt49c_pt49  
            Tt49c_Tb = Tt49_Tb * Tt49c_Tt49
            ht49c_Tb = ht49_Tb * ht49c_ht49
            st49c_Tb = st49_Tb * st49c_st49

            pt49c_Pc = pt49_Pc * pt49c_pt49  
            Tt49c_Pc = Tt49_Pc * Tt49c_Tt49
            ht49c_Pc = ht49_Pc * ht49c_ht49
            st49c_Pc = st49_Pc * st49c_st49

            pt49c_Mi = pt49_Mi * pt49c_pt49  
            Tt49c_Mi = Tt49_Mi * Tt49c_Tt49
            ht49c_Mi = ht49_Mi * ht49c_ht49
            st49c_Mi = st49_Mi * st49c_st49

            #---- apply core nozzle loss
            pt5 = pt49c * pitn
            Tt5 = Tt49c
            ht5 = ht49c
            st5 = st49c
            cpt5 = cpt49c
            Rt5 = Rt49c

            pt5_pf = pt49c_pf * pitn
            Tt5_pf = Tt49c_pf
            ht5_pf = ht49c_pf
            st5_pf = st49c_pf

            pt5_pl = pt49c_pl * pitn
            Tt5_pl = Tt49c_pl
            ht5_pl = ht49c_pl
            st5_pl = st49c_pl

            pt5_ph = pt49c_ph * pitn
            Tt5_ph = Tt49c_ph
            ht5_ph = ht49c_ph
            st5_ph = st49c_ph

            pt5_mf = pt49c_mf * pitn
            Tt5_mf = Tt49c_mf
            ht5_mf = ht49c_mf
            st5_mf = st49c_mf

            pt5_ml = pt49c_ml * pitn
            Tt5_ml = Tt49c_ml
            ht5_ml = ht49c_ml
            st5_ml = st49c_ml

            pt5_mh = pt49c_mh * pitn
            Tt5_mh = Tt49c_mh
            ht5_mh = ht49c_mh
            st5_mh = st49c_mh

            pt5_Tb = pt49c_Tb * pitn
            Tt5_Tb = Tt49c_Tb
            ht5_Tb = ht49c_Tb
            st5_Tb = st49c_Tb

            pt5_Pc = pt49c_Pc * pitn
            Tt5_Pc = Tt49c_Pc
            ht5_Pc = ht49c_Pc
            st5_Pc = st49c_Pc

            pt5_Mi = pt49c_Mi * pitn
            Tt5_Mi = Tt49c_Mi
            ht5_Mi = ht49c_Mi
            st5_Mi = st49c_Mi

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

            if (M7 < 1.0)
                  #----- fan nozzle is not choked... assumed p7 = p0 is correct
                  p7_pt7 = p7_pfn * pfn_pt7 + p7_pt7
                  T7_pt7 = T7_pfn * pfn_pt7
                  h7_pt7 = h7_pfn * pfn_pt7
                  s7_pt7 = s7_pfn * pfn_pt7

                  p7_pf = p7_pt7 * pt7_pf + p7_st7 * st7_pf
                  T7_pf = T7_pt7 * pt7_pf + T7_st7 * st7_pf
                  h7_pf = h7_pt7 * pt7_pf + h7_st7 * st7_pf
                  s7_pf = s7_pt7 * pt7_pf + s7_st7 * st7_pf

                  p7_mf = p7_pt7 * pt7_mf + p7_st7 * st7_mf
                  T7_mf = T7_pt7 * pt7_mf + T7_st7 * st7_mf
                  h7_mf = h7_pt7 * pt7_mf + h7_st7 * st7_mf
                  s7_mf = s7_pt7 * pt7_mf + s7_st7 * st7_mf

                  p7_Mi = p7_pt7 * pt7_Mi
                  T7_Mi = T7_pt7 * pt7_Mi
                  h7_Mi = h7_pt7 * pt7_Mi
                  s7_Mi = s7_pt7 * pt7_Mi

            else
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
                  p7_pf = p7_st7 * st7_pf +
                          p7_pt7 * pt7_pf +
                          p7_ht7 * ht7_pf + p7_Tt7 * Tt7_pf
                  T7_pf = T7_ht7 * ht7_pf + T7_Tt7 * Tt7_pf
                  h7_pf = h7_ht7 * ht7_pf + h7_Tt7 * Tt7_pf
                  s7_pf = s7_ht7 * ht7_pf + s7_Tt7 * Tt7_pf

                  p7_mf = p7_st7 * st7_mf +
                          p7_pt7 * pt7_mf +
                          p7_ht7 * ht7_mf + p7_Tt7 * Tt7_mf
                  T7_mf = T7_ht7 * ht7_mf + T7_Tt7 * Tt7_mf
                  h7_mf = h7_ht7 * ht7_mf + h7_Tt7 * Tt7_mf
                  s7_mf = s7_ht7 * ht7_mf + s7_Tt7 * Tt7_mf

                  p7_Mi = p7_pt7 * pt7_Mi
                  T7_Mi = 0.0

            end

            if (ht7 > h7)
                  u7 = sqrt(2.0 * (ht7 - h7))
                  u7_pf = (ht7_pf - h7_pf) / u7
                  u7_mf = (ht7_mf - h7_mf) / u7
                  u7_Mi = ht7_Mi / u7
            else
                  u7 = 0.0
                  u7tmp = max(0.2 * u0, 0.02 * sqrt(Rt7 * Tt7))
                  u7_pf = (ht7_pf - h7_pf) / u7tmp
                  u7_mf = (ht7_mf - h7_mf) / u7tmp
                  u7_Mi = ht7_Mi / u7tmp
            end

            rho7 = p7 / (R7 * T7)
            rho7_pf = p7_pf / (R7 * T7) - (rho7 / T7) * T7_pf
            rho7_mf = p7_mf / (R7 * T7) - (rho7 / T7) * T7_mf
            rho7_Mi = p7_Mi / (R7 * T7) - (rho7 / T7) * T7_Mi

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

            if (M5 < 1.0)
                  #----- core nozzle is not choked... assumed p5 = p0 is correct
                  p5_pt5 = p5_pcn * pcn_pt5 + p5_pt5
                  T5_pt5 = T5_pcn * pcn_pt5
                  h5_pt5 = h5_pcn * pcn_pt5
                  s5_pt5 = s5_pcn * pcn_pt5

                  p5_pf = p5_pt5 * pt5_pf + p5_st5 * st5_pf
                  T5_pf = T5_pt5 * pt5_pf + T5_st5 * st5_pf
                  h5_pf = h5_pt5 * pt5_pf + h5_st5 * st5_pf
                  s5_pf = s5_pt5 * pt5_pf + s5_st5 * st5_pf

                  p5_pl = p5_pt5 * pt5_pl + p5_st5 * st5_pl
                  T5_pl = T5_pt5 * pt5_pl + T5_st5 * st5_pl
                  h5_pl = h5_pt5 * pt5_pl + h5_st5 * st5_pl
                  s5_pl = s5_pt5 * pt5_pl + s5_st5 * st5_pl

                  p5_ph = p5_pt5 * pt5_ph + p5_st5 * st5_ph
                  T5_ph = T5_pt5 * pt5_ph + T5_st5 * st5_ph
                  h5_ph = h5_pt5 * pt5_ph + h5_st5 * st5_ph
                  s5_ph = s5_pt5 * pt5_ph + s5_st5 * st5_ph

                  p5_mf = p5_pt5 * pt5_mf + p5_st5 * st5_mf
                  T5_mf = T5_pt5 * pt5_mf + T5_st5 * st5_mf
                  h5_mf = h5_pt5 * pt5_mf + h5_st5 * st5_mf
                  s5_mf = s5_pt5 * pt5_mf + s5_st5 * st5_mf

                  p5_ml = p5_pt5 * pt5_ml + p5_st5 * st5_ml
                  T5_ml = T5_pt5 * pt5_ml + T5_st5 * st5_ml
                  h5_ml = h5_pt5 * pt5_ml + h5_st5 * st5_ml
                  s5_ml = s5_pt5 * pt5_ml + s5_st5 * st5_ml

                  p5_mh = p5_pt5 * pt5_mh + p5_st5 * st5_mh
                  T5_mh = T5_pt5 * pt5_mh + T5_st5 * st5_mh
                  h5_mh = h5_pt5 * pt5_mh + h5_st5 * st5_mh
                  s5_mh = s5_pt5 * pt5_mh + s5_st5 * st5_mh

                  p5_Tb = p5_pt5 * pt5_Tb + p5_st5 * st5_Tb
                  T5_Tb = T5_pt5 * pt5_Tb + T5_st5 * st5_Tb
                  h5_Tb = h5_pt5 * pt5_Tb + h5_st5 * st5_Tb
                  s5_Tb = s5_pt5 * pt5_Tb + s5_st5 * st5_Tb

                  p5_Pc = p5_pt5 * pt5_Pc + p5_st5 * st5_Pc
                  T5_Pc = T5_pt5 * pt5_Pc + T5_st5 * st5_Pc
                  h5_Pc = h5_pt5 * pt5_Pc + h5_st5 * st5_Pc
                  s5_Pc = s5_pt5 * pt5_Pc + s5_st5 * st5_Pc

                  p5_Mi = p5_pt5 * pt5_Mi + p5_st5 * st5_Mi
                  T5_Mi = T5_pt5 * pt5_Mi + T5_st5 * st5_Mi
                  h5_Mi = h5_pt5 * pt5_Mi + h5_st5 * st5_Mi
                  s5_Mi = s5_pt5 * pt5_Mi + s5_st5 * st5_Mi

            else
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

                  p5_pf = p5_st5 * st5_pf +
                          p5_pt5 * pt5_pf +
                          p5_Tt5 * Tt5_pf + p5_ht5 * ht5_pf
                  T5_pf = T5_Tt5 * Tt5_pf + T5_ht5 * ht5_pf
                  h5_pf = h5_Tt5 * Tt5_pf + h5_ht5 * ht5_pf
                  s5_pf = s5_Tt5 * Tt5_pf + s5_ht5 * ht5_pf

                  p5_pl = p5_st5 * st5_pl +
                          p5_pt5 * pt5_pl +
                          p5_Tt5 * Tt5_pl + p5_ht5 * ht5_pl
                  T5_pl = T5_Tt5 * Tt5_pl + T5_ht5 * ht5_pl
                  h5_pl = h5_Tt5 * Tt5_pl + h5_ht5 * ht5_pl
                  s5_pl = s5_Tt5 * Tt5_pl + s5_ht5 * ht5_pl

                  p5_ph = p5_st5 * st5_ph +
                          p5_pt5 * pt5_ph +
                          p5_Tt5 * Tt5_ph + p5_ht5 * ht5_ph
                  T5_ph = T5_Tt5 * Tt5_ph + T5_ht5 * ht5_ph
                  h5_ph = h5_Tt5 * Tt5_ph + h5_ht5 * ht5_ph
                  s5_ph = s5_Tt5 * Tt5_ph + s5_ht5 * ht5_ph

                  p5_mf = p5_st5 * st5_mf +
                          p5_pt5 * pt5_mf +
                          p5_Tt5 * Tt5_mf + p5_ht5 * ht5_mf
                  T5_mf = T5_Tt5 * Tt5_mf + T5_ht5 * ht5_mf
                  h5_mf = h5_Tt5 * Tt5_mf + h5_ht5 * ht5_mf
                  s5_mf = s5_Tt5 * Tt5_mf + s5_ht5 * ht5_mf

                  p5_ml = p5_st5 * st5_ml +
                          p5_pt5 * pt5_ml +
                          p5_Tt5 * Tt5_ml + p5_ht5 * ht5_ml
                  T5_ml = T5_Tt5 * Tt5_ml + T5_ht5 * ht5_ml
                  h5_ml = h5_Tt5 * Tt5_ml + h5_ht5 * ht5_ml
                  s5_ml = s5_Tt5 * Tt5_ml + s5_ht5 * ht5_ml

                  p5_mh = p5_st5 * st5_mh +
                          p5_pt5 * pt5_mh +
                          p5_Tt5 * Tt5_mh + p5_ht5 * ht5_mh
                  T5_mh = T5_Tt5 * Tt5_mh + T5_ht5 * ht5_mh
                  h5_mh = h5_Tt5 * Tt5_mh + h5_ht5 * ht5_mh
                  s5_mh = s5_Tt5 * Tt5_mh + s5_ht5 * ht5_mh

                  p5_Tb = p5_st5 * st5_Tb +
                          p5_pt5 * pt5_Tb +
                          p5_Tt5 * Tt5_Tb + p5_ht5 * ht5_Tb
                  T5_Tb = T5_Tt5 * Tt5_Tb + T5_ht5 * ht5_Tb
                  h5_Tb = h5_Tt5 * Tt5_Tb + h5_ht5 * ht5_Tb
                  s5_Tb = s5_Tt5 * Tt5_Tb + s5_ht5 * ht5_Tb

                  p5_Pc = p5_st5 * st5_Pc +
                          p5_pt5 * pt5_Pc +
                          p5_Tt5 * Tt5_Pc + p5_ht5 * ht5_Pc
                  T5_Pc = T5_Tt5 * Tt5_Pc + T5_ht5 * ht5_Pc
                  h5_Pc = h5_Tt5 * Tt5_Pc + h5_ht5 * ht5_Pc
                  s5_Pc = s5_Tt5 * Tt5_Pc + s5_ht5 * ht5_Pc

                  p5_Mi = p5_st5 * st5_Mi +
                          p5_pt5 * pt5_Mi +
                          p5_Tt5 * Tt5_Mi + p5_ht5 * ht5_Mi
                  T5_Mi = T5_Tt5 * Tt5_Mi + T5_ht5 * ht5_Mi
                  h5_Mi = h5_Tt5 * Tt5_Mi + h5_ht5 * ht5_Mi
                  s5_Mi = s5_Tt5 * Tt5_Mi + s5_ht5 * ht5_Mi

            end

            R5_pl = zeros(T, 1)[1]
            R5_ph = zeros(T, 1)[1]
            R5_mf = zeros(T, 1)[1]
            R5_ml = zeros(T, 1)[1]
            R5_mh = zeros(T, 1)[1]
            R5_Tb = zeros(T, 1)[1]
            R5_Pc = zeros(T, 1)[1]
            R5_Mi = zeros(T, 1)[1]

            for i = 1:nair

                  p5_pl = p5_pl + p_al[i] * lamp_pl[i]
                  T5_pl = T5_pl + T_al[i] * lamp_pl[i]
                  h5_pl = h5_pl + h_al[i] * lamp_pl[i]
                  s5_pl = s5_pl + s_al[i] * lamp_pl[i]
                  R5_pl = R5_pl + R_al[i] * lamp_pl[i]

                  p5_ph = p5_ph + p_al[i] * lamp_ph[i]
                  T5_ph = T5_ph + T_al[i] * lamp_ph[i]
                  h5_ph = h5_ph + h_al[i] * lamp_ph[i]
                  s5_ph = s5_ph + s_al[i] * lamp_ph[i]
                  R5_ph = R5_ph + R_al[i] * lamp_ph[i]

                  p5_mf = p5_mf + p_al[i] * lamp_mf[i]
                  T5_mf = T5_mf + T_al[i] * lamp_mf[i]
                  h5_mf = h5_mf + h_al[i] * lamp_mf[i]
                  s5_mf = s5_mf + s_al[i] * lamp_mf[i]
                  R5_mf = R5_mf + R_al[i] * lamp_mf[i]

                  p5_ml = p5_ml + p_al[i] * lamp_ml[i]
                  T5_ml = T5_ml + T_al[i] * lamp_ml[i]
                  h5_ml = h5_ml + h_al[i] * lamp_ml[i]
                  s5_ml = s5_ml + s_al[i] * lamp_ml[i]
                  R5_ml = R5_ml + R_al[i] * lamp_ml[i]

                  p5_mh = p5_mh + p_al[i] * lamp_mh[i]
                  T5_mh = T5_mh + T_al[i] * lamp_mh[i]
                  h5_mh = h5_mh + h_al[i] * lamp_mh[i]
                  s5_mh = s5_mh + s_al[i] * lamp_mh[i]
                  R5_mh = R5_mh + R_al[i] * lamp_mh[i]

                  p5_Tb = p5_Tb + p_al[i] * lamp_Tb[i]
                  T5_Tb = T5_Tb + T_al[i] * lamp_Tb[i]
                  h5_Tb = h5_Tb + h_al[i] * lamp_Tb[i]
                  s5_Tb = s5_Tb + s_al[i] * lamp_Tb[i]
                  R5_Tb = R5_Tb + R_al[i] * lamp_Tb[i]

                  p5_Mi = p5_Mi + p_al[i] * lamp_Mi[i]
                  T5_Mi = T5_Mi + T_al[i] * lamp_Mi[i]
                  h5_Mi = h5_Mi + h_al[i] * lamp_Mi[i]
                  s5_Mi = s5_Mi + s_al[i] * lamp_Mi[i]
                  R5_Mi = R5_Mi + R_al[i] * lamp_Mi[i]
            end

            if (ht5 > h5)
                  u5 = sqrt(2.0 * (ht5 - h5))
                  u5_pf = (ht5_pf - h5_pf) / u5
                  u5_pl = (ht5_pl - h5_pl) / u5
                  u5_ph = (ht5_ph - h5_ph) / u5
                  u5_mf = (ht5_mf - h5_mf) / u5
                  u5_ml = (ht5_ml - h5_ml) / u5
                  u5_mh = (ht5_mh - h5_mh) / u5
                  u5_Tb = (ht5_Tb - h5_Tb) / u5
                  u5_Pc = (ht5_Pc - h5_Pc) / u5
                  u5_Mi = (ht5_Mi - h5_Mi) / u5
            else
                  u5 = 0.0
                  u5tmp = 0.05 * sqrt(Rt5 * Tt5)
                  u5_pf = (ht5_pf - h5_pf) / u5tmp
                  u5_pl = (ht5_pl - h5_pl) / u5tmp
                  u5_ph = (ht5_ph - h5_ph) / u5tmp
                  u5_mf = (ht5_mf - h5_mf) / u5tmp
                  u5_ml = (ht5_ml - h5_ml) / u5tmp
                  u5_mh = (ht5_mh - h5_mh) / u5tmp
                  u5_Tb = (ht5_Tb - h5_Tb) / u5tmp
                  u5_Pc = (ht5_Pc - h5_Pc) / u5tmp
                  u5_Mi = (ht5_Mi - h5_Mi) / u5tmp
            end

            rho5 = p5 / (R5 * T5)
            rho5_pf = p5_pf / (R5 * T5) - (rho5 / T5) * T5_pf
            rho5_pl = p5_pl / (R5 * T5) - (rho5 / T5) * T5_pl - (rho5 / R5) * R5_pl
            rho5_ph = p5_ph / (R5 * T5) - (rho5 / T5) * T5_ph - (rho5 / R5) * R5_ph
            rho5_mf = p5_mf / (R5 * T5) - (rho5 / T5) * T5_mf - (rho5 / R5) * R5_mf
            rho5_ml = p5_ml / (R5 * T5) - (rho5 / T5) * T5_ml - (rho5 / R5) * R5_ml
            rho5_mh = p5_mh / (R5 * T5) - (rho5 / T5) * T5_mh - (rho5 / R5) * R5_mh
            rho5_Tb = p5_Tb / (R5 * T5) - (rho5 / T5) * T5_Tb - (rho5 / R5) * R5_Tb
            rho5_Pc = p5_Pc / (R5 * T5) - (rho5 / T5) * T5_Pc - (rho5 / R5) * R5_Pc
            rho5_Mi = p5_Mi / (R5 * T5) - (rho5 / T5) * T5_Mi - (rho5 / R5) * R5_Mi

            # ===========================================================================
            #---- set up Newton system

            #---- fan/LPC speed constraint
            trf = sqrt(Tt2 / Tref)
            trl = sqrt(Tt19c / Tref)
            res[1, 1] = Gearf * trf * Nf - trl * Nl
            a[1, 1] = Gearf * trf * Nf_pf
            a[1, 2] = -trl * Nl_pl
            a[1, 3] = 0.0
            a[1, 4] = Gearf * trf * Nf_mf
            a[1, 5] = -trl * Nl_ml
            a[1, 6] = 0.0
            a[1, 7] = 0.0
            a[1, 8] = 0.0
            a[1, 9] = 0.0
            rrel[1] = res[1] / Nl

            #-------------------------------------------------------------------------
            #---- fixed corrected mass flow at LPT IGV (vertical-line LPT map)
            res[2, 1] = mblt - mbltD
            a[2, 1] = 0.0
            a[2, 2] = mblt_pl
            a[2, 3] = mblt_ph
            a[2, 4] = mblt_mf
            a[2, 5] = mblt_ml
            a[2, 6] = mblt_mh
            a[2, 7] = mblt_Tb
            a[2, 8] = 0.0
            a[2, 9] = mblt_Mi
            rrel[2] = res[2] / mbltD

            #-------------------------------------------------------------------------
            #---- fixed corrected mass flow at HPT IGV (vertical-line HPT map)
            res[3, 1] = mbht - mbhtD
            a[3, 1] = 0.0
            a[3, 2] = mbht_pl
            a[3, 3] = mbht_ph
            a[3, 4] = mbht_mf
            a[3, 5] = mbht_ml
            a[3, 6] = mbht_mh
            a[3, 7] = mbht_Tb
            a[3, 8] = 0.0
            a[3, 9] = mbht_Mi
            rrel[3] = res[3, 1] / mbhtD

            #-------------------------------------------------------------------------
            #---- fan nozzle mass flow, choked or unchoked
            mdotf = mf * sqrt(Tref / Tt2) * pt2 / pref
            mdotf_mf = sqrt(Tref / Tt2) * pt2 / pref
            mdotf_pt2 = mf * sqrt(Tref / Tt2) / pref
            res[4, 1] = mdotf - rho7 * u7 * A7
            a[4, 1] = -(rho7_pf * u7 + rho7 * u7_pf) * A7
            a[4, 2] = 0.0
            a[4, 3] = 0.0
            a[4, 4] = mdotf_mf + mdotf_pt2 * pt2_mf -
                      (rho7_mf * u7 + rho7 * u7_mf) * A7
            a[4, 5] = mdotf_pt2 * pt2_ml
            a[4, 6] = 0.0
            a[4, 7] = 0.0
            a[4, 8] = 0.0
            a[4, 9] = mdotf_pt2 * pt2_Mi -
                      (rho7_Mi * u7 + rho7 * u7_Mi) * A7
            rrel[4] = res[4, 1] / mdotf

            #-------------------------------------------------------------------------
            #---- core nozzle mass flow, choked or unchoked
            mdotc = (1.0 - fo + ff) * ml * sqrt(Tref / Tt19c) * pt19c / pref
            mdotc_ml = (1.0 - fo + ff) * sqrt(Tref / Tt19c) * pt19c / pref
            mdotc_fo = -ml * sqrt(Tref / Tt19c) * pt19c / pref
            mdotc_ff = ml * sqrt(Tref / Tt19c) * pt19c / pref
            mdotc_pt19c = (1.0 - fo + ff) * ml * sqrt(Tref / Tt19c) / pref

            mdotc_pl = mdotc_ff * ff_pl
            mdotc_ph = mdotc_ff * ff_ph
            mdotc_mf = mdotc_ff * ff_mf + mdotc_fo * fo_mf + mdotc_pt19c * pt19c_mf
            mdotc_ml = mdotc_ff * ff_ml + mdotc_fo * fo_ml + mdotc_pt19c * pt19c_ml +
                       mdotc_ml
            mdotc_mh = mdotc_ff * ff_mh
            mdotc_Tb = mdotc_ff * ff_Tb
            mdotc_Mi = mdotc_ff * ff_Mi + mdotc_fo * fo_Mi + mdotc_pt19c * pt19c_Mi

            res[5, 1] = mdotc - rho5 * u5 * A5
            a[5, 1] = -(rho5_pf * u5 + rho5 * u5_pf) * A5
            a[5, 2] = mdotc_pl - (rho5_pl * u5 + rho5 * u5_pl) * A5
            a[5, 3] = mdotc_ph - (rho5_ph * u5 + rho5 * u5_ph) * A5
            a[5, 4] = mdotc_mf - (rho5_mf * u5 + rho5 * u5_mf) * A5
            a[5, 5] = mdotc_ml - (rho5_ml * u5 + rho5 * u5_ml) * A5
            a[5, 6] = mdotc_mh - (rho5_mh * u5 + rho5 * u5_mh) * A5
            a[5, 7] = mdotc_Tb - (rho5_Tb * u5 + rho5 * u5_Tb) * A5
            a[5, 8] = -(rho5_Pc * u5 + rho5 * u5_Pc) * A5
            a[5, 9] = mdotc_Mi - (rho5_Mi * u5 + rho5 * u5_Mi) * A5
            rrel[5] = res[5, 1] / mdotc

            #-------------------------------------------------------------------------
            #---- LPC mass flow = HPC mass flow + mass offtake
            mdotl = ml * sqrt(Tref / Tt19c) * pt19c / pref
            mdotl_ml = sqrt(Tref / Tt19c) * pt19c / pref
            mdotl_pt19c = mdotl / pt19c

            mdotl_mf = mdotl_pt19c * pt19c_mf
            mdotl_ml = mdotl_pt19c * pt19c_ml + mdotl_ml
            mdotl_Mi = mdotl_pt19c * pt19c_Mi


            mdoth = mh * sqrt(Tref / Tt25c) * pt25c / pref
            mdoth_mh = sqrt(Tref / Tt25c) * pt25c / pref
            mdoth_Tt25c = -0.5 * mdoth / Tt25c
            mdoth_pt25c = mdoth / pt25c

            mdoth_pl = mdoth_Tt25c * Tt25c_pl + mdoth_pt25c * pt25c_pl
            mdoth_mf = mdoth_pt25c * pt25c_mf
            mdoth_ml = mdoth_Tt25c * Tt25c_ml + mdoth_pt25c * pt25c_ml
            mdoth_Mi = mdoth_pt25c * pt25c_Mi

            res[6, 1] = mdoth - mdotl + mofft
            a[6, 1] = 0.0
            a[6, 2] = mdoth_pl
            a[6, 3] = 0.0
            a[6, 4] = mdoth_mf - mdotl_mf
            a[6, 5] = mdoth_ml - mdotl_ml
            a[6, 6] = mdoth_mh
            a[6, 7] = 0.0
            a[6, 8] = 0.0
            a[6, 9] = mdoth_Mi - mdotl_Mi
            rrel[6] = res[6, 1] / mdoth

            #-------------------------------------------------------------------------

            if compare_strings(opt_calc_call, "oper_fixedTt4")
                  #----- specified Tt4 constraint
                  res[7, 1] = Tt4 - Tt4spec
                  a[7, 1] = 0.0
                  a[7, 2] = 0.0
                  a[7, 3] = 0.0
                  a[7, 4] = 0.0
                  a[7, 5] = 0.0
                  a[7, 6] = 0.0
                  a[7, 7] = Tt4_Tb
                  a[7, 8] = 0.0
                  a[7, 9] = 0.0
                  rrel[7] = res[7, 1] / Tt4spec

                  #--------------------------------------------------------
            else
                  #----- specified thrust constraint

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  #----- fan nozzle flow 7-8, use alpha mass fraction (air)
                  pfn = p0 / pt7
                  pfn_pt7 = -pfn / pt7
                  epi = 1.0
                  p8, T8, h8, s8, cp8, R8,
                  p8_pt7,
                  p8_st7, T8_st7, h8_st7, s8_st7,
                  p8_pfn, T8_pfn, h8_pfn, s8_pfn,
                  p8_epi, T8_epi, h8_epi, s8_epi,
                  p_al, T_al, h_al, s_al,
                  cp_al, R_al = gas_pratd(alpha, nair,
                        pt7, Tt7, ht7, st7, cpt7, Rt7, pfn, epi)
                  h8_pt7 = h8_pfn * pfn_pt7

                  h8_pf = h8_pt7 * pt7_pf + h8_st7 * st7_pf
                  h8_mf = h8_pt7 * pt7_mf + h8_st7 * st7_mf
                  h8_ml = h8_pt7 * pt7_ml
                  h8_Mi = h8_pt7 * pt7_Mi

                  if (ht7 > h8)
                        u8 = sqrt(2.0 * (ht7 - h8))
                        u8_pf = (ht7_pf - h8_pf) / u8
                        u8_mf = (ht7_mf - h8_mf) / u8
                        u8_ml = (ht7_ml - h8_ml) / u8
                        u8_Mi = (ht7_Mi - h8_Mi) / u8
                  else
                        u8 = 0.0
                        u8tmp = max(0.2 * u0, 0.01 * sqrt(Rt7 * Tt7))
                        u8_pf = (ht7_pf - h8_pf) / u8tmp
                        u8_mf = (ht7_mf - h8_mf) / u8tmp
                        u8_ml = (ht7_ml - h8_ml) / u8tmp
                        u8_Mi = (ht7_Mi - h8_Mi) / u8tmp
                  end
                  rho8 = p8 / (R8 * T8)

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  #----- core nozzle flow 5-6, use lambda mass fraction (combustion products)
                  pcn = p0 / pt5
                  pcn_pt5 = -pcn / pt5
                  epi = 1.0
                  p6, T6, h6, s6, cp6, R6,
                  p6_pt5,
                  p6_st5, T6_st5, h6_st5, s6_st5,
                  p6_pcn, T6_pcn, h6_pcn, s6_pcn,
                  p6_epi, T6_epi, h6_epi, s6_epi,
                  p_al, T_al, h_al, s_al,
                  cp_al, R_al = gas_pratd(lambdap, nair,
                        pt5, Tt5, ht5, st5, cpt5, Rt5, pcn, epi)
                  h6_pt5 = h6_pcn * pcn_pt5

                  h6_pf = h6_pt5 * pt5_pf + h6_st5 * st5_pf
                  h6_pl = h6_pt5 * pt5_pl + h6_st5 * st5_pl
                  h6_ph = h6_pt5 * pt5_ph + h6_st5 * st5_ph
                  h6_mf = h6_pt5 * pt5_mf + h6_st5 * st5_mf
                  h6_ml = h6_pt5 * pt5_ml + h6_st5 * st5_ml
                  h6_mh = h6_pt5 * pt5_mh + h6_st5 * st5_mh
                  h6_Tb = h6_pt5 * pt5_Tb + h6_st5 * st5_Tb
                  h6_Pc = h6_pt5 * pt5_Pc + h6_st5 * st5_Pc
                  h6_Mi = h6_pt5 * pt5_Mi + h6_st5 * st5_Mi
                  for i = 1:nair
                        h6_pl = h6_pl + h_al[i] * lamp_pl[i]
                        h6_ph = h6_ph + h_al[i] * lamp_ph[i]
                        h6_mf = h6_mf + h_al[i] * lamp_mf[i]
                        h6_ml = h6_ml + h_al[i] * lamp_ml[i]
                        h6_mh = h6_mh + h_al[i] * lamp_mh[i]
                        h6_Tb = h6_Tb + h_al[i] * lamp_Tb[i]
                        h6_Mi = h6_Mi + h_al[i] * lamp_Mi[i]
                  end

                  if (ht5 > h6)
                        u6 = sqrt(2.0 * (ht5 - h6))
                        u6_pf = (ht5_pf - h6_pf) / u6
                        u6_pl = (ht5_pl - h6_pl) / u6
                        u6_ph = (ht5_ph - h6_ph) / u6
                        u6_mf = (ht5_mf - h6_mf) / u6
                        u6_ml = (ht5_ml - h6_ml) / u6
                        u6_mh = (ht5_mh - h6_mh) / u6
                        u6_Tb = (ht5_Tb - h6_Tb) / u6
                        u6_Pc = (ht5_Pc - h6_Pc) / u6
                        u6_Mi = (ht5_Mi - h6_Mi) / u6
                  else
                        u6 = 0.0
                        u6tmp = max(0.2 * u0, 0.01 * sqrt(Rt5 * Tt5))
                        u6_pf = (ht5_pf - h6_pf) / u6tmp
                        u6_pl = (ht5_pl - h6_pl) / u6tmp
                        u6_ph = (ht5_ph - h6_ph) / u6tmp
                        u6_mf = (ht5_mf - h6_mf) / u6tmp
                        u6_ml = (ht5_ml - h6_ml) / u6tmp
                        u6_mh = (ht5_mh - h6_mh) / u6tmp
                        u6_Tb = (ht5_Tb - h6_Tb) / u6tmp
                        u6_Mi = (ht5_Mi - h6_Mi) / u6tmp
                  end

                  BPR = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c
                  BPR_mf = 1.0 / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c +
                           BPR / pt2 * pt2_mf - BPR / pt19c * pt19c_mf
                  BPR_ml = -BPR / ml +
                           BPR / pt2 * pt2_ml - BPR / pt19c * pt19c_ml
                  BPR_Mi = BPR / pt2 * pt2_Mi - BPR / pt19c * pt19c_Mi

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                  #----- overall thrust
                  if (u0 == 0.0)
                        Finl = 0.0
                  else
                        Finl = Phiinl / u0
                  end
                  F = ((1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9) * mdotl + Finl
                  F_ff = u6 * mdotl
                  F_u6 = (1.0 - fo + ff) * mdotl
                  F_BPR = (u8 - u0) * mdotl
                  F_u8 = BPR * mdotl
                  F_fo = (-u6 + u9) * mdotl
                  F_mdotl = (1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9

                  F_pf = F_u6 * u6_pf + F_u8 * u8_pf
                  F_pl = F_ff * ff_pl + F_u6 * u6_pl
                  F_ph = F_ff * ff_ph + F_u6 * u6_ph
                  F_mf = F_ff * ff_mf + F_u6 * u6_mf + F_u8 * u8_mf + F_BPR * BPR_mf
                  F_ml = F_ff * ff_ml + F_u6 * u6_ml + F_u8 * u8_ml + F_BPR * BPR_ml
                  F_mh = F_ff * ff_mh + F_u6 * u6_mh
                  F_Tb = F_ff * ff_Tb + F_u6 * u6_Tb
                  F_Pc = F_u6 * u6_Pc
                  F_Mi = F_ff * ff_Mi + F_u6 * u6_Mi + F_u8 * u8_Mi + F_BPR * BPR_Mi

                  F_mf = F_mf + F_fo * fo_mf + F_mdotl * mdotl_mf
                  F_ml = F_ml + F_fo * fo_ml + F_mdotl * mdotl_ml
                  F_Mi = F_Mi + F_fo * fo_Mi + F_mdotl * mdotl_Mi

                  #----- specified-F constraint
                  res[7, 1] = F - Fspec
                  a[7, 1] = F_pf
                  a[7, 2] = F_pl
                  a[7, 3] = F_ph
                  a[7, 4] = F_mf
                  a[7, 5] = F_ml
                  a[7, 6] = F_mh
                  a[7, 7] = F_Tb
                  a[7, 8] = F_Pc
                  a[7, 9] = F_Mi
                  rrel[7] = res[7, 1] / max(Fspec, F, 1e-6)


                  #      res[7] = Tt5
                  #      a(7,1) = Tt5_pf
                  #      a(7,2) = Tt5_pl
                  #      a(7,3) = Tt5_ph
                  #      a(7,4) = Tt5_mf
                  #      a(7,5) = Tt5_ml
                  #      a(7,6) = Tt5_mh
                  #      a(7,7) = Tt5_Tb
                  #      a(7,8) = Tt5_Pc
                  #      a(7,9) = Tt5_Mi


            end


            #------------------------------------------------------
            #---- LPT power matched to fan+LPC power
            epi = 1.0 / eplt
            epi_eplt = -epi / eplt
            pt5h, Tt5h, ht5h, st5h, cpt5h, Rt5h,
            pt5h_st45,
            pt5h_pt45,
            pt5h_epi,
            pt5h_ht45, Tt5h_ht45, ht5h_ht45, st5h_ht45,
            pt5h_dhlt, Tt5h_dhlt, ht5h_dhlt, st5h_dhlt,
            p_al, T_al, h_al, s_al,
            cp_al, R_al = gas_delhd(lambdap, nair,
                  pt45, Tt45, ht45, st45, cpt45, Rt45, dhlt, epi)

            pt5h_eplt = pt5h_epi * epi_eplt

            pt5h_pf = pt5h_eplt * eplt_pf +
                      pt5h_dhlt * dhlt_pf
            pt5h_pl = pt5h_st45 * st45_pl +
                      pt5h_pt45 * pt45_pl + pt5h_eplt * eplt_pl +
                      pt5h_ht45 * ht45_pl + pt5h_dhlt * dhlt_pl
            pt5h_ph = pt5h_st45 * st45_ph +
                      pt5h_pt45 * pt45_ph + pt5h_eplt * eplt_ph +
                      pt5h_ht45 * ht45_ph + pt5h_dhlt * dhlt_ph
            pt5h_mf = pt5h_st45 * st45_mf +
                      pt5h_pt45 * pt45_mf + pt5h_eplt * eplt_mf +
                      pt5h_ht45 * ht45_mf + pt5h_dhlt * dhlt_mf
            pt5h_ml = pt5h_st45 * st45_ml +
                      pt5h_pt45 * pt45_ml + pt5h_eplt * eplt_ml +
                      pt5h_ht45 * ht45_ml + pt5h_dhlt * dhlt_ml
            pt5h_mh = pt5h_st45 * st45_mh +
                      pt5h_pt45 * pt45_mh + pt5h_eplt * eplt_mh +
                      pt5h_ht45 * ht45_mh + pt5h_dhlt * dhlt_mh
            pt5h_Tb = pt5h_st45 * st45_Tb +
                      pt5h_pt45 * pt45_Tb + pt5h_eplt * eplt_Tb +
                      pt5h_ht45 * ht45_Tb + pt5h_dhlt * dhlt_Tb
            pt5h_Mi = pt5h_st45 * st45_Mi +
                      pt5h_pt45 * pt45_Mi + pt5h_eplt * eplt_Mi +
                      pt5h_ht45 * ht45_Mi + pt5h_dhlt * dhlt_Mi

            for i = 1:nair
                  pt5h_pl = pt5h_pl + p_al[i] * lamp_pl[i]
                  pt5h_ph = pt5h_ph + p_al[i] * lamp_ph[i]
                  pt5h_mf = pt5h_mf + p_al[i] * lamp_mf[i]
                  pt5h_ml = pt5h_ml + p_al[i] * lamp_ml[i]
                  pt5h_mh = pt5h_mh + p_al[i] * lamp_mh[i]
                  pt5h_Tb = pt5h_Tb + p_al[i] * lamp_Tb[i]
                  pt5h_Mi = pt5h_Mi + p_al[i] * lamp_Mi[i]
            end

            res[8, 1] = pt49 - pt5h
            a[8, 1] = pt49_pf - pt5h_pf
            a[8, 2] = pt49_pl - pt5h_pl
            a[8, 3] = pt49_ph - pt5h_ph
            a[8, 4] = pt49_mf - pt5h_mf
            a[8, 5] = pt49_ml - pt5h_ml
            a[8, 6] = pt49_mh - pt5h_mh
            a[8, 7] = pt49_Tb - pt5h_Tb
            a[8, 8] = pt49_Pc
            a[8, 9] = pt49_Mi - pt5h_Mi
            rrel[8] = res[8, 1] / pt49
            # ===========================================================================
            #---- inlet Mach - area constraint
            mfA = mf * sqrt(Tref / Tt2) * pt2 / pref / (rho2 * u2)
            mfA_mf = mfA * (pt2_mf / pt2 - rho2_mf / rho2) +
                     sqrt(Tref / Tt2) * pt2 / pref / (rho2 * u2)
            mfA_ml = mfA * (pt2_ml / pt2 - rho2_ml / rho2)
            mfA_Mi = mfA * (pt2_Mi / pt2 - rho2_Mi / rho2 - u2_Mi / u2)

            mlA = ml * sqrt(Tref / Tt19c) * pt19c / pref / (rho19c * u19c)
            mlA_mf = mlA * (pt19c_mf / pt19c - rho19c_mf / rho19c)
            mlA_ml = mlA * (pt19c_ml / pt19c - rho19c_ml / rho19c) +
                     sqrt(Tt19c / Tref) * pref / pt19c / (rho19c * u19c)
            mlA_Mi = mlA * (pt19c_Mi / pt19c - rho19c_Mi / rho19c - u19c_Mi / u19c)

            res[9, 1] = mfA + mlA - A2
            a[9, 1] = 0.0
            a[9, 2] = 0.0
            a[9, 3] = 0.0
            a[9, 4] = mfA_mf + mlA_mf
            a[9, 5] = mfA_ml + mlA_ml
            a[9, 6] = 0.0
            a[9, 7] = 0.0
            a[9, 8] = 0.0
            a[9, 9] = mfA_Mi + mlA_Mi
            rrel[9] = res[9, 1]


            # ===========================================================================

            for i = 1:9
                  rsav[i] = res[i, 1]
                  for k = 1:9
                        asav[i, k] = a[i, k]
                  end
            end


            if (iter == -2)
                  for i = 1:9
                        a_dlls[i, 1] = a[i, 1]
                        a_dlls[i, 2] = a[i, 2]
                        a_dlls[i, 3] = a[i, 3]
                        a_dlls[i, 4] = a[i, 4]
                        a_dlls[i, 5] = a[i, 5]
                        a_dlls[i, 6] = a[i, 6]
                        a_dlls[i, 7] = a[i, 7]
                        a_dlls[i, 8] = a[i, 8]
                        a_dlls[i, 9] = a[i, 9]
                        res_dlls[i] = res[i, 1]
                  end
            elseif (iter == -1)

                  for i = 1:9
                        dd = (res[i, 1] - res_dlls[i]) / eps
                        aa = (a[i, j] + a_dlls[i, j]) * 0.5
                        compare(ss, aa, dd)
                        println("Ja", i, j, aa, ss, "Jd", i, j, dd, ss)
                  end

                  error("tfoper.jl error due to iter == -1")

            end


            if (iter < 0)

                  println("TFOPER: Convergence failed.  opt_call_calc=", opt_calc_call)

                  Tt4 = Tb
                  pt5 = Pc
                  BPR = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c

                  println("pt18 Tt18 =", pt18, Tt18)
                  println("pt2  Tt2  =", pt2, Tt2)
                  println("pt25 Tt25 =", pt25, Tt25)
                  println("pt3  Tt3  =", pt3, Tt3)
                  println("pt4  Tt4  =", pt4, Tt4)
                  println("pt41 Tt41 =", pt41, Tt41)
                  println("pt45 Tt45 =", pt45, Tt45)
                  println("pt5  Tt5  =", pt5, Tt5)
                  println("p5        =", p5)
                  println("FPR  BPR  =", pf, BPR)

                  Lconv = false
                  return TSFC, Fsp, hfuel, ff,
                  Feng, mcore,
                  pif, pilc, pihc,
                  mbf, mblc, mbhc,
                  Nbf, Nblc, Nbhc,
                  Tt0, ht0, pt0, cpt0, Rt0,
                  Tt18, ht18, pt18, cpt18, Rt18,
                  Tt19, ht19, pt19, cpt19, Rt19,
                  Tt2, ht2, pt2, cpt2, Rt2,
                  Tt21, ht21, pt21, cpt21, Rt21,
                  Tt25, ht25, pt25, cpt25, Rt25,
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

            #---- solve Newton system and set Newton deltas
            res = gaussn(9, 9, a, res, 1)

            dMi = -res[9, 1]
            if (Mi >= Mimax && dMi > 0.0)
                  #----- Fan face is approaching choking, possibly due to an iteration transient
                  #-     Artificially limit it to Mimax
                  for i = 1:9
                        res[i, 1] = rsav[i]
                        for k = 1:9
                              a[i, k] = asav[i, k]
                        end
                  end

                  i = 9
                  for k = 1:9
                        a[i, k] = 0.0
                  end
                  res[i, 1] = Mi - Mimax
                  a[i, i] = 1.0

                  res = gaussn(9, 9, a, res, 1)
            end


            dpf = -res[1, 1]
            dpl = -res[2, 1]
            dph = -res[3, 1]
            dmf = -res[4, 1]
            dml = -res[5, 1]
            dmh = -res[6, 1]
            dTb = -res[7, 1]
            dPc = -res[8, 1]
            dMi = -res[9, 1]

            #---- max relative change
            dmax = max(abs(dpf) / pf,
                  abs(dpl) / pl,
                  abs(dph) / ph,
                  abs(dmf) / mf,
                  abs(dml) / ml,
                  abs(dmh) / mh,
                  abs(dTb) / Tb,
                  abs(dPc) / Pc,
                  abs(dMi) / Mi)

            #---- max,min allowed changes 
            pf0 = 1.0
            pl0 = 1.0

            dpfmax = 0.30 * (pf - pf0)
            dplmax = 0.25 * (pl - pl0)
            dphmax = 0.25 * (ph - 1.0)
            dmfmax = 0.20 * mf
            dmlmax = 0.20 * ml
            dmhmax = 0.20 * mh
            dTbmax = 0.50 * (Tb - Tt3)
            dPcmax = 1.00 * (Pc - p0)
            dMimax = 1.001 * (Mimax - Mi)

            dpfmin = -0.30 * (pf - pf0)
            dplmin = -0.25 * (pl - pl0)
            dphmin = -0.25 * (ph - 1.0)
            dmfmin = -0.20 * mf
            dmlmin = -0.20 * ml
            dmhmin = -0.20 * mh
            dTbmin = -0.30 * (Tb - Tt3)
            dPcmin = -0.50 * (Pc - p0)
            dMimin = -0.20 * Mi

            #---- set underrelaxation factor, 
            #-    if needed to limit any one change to its max value
            rlx = 1.0
            
            if (rlx * dpf > dpfmax)
                  rlx = dpfmax / dpf
            end
            if (rlx * dpl > dplmax)
                  rlx = dplmax / dpl
            end
            if (rlx * dph > dphmax)
                  rlx = dphmax / dph
            end
            if (rlx * dmf > dmfmax)
                  rlx = dmfmax / dmf
            end
            if (rlx * dml > dmlmax)
                  rlx = dmlmax / dml
            end
            if (rlx * dmh > dmhmax)
                  rlx = dmhmax / dmh
            end
            if (rlx * dTb > dTbmax)
                  rlx = dTbmax / dTb
            end
            if (rlx * dPc > dPcmax)
                  rlx = dPcmax / dPc
            end
            if (rlx * dMi > dMimax)
                  rlx = dMimax / dMi
            end

            if (rlx * dpf < dpfmin)
                  rlx = dpfmin / dpf
            end
            if (rlx * dpl < dplmin)
                  rlx = dplmin / dpl
            end
            if (rlx * dph < dphmin)
                  rlx = dphmin / dph
            end
            if (rlx * dmf < dmfmin)
                  rlx = dmfmin / dmf
            end
            if (rlx * dml < dmlmin)
                  rlx = dmlmin / dml
            end
            if (rlx * dmh < dmhmin)
                  rlx = dmhmin / dmh
            end
            if (rlx * dTb < dTbmin)
                  rlx = dTbmin / dTb
            end
            if (rlx * dPc < dPcmin)
                  rlx = dPcmin / dPc
            end
            if (rlx * dMi < dMimin)
                  rlx = dMimin / dMi
            end



            rlx = 1.0
            vrlx = " "

            if (rlx * dpf > dpfmax)
                  rlx = dpfmax / dpf
                  vrlx = "pf"
            end
            if (rlx * dpl > dplmax)
                  rlx = dplmax / dpl
                  vrlx = "pl"
            end
            if (rlx * dph > dphmax)
                  rlx = dphmax / dph
                  vrlx = "ph"
            end
            if (rlx * dmf > dmfmax)
                  rlx = dmfmax / dmf
                  vrlx = "mf"
            end
            if (rlx * dml > dmlmax)
                  rlx = dmlmax / dml
                  vrlx = "ml"
            end
            if (rlx * dmh > dmhmax)
                  rlx = dmhmax / dmh
                  vrlx = "mh"
            end
            if (rlx * dTb > dTbmax)
                  rlx = dTbmax / dTb
                  vrlx = "Tb"
            end
            if (rlx * dPc > dPcmax)
                  rlx = dPcmax / dPc
                  vrlx = "Pc"
            end
            if (rlx * dMi > dMimax)
                  rlx = dMimax / dMi
                  vrlx = "Mi"
            end

            if (rlx * dpf < dpfmin)
                  rlx = dpfmin / dpf
                  vrlx = "pf"
            end
            if (rlx * dpl < dplmin)
                  rlx = dplmin / dpl
                  vrlx = "pl"
            end
            if (rlx * dph < dphmin)
                  rlx = dphmin / dph
                  vrlx = "ph"
            end
            if (rlx * dmf < dmfmin)
                  rlx = dmfmin / dmf
                  vrlx = "mf"
            end
            if (rlx * dml < dmlmin)
                  rlx = dmlmin / dml
                  vrlx = "ml"
            end
            if (rlx * dmh < dmhmin)
                  rlx = dmhmin / dmh
                  vrlx = "mh"
            end
            if (rlx * dTb < dTbmin)
                  rlx = dTbmin / dTb
                  vrlx = "Tb"
            end
            if (rlx * dPc < dPcmin)
                  rlx = dPcmin / dPc
                  vrlx = "Pc"
            end
            if (rlx * dMi < dMimin)
                  rlx = dMimin / dMi
                  vrlx = "Mi"
            end

            #If iter>10, a limit cycle may have been reached
            #Apply a relaxation factor that does not oscillate
            rlx_it = 1.0 - 0.6*rand() * (iter > 10) #Relaxation based on RNG
            rlx = rlx * rlx_it

            #---- exit if convergence test is met or if max iterations reached
            if (dmax < toler) | (iter == itmax)

                  # ===============================================================
                  #---- pick up here if converged normally

                  #---- fan nozzle flow 7-8, use alpha mass fraction (air)
                  pt8 = pt7
                  ht8 = ht7
                  Tt8 = Tt7
                  st8 = st7
                  cpt8 = cpt7
                  Rt8 = Rt7

                  pfn = p0 / pt8
                  p8, T8, h8, s8, cp8, R8 = gas_prat(alpha, nair,
                        pt8, Tt8, ht8, st8, cpt8, Rt8, pfn, 1.0)
                  if (ht8 > h8)
                        u8 = sqrt(2.0 * (ht8 - h8))
                  else
                        u8 = 0.0
                  end
                  rho8 = p8 / (R8 * T8)
                  #
                  #---------------------------------------------------------------
                  #---- core nozzle flow 5-6, use lambda mass fraction (combustion products)
                  pt6 = pt5
                  ht6 = ht5
                  Tt6 = Tt5
                  st6 = st5
                  cpt6 = cpt5
                  Rt6 = Rt5

                  pcn = p0 / pt6
                  p6, T6, h6, s6, cp6, R6 = gas_prat(lambdap, nair,
                        pt6, Tt6, ht6, st6, cpt6, Rt6, pcn, 1.0)
                  if (ht6 > h6)
                        u6 = sqrt(2.0 * (ht6 - h6))
                  else
                        u6 = 0.0
                  end
                  rho6 = p6 / (R6 * T6)

                  # ===============================================================
                  #---- set final corrected mass flows, pressure ratios, corrected speeds
                  mbf = mf
                  mblc = ml
                  mbhc = mh

                  pif = pf
                  pilc = pl
                  pihc = ph

                  Nbf = Nl / Gearf * sqrt(Tt19c / Tt2)
                  Nblc = Nl
                  Nbhc = Nh

                  Tt4 = Tb
                  pt5 = Pc

                  #---- final core mass flow and bypass ratio
                  mcore = ml * sqrt(Tref / Tt19c) * pt19c / pref
                  BPR = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c

                  #---- offtake mass ratio
                  fo = mofft / mcore

                  #- - - - - - - - - - - - - - - - - - - - - - - - - - -
                  #---- effective fuel heating value, from states 3, 4  (just for info)
                  cpa = 0.5 * (cpt3 + cpt4)
                  hfuel = cpa * (Tt4 - Tt3 + ffb * (Tt4 - Ttf)) / (etab * ffb)

                  #---- effective fuel heating value, from states 3, 4.1  (just for info)
                  #      cpa = 0.5*(cpt3+cpt41)
                  #      hfuel = cpa*((1.0-fo)*(Tt41-Tt3) + ff*(Tt41-Ttf)) / (etab*ff)
                  #- - - - - - - - - - - - - - - - - - - - - - - - - - -

                  #---- overall effective thrust, including inlet defect credit
                  if (u0 == 0.0)
                        Finl = 0.0
                  else
                        Finl = Phiinl / u0
                  end
                  Feng = ((1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9) * mcore + Finl

                  #---- overall Fsp and TSFC
                  if (u0 == 0.0)
                        Fsp = 0.0
                  else
                        Fsp = Feng / (u0 * mcore * (1.0 + BPR))
                  end

                  if (Feng == 0.0)
                        TSFC = 0.0
                  else
                        TSFC = (gee * ff * mcore) / Feng
                  end

                  #-------------------------------------------------------------
                  #---- plume areas
                  A8 = BPR * mcore / (rho8 * u8)
                  A6 = (1.0 - fo + ff) * mcore / (rho6 * u6)

                  if (u9 == 0.0)
                        A9 = 0.0
                  else
                        A9 = fo * mcore / (rho9 * u9)
                  end

                  #---- plume Mach numbers
                  M8 = u8 / sqrt(T8 * R8 * cp8 / (cp8 - R8))
                  M6 = u6 / sqrt(T6 * R6 * cp6 / (cp6 - R6))

                  #--------------------------------------------------------------
                  #---- calculate component efficiencies  (informative only -- not needed here)
                  etaf = 0.0
                  etalc = 0.0
                  etahc = 0.0
                  etaht = 0.0
                  etalt = 0.0

                  #---- fan
                  pt21i, Tt21i, ht21i, st21i, cpt21i, Rt21i = gas_prat(alpha, nair,
                        pt2, Tt2, ht2, st2, cpt2, Rt2, pif, 1.0)
                  etaf = (ht21i - ht2) / (ht21 - ht2)

                  #---- LP compressor
                  pt25i, Tt25i, ht25i, st25i, cpt25i, Rt25i = gas_prat(alpha, nair,
                        pt19c, Tt19c, ht19c, st19c, cpt19c, Rt19c, pilc, 1.0)
                  etalc = (ht25i - ht19c) / (ht25 - ht19c)

                  #---- HP compressor

                  pt3i, Tt3i, ht3i, st3i, cpt3i, Rt3i = gas_prat(alpha, nair,
                        pt25c, Tt25c, ht25c, st25c, cpt25c, Rt25c, pihc, 1.0)
                  etahc = (ht3i - ht25c) / (ht3 - ht25c)

                  #---- HP turbine
                  piht = pt45 / pt41
                  pt45i, Tt45i, ht45i, st45i, cpt45i, Rt45i = gas_prat(lambdap, nair,
                        pt41, Tt41, ht41, st41, cpt41, Rt41, piht, 1.0)
                  etaht = (ht45 - ht41) / (ht45i - ht41)

                  #---- LP turbine
                  pilt = pt49 / pt45
                  pt49i, Tt49i, ht49i, st49i, cpt49i, Rt49i = gas_prat(lambdap, nair,
                        pt45, Tt45, ht45, st45, cpt45, Rt45, pilt, 1.0)
                  etalt = (ht5 - ht45) / (ht49i - ht45)

                  #--------------------------------------------------------------
                  #---- calculate fan-face static quantities from Mach number variable Mi
                  M2 = Mi
                  p2, T2, h2, s2, cp2, R2 = gas_mach(alpha, nair,
                        pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, M2, 1.0)
                  u2 = sqrt(2.0 * (ht2 - h2))

                  #---- calculate HPC-face static quantities from core mass flow
                  #-    (first check for choking)
                  M25s = 0.99
                  p25s, T25s, h25s, s25s, cp25s, R25s = gas_mach(alpha, nair,
                        pt25c, Tt25c, ht25c, st25c, cpt25c, Rt25c, 0.0, M25s, 1.0)
                  u25s = sqrt(2.0 * (ht25c - h25s))
                  mdot25s = p25s / (R25s * T25s) * u25s * A25

                  if ((1.0 - fo) * mcore >= mdot25s)
                        #----- HPC face is choked... artificially use sonic quantities
                        u25c = u25s
                        p25c = p25s
                        T25c = T25s
                        cp25c = cp25s
                        R25c = R25s
                        M25c = M25s
                  else
                        #----- normal case... set static quantities from prescribed mass flux
                        Mguess = min(M25, 0.90)
                        mp25 = (1.0 - fo) * mcore / A25
                        p25c, T25c, h25c, s25c, cp25c, R25c = gas_mass(alpha, nair,
                              pt25c, Tt25c, ht25c, st25c, cpt25c, Rt25c, mp25, Mguess)
                        u25c = sqrt(2.0 * (ht25c - h25c))
                        M25c = u25c / sqrt(T25c * R25c * cp25c / (cp25c - R25c))
                  end

                  if (ht5 < h5) #error if total enthalpy at nozzle inlet is lower than static enthalpy
                        error("ht5 < h5 : ", ht5, h5)
                  end
                  
                  if (dmax < toler) #Convergence achieved
                        Lconv = true

                  else #Finish run as max iterations reached
                        Lconv = false
                        
                        # println("TFOPER: Convergence failed.  iTFspec=", iTFspec)

                        # Tt4 = Tb
                        # pt5 = Pc
                        # BPR = mf / ml * sqrt(Tt19c / Tt2) * pt2 / pt19c

                        # println("pt18 Tt18 =", pt18, Tt18)
                        # println("pt2  Tt2  =", pt2, Tt2)
                        # println("pt25 Tt25 =", pt25, Tt25)
                        # println("pt3  Tt3  =", pt3, Tt3)
                        # println("pt4  Tt4  =", pt4, Tt4)
                        # println("pt41 Tt41 =", pt41, Tt41)
                        # println("pt45 Tt45 =", pt45, Tt45)
                        # println("pt5  Tt5  =", pt5, Tt5)
                        # println("p5        =", p5)
                        # println("FPR  BPR  =", pf, BPR)
                  end
                 
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
            

            #---- Newton update
            pf = pf + rlx * dpf
            pl = pl + rlx * dpl
            ph = ph + rlx * dph
            mf = mf + rlx * dmf
            ml = ml + rlx * dml
            mh = mh + rlx * dmh
            Tb = Tb + rlx * dTb
            Pc = Pc + rlx * dPc
            Mi = Mi + rlx * dMi

            Mi = min(Mi, Mimax)

      end # with next Newton iteration

end # tfoper


