"""
    tfsize!(gee, M0, T0, p0, a0, M2, M25,
      Feng, Phiinl, Kinl, eng_has_BLI_cores,
      BPR, pif, pilc, pihc,
      pid, pib, pifn, pitn,
      Ttf, ifuel, etab,
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

Turbofan performance and sizing routine.

Calculation procedure follows that of Kerrebrock, but the usual gas property formulas
are replaced by function calls (described in [Gas Calculations](@ref)), which can
therefore implement more general gas models. In addition, a turbine cooling model is
added.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gee`:     gravity acceleration
    - `M0`:      freestream Mach
    - `T0`:      freestream temperature  [K]
    - `p0`:      freestream pressure  [Pa]
    - `M2`:      fan-face Mach number
    - `M25`:     HPC-face Mach number
    - `Feng`:    required net thrust  (PK_inl+PK_out-Phi_jet)/u0  =  sum( mdot u)
    - `Phiinl`:  inlet ingested dissipation
    - `eng_has_BLI_cores`:   false=core in clear flow, true=core sees `Phiinl`
    - `BPR`:     bypass ratio  = `mdot_fan/mdot_core`
    - `pif`:     fan      pressure ratio  ( = pt7 /pt2)
    - `pilc`:    LP comp  pressure ratio  ( = pt25/pt2)
    - `pihc`:    HP comp  pressure ratio  ( = pt3 /pt25)
    - `pid`:     diffuser pressure ratio  ( = pt2 /pt0)
    - `pib`:     burner   pressure ratio  ( = pt4 /pt3)
    - `pifn`:    fan     nozzle pressure ratio  ( = pt7/pt2.1)
    - `pitn`:    turbine nozzle pressure ratio  ( = pt5/pt4.9)
    - `Ttf`:     fuel temperature entering combustor
    - `ifuel`:   fuel index, see function [`gasfun`](@ref)
    - `hvap`:    fuel enthalpy of vaporization (J/kg)
    - `etab`:    combustor efficiency (fraction of fuel burned)
    - `epf0`:    fan max polytropic efficiency
    - `eplc0`:   LPC max polytropic efficiency
    - `ephc0`:   HPC max polytropic efficiency
    - `epht0`:   HPT max polytropic efficiency
    - `eplt0`:   LPT max polytropic efficiency
      
    - `mofft`:   mass flow offtake at LPC discharge station 2.5
    - `Pofft`:   low spool power offtake
    - `Tt9`:     offtake air discharge total temperature
    - `pt9`:     offtake air discharge total pressure
    - `epsl`:    low  spool power loss fraction
    - `epsh`:    high spool power loss fraction
      
    - `opt_cooling`: turbine cooling flag
      - `"none"` = no cooling, ignore all cooling parameters below
      - `"fixed_coolingflowratio"` = usual cooling, using passed-in fcool
      - `"fixed_Tmetal"` = usual cooling, but set (and return) `fcool` from `Tmetal`
    - `Mtexit`:   turbine blade-row exit Mach, for setting temperature drops
    - `dTstrk`:   hot-streak temperature delta [K], used only if `opt_cooling="fixed_Tmetal"`
    - `StA`:      area-weighted Stanton number    , used only if `opt_cooling="fixed_Tmetal"`
    - `M4a`:      effective Mach at cooling-flow outlet (start of mixing)
    - `ruc`:      cooling-flow outlet velocity ratio, u/ue
    - `ncrowx`:      dimension of epsrow array
    - `ncrow`:       number of blade rows requiring cooling
    - `epsrow(.)`:   specified cooling-flow bypass ratio if `opt_cooling="fixed_coolingflowratio"`
    - `Tmrow(.)`:    specified metal temperature [K] if `opt_cooling="fixed_Tmetal"`


      **Outputs:**
    - `epsrow(.)`:   resulting cooling-flow bypass ratio if `opt_cooling="fixed_Tmetal"`
    - `Tmrow(.)`:    resulting metal temperature [K] if `opt_cooling="fixed_coolingflowratio"`
    - `TSFC`:    thrust specific fuel consumption = `mdot_fuel g / F`   [1/s]
    - `Fsp`:     specific thrust  = `F / (mdot u0) = F / ((1+BPR) mdot_core u0)`
    - `hfuel`:   fuel heating value   [J / kg K]
    - `ff`:      fuel mass flow fraction  =  `mdot_fuel / mdot_core`
    - `mcore`:   core mass flow = `mdot_core`  [kg/s]
    - `A2`:      fan-face area [m^2]
    - `A25`:     HPC-face area [m^2]
    - `A5`:      core nozzle area [m^2]
    - `A7`:      fan  nozzle area [m^2]
    - `A6`:      core plume  area [m^2]
    - `A8`:      fan  plume  area [m^2]
    - `Tt?`:     total temperature
    - `ht?`:     total complete enthalpy (includes heat of formation)
    - `pt?`:     total pressure
    - `cpt?`:    specific heat at stagnation temperature  (= `dh/dT`)
    - `Rt?`:     gas constant  at stagnation conditions
    - `T?`:      static temperature
    - `u?`:      velocity
    - `epf`:     fan polytropic efficiency
    - `eplc`:    LPC polytropic efficiency
    - `ephc`:    HPC polytropic efficiency
    - `epht`:    HPT polytropic efficiency
    - `eplt`:    LPT polytropic efficiency
    - `etaf`:    fan overall efficiency
    - `etalc`:   LPC overall efficiency
    - `etahc`:   HPC overall efficiency
    - `etaht`:   HPT overall efficiency
    - `etalt`:   LPT overall efficiency
    - `Lconv`:   `true` if convergence was successful, `false` otherwise

    The "?" symbol denotes the station index:
    - 0: freestream
    - 18: fan face outside of casing BLs
    - 19: fan face over LPC portion
    - 2: fan face over fan portion
    - 21: fan exit
    - 25: LPC exit, HPC inlet
    - 3: compressor exit
    - 4: combustor exit before cooling air addition
    - 41: turbine inlet after cooling air addition
    - 45: HPT exit, LPT inlet
    - 49: LPT exit
    - 5: core nozzle
    - 6: core flow downstream
    - 7: fan nozzle
    - 8: fan flow downstream
"""
function tfsize!(gee, M0, T0, p0, a0, M2, M25,
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

      n = 6

      # from 'airfrac.inc'
      # air fractions  
      #        N2      O2      CO2    H2O      Ar       fuel
      alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127, 0.0]

      # fuel fractions
      beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

      #---- fractional core mass flow convergence tolerance
      toler = 1.0e-12

      #---- mass offtake fraction update under-relaxation factor
      rlxfo = 0.8

      mcore = 0.0
      fo = 0.0
      Pom = 0.0

      #---- overall pressure ratio
      pic = pilc * pihc

      #---- number of air constitutents (all but fuel)
      nair = n - 1

      # ===============================================================
      #---- set combustion-change mass fractions gamma[i] for specified fuel
      gamma = gasfuel(ifuel, n)

      # Create buffer
      buf = Zygote.Buffer(gamma, length(gamma))
      for i = 1:length(gamma)
            buf[i] = gamma[i]
      end

      #---- apply combustor efficiency
      # Zygote.jl can not handle this...
      # for i = 1:nair
      #       gamma[i] = etab * gamma[i]
      # end
      # gamma[n] = 1.0 - etab

      # Zygote can handle this
      for i = 1:nair
            buf[i] = etab * buf[i]
      end
      buf[n] = 1.0 - etab

      gamma = copy(buf)
      #
      # ===============================================================
      #---- freestream static quantities
      s0, dsdt, h0, dhdt, cp0, R0 = gassum(alpha, nair, T0)
      gam0 = cp0 / (cp0 - R0)
      u0 = M0 * a0

      # ===============================================================
      #---- freestream total quantities
      hspec = h0 + 0.5 * u0^2
      Tguess = T0 * (1.0 + 0.5 * (gam0 - 1.0) * M0^2)
      Tt0 = gas_tset(alpha, nair, hspec, Tguess)

      st0, dsdt, ht0, dhdt, cpt0, Rt0 = gassum(alpha, nair, Tt0)
      pt0 = p0 * exp((st0 - s0) / Rt0)
      at0 = sqrt(Tt0 * Rt0 * cpt0 / (cpt0 - Rt0))

      # ===============================================================
      #---- offtake plume flow 9
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
      
      #---- initial guesses for station 2, 1.9 and 1.9c
      pt2 = pt18
      Tt2 = Tt18
      pt19 = pt18
      Tt19 = Tt18
      pt19c = pt19
      Tt19c = Tt19


      if (Kinl == 0.0 && mofft == 0.0 && Pofft == 0.0)
            #----- single design pass will be sufficient
            npass = 1
      else
            #----- must use multiple passes to converge pt2,pt19 from inlet defect Kinl,
            #-     and to converge on offtake fractions fo, Pom
            npass = 60
      end

      sbfan = 0.0
      sbcore = 0.0
      for ipass = 1:npass

            # ===============================================================
            #---- set fan inlet conditions corrected for BLI
            if (ipass == 1)
                  #c      if(mcore == 0.0)
                  #----- don't know engine mass flow yet, so ignore any BLI mixing
                  if (eng_has_BLI_cores)
                        sbfan = 0.0
                        sbcore = 0.0
                  else #clean flow
                        sbfan = 0.0
                        sbcore = 0.0
                  end

            else
                  #----- account for inlet BLI defect via mass-averaged entropy
                  a2sq = at0^2 / (1.0 + 0.5 * (gam0 - 1.0) * M2^2)

                  if eng_has_BLI_cores
                        #------ BL mixes with fan + core flow
                        mmix = BPR * mcore * sqrt(Tt2 / Tt0) * pt0 / pt2 +
                                    mcore * sqrt(Tt19c / Tt0) * pt0 / pt19c
                        sbfan2 = Kinl * gam0 / (mmix * a2sq)
                        sbcore2 = sbfan
                  else
                        #------ BL mixes with fan flow only
                        mmix = BPR * mcore * sqrt(Tt2 / Tt0) * pt0 / pt2
                        sbfan2 = Kinl * gam0 / (mmix * a2sq)
                        sbcore2 = 0.0
                  end

                  #----- update mixed-out entropies, with some underrelaxation       
                  rlxs = 0.85
                  sbfan = sbfan + rlxs * (sbfan2 - sbfan)
                  sbcore = sbcore + rlxs * (sbcore2 - sbcore)

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

            # ===============================================================
            #---- fan flow 2-7
            pifD = pif
            mbfD = 1.0
            mf = 1.0

            _, epf, _, _, _, _, _, _ = 
                  calculate_compressor_speed_and_efficiency(FanMap, pif, mf, pifD, mbfD, 1.0, epf0, Ng = 1.0, Rg = 2.0)

            pt21, Tt21, ht21, st21, cpt21, Rt21 = gas_prat(alpha, nair,
                  pt2, Tt2, ht2, st2, cpt2, Rt2, pif, epf)

            #---- fan duct nozzle total quantities
            pt7 = pt21 * pifn
            Tt7 = Tt21
            ht7 = ht21
            st7 = st21
            cpt7 = cpt21
            Rt7 = Rt21

            # ===============================================================
            #---- Compressor precooler 19 - 19c
            pt19c = pt19 - Δp_PreC
            ht19c = ht19 + Δh_PreC
            Tt19c = gas_tset(alpha, nair, ht19c, Tt19)
            st19c, _, ht19c, _, cpt19c, Rt19c = gassum(alpha, nair, Tt19c)


            #c      write(*,*) 'epf0 epf ', epf0, epf
            #
            # ===============================================================
            #---- LP compressor flow 19c - 25
            pilcD = pilc
            mblcD = 1.0
            ml = 1.0
            
            _, eplc, _, _, _, _, _, _ = 
                  calculate_compressor_speed_and_efficiency(LPCMap, pilc, ml, pilcD, mblcD, 1.0, eplc0, Ng = 1.0, Rg = 2.0)

            pt25, Tt25, ht25, st25, cpt25, Rt25 = gas_prat(alpha, nair,
                  pt19c, Tt19c, ht19c, st19c, cpt19c, Rt19c, pilc, eplc)

            # ===============================================================
            #---- Compressor intercooler 25 - 25c
            pt25c = pt25 - Δp_InterC
            ht25c = ht25 + Δh_InterC
            Tt25c = gas_tset(alpha, nair, ht25c, Tt25)
            st25c, _, ht25c, _, cpt25c, Rt25c = gassum(alpha, nair, Tt25c)

            # ===============================================================
            #---- HP compressor flow 25c - 3
            pihcD = pihc
            mbhcD = 1.0
            mh = 1.0

            _, ephc, _, _, _, _, _, _ = 
                  calculate_compressor_speed_and_efficiency(HPCMap, pihc, mh, pihcD, mbhcD, 1.0, ephc0, Ng = 1.0, Rg = 2.0)

            pt3, Tt3, ht3, st3, cpt3, Rt3 = gas_prat(alpha, nair,
                  pt25c, Tt25c, ht25c, st25c, cpt25c, Rt25c, pihc, ephc)

            # ===============================================================
            #---- combustor flow 3-4   (ffb = mdot_fuel/mdot_burner)
            ffb, lambda = gas_burn(alpha, beta, gamma, n, ifuel, Tt3, Ttf, Tt4, hvap)
            st4, dsdt, ht4, dhdt, cpt4, Rt4 = gassum(lambda, nair, Tt4)
            pt4 = pt3 * pib
            gam4 = cpt4 / (cpt4 - Rt4)

            # ===============================================================

            lambdap = zeros(nair)
            if compare_strings(opt_cooling, "none")
                  #----- no cooling air present... station 41 is same as 4
                  ff = ffb * (1.0 - fo)

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
                  #----------------------------------------------------------------
            else
                  #----- cooling air is present... calculate station 41

                  #----- hot-section temperature ratio for each blade row (for cooling model)
                  gmi4 = Rt4 / (cpt4 - Rt4)
                  Trrat = 1.0 / (1.0 + 0.5 * gmi4 * Mtexit^2)

                  # Heat exchanger to cool turbine cooling air
                  ht_tc = ht3 + Δh_TurbC #Specific enthalpy of turbine cooling air
                  Tt_tc = gas_tset(alpha, nair, ht_tc, Tt3) #Temperature of turbine cooling air

                  if compare_strings(opt_cooling, "fixed_coolingflowratio")
                        #------ epsrow(.) is assumed to be passed in.. calculate Tmrow(.)
                        Tmrow = Tmcalc(ncrowx, ncrow,
                        Tt_tc, Tt4, dTstrk, Trrat,
                              efilm, tfilm, StA, epsrow)
                  elseif compare_strings(opt_cooling, "fixed_Tmetal")
                        #------ calculate cooling mass flow ratios epsrow(.) to get specified Tmrow(.)
                        ncrow, epsrow, epsrow_Tt3, epsrow_Tt4, epsrow_Trr = mcool(ncrowx,
                              Tmrow, Tt_tc, Tt4, dTstrk, Trrat,
                              efilm, tfilm, StA)
                  end

                  #----- total cooling-flow fraction
                  fc = 0.0
                  for icrow = 1:ncrow
                        fc = fc + (1.0 - fo) * epsrow[icrow]
                  end

                  if (fc >= 0.99)
                        error("TFSIZE: Excessive cooling flow", 
                              "\n\tmcool/mcore = ", fc, 
                              "\n\tTt3 Tt4 Tmetal ", Tt3, " K, ", Tt4, " K, ", Tmrow[1], " K")
                  end

                  #----- set ff = mdot_fuel/mdot_core = ffb * mdot_burner/mdot_core
                  ff = (1.0 - fo - fc) * ffb

                  pt4a = pt4
                  Tt4a = Tt4
                  ht4a = ht4
                  st4a = st4
                  cpt4a = cpt4
                  Rt4a = Rt4

                  #----- speed at start-of-mixing station 4a
                  p4a, T4a, h4a, s4a, cp4a, R4a = gas_mach(lambda, nair,
                        pt4a, Tt4a, ht4a, st4a, cpt4a, Rt4a, 0.0, M4a, 1.0)
                  u4sq = max(2.0 * (ht4a - h4a), 0.0)
                  u4a = sqrt(u4sq)

                  #----- exit speed of cooling air at station 4a
                  uc = ruc * u4a

                  #----- IGV exit mixing
                  frac4 = (1.0 - fo - fc + ff) / (1.0 - fo + ff)
                  fracm = fc / (1.0 - fo + ff)

                  #----- mixed constituent fraction vector from mass equation
                  # for i = 1:nair
                  #       lambdap[i] = frac4 * lambda[i] + fracm * alpha[i]
                  # end

                  buf = Zygote.Buffer(lambdap, length(lambdap))
                  for i = 1:nair
                        buf[i] = frac4 * lambda[i] + fracm * alpha[i]
                  end

                  lambdap = copy(buf)

                  #----- mixed total enthalpy from enthalpy equation
                  ht41 = frac4 * ht4 + fracm * ht_tc

                  #----- total temperature from total enthalpy
                  Tguess = frac4 * Tt4 + fracm * Tt_tc
                  Tt41 = gas_tset(lambdap, nair, ht41, Tguess)

                  #----- all total quantities (except for total pressure), from total temperature
                  st41, dsdt, ht41, dhdt, cpt41, Rt41 = gassum(lambdap, nair, Tt41)

                  #----- mixed velocity from momentum equation, assuming constant static pressure
                  p41 = p4a
                  u41 = frac4 * u4a + fracm * uc

                  #----- static temperature from static enthalpy
                  h41 = ht41 - 0.5 * u41^2
                  Tguess = T4a + (h41 - h4a) / cp4a
                  T41 = gas_tset(lambdap, nair, h41, Tguess)

                  #----- all static quantities, from static temperature
                  s41, dsdt, h41, dhdt, cp41, R41 = gassum(lambdap, nair, T41)

                  #----- all stagnation quantities, from total-static enthalpy difference
                  dhb = ht41 - h41
                  epi = 1.0
                  pt41, Tt41, ht41, st41, cpt41, Rt41 = gas_delh(lambdap, nair,
                        p41, T41, h41, s41, cp41, R41, dhb, epi)

            end

            # ===============================================================
            #---- LPT and HPT work, per unit mass flow
            dhfac = -(1.0 - fo) / (1.0 - fo + ff) / (1.0 - epsh)
            dlfac = -1.0 / (1.0 - fo + ff) / (1.0 - epsl)

            dhht = (ht3 - ht25c) * dhfac
            dhlt = (ht25 - ht19c + BPR * (ht21 - ht2) + Pom) * dlfac

            #---- HPT flow
            #     Trh =  Tt41/(Tt41 + dhht/cpt41)
            #     gexh = cpt41/(Rt41*epht0)
            #     pihtD = Trh^gexh
            epht1 = epht0 #Assume same as design point prior to cooling

            #Find cooled HPT efficiency
            epht = find_cooled_hpt_efficiency(epht1, epht_fc, fc0, fc)
            
            epi = 1.0 / epht
            pt45, Tt45, ht45, st45, cpt45, Rt45 = gas_delh(lambdap, nair,
                  pt41, Tt41, ht41, st41, cpt41, Rt41, dhht, epi)


            eplt = eplt0 #Assume same as design turbine efficiency
            epi = 1.0 / eplt
            pt49, Tt49, ht49, st49, cpt49, Rt49 = gas_delh(lambdap, nair,
                  pt45, Tt45, ht45, st45, cpt45, Rt45, dhlt, epi)

            # ===============================================================
            #---- Regenerative cooling heat exchanger 49 - 49c
            pt49c = pt49 - Δp_Regen
            ht49c = ht49 + Δh_Regen

            Tt49c = gas_tset(lambdap, nair, ht49c, Tt49)
            st49c, _, ht49c, _, cpt49c, Rt49c = gassum(lambdap, nair, Tt49c)

            # ===============================================================
            #---- Turbine nozzle 49c - 5

            pt5 = pt49c * pitn
            Tt5 = Tt49c
            ht5 = ht49c
            st5 = st49c
            cpt5 = cpt49c
            Rt5 = Rt49c


            #
            # ===============================================================
            #---- fan plume flow 7-8, use alpha mass fraction (air)
            pt8 = pt7
            ht8 = ht7
            Tt8 = Tt7
            st8 = st7
            cpt8 = cpt7
            Rt8 = Rt7
            pratfn = p0 / pt8
            p8, T8, h8, s8, cp8, R8 = gas_prat(alpha, nair, pt8, Tt8, ht8, st8, cpt8, Rt8, pratfn, 1.0)
            if (h8 >= ht8)

                  Lconv = false

                  u8 = 0.001 * sqrt(R8 * T8)
                  error("TFSIZE: Negative fan plume velocity", 
                        "\n\tpt2,  Tt2  = ", pt2, " Pa, ",  Tt2, " K",
                        "\n\tpt8,  Tt8  = ", pt8, " Pa, ",  Tt8, " K", 
                        "\n\tp8,  T8  = "  , p8,  " Pa, ",  T8,  " K", 
                        "\n\tpif,  BPR  = ", pif, " Pa, ",  BPR)
            else
                  u8 = sqrt(2.0 * (ht8 - h8))
            end
            rho8 = p8 / (R8 * T8)

            # ===============================================================
            #---- core plume flow 5-6, use lambdap mass fraction (combustion products)
            pt6 = pt5
            ht6 = ht5
            Tt6 = Tt5
            st6 = st5
            cpt6 = cpt5
            Rt6 = Rt5
            prattn = p0 / pt6
            p6, T6, h6, s6, cp6, R6 = gas_prat(lambdap, nair, pt6, Tt6, ht6, st6, cpt6, Rt6, prattn, 1.0)
            if (h6 >= ht6)
                  
                  Lconv = false
                  u6 = 0.001 * sqrt(R6 * T6)
                  error("TFSIZE: Negative core plume velocity", 
                        "\n\tpt2,  Tt2  = ", pt2, " Pa, ",  Tt2, " K",
                        "\n\tpt3,  Tt3  = ", pt3, " Pa, ",  Tt3, " K", 
                        "\n\tpt4,  Tt4  = ", pt4, " Pa, ",  Tt4, " K", 
                      "\n\tpt41,  Tt41  = ", pt41," Pa, ",  Tt41," K", 
                        "\n\tpt6,  Tt6  = ", pt6, " Pa, ",  Tt6, " K", 
                        "\n\tp6,  T6  = "  , p6,  " Pa, ",  T6,  " K", 
                        "\n\tpif,  BPR  = ", pif, " Pa, ",  BPR)
            else
                  u6 = sqrt(2.0 * (ht6 - h6))
            end

            rho6 = p6 / (R6 * T6)

            #      write(*,*) 'Pt6   u6 ', pt6, u6
            #
            # ===============================================================
            #---- effective fuel heating value, over states 3, 4  (just for info)
            cpa = 0.5 * (cpt3 + cpt4)
            hfuel = cpa * (Tt4 - Tt3 + ffb * (Tt4 - Ttf)) / (etab * ffb)

            #---- effective fuel heating value, over states 3, 4.1  (just for info)
            #      cpa = 0.5*(cpt3+cpt41)
            #      hfuel = cpa*((1.0-fo)*(Tt41-Tt3) + ff*(Tt41-Ttf)) / (etab*ff)

            # ===============================================================
            #---- size core mass flow

            #---- store current values for better update, convergence checks
            mcold = mcore
            foold = fo

            #---- added effective net thrust from dissipation in ingested streamtube
            if (u0 == 0.0)
                  Finl = 0.0
            else
                  Finl = Phiinl / u0
            end

            #---- set core mass flow from specified effective net thrust
            mcore = (Feng - Finl) /
                    ((1.0 - fo + ff) * u6 - u0 + BPR * (u8 - u0) + fo * u9)

            #---- corresponding new offtake mass flow fraction
            fonew = mofft / mcore
            dfo = fonew - foold

            fo = fo + rlxfo * dfo

            #---- estimate better new mass flow, compensating for change in mass offtake
            mfac = min(2.0, 1.0 / (1.0 - dfo))
            mcore = mcore * mfac

            #---- power offtake per mass flow
            Pom = Pofft / mcore

            #---- overall Fsp and TSFC
            Fsp = Feng / (u0 * mcore * (1.0 + BPR))
            if (Feng <= 0.0)
                  TSFC = 0.0
            else
                  TSFC = (gee * ff * mcore) / Feng
            end

            # ===============================================================

            M8 = u8 / sqrt(cp8 * R8 / (cp8 - R8) * T8)
            if (M8 < 1.0)
                  #----- subsonic fan plume... fan nozzle flow is same as plume
                  p7 = p8
                  T7 = T8
                  h7 = h8
                  s7 = s8
                  cp7 = cp8
                  R7 = R8
                  u7 = u8
            else
                  #----- supersonic fan plume... fan nozzle is choked
                  M7 = 1.0
                  p7, T7, h7, s7, cp7, R7 = gas_mach(alpha, nair,
                        pt7, Tt7, ht7, st7, cpt7, Rt7, 0.0, M7, 1.0)
                  u7 = sqrt(2.0 * (ht7 - h7))
            end
            rho7 = p7 / (R7 * T7)

            #---- size fan  nozzle and plume areas
            A7 = BPR * mcore / (rho7 * u7)
            A8 = BPR * mcore / (rho8 * u8)

            # ===============================================================
            M6 = u6 / sqrt(cp6 * R6 / (cp6 - R6) * T6)
            #      write(*,*) 'u6,M6', u6,M6
            if (M6 < 1.0)
                  #----- subsonic core plume... core nozzle flow is same as plume
                  p5 = p6
                  T5 = T6
                  h5 = h6
                  s5 = s6
                  cp5 = cp6
                  R5 = R6
                  u5 = u6
            else
                  #----- supersonic core plume... core nozzle is choked
                  M5 = 1.0
                  p5, T5, h5, s5, cp5, R5 = gas_mach(lambdap, nair,
                        pt5, Tt5, ht5, st5, cpt5, Rt5, 0.0, M5, 1.0)
                  u5 = sqrt(2.0 * (ht5 - h5))
            end

            rho5 = p5 / (R5 * T5)

            #---- size core nozzle and plume areas
            A5 = (1.0 - fo + ff) * mcore / (rho5 * u5)
            A6 = (1.0 - fo + ff) * mcore / (rho6 * u6)

            if (u9 == 0.0)
                  A9 = 0.0
            else
                  A9 = fo * mcore / (rho9 * u9)
            end

            # ===============================================================
            #---- size fan and compressor areas
            p2, T2, h2, s2, cp2, R2 = gas_mach(alpha, nair,
                  pt2, Tt2, ht2, st2, cpt2, Rt2, 0.0, M2, 1.0)
            u2 = sqrt(2.0 * (ht2 - h2))
            rho2 = p2 / (R2 * T2)

            p19c, T19c, h19c, s19c, cp19c, R19c = gas_mach(alpha, nair,
                  pt19c, Tt19c, ht19c, st19c, cpt19c, Rt19c, 0.0, M2, 1.0)
            u19c = sqrt(2.0 * (ht19c - h19c))
            rho19c = p19c / (R19c * T19c)
            A2 = BPR * mcore / (rho2 * u2) + mcore / (rho19c * u19c)

            p25c, T25c, h25c, s25c, cp25c, R25c = gas_mach(alpha, nair,
                  pt25c, Tt25c, ht25c, st25c, cpt25c, Rt25c, 0.0, M25, 1.0)
            u25c = sqrt(2.0 * (ht25c - h25c))
            rho25c = p25c / (R25c * T25c)
            A25 = (1.0 - fo) * mcore / (rho25c * u25c)

            if (ipass >= 2) || (npass == 1)
                  dmfrac = 1.0 - mcold / mcore

                  if (abs(dmfrac) < toler) || (npass == 1)

                        # ===============================================================
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
                        etalt = (ht49 - ht45) / (ht49i - ht45)
                        
                        Lconv = true
                        return epsrow, Tmrow,
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
                        Lconv
                  end

            end


      end

      if (npass > 1)
            println("TFSIZE: Convergence failed.  dm/m = ", dmfrac)
      end

end # tfsize