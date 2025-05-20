"""
    ductedfanoper!(M0, T0, p0, a0, Tref, pref,
        Phiinl, Kinl, iBLIc,
        pid, pifn, 
        pifD, 
        mbfD, 
        A2, A7,
        epf0,
        Feng, Peng,
        M2, pif, mbf, 
        iPspec)

Ducted fan operation routine

    This code follows a similar strategy as that in `tfoper()`. However, as the problem has fewer unknowns, the 
    Jacobian is not computed manually and a black-box non-linear solver is used instead.

    The gas routines reside in the following source files:
     gascalc.f  Routines for various processes 
                (compressor, turbine, combustor, etc)
     gasfun.f   Routines for computing cp[T], h[t], sigma[T], R,
                called by the routines in gascalc.f

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `M0`:      freestream Mach
    - `T0`:      freestream temperature  [K]
    - `p0`:      freestream pressure  [Pa]
    - `Tref`:    reference temperature for corrected mass flow and speed
    - `pref`:    reference pressure for corrected mass flow
    - `Phiinl`:  inlet ingested dissipation  Phi_inl
    - `iBLIc`:   0=core in clear flow, 1=core sees Phiinl
    - `pid`:     diffuser pressure ratio  ( = pt2/pt0)
    - `pifn`:    fan     nozzle pressure ratio  ( = pt7/pt6.9)
    - `pifD`:    design fan pressure ratio  ( = pt21/pt2 )
    - `mbfD`:    design corrected fan mass flow ( = mf*sqrt(Tt2 /Tref)/(pt2 /pref) )
    - `A2`:      fan-face area [m^2]                mf = mc*BPR, mt = mc*(1+ff)
    - `A7`:      fan  nozzle area [m^2]
    - `epf0`:    max fan polytropic efficiency
    - `Feng`:    required net thrust [N]
    - `Peng`:    power required to drive fan [W]
    - `M2`:      guess Mach number at station 2
    - `pif`:     guess fan pressure ratio
    - `mbf`:     guess corrected fan mass flow rate
    - `iPspec`:  true if power is specified, false if thrust is specified

      **Output:**
    ------
    - `pif`:     fan pressure ratio
    - `mbf`:     corrected fan mass flow rate
    - `M2`:      Mach number at station 2

    The "?" symbol denotes the station index:
      0   freestream
      18  fan face outside of casing BLs
      2   fan face over fan portion
      7   fan nozzle
      8   fan flow downstream
"""
function ductedfanoper!(M0, T0, p0, a0, Tref, pref,
      Phiinl, Kinl, iBLIc,
      pid, pifn, 
      pifD, 
      mbfD,
      NbfD, 
      A2, A7,
      epf0,
      Feng, Peng,
      M2, pif, mbf, 
      Δh_radiator, Δp_radiator,
      iPspec)
    
    tol = 1e-10 #convergence tolerance

    pf = pif
    mf = mbf
    Mi = M2
    if (pf == 0.0)
        pf = pifD
    end
    if (mf == 0.0)
        mf = mbfD
    end
    if (Mi == 0.0)
        Mi = 0.6
    end
    if (Mi == 1.0)
        Mi = 0.6
    end

    guess = zeros(3)
    guess[1] = pf
    guess[2] = mf
    guess[3] = Mi

    #This function returns the residual of the non-linear engine problem. it
    #can also return the engine performance results.
    function DuctedFanOffDesign(x; iPspec = false, store_data = false)
        #Extract unknowns
        pf = x[1]
        mf = x[2]
        Mi = x[3]
        #Note that this function can access the scope of the outer function, so that 
        #engine variables do not need to be loaded again.

        # Constants
        alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127, 0.0]
        nair = 5

        #---- minimum allowable fan efficiency
        epfmin = 0.60

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

        # ===============================================================
        #---- diffuser flow 0-2
        Tt18 = Tt0
        st18 = st0
        ht18 = ht0
        cpt18 = cpt0
        Rt18 = Rt0
        pt18 = pt0 * pid

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
        # ===============================================================
        #---- Fan flow flow 2-21
        Tt18 = Tt0
        st18 = st0
        ht18 = ht0
        cpt18 = cpt0
        Rt18 = Rt0
        pt18 = pt0 * pid

        _, epf, _, _, _, _, _, _ = 
            calculate_compressor_speed_and_efficiency(FanMap, pf, mf, pifD, mbfD, 1.0, epf0)
  
        if (pf < 1.0)
                epf = 1.0 / epf
        end

        pt21, Tt21, ht21, st21, cpt21, Rt21 = gas_prat(alpha, nair,
                pt2, Tt2, ht2, st2, cpt2, Rt2, pf, epf)

        # ===============================================================
        #---- Radiator heat exchanger
        pt7 = pt21 * pifn - Δp_radiator
        ht7 = ht21 + Δh_radiator

        Tt7 = gas_tset(alpha, nair, ht7, Tt21)
        st7, _, ht7, _, cpt7, Rt7 = gassum(alpha, nair, Tt7)

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

        #Calculate thrust
        epi = 1.0
        p8, T8, h8, s8, cp8, R8 = gas_prat(alpha, nair,
                        pt7, Tt7, ht7, st7, cpt7, Rt7, pfn, epi)

        if (ht7 > h8)
            u8 = sqrt(2.0 * (ht7 - h8))
        else
            u8 = 0.0
        end
        M8 = u8 / sqrt(T8 * R8 * cp8 / (cp8 - R8))
        rho8 = p8 / (R8 * T8)

        #----- overall thrust and power
        if (u0 == 0.0)
            Finl = 0.0
        else
            Finl = Phiinl / u0
        end
        mfan = mf * sqrt(Tref / Tt2) * pt2 / pref #Fan mass flow rate
        F = mfan * (u7 - u0) + A7 * (p7 - p0) + Finl #Total thrust
        P = mfan * (ht21 - ht2) #Fan power

        A8 = mfan / (rho8 * u8) #Plume area

    # ===============================================================
        #---- #Set up residuals
        res = zeros(3)

        #Fan nozzle mass flow, choked or unchoked
        res[1] = (mfan - rho7 * A7 * u7)/mfan

        #Front mass flow
        res[2] = (mfan - rho2 * A2 * u2)/mfan

        if iPspec #Specified power constraint
            res[3] = (P - Peng)/Peng

        else #Specified thrust constraint
            res[3] = (F- Feng)/Feng
        end

        if store_data
            pif = pf
            mbf = mf
            M2 = Mi
            #---- fan isentropic efficiency
            pt21i, Tt21i, ht21i, st21i, cpt21i, Rt21i = gas_prat(alpha, nair,
            pt2, Tt2, ht2, st2, cpt2, Rt2, pif, 1.0)
            etaf = (ht21i - ht2) / (ht21 - ht2)

            Nbf, _, _ = Ncmap(pf, mf, pifD, mbfD, NbfD, Cmapf)

            #---- overall Fsp
            if (u0 == 0.0)
                Fsp = 0.0
            else
                Fsp = F / (u0 * mfan)
            end

            #---- thrust-specific energy consumption
            TSEC = P/F

            outputs = (res, TSEC, Fsp, F, P, mfan, 
                        pif, mbf, Nbf, 
                        Tt0, ht0, pt0, cpt0, Rt0,
                        Tt18, ht18, pt18, cpt18, Rt18,
                        Tt2, ht2, pt2, cpt2, Rt2,
                        Tt21, ht21, pt21, cpt21, Rt21,
                        Tt7, ht7, pt7, cpt7, Rt7,
                        u0, T2, u2, p2, cp2, R2, M2,
                        T7, u7, p7, cp7, R7, M7,
                        T8, u8, p8, cp8, R8, M8, A8,
                        epf, etaf)
        else 
            outputs = res
        end
        return outputs
    end

    residual(x) = DuctedFanOffDesign(x, iPspec = iPspec)
    sol = nlsolve(residual, guess, ftol = tol) #Use NLsolve.jl to solve for ducted fan state

    #Evaluate residual once more, storing parameters
    outputs = DuctedFanOffDesign(sol.zero, iPspec = iPspec, store_data = true)
    res = outputs[1]
    if maximum(abs.(res)) > tol #If infinity norm is above tolerance
        @warn "DUCTEDFANOPER: convergence failed, iPspec = $iPspec"
    end

    return outputs[2:end]
end