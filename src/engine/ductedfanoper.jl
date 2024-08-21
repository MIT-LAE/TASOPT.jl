mutable struct DuctedFanData
    Fsp::Float64
    Pfan::Float64
    Feng::Float64
    Peng::Float64
    mfan::Float64
    pif::Float64
    mbf::Float64
    Nbf::Float64
    Tt0::Float64
    ht0::Float64
    pt0::Float64
    cpt0::Float64
    Rt0::Float64
    Tt18::Float64
    ht18::Float64
    pt18::Float64
    cpt18::Float64
    Rt18::Float64
    Tt2::Float64
    ht2::Float64
    pt2::Float64
    cpt2::Float64
    Rt2::Float64
    Tt21::Float64
    ht21::Float64
    pt21::Float64
    cpt21::Float64
    Rt21::Float64
    Tt7::Float64
    ht7::Float64
    pt7::Float64
    cpt7::Float64
    Rt7::Float64
    u0::Float64
    T2::Float64
    u2::Float64
    p2::Float64
    cp2::Float64
    R2::Float64
    M2::Float64
    A2::Float64
    T7::Float64
    u7::Float64
    p7::Float64
    cp7::Float64
    R7::Float64
    M7::Float64
    A7::Float64
    T8::Float64
    u8::Float64
    p8::Float64
    cp8::Float64
    R8::Float64
    M8::Float64
    A8::Float64
    epf::Float64
    etaf::Float64
    M0::Float64
    T0::Float64
    p0::Float64
    a0::Float64
    Tref::Float64
    pref::Float64
    Phiinl::Float64
    Kinl::Float64
    iBLIc::Float64
    pid::Float64
    pifn::Float64
    pifD::Float64
    mbfD::Float64
    epf0::Float64
    pifK::Float64
    epfK::Float64
    DuctedFanData() = new() 
end
"""
    ductedfanoper!(M0, T0, p0, a0, Tref, pref,
        Phiinl, Kinl, iBLIc,
        pid, pifn, 
        pifD, 
        mbfD, 
        A2, A7,
        epf0,
        pifK, epfK,
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
    - `pifK`:    fan efficiency FPR offset:    epolf = epf0 + epfK*(pif-pifK)
    - `epfK`:    fan efficiency pif derivative
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
      A2, A7,
      epf0,
      pifK, epfK,
      Feng, Peng,
      M2, pif, mbf, 
      iPspec)
    
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
    engdata = DuctedFanData()
    inputs = (
        M0, T0, p0, a0, Tref, pref, Phiinl, Kinl, iBLIc, pid, pifn, 
        pifD, mbfD, A2, A7, epf0, pifK, epfK, Feng, Peng
    ) 
    update_engine_data!(engdata, inputs)

    guess = zeros(3)
    guess[1] = pf
    guess[2] = mf
    guess[3] = Mi

    residual(x) = res_df(x, engdata, iPspec = iPspec)
    sol = nlsolve(residual, guess, ftol = 1e-8) 

    return sol.zero
end

function unpack_input_data(data::DuctedFanData)
    M0 = data.M0
    T0 = data.T0
    p0 = data.p0
    a0 = data.a0
    Tref = data.Tref
    pref = data.pref
    Phiinl = data.Phiinl
    Kinl = data.Kinl
    iBLIc = data.Phiinl
    pid = data.pid
    pifn = data.pifn
    pifD = data.pifD
    mbfD = data.mbfD
    A2 = data.A2
    A7 = data.A7
    epf0 = data.epf0
    pifK = data.pifK
    epfK = data.epfK
    Feng = data.Feng
    Peng = data.Peng
    
    return (M0, T0, p0, a0, Tref, pref, Phiinl, Kinl, iBLIc, pid, pifn, pifD, mbfD, A2, A7, epf0, pifK, epfK, Feng, Peng)
end

function update_engine_data!(data::DuctedFanData, inputs)
    (
        M0, T0, p0, a0, Tref, pref, Phiinl, Kinl, iBLIc, pid, pifn, 
        pifD, mbfD, A2, A7, epf0, pifK, epfK, Feng, Peng
    ) = inputs

    data.M0 = M0
    data.T0 = T0
    data.p0 = p0
    data.a0 = a0
    data.Tref = Tref
    data.pref = pref
    data.Phiinl = Phiinl
    data.Kinl = Kinl
    data.iBLIc = iBLIc
    data.pid = pid
    data.pifn = pifn
    data.pifD = pifD
    data.mbfD = mbfD
    data.A2 = A2
    data.A7 = A7
    data.epf0 = epf0
    data.pifK = pifK
    data.epfK = epfK
    data.Feng = Feng
    data.Peng = Peng
end


function res_df(x, engdata; iPspec = false, store_data = false)
    #Extract unknowns
    pf = x[1]
    mf = x[2]
    Mi = x[3]

    #Extract Inputs
    (M0, T0, p0, a0, Tref, pref, Phiinl, Kinl, iBLIc, pid, pifn, pifD, mbfD, A2, A7, epf0, pifK, epfK, Fspec, Pspec) =
        unpack_input_data(engdata)

    # Constants
    alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127, 0.0]
    nair = 5
    #        a     b     k     mo     da    c    d     C    D
    Cmapf = [3.50, 0.80, 0.03, 0.95, -0.50, 3.0, 6.0, 0.0, 0.0]
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

    epf, _, _ = ecmap(pf, mf, pifD, mbfD, Cmapf, epf0, pifK, epfK)

    if (epf < epfmin)
            epf = epfmin

    end
    if (pf < 1.0)
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

    mfan = mf * sqrt(Tref / Tt2) * pt2 / pref #Fan mass flow ref
    P = mfan * (ht21 - ht2) #Fan power
    #Fan nozzle mass flow, choked or unchoked
    res[1] = mfan - rho7 * A7 * u7

    #Inlet Mach area constraint
    mfA = mf * sqrt(Tref / Tt2) * pt2 / pref / (rho2 * u2)
    res[2] = mfA - A2

    if iPspec #Specified power constraint
       
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
