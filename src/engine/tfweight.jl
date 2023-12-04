"""
      tfweight(iengwgt, Gearf, OPR, BPR, mdotc, dfan, rSnace,
      dlcomp, neng, feadd, fpylon)

Engine weight estimation function using Giulia Pantalone, Drela, or Fitzgerald model.
      
!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `iengwgt`: Engine model index, Drela=0, Fitzgerald=1, and Pantalone>=3, 
    - `OPR`: Overall pressure ratio.
    - `BPR`: By-pass ratio.
    - `mdotc`: Engine core mass flow rate.
    - `dfan`: Fan diameter.
    - `rSnace`: 
    - `dlcomp`:
    - `neng`: Number of engines.
    - `feadd`: Fuel system weight ratio.
    - `fpylon`: Pylon weight fraction.

    **Output:**
    - `Weng`: Total engine weight.
    - `Wnac`: Nacelle weight. 
    - `Webare`: Bare engine weight.
    - `Snace1`: Nacelle area.
"""
function tfweight(iengwgt, Gearf, OPR, BPR, mdotc, dfan, rSnace,
    dlcomp, neng, feadd, fpylon)

    # include("constants.inc")

    if (iengwgt > 2)

        # Giulia Pantalone's weight model (iengwgt = 3 or 4)
        if (abs(Gearf - 1.0) < 0.001)

            # ungeared fan
            if (iengwgt == 3)

                # basic technology
                Wcore, Wfan, Wcomb, Wnozz, Wnace = ddct(mdotc / 0.4535, OPR, BPR)

            else

                # advanced technology
                Wcore, Wfan, Wcomb, Wnozz, Wnace = ddat(mdotc / 0.4535, OPR, BPR)

            end

        else

            # geared fan
            if (iengwgt == 3)

                # basic technology
                Wcore, Wfan, Wcomb, Wnozz, Wnace = gct(mdotc / 0.4535, OPR, BPR)

            else

                # advanced technology
                Wcore, Wfan, Wcomb, Wnozz, Wnace = gat(mdotc / 0.4535, OPR, BPR)

            end
        end

        # total engine weight
        Weng1 = (Wcore + Wfan + Wcomb + Wnace + Wnozz) / lb_N
        Wpylon = Weng1 * fpylon
        We1 = (Wcore + Wcomb + Wfan) / lb_N
        Weng = (Weng1 + Wpylon) * neng
        Webare = We1 * neng
        Wnac = neng * Wnace / lb_N

    else

        # Use Drela's or Fitzgerald's engine model
        if (iengwgt == 0)

            #  Drela's original weight model
            We1 = (mdotc / 45.35) *
                  (1684.5 +
                   17.7 * (OPR / 30.0) +
                   1662.2 * (BPR / 5.0)^1.2) / lb_N

        else

            if (abs(Gearf - 1.0) < 0.001)

                # Nate Fitzgerald's weight model, ungeared fan
                if (iengwgt == 1)

                    # basic technology
                    acon = 18.09 * BPR^2 + 476.9 * BPR + 701.3
                    bcon = 0.001077 * BPR^2 - 0.03716 * BPR + 1.190
                    ccon = -0.01058 * BPR + 0.3259

                else

                    # advanced technology
                    acon = 15.38 * BPR^2 + 401.1 * BPR + 631.5
                    bcon = 0.001057 * BPR^2 - 0.03693 * BPR + 1.171
                    ccon = -0.01022 * BPR + 0.2321
                end

            else

                # Nate Fitzgerald's weight model, geared fan
                if (iengwgt == 1)

                    # basic technology
                    acon = -0.6590 * BPR^2 + 292.8 * BPR + 1915.0
                    bcon = 0.00006784 * BPR^2 - 0.006488 * BPR + 1.061
                    ccon = -0.001969 * BPR + 0.07107

                else
                    # advanced technology
                    acon = -0.6204 * BPR^2 + 237.3 * BPR + 1702.0
                    bcon = 0.00005845 * BPR^2 - 0.005866 * BPR + 1.045
                    ccon = -0.001918 * BPR + 0.06765

                end
            end

            We1 = (acon / lb_N) * (mdotc / 45.35)^bcon * (OPR / 40.0)^ccon
        end

        # bare engine weight
        Webare = We1 * neng
        # nacelle weight model from NASA CR 151970
        Snace1 = rSnace * 0.25 * pi * dfan^2
        Ainlet = 0.4 * Snace1
        Acowl = 0.2 * Snace1
        Aexh = 0.4 * Snace1
        Acore = pi * dlcomp * (3.0 * dlcomp)
        Wnace1 = 4.45 * (Ainlet / 0.3048^2) * (2.5 + 0.0238 * dfan / 0.0254) +
                 4.45 * (Acowl / 0.3048^2) * 1.9 +
                 4.45 * (Aexh / 0.3048^2) * (2.5 + 0.0363 * dfan / 0.0254) +
                 4.45 * (Acore / 0.3048^2) * 1.9
        Wnac = Wnace1 * neng

        # engine accessories + fuel system
        Weadd = Webare * feadd

        # engine pylons
        Wpylon = (Webare + Weadd + Wnac) * fpylon

        # total engine weight
        Weng = Webare + Weadd + Wnac + Wpylon

    end


    return Weng, Wnac, Webare, Snace1
end


"""
      ddct(mdotc,OPR,BPR)

      Calculates engine component weights for Direct-drive turbofan with current technology
      using Pantalone model.
            
      !!! details "🔃 Inputs and Outputs"
      **Input:**
      - `mdotc:`:  Core mass flow rate.
      - `OPR:`:  Overall compression ratio.
      - `BPR:`:  By-pass ratio.
      **Output:**
      - `Wcore:`:  Core weight.
      - `Wfan:`:  Fan weight.
      - `Wcomb:`:  Combustor weight.
      - `Wnozz:`:  Nozzle weight.
      - `Wnace:`:  Nacelle weight.
"""
function ddct(mdotc, OPR, BPR)

    #     this file provides coefficients for the LS models
    include("dir_cur.inc")
    #     this file provides data for GP model for the core weight: 
    #       alfac, OPRc, BPRc, M25c, sfc, Pc
    include("crddc.inc")
    xstar = zeros(3)
    Xc = zeros(2000, 3)

    mdoti = mdotc * (1 + BPR) # inlet mass flow
    xstar[1] = OPR
    xstar[2] = BPR
    xstar[3] = mdotc

    #-----------------------------------------------------------
    #     Populate X vectors
    Xf = cubprm(OPR, BPR, mdoti)
    Xb = qdprm(OPR, BPR, mdotc)
    Xnz = qrtprm(OPR, BPR, mdotc)
    Xnc = cubprm(OPR, BPR, mdotc)

    #-----------------------------------------------------------
    #     Core weight  
    for i = 1:2000
        Xc[i, 1] = OPRc[i]
        Xc[i, 2] = BPRc[i]
        Xc[i, 3] = M25c[i]
    end

    Wcore = gppre(Xc, alfac, xstar, 2000, 3, Pc, sfc)
    #-----------------------------------------------------------
    #     Fan weight - cubic f(OPR,BPR,mdoti)
    Wfan = 0
    for i = 1:20
        Wfan = Wfan + Xf[i] * Afan[i]
    end
    #-----------------------------------------------------------
    #     Combustor weight - quadratic f(OPR,BPR,mdotc)
    Wcomb = 0
    for i = 1:10
        Wcomb = Wcomb + Xb[i] * Acomb[i]
    end
    #-----------------------------------------------------------
    #     Nozzle weight - quartic f(OPR,BPR,mdotc)
    Wnozz = 0
    for i = 1:35
        Wnozz = Wnozz + Xnz[i] * Anozz[i]
    end
    #-----------------------------------------------------------
    #     Nacelle weight - cubic f(OPR,BPR,mdotc)
    Wnace = 0
    for i = 1:20
        Wnace = Wnace + Xnc[i] * Anace[i]
    end

    return Wcore, Wfan, Wcomb, Wnozz, Wnace
end

"""
      ddat(mdotc,OPR,BPR)
      
      Calculates engine component weights for Direct-drive turbofan with advanced technology 
      using Pantalone's model.

      !!! details "🔃 Inputs and Outputs"
      **Input:**
      - `mdotc:`:  Core mass flow rate.
      - `OPR:`:  Overall compression ratio.
      - `BPR:`:  By-pass ratio.
      **Output:**
      - `Wcore:`:  Core weight.
      - `Wfan:`:  Fan weight.
      - `Wcomb:`:  Combustor weight.
      - `Wnozz:`:  Nozzle weight.
      - `Wnace:`:  Nacelle weight.
"""
function ddat(mdotc, OPR, BPR)

    #     this file provides coefficients for the LS models
    include("dir_adv.inc")
    #     this file provides data for GP model for the core weight: 
    #       alfac, OPRc, BPRc, M25c, sfc, Pc
    include("crdda.inc")

    xstar = zeros(3)
    Xc = zeros(2000, 3)

    mdoti = mdotc * (1 + BPR) # inlet mass flow
    xstar[1] = OPR
    xstar[2] = BPR
    xstar[3] = mdotc

    #-----------------------------------------------------------
    #     Populate X vectors
    Xf = cubprm(OPR, BPR, mdoti)
    Xb = qdprm(OPR, BPR, mdotc)
    Xnz = qrtprm(OPR, BPR, mdotc)
    Xnc = cubprm(OPR, BPR, mdotc)

    #-----------------------------------------------------------
    #     Core weight  
    for i = 1:2000
        Xc[i, 1] = OPRc[i]
        Xc[i, 2] = BPRc[i]
        Xc[i, 3] = M25c[i]
    end

    Wcore = gppre(Xc, alfac, xstar, 2000, 3, Pc, sfc)
    #-----------------------------------------------------------
    #     Fan weight - cubic f(OPR,BPR,mdoti)
    Wfan = 0
    for i = 1:20
        Wfan = Wfan + Xf[i] * Afan[i]
    end
    #-----------------------------------------------------------
    #     Combustor weight - quadratic f(OPR,BPR,mdotc)
    Wcomb = 0
    for i = 1:10
        Wcomb = Wcomb + Xb[i] * Acomb[i]
    end
    #-----------------------------------------------------------
    #     Nozzle weight - quartic f(OPR,BPR,mdotc)
    Wnozz = 0
    for i = 1:35
        Wnozz = Wnozz + Xnz[i] * Anozz[i]
    end
    #-----------------------------------------------------------
    #     Nacelle weight - cubic f(OPR,BPR,mdotc)
    Wnace = 0
    for i = 1:20
        Wnace = Wnace + Xnc[i] * Anace[i]
    end

    return Wcore, Wfan, Wcomb, Wnozz, Wnace
end # ddat

"""
      gct(mdotc,OPR,BPR)

      Calculates engine component weights for geared turbofan with current technology 
      using Pantalone's model.

      !!! details "🔃 Inputs and Outputs"
      **Input:**
      - `mdotc: Core mass flow rate.
      - `OPR: Overall compression ratio.
      - `BPR: By-pass ratio.
      **Output:**
      - `Wcore: Core weight.
      - `Wfan: Fan weight.
      - `Wcomb: Combustor weight.
      - `Wnozz: Nozzle weight.
      - `Wnace: Nacelle weight.
"""
function gct(mdotc, OPR, BPR)

    #     this file provides coefficients for the LS models
    include("g_cur.inc")
    #     these files provide data for GP model for the core and nozzlec      weights respectively: 
    #       alfac, OPRc, BPRc, M25c, sfc, Pc
    #       alfan, OPRn, BPRn, M25n, sfn, Pn
    include("nzgc.inc")
    include("crgc.inc")

    xstar = zeros(3)
    Xc = zeros(2000, 3)
    Xnz = zeros(2000, 3)

    mdoti = mdotc * (1 + BPR) # inlet mass flow
    xstar[1] = OPR
    xstar[2] = BPR
    xstar[3] = mdotc

    #     Populate X vectors
    Xf = cubprm(OPR, BPR, mdoti)
    Xb = cubprm(OPR, BPR, mdotc)
    Xnc = cubprm(OPR, BPR, mdotc)

    #-----------------------------------------------------------
    #     Core weight  
    for i = 1:2000
        Xc[i, 1] = OPRc[i]
        Xc[i, 2] = BPRc[i]
        Xc[i, 3] = M25c[i]
    end

    Wcore = gppre(Xc, alfac, xstar, 2000, 3, Pc, sfc)
    #-----------------------------------------------------------
    #     Nozzle weight  
    for i = 1:2000
        Xnz[i, 1] = OPRn[i]
        Xnz[i, 2] = BPRn[i]
        Xnz[i, 3] = M25n[i]
    end

    Wnozz = gppre(Xnz, alfan, xstar, 2000, 3, Pn, sfn)
    #-----------------------------------------------------------
    #     Fan weight - cubic f(OPR,BPR,mdoti)
    Wfan = 0
    for i = 1:20
        Wfan = Wfan + Xf[i] * Afan[i]
    end
    #-----------------------------------------------------------
    #     Combustor weight - cubic f(OPR,BPR,mdotc)
    Wcomb = 0
    for i = 1:20
        Wcomb = Wcomb + Xb[i] * Acomb[i]
    end
    #-----------------------------------------------------------
    #     Nacelle weight - cubic f(OPR,BPR,mdotc)
    Wnace = 0
    for i = 1:20
        Wnace = Wnace + Xnc[i] * Anace[i]
    end

    return Wcore, Wfan, Wcomb, Wnozz, Wnace
end # gct

"""
      gat(mdotc,OPR,BPR)

      Calculates engine component weights for geared turbofan with advanced technology
      using Pantalone's model.

      !!! details "🔃 Inputs and Outputs"
      **Input:**
      - `mdotc:`:  Core mass flow rate.
      - `OPR:`:  Overall compression ratio.
      - `BPR:`:  By-pass ratio.
      **Output:**
      - `Wcore:`:  Core weight.
      - `Wfan:`:  Fan weight.
      - `Wcomb:`:  Combustor weight.
      - `Wnozz:`:  Nozzle weight.
      - `Wnace:`:  Nacelle weight.
"""
function gat(mdotc, OPR, BPR)

    #     this file provides coefficients for the LS models
    include("g_adv.inc")
    #     these files provide data for GP model for the core and nozzlec      weights respectively: 
    #       alfac, OPRc, BPRc, M25c, sfc, Pc
    #       alfan, OPRn, BPRn, M25n, sfn, Pn
    include("nzga.inc")
    include("crga.inc")

    xstar = zeros(3)
    Xc = zeros(2000, 3)
    Xnz = zeros(2000, 3)


    mdoti = mdotc * (1 + BPR) # inlet mass flow
    xstar[1] = OPR
    xstar[2] = BPR
    xstar[3] = mdotc
    #-----------------------------------------------------------
    #     Populate X vectors
    Xf = cubprm(OPR, BPR, mdoti)
    Xb = cubprm(OPR, BPR, mdotc)
    Xnc = cubprm(OPR, BPR, mdotc)

    #-----------------------------------------------------------
    #     Core weight  
    for i = 1:2000
        Xc[i, 1] = OPRc[i]
        Xc[i, 2] = BPRc[i]
        Xc[i, 3] = M25c[i]
    end

    Wcore = gppre(Xc, alfac, xstar, 2000, 3, Pc, sfc)
    #-----------------------------------------------------------
    #     Nozzle weight  
    for i = 1:2000
        Xnz[i, 1] = OPRn[i]
        Xnz[i, 2] = BPRn[i]
        Xnz[i, 3] = M25n[i]
    end

    Wnozz = gppre(Xnz, alfan, xstar, 2000, 3, Pn, sfn)
    #-----------------------------------------------------------
    #     Fan weight - cubic f(OPR,BPR,mdoti)
    Wfan = 0
    for i = 1:20
        Wfan = Wfan + Xf[i] * Afan[i]
    end
    #-----------------------------------------------------------
    #     Combustor weight - cubic f(OPR,BPR,mdotc)
    Wcomb = 0
    for i = 1:20
        Wcomb = Wcomb + Xb[i] * Acomb[i]
    end
    #-----------------------------------------------------------
    #     Nacelle weight - cubic f(OPR,BPR,mdotc)
    Wnace = 0
    for i = 1:20
        Wnace = Wnace + Xnc[i] * Anace[i]
    end

    return Wcore, Wfan, Wcomb, Wnozz, Wnace
end

"""
      X = qdprm(x1,x2,x3)      

      Populates vector for quadratic fit.
"""
function qdprm(x1, x2, x3)

    X = zeros(10)

    X[1] = 1
    X[2] = x1
    X[3] = x2
    X[4] = x3
    X[5] = x1^2
    X[6] = x2^2
    X[7] = x3^2
    X[8] = x1 * x2
    X[9] = x2 * x3
    X[10] = x1 * x3

    return X
end

"""
      cubprm(x1,x2,x3)

      Populates vector for cubic fit.
"""
function cubprm(x1, x2, x3)

    X = zeros(20)

    X[1] = 1
    X[2] = x1
    X[3] = x2
    X[4] = x3
    X[5] = x1^2
    X[6] = x2^2
    X[7] = x3^2
    X[8] = x1 * x2
    X[9] = x2 * x3
    X[10] = x1 * x3
    X[11] = x1^3
    X[12] = x2^3
    X[13] = x3^3
    X[14] = x1^2 * x2
    X[15] = x1^2 * x3
    X[16] = x2^2 * x1
    X[17] = x2^2 * x3
    X[18] = x3^2 * x1
    X[19] = x3^2 * x2
    X[20] = x1 * x2 * x3

    return X
end

"""
      qrtprm(x1,x2,x3)

      Populates vector for quartic fit.
"""

function qrtprm(x1, x2, x3)

    X = zeros(35)

    X[1] = 1
    X[2] = x1
    X[3] = x2
    X[4] = x3
    X[5] = x1^2
    X[6] = x2^2
    X[7] = x3^2
    X[8] = x1 * x2
    X[9] = x2 * x3
    X[10] = x1 * x3
    X[11] = x1^3
    X[12] = x2^3
    X[13] = x3^3
    X[14] = x1^2 * x2
    X[15] = x1^2 * x3
    X[16] = x2^2 * x1
    X[17] = x2^2 * x3
    X[18] = x3^2 * x1
    X[19] = x3^2 * x2
    X[20] = x1 * x2 * x3
    X[21] = x1^4
    X[22] = x2^4
    X[23] = x3^4
    X[24] = x1^3 * x2
    X[25] = x1^3 * x3
    X[26] = x2^3 * x1
    X[27] = x2^3 * x3
    X[28] = x3^3 * x1
    X[29] = x3^3 * x2
    X[30] = x1^2 * x2^2
    X[31] = x2^2 * x3^2
    X[32] = x1^2 * x3^2
    X[33] = x1^2 * x2 * x3
    X[34] = x2^2 * x1 * x3
    X[35] = x3^2 * x1 * x2

    return X
end
