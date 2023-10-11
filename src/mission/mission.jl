"""
Runs aircraft through mission, calculating fuel burn
and other mission variables.

Input:
 pari[.]   integer flags
 parg[.]   geometry parameters
 parm[.]   mission parameters
 iairf     index of airfoil database to use
 initeng    0 = engine state will be initialized for all points
            1 = engine state is assumed to be initialized
 ipc1       0 = ipcruise1 aero and engine point needs to be calculated
            1 = ipcruise1 aero and engine point assumed calculated

Input/Output:
 para[.p]  aero     parameters for points p=1..iptotal
 pare[.p]  engine   parameters for points p=1..iptotal


NOTE: 
 This routine assumes that estimates of the climb-leg flight path 
 gamma angles are passed in via para[iagamV,ipclimb1:ipclimbn].
 These appear as cos(gamma) factors in the climb equations,
 and can be passed in as zero with only a minor error.
 They are updated and returned in the same para[iagamV,ip] array.
"""
function mission!(pari, parg, parm, para, pare, Ldebug, NPSS_PT, NPSS, ipc1)#, iairf, initeng, ipc1)

      t_prop = 0.0
      calc_ipc1 = true
      ifirst = true

      if pari[iiengmodel] == 0
            # Drela engine model
            use_NPSS = false
      else
            # NPSS
            use_NPSS = true
      end

      # HACK TODO add the para back
      # iairf
      initeng = 0
      # ipc1 = 0
      ipc1 = 0 # HACK

      itergmax::Int64 = 15
      gamVtol = 1.0e-12
      #     gamVtol  = 1.0e-10

      # unpack flags
      iengloc = pari[iiengloc]
      ifclose = pari[iifclose]
      ifuel = pari[iifuel]

      if use_NPSS
            mofWpay = parg[igmofWpay]
            mofWMTO = parg[igmofWMTO]
            PofWpay = parg[igPofWpay]
            PofWMTO = parg[igPofWMTO]


            # BLI
            fBLIf = parg[igfBLIf]
            DAfsurf = para[iaDAfsurf, ipcruise1]
            KAfTE = para[iaKAfTE, ipcruise1]
      end

      # Mission range
      Rangetot = parm[imRange]
      para[iaRange, ipdescentn] = Rangetot

      # MTOW
      WMTO = parg[igWMTO]

      # Payload fraction for this mission
      rpay = parm[imWpay] / parg[igWpay]
      ξpay = 0.0

      # Zero-fuel weight for this mission
      Wzero = WMTO - parg[igWfuel] - parg[igWpay] + parm[imWpay]

      # mission TO weight
      WTO = parm[imWTO]

      # mission TO fuel weight
      WfTO = WTO - Wzero

      # set known operating conditions

      # takeoff altitude conditions
      ip = ipstatic
      altkm = para[iaalt, ip] / 1000.0
      T_std, p_std, ρ_std, a_std, μ_std = atmos(altkm)
      T0 = parm[imT0TO]
      p0 = p_std
      ρ0 = ρ_std * (T_std / T0)
      a0 = a_std * sqrt(T0 / T_std)
      μ0 = μ_std * (T0 / T_std)^0.8

      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = ρ0
      pare[iemu0, ip] = μ0
      pare[ieM0, ip] = 0.0
      pare[ieu0, ip] = 0.0

      para[iaMach, ip] = 0.0
      para[iaCL, ip] = 0.0
      para[iaReunit, ip] = 0.0
      para[iaWbuoy, ip] = 0.0

      # ip ∈ iprotate:ipclimb1
      pare[iep0, iprotate:ipclimb1] .= p0
      pare[ieT0, iprotate:ipclimb1] .= T0
      pare[iea0, iprotate:ipclimb1] .= a0
      pare[ierho0, iprotate:ipclimb1] .= ρ0
      pare[iemu0, iprotate:ipclimb1] .= μ0
      para[iaWbuoy, iprotate:ipclimb1] .= 0.0


      # Start-of-cruise altitude conditions
      ip = ipcruise1
      Mach = para[iaMach, ip]
      altkm = para[iaalt, ip] / 1000.0

      T0, p0, rho0, a0, mu0 = atmos(altkm)
      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = rho0
      pare[iemu0, ip] = mu0
      pare[ieM0, ip] = Mach
      pare[ieu0, ip] = Mach * a0
      para[iaReunit, ip] = Mach * a0 * rho0 / mu0

      # End-of-descent altitude conditions
      ip = ipdescentn
      altkm = para[iaalt, ip] / 1000.0
      T_std, p_std, rho_std, a_std, mu_std = atmos(altkm)
      T0 = parm[imT0TO]
      p0 = p_std
      rho0 = rho_std * (T_std / T0)
      a0 = a_std * sqrt(T0 / T_std)
      mu0 = mu_std * (T0 / T_std)^0.8
      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = rho0
      pare[iemu0, ip] = mu0
      para[iaWbuoy, ip] = 0.0

      # interpolate CL over climb points, 
      #  between specified ipclimb1+1, ipclimbn values
      CLa = para[iaCL, ipclimb1+1]
      CLb = para[iaCL, ipcruise1]

      @inbounds for ip = ipclimb1+1:ipclimbn
            frac = float(ip - (ipclimb1 + 1)) /
                   float(ipclimbn - (ipclimb1 + 1))
            para[iaCL, ip] = CLa * (1.0 - frac^2) +
                             CLb * frac^2
      end

      #---- interpolate CL over descent points, 
      #-     between specified ipdescent1, ipdescentn-1 values
      #      CLd = para[iaCL,ipcruisen]
      #      CLe = para[iaCL,ipcruisen]
      CLd = para[iaCL, ipcruisen] * 0.96
      CLe = para[iaCL, ipcruisen] * 0.50

      @inbounds for ip = ipdescent1:ipdescentn-1
            frac = float(ip - ipdescent1) /
                   float((ipdescentn - 1) - ipdescent1)
            fb = 1.0 - frac
            para[iaCL, ip] = CLd * fb^2 +
                             CLe * (1.0 - fb^2)
      end


      #---- estimate takeoff speed and set V,Re over climb and descent
      cosL = cos(parg[igsweep] * π / 180.0)
      CLTO = para[iaclpmax, iptakeoff] * cosL^2
      # [prash] I think the assumption here is that Wcruise/WTO~1 and
      # just scaling it with rho and CL to get an initial estimate for V
      VTO = pare[ieu0, ipcruise1] *
            sqrt(pare[ierho0, ipcruise1] / pare[ierho0, iptakeoff]) *
            sqrt(para[iaCL, ipcruise1] / CLTO)
      ReTO = VTO * pare[ierho0, iptakeoff] / pare[iemu0, iptakeoff]
      pare[ieu0, iprotate:ipclimb1] .= VTO
      para[iaReunit, iprotate:ipclimb1] .= ReTO
      @inbounds for ip = ipclimb1+1:ipclimbn
            frac = float(ip - ipclimb1) / float(ipclimbn - ipclimb1)
            V = VTO * (1.0 - frac) + pare[ieu0, ipcruise1] * frac
            Re = ReTO * (1.0 - frac) + para[iaReunit, ipcruise1] * frac
            pare[ieu0, ip] = V
            para[iaReunit, ip] = Re
      end
      @inbounds for ip = ipdescent1:ipdescentn
            frac = float(ip - ipdescent1) / float(ipdescentn - ipdescent1)
            V = VTO * frac + pare[ieu0, ipcruisen] * (1.0 - frac)
            Re = ReTO * frac + para[iaReunit, ipcruisen] * (1.0 - frac)
            pare[ieu0, ip] = V
            para[iaReunit, ip] = Re
      end



      #---- takeoff CLmax via section clmax and sweep correction
      ip = iprotate
      clpmax = para[iaclpmax, ip]
      sweep = parg[igsweep]
      cosL = cos(sweep * π / 180.0)
      CLmax = clpmax * cosL^2

      #---- Vs stall speed (takeoff condition)
      rho0 = pare[ierho0, ip]
      a0 = pare[iea0, ip]
      S = parg[igS]
      Vstall = sqrt(2.0 * WTO / (rho0 * S * CLmax))
      Mstall = Vstall / a0
      pare[ieu0, ip] = Vstall
      pare[ieM0, ip] = Mstall
      para[iaMach, ip] = Mstall
      para[iaReunit, ip] = Vstall * pare[ierho0, ip] / pare[iemu0, ip]

      #---- V2 speed per FAR-25  (takeoff,cutback,climb1 condition)
      V2 = Vstall * 1.2
      M2 = Mstall * 1.2
      CL2 = CLmax / 1.2^2

      ip = iptakeoff
      pare[ieu0, ip] = V2
      pare[ieM0, ip] = M2
      para[iaMach, ip] = M2
      para[iaReunit, ip] = V2 * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2

      #---- set pitch trim by adjusting CLh
      Wf = WTO - Wzero
      rfuel = Wf / parg[igWfuel]
      itrim = 1
      balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)

      CLh2 = para[iaCLh, ip]
      xCG2 = para[iaxCG, ip]
      xCP2 = para[iaxCP, ip]
      xNP2 = para[iaxNP, ip]

      para[iaxCG, 1:iprotate] .= xCG2
      para[iaxCP, 1:iprotate] .= xCP2
      para[iaxNP, 1:iprotate] .= xNP2
      para[iaCLh, 1:iprotate] .= 0.0

      #---- initial guesses for climb points
      ip = ipcutback
      pare[ieu0, ip] = V2
      pare[ieM0, ip] = M2
      para[iaMach, ip] = M2
      para[iaReunit, ip] = V2 * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2
      para[iaxCG, ip] = xCG2
      para[iaxCP, ip] = xCP2
      para[iaxNP, ip] = xNP2
      para[iaCLh, ip] = CLh2

      ip = ipclimb1
      pare[ieu0, ip] = V2
      pare[ieM0, ip] = M2
      para[iaMach, ip] = M2
      para[iaReunit, ip] = V2 * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2
      para[iaxCG, ip] = xCG2
      para[iaxCP, ip] = xCP2
      para[iaxNP, ip] = xNP2
      para[iaCLh, ip] = CLh2

      #---- also use V2 speed for end-of-descent condition, with weight correction
      ip = ipdescentn
      Vrat = sqrt(para[iafracW, ip] / para[iafracW, ipclimb1])
      pare[ieu0, ip] = V2 * Vrat
      pare[ieM0, ip] = M2 * Vrat
      para[iaMach, ip] = M2 * Vrat
      para[iaReunit, ip] = V2 * Vrat * pare[ierho0, ip] / pare[iemu0, ip]
      para[iaCL, ip] = CL2
      #
      # ============================================================================
      #---- set up climb points at equal altitude intervals, from altb to altc
      #-    (takeoff ground temperature is neglected here -- std atmosphere is used)
      altb = para[iaalt, iptakeoff]
      altc = para[iaalt, ipcruise1]
      @inbounds for ip = ipclimb1+1:ipclimbn
            frac = float(ip - ipclimb1) / float(ipclimbn - ipclimb1)
            para[iaalt, ip] = altb * (1.0 - frac) + altc * frac

            altkm = para[iaalt, ip] / 1000.0

            T0, p0, rho0, a0, mu0 = atmos(altkm)
            pare[iep0, ip] = p0
            pare[ieT0, ip] = T0
            pare[iea0, ip] = a0
            pare[ierho0, ip] = rho0
            pare[iemu0, ip] = mu0

            rhocab = max(parg[igpcabin], p0) / (RSL * TSL)
            para[iaWbuoy, ip] = (rhocab - rho0) * gee * parg[igcabVol]
      end

      #---- set climb Tt4's from fractions
      fT1 = parg[igfTt4CL1]
      fTn = parg[igfTt4CLn]
      Tt4TO = pare[ieTt4, iptakeoff]
      Tt4CR = pare[ieTt4, ipcruise1]
      @inbounds for ip = ipclimb1:ipclimbn
            frac = float(ip - ipclimb1) /
                   float(ipclimbn - ipclimb1)
            Tfrac = fT1 * (1.0 - frac) + fTn * frac
            pare[ieTt4, ip] = Tt4TO * (1.0 - Tfrac) + Tt4CR * Tfrac
      end
      #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #---- initial values for range, time, weight fraction
      para[iaRange, iprotate:ipclimb1] .= 0.0
      para[iatime, iprotate:ipclimb1] .= 0.0
      para[iafracW, iprotate:ipclimb1] .= WTO / WMTO
      para[iaWbuoy, iprotate:ipclimb1] .= 0.0

      # Initialize climb integrands
      FoW = zeros(Float64, iptotal)
      FFC = zeros(Float64, iptotal)
      Vgi = zeros(Float64, iptotal)
      if use_NPSS
            dPdt = zeros(Float64, iptotal)
            dygdt = zeros(Float64, iptotal)


            pare[iePLH2, ipclimb1] = 1.5
            pare[ieyg, ipclimb1] = 0.1
      end



      # integrate trajectory over climb
      @inbounds for ip = ipclimb1:ipclimbn
            if (Ldebug)
                  printstyled("Climb angle integration - ip = ", ip - ipclimb1 + 1, "\n"; color=:light_green)
            end

            # velocity calculation from CL, Weight, altitude
            W = para[iafracW, ip] * WMTO
            CL = para[iaCL, ip]
            ρ = pare[ierho0, ip]
            μ = pare[iemu0, ip]
            Vsound = pare[iea0, ip]
            cosg = cos(para[iagamV, ip])
            BW = W + para[iaWbuoy, ip]
            Wpay = parg[igWpay]
            neng = parg[igneng]

            if use_NPSS

                  if pari[iiengtype] == 1
                        initfanPCT = 80.0
                        finalfanPCT = parg[igfanPCT]
                        frac = float(ip - ipclimb1) / float(ipclimbn - ipclimb1)
                        fanPCT = initfanPCT * (1.0 - frac) + finalfanPCT * frac

                        mofft = (mofWpay * Wpay + mofWMTO * W) / neng
                        Pofft = (PofWpay * Wpay + PofWMTO * W) / neng
                  end

                  dgamV = 1.0
                  Ftotal = 0.0
                  DoL = 0.0
                  mdotf = 0.0
                  V = 0.0
                  if (Ldebug)
                        printstyled(@sprintf("\t%5s  %10s  %10s  %10s  %10s  %10s  %10s  %10s \n",
                                    "iterg", "dgamV", "gamV", "cosg", "BW", "Ftotal", "DoL", "V"); color=:light_green)
                  end
            end

            @inbounds for iterg = 1:itergmax
                  V = sqrt(2.0 * BW * cosg / (ρ * S * CL))
                  Mach = V / Vsound

                  para[iaMach, ip] = Mach
                  para[iaReunit, ip] = V * ρ / μ

                  pare[ieu0, ip] = V
                  pare[ieM0, ip] = Mach

                  # Set pitch trim by adjusting CLh
                  Wf = W - Wzero
                  rfuel = Wf / parg[igWfuel]
                  itrim = 1
                  balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)

                  if (ip == ipclimb1)
                        icdfun = 0 #use explicitly specified wing cdf, cdp
                  else
                        icdfun = 1 #use airfoil database
                  end
                  cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)


                  if use_NPSS
                        ρ0 = pare[ierho0, ip]
                        u0 = pare[ieu0, ip]
                        Φinl = 0.5 * ρ0 * u0^3 * (DAfsurf * fBLIf) / 2.0
                        Kinl = 0.5 * ρ0 * u0^3 * (KAfTE * fBLIf) / 2.0 # Assume 2 engines

                        if pari[iiengtype] == 0
                              NPSS_success, Ftotal, η, P, Hrej, heatexcess,
                              mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, para[iaalt, ip], Mach, 0.0, pare[ieTt4, ip], Kinl, Φinl, 0.0, 0.0, ifirst, parg, parpt, pare, ip)

                              pare[iedeNOx, ip] = deNOx
                              pare[ieEINOx1, ip] = EINOx1
                              pare[ieEINOx2, ip] = EINOx2
                              pare[ieemot:ieethermal, ip] .= η
                              pare[ieHrejmot:ieHrejtot, ip] .= Hrej
                              pare[ieHexcess, ip] = heatexcess
                        else
                              NPSS_success, Ftotal, heatexcess,
                              mdotf, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, para[iaalt, ip], Mach, 0.0, pare[ieTt4, ip], ifirst, parg, parpt, pare, ip, fanPCT, mofft, Pofft)

                              pare[ieEINOx1, ip] = EINOx1

                              pare[ieemot:ieethermal, ip] .= 0.0
                              pare[ieHrejmot:ieHrejtot, ip] .= 0.0
                              pare[ieHexcess, ip] = heatexcess
                        end

                        ifirst = false
                        pare[ieOPR, ip] = OPR
                        pare[ieTt3, ip] = Tt3
                        pare[ieWc3, ip] = Wc3
                        pare[iemdotf, ip] = mdotf

                        DoL = para[iaCD, ip] / para[iaCL, ip]
                  else
                        icall = 1
                        icool = 1

                        ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, initeng)

                        Ftotal = pare[ieFe, ip] * parg[igneng]
                        TSFC = pare[ieTSFC, ip]
                        DoL = para[iaCD, ip] / para[iaCL, ip]

                  end

                  # Calculate improved flight angle
                  ϕ = Ftotal / BW
                  sing = (ϕ - DoL * sqrt(1.0 - ϕ^2 + DoL^2)) / (1.0 + DoL^2)
                  gamV = asin(sing)
                  cosg = sqrt(1.0 - sing^2)
                  # gamV = atan(sing, cosg)

                  dgamV = gamV - para[iagamV, ip]

                  para[iagamV, ip] = gamV
                  para[iaROC, ip] = sing * V * 60 / ft_to_m #ft per min
                  if (Ldebug)
                        printstyled(@sprintf("\t%5d  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e\n",
                                    iterg, abs(dgamV), gamV * 180 / π, cosg, BW, Ftotal, DoL, V); color=:light_green)
                  end

                  if (abs(dgamV) < gamVtol)
                        break
                  end
            end

            if (abs(dgamV) > gamVtol)
                  println("Climb gamV not converged")
            end

            pare[ieFe, ip] = Ftotal
            # Store integrands for range and weight integration using a predictor-corrector scheme
            FoW[ip] = Ftotal / (BW * cosg) - DoL
            if (!use_NPSS)
                  FFC[ip] = Ftotal * TSFC / (W * V * cosg)
            else
                  FFC[ip] = mdotf * gee / (W * V * cosg)
            end
            Vgi[ip] = 1.0 / (V * cosg)

            if (use_NPSS)
                  ΔT = pare[ieT0, ip] - 20.0
                  P = pare[iePLH2, ip]
                  yg = pare[ieyg, ip]
                  dygdt[ip] = mdotf / ρmix(yg, P)
                  dPdt[ip] = dPdt_LH2(ΔT, mdotf, P, yg)
            end

            Mach = para[iaMach, ip]
            CL = para[iaCL, ip]
            CD = para[iaCD, ip]
            gamV = para[iagamV, ip]

            if (ip > ipclimb1)
                  # Corrector step
                  dh = para[iaalt, ip] - para[iaalt, ip-1]
                  dVsq = pare[ieu0, ip]^2 - pare[ieu0, ip-1]^2

                  FoWavg = 0.5 * (FoW[ip] + FoW[ip-1])
                  FFCavg = 0.5 * (FFC[ip] + FFC[ip-1])
                  Vgiavg = 0.5 * (Vgi[ip] + Vgi[ip-1])

                  if (use_NPSS)
                        dPdtavg = 0.5 * (dPdt[ip] + dPdt[ip-1])
                        dygdtavg = 0.5 * (dygdt[ip] + dygdt[ip-1])
                  end

                  dR = (dh + 0.5 * dVsq / gee) / FoWavg
                  dt = dR * Vgiavg
                  rW = exp(-dR * FFCavg)    #Ratio of weights Wᵢ₊₁/Wᵢ

                  if (use_NPSS)
                        dP = dPdtavg * Vgiavg * dR
                        dyg = dygdtavg * Vgiavg * dR
                  end

                  para[iaRange, ip] = para[iaRange, ip-1] + dR
                  para[iatime, ip] = para[iatime, ip-1] + dt
                  para[iafracW, ip] = para[iafracW, ip-1] * rW

                  if (use_NPSS)
                        pare[iePLH2, ip] = pare[iePLH2, ip-1] + dP
                        pare[ieyg, ip] = pare[ieyg, ip-1] + dyg
                  end

                  ifirst = false

            end
            if (ip < ipclimbn)
                  # Predictor integration step, forward Euler
                  if (para[iagamV, ip+1] ≤ 0.0)
                        # if gamV guess is not passed in, use previous point as the guess
                        para[iagamV, ip+1] = para[iagamV, ip]
                  end

                  W = para[iafracW, ip+1] * WMTO   # Initial weight fractions have been set in wsize.jl
                  CL = para[iaCL, ip+1]
                  ρ = pare[ierho0, ip+1]
                  cosg = cos(para[iagamV, ip+1])

                  BW = W + para[iaWbuoy, ip+1]

                  V = sqrt(2 * BW * cosg / (ρ * S * CL))
                  pare[ieu0, ip+1] = V
                  dh = para[iaalt, ip+1] - para[iaalt, ip]
                  dVsq = pare[ieu0, ip+1]^2 - pare[ieu0, ip]^2

                  dR = (dh + 0.5 * dVsq / gee) / FoW[ip]
                  dt = dR * Vgi[ip]
                  rW = exp(-dR * FFC[ip])

                  if (use_NPSS)
                        dP = dPdt[ip] * Vgi[ip] * dR
                        dyg = dygdt[ip] * Vgi[ip] * dR
                  end

                  para[iaRange, ip+1] = para[iaRange, ip] + dR
                  para[iatime, ip+1] = para[iatime, ip] + dt
                  para[iafracW, ip+1] = para[iafracW, ip] * rW

                  if (use_NPSS)
                        pare[iePLH2, ip+1] = pare[iePLH2, ip] + dP
                        pare[ieyg, ip+1] = pare[ieyg, ip] + dyg
                  end
            end

      end # done integrating climb

      # First cruise point is last climb point
      para[iaRange, ipcruise1] = para[iaRange, ipclimbn]
      para[iatime, ipcruise1] = para[iatime, ipclimbn]
      para[iafracW, ipcruise1] = para[iafracW, ipclimbn]
      para[iaWbuoy, ipcruise1] = para[iaWbuoy, ipclimbn]

      pare[iePLH2, ipcruise1] = pare[iePLH2, ipclimbn]
      pare[ieyg, ipcruise1] = pare[ieyg, ipclimbn]

      # Cruise start
      ip = ipcruise1
      # Set pitch trim by adjusting CLh
      Wf = para[iafracW, ip] * WMTO - Wzero
      rfuel = Wf / parg[igWfuel]
      itrim = 1
      balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)

      # if (calc_ipc1)
      if (ipc1 == 0)
            # println("Calculating cruise point")
            # println(pare[ieFe, ip])
            # Calculate only if requested since for design mission start of cruise is the des point and ∴ already calcualted 
            # Calculate drag
            icdfun = 1
            cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)
            DoL = para[iaCD, ip] / para[iaCL, ip]
            W = para[iafracW, ip] * WMTO
            BW = W + para[iaWbuoy, ip]
            F = BW * (DoL + para[iagamV, ip])
            Wpay = parg[igWpay]

            if use_NPSS
                  if pari[iiengtype] == 1
                        mofWpay = parg[igmofWpay]
                        mofWMTO = parg[igmofWMTO]
                        PofWpay = parg[igPofWpay]
                        PofWMTO = parg[igPofWMTO]
                        neng = parg[igneng]

                        mofft = (mofWpay * Wpay + mofWMTO * W) / neng
                        Pofft = (PofWpay * Wpay + PofWMTO * W) / neng
                  end

                  Mach = pare[ieM0, ip]
                  # Run powertrain [TODO] actually should run this to a required thrust not Tt4
                  ρ0 = pare[ierho0, ip]
                  u0 = pare[ieu0, ip]
                  Φinl = 0.5 * ρ0 * u0^3 * (DAfsurf * fBLIf) / 2.0
                  Kinl = 0.5 * ρ0 * u0^3 * (KAfTE * fBLIf) / 2.0 # Assume 2 engines

                  if pari[iiengtype] == 0
                        NPSS_success, Ftotal, η, P, Hrej, heatexcess,
                        mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, para[iaalt, ip], Mach, F, 0 * pare[ieTt4, ip], Kinl, Φinl, 0.0, 0.0, ifirst, parg, parpt, pare, ip)

                        pare[ieOPR, ip] = OPR
                        pare[ieTt3, ip] = Tt3
                        pare[ieWc3, ip] = Wc3
                        pare[iedeNOx, ip] = deNOx
                        pare[ieEINOx1, ip] = EINOx1
                        pare[ieEINOx2, ip] = EINOx2
                        pare[iemdotf, ip] = mdotf
                        pare[ieemot:ieethermal, ip] .= η
                        pare[ieHrejmot:ieHrejtot, ip] .= Hrej
                        pare[ieHexcess, ip] = heatexcess
                        pare[ieFe, ip] = Ftotal

                  else
                        NPSS_success, Ftotal, heatexcess,
                        mdotf, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, para[iaalt, ip], Mach, F, pare[ieTt4, ip], ifirst, parg, parpt, pare, ip, 0.0, mofft, Pofft)

                        ifirst = false

                        pare[ieOPR, ip] = OPR
                        pare[ieTt3, ip] = Tt3
                        pare[ieWc3, ip] = Wc3
                        pare[ieEINOx1, ip] = EINOx1
                        pare[iemdotf, ip] = mdotf
                        pare[ieemot:ieethermal, ip] .= 0.0
                        pare[ieHrejmot:ieHrejtot, ip] .= 0.0
                        pare[ieHexcess, ip] = heatexcess
                  end
            else
                  icall = 2
                  icool = 1
                  ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, initeng)
            end
      end

      # set cruise-climb climb angle, from fuel burn rate and atmospheric dp/dz
      if (!use_NPSS)
            TSFC = pare[ieTSFC, ip]
      end
      V = pare[ieu0, ip]
      p0 = pare[iep0, ip]
      ρ0 = pare[ierho0, ip]
      DoL = para[iaCD, ip] / para[iaCL, ip]
      W = para[iafracW, ip] * WMTO
      BW = W + para[iaWbuoy, ip]

      # Calculate Cruise-climb angle:
      # ------
      # γ ≈ dh/dR = dh/dp × dp/dW × dW/dR
      #           = -1/(ρg) × p/W × -̇mf g/(Vcosγ)
      #        ∴γ ≈ ̇mf*p/(ρWV)
      # For conventional aircraft can represent this as function of TSFC: (see TASOPT docs Eq 438)
      # gamVcr1 = DoL*p0*TSFC/(ρ0*gee*V - p0*TSFC)
      if use_NPSS
            gamVcr1 = mdotf * p0 / (ρ0 * W * V) # TSFC not the most useful for turbo-electric systems. Buoyancy weight has been neglected.
      else
            gamVcr1 = DoL * p0 * TSFC / (ρ0 * gee * V - p0 * TSFC)
      end
      para[iagamV, ip] = gamVcr1

      if use_NPSS
            Ftotal = pare[ieFe, ip]
      else
            Ftotal = pare[ieFe, ip] * parg[igneng]
      end
      cosg = cos(gamVcr1)

      FoW[ip] = Ftotal / (BW * cosg) - DoL
      if use_NPSS
            FFC[ip] = mdotf * gee / (W * V * cosg)
      else
            FFC[ip] = Ftotal * TSFC / (W * V * cosg)
      end
      Vgi[ip] = 1.0 / (V * cosg)

      if use_NPSS
            ΔT = pare[ieT0, ip] - 20.0
            P = pare[iePLH2, ip]
            yg = pare[ieyg, ip]
            dygdt[ip] = mdotf / ρmix(yg, P)
            dPdt[ip] = dPdt_LH2(ΔT, mdotf, P, yg)
      end

      #---- set end-of-cruise point "d" using cruise and descent angles, 
      #-     and remaining range legs
      gamVde1 = parm[imgamVDE1]
      gamVden = parm[imgamVDEn]
      gamVdeb = 0.5 * (gamVde1 + gamVden)
      alte = para[iaalt, ipdescentn]
      dRclimb = para[iaRange, ipclimbn] - para[iaRange, ipclimb1]
      dRcruise = (alte - altc - gamVdeb * (Rangetot - dRclimb)) / (gamVcr1 - gamVdeb)

      altd = altc + gamVcr1 * dRcruise
      # println("Cruise Fe inside mission", pare[ieFe, ipcruise1])

      # Final cruise point
      ip = ipcruisen
      Mach = para[iaMach, ip]
      altkm = altd / 1000.0

      T0, p0, ρ0, a0, μ0 = atmos(altkm)
      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = ρ0
      pare[iemu0, ip] = μ0
      pare[ieM0, ip] = Mach
      pare[ieu0, ip] = Mach * a0
      para[iaReunit, ip] = Mach * a0 * ρ0 / μ0
      para[iaalt, ip] = altd

      ρcab = max(parg[igpcabin], p0) / (RSL * TSL)
      para[iaWbuoy, ip] = (ρcab - ρ0) * gee * parg[igcabVol]

      # Set pitch trim by adjusting CLh
      Wf = para[iafracW, ip] * WMTO - Wzero
      rfuel = Wf / parg[igWfuel]
      itrim = 1
      balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)

      # Calc Drag
      icdfun = 1
      cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)
      DoL = para[iaCD, ip] / para[iaCL, ip]
      W = para[iafracW, ip] * WMTO
      BW = W + para[iaWbuoy, ip]
      Ftotal = BW * (DoL + para[iagamV, ip])
      if !use_NPSS
            pare[ieFe, ip] = Ftotal / parg[igneng]
      end

      if use_NPSS
            Mach = pare[ieM0, ip]
            # Run powertrain [TODO] actually should run this to a required thrust not Tt4
            ρ0 = pare[ierho0, ip]
            u0 = pare[ieu0, ip]
            Φinl = 0.5 * ρ0 * u0^3 * (DAfsurf * fBLIf) / 2.0
            Kinl = 0.5 * ρ0 * u0^3 * (KAfTE * fBLIf) / 2.0 # Assume 2 engines

            if pari[iiengtype] == 0
                  NPSS_success, Ftotal, η, P, Hrej, heatexcess,
                  mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, para[iaalt, ip], Mach, F, pare[ieTt4, ip], Kinl, Φinl, 0.0, 0.0, ifirst, parg, parpt, pare, ip)

                  pare[ieOPR, ip] = OPR
                  pare[ieTt3, ip] = Tt3
                  pare[ieWc3, ip] = Wc3
                  pare[iedeNOx, ip] = deNOx
                  pare[ieEINOx1, ip] = EINOx1
                  pare[ieEINOx2, ip] = EINOx2
                  pare[iemdotf, ip] = mdotf
                  pare[ieemot:ieethermal, ip] .= η
                  pare[ieHrejmot:ieHrejtot, ip] .= Hrej
                  pare[ieHexcess, ip] = heatexcess
            else
                  NPSS_success, Ftotal, heatexcess,
                  mdotf, EINOx1, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, para[iaalt, ip], Mach, F, pare[ieTt4, ip], ifirst, parg, parpt, pare, ip, 0.0, mofft, Pofft)

                  pare[ieOPR, ip] = OPR
                  pare[ieTt3, ip] = Tt3
                  pare[ieWc3, ip] = Wc3
                  pare[ieEINOx1, ip] = EINOx1
                  pare[iemdotf, ip] = mdotf
                  pare[ieemot:ieethermal, ip] .= 0.0
                  pare[ieHrejmot:ieHrejtot, ip] .= 0.0
                  pare[ieHexcess, ip] = heatexcess
            end
      else
            icall = 2
            icool = 1

            ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, initeng)

      end

      if (!use_NPSS)
            TSFC = pare[ieTSFC, ip]
      end
      V = pare[ieu0, ip]
      p0 = pare[iep0, ip]
      ρ0 = pare[ierho0, ip]
      DoL = para[iaCD, ip] / para[iaCL, ip]
      if use_NPSS
            gamVcr2 = mdotf * p0 / (ρ0 * W * V)
      else
            gamVcr2 = DoL * p0 * TSFC / (ρ0 * gee * V - p0 * TSFC)
      end
      para[iagamV, ip] = gamVcr2

      cosg = cos(gamVcr1)

      FoW[ip] = Ftotal / (BW * cosg) - DoL
      if (!use_NPSS)
            FFC[ip] = Ftotal * TSFC / (W * V * cosg)
      else
            FFC[ip] = mdotf * gee / (W * V * cosg)
      end
      Vgi[ip] = 1.0 / (V * cosg)

      if use_NPSS
            ΔT = pare[ieT0, ip] - 20.0
            P = pare[iePLH2, ip]
            dygdt[ip] = mdotf / ρmix(yg, P)
            dPdt[ip] = dPdt_LH2(ΔT, mdotf, P, yg)
      end

      ip1 = ipcruise1
      ipn = ipcruisen

      FoWavg = 0.5 * (FoW[ipn] + FoW[ip1])
      FFCavg = 0.5 * (FFC[ipn] + FFC[ip1])
      Vgiavg = 0.5 * (Vgi[ipn] + Vgi[ip1])

      if use_NPSS
            dPdtavg = 0.5 * (dPdt[ipn] + dPdt[ip1])
            dygdtavg = 0.5 * (dygdt[ipn] + dygdt[ip1])
            dP = dPdtavg * Vgiavg * dRcruise
            dyg = dygdtavg * Vgiavg * dRcruise
      end

      dtcruise = dRcruise * Vgiavg
      rWcruise = exp(-dRcruise * FFCavg)

      para[iaRange, ipn] = para[iaRange, ip1] + dRcruise
      para[iatime, ipn] = para[iatime, ip1] + dtcruise
      para[iafracW, ipn] = para[iafracW, ip1] * rWcruise

      if use_NPSS
            pare[iePLH2, ipn] = pare[iePLH2, ip1] + dP
            pare[ieyg, ipn] = pare[ieyg, ip1] + dyg
      end

      # set intermediate points over cruise, if any, just by interpolating
      for ip = ipcruise1+1:ipcruisen-1
            frac = float(ip - ipcruise1) / float(ipcruisen - ipcruise1)
            Mach = para[iaMach, ip]
            para[iaalt, ip] = altc * (1.0 - frac) + altd * frac
            altkm = para[iaalt, ip] / 1000.0
            T0, p0, rho0, a0, mu0 = atmos(altkm)
            pare[iep0, ip] = p0
            pare[ieT0, ip] = T0
            pare[iea0, ip] = a0
            pare[ierho0, ip] = rho0
            pare[iemu0, ip] = mu0
            pare[ieM0, ip] = Mach
            pare[ieu0, ip] = Mach * a0
            para[iaReunit, ip] = Mach * a0 * rho0 / mu0

            rhocab = max(parg[igpcabin], p0) / (RSL * TSL)
            para[iaWbuoy, ip] = (rhocab - rho0) * gee * parg[igcabVol]

            para[iaRange, ip] = para[iaRange, ipcruise1] + dRcruise * frac
            para[iatime, ip] = para[iatime, ipcruise1] + dtcruise * frac
            para[iafracW, ip] = para[iafracW, ipcruise1] * rWcruise^frac

      end

      # Descent
      #TODO descent has not been properly implemented yet - below it just performs an approximate scaling based
      #     on TASOPT to estimate the descent fuel burn
      ip = ipdescent1
      pare[iep0, ip] = pare[iep0, ipcruisen]
      pare[ieT0, ip] = pare[ieT0, ipcruisen]
      pare[iea0, ip] = pare[iea0, ipcruisen]
      pare[ierho0, ip] = pare[ierho0, ipcruisen]
      pare[iemu0, ip] = pare[iemu0, ipcruisen]
      pare[ieM0, ip] = pare[ieM0, ipcruisen]
      pare[ieu0, ip] = pare[ieu0, ipcruisen]

      para[iaMach, ip] = para[iaMach, ipcruisen]
      para[iaReunit, ip] = para[iaReunit, ipcruisen]
      para[iaalt, ip] = para[iaalt, ipcruisen]

      para[iaRange, ip] = para[iaRange, ipcruisen]
      para[iatime, ip] = para[iatime, ipcruisen]
      para[iafracW, ip] = para[iafracW, ipcruisen]
      para[iaWbuoy, ip] = para[iaWbuoy, ipcruisen]


      if use_NPSS
            # Mission fuel fractions and weights
            fracWa = para[iafracW, ipclimb1]
            # fracWe = para[iafracW, ipdescentn] 
            fracWe = para[iafracW, ipdescent1] * (1 - (1 - 0.993) * ((para[iaalt, ipdescent1] / ft_to_m) / 39858.0))# Pp temp set to 95% of ipdescent1 instead of ipdescentn since descent has not been calculated yet 
            freserve = parg[igfreserve]
            fburn = fracWa - fracWe            # Burnt fuel fraction
            ffuel = fburn * (1.0 + freserve)
            Wfuel = WMTO * ffuel
            WTO = Wzero + Wfuel

            parm[imWTO] = WTO
            parm[imWfuel] = Wfuel
            # printstyled("Wfuel = $Wfuel \n fburn = $fburn \n", color=:red)

            Wburn = WMTO * fburn
            parg[igWfburn] = Wburn
            parm[imPFEI] = Wburn / gee * parg[igLHVfuel] * 1e6 / (parm[imWpay] * parm[imRange])
      end

      if (!use_NPSS)
            Rd = para[iaRange, ipdescent1]
            Re = para[iaRange, ipdescentn]
            altd = para[iaalt, ipdescent1]
            alte = para[iaalt, ipdescentn]
            para[iagamV, ipdescent1] = gamVde1
            para[iagamV, ipdescentn] = gamVden

            for ip = ipdescent1+1:ipdescentn-1
                  frac = float(ip - ipdescent1) / float(ipdescentn - ipdescent1)
                  R = Rd * (1.0 - frac) + Re * frac

                  alt = altd + (Re - Rd) * (gamVde1 * (frac - 0.5 * frac^2) + gamVden * 0.5 * frac^2)
                  gamVde = gamVde1 * (1.0 - frac) + gamVden * frac
                  para[iagamV, ip] = gamVde

                  altkm = alt / 1000.0
                  T0, p0, ρ0, a0, μ0 = atmos(altkm)
                  pare[iep0, ip] = p0
                  pare[ieT0, ip] = T0
                  pare[iea0, ip] = a0
                  pare[ierho0, ip] = ρ0
                  pare[iemu0, ip] = μ0

                  para[iaRange, ip] = R
                  para[iaalt, ip] = alt

                  rhocab = max(parg[igpcabin], p0) / (RSL * TSL)
                  para[iaWbuoy, ip] = (rhocab - rho0) * gee * parg[igcabVol]
            end
            para[iaWbuoy, ipdescentn] = 0.0

            # integrate time and weight over descent
            for ip = ipdescent1:ipdescentn

                  # velocity calculation from CL, Weight, altitude
                  gamVde = para[iagamV, ip]
                  cosg = cos(gamVde)
                  W = para[iafracW, ip] * WMTO
                  BW = W + para[iaWbuoy, ip]
                  CL = para[iaCL, ip]
                  rho = pare[ierho0, ip]
                  V = sqrt(2.0 * BW * cosg / (rho * S * CL))
                  Mach = V / pare[iea0, ip]

                  para[iaMach, ip] = Mach
                  para[iaReunit, ip] = V * rho / pare[iemu0, ip]

                  pare[ieu0, ip] = V
                  pare[ieM0, ip] = Mach

                  # set pitch trim by adjusting CLh
                  Wf = W - Wzero
                  rfuel = Wf / parg[igWfuel]
                  itrim = 1
                  balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)

                  if (ip == ipdescentn)
                        # use explicitly specified wing cdf,cdp
                        icdfun = 0
                  else
                        # use airfoil database for wing cdf,cdp
                        icdfun = 1
                  end
                  cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)

                  # set up for engine calculation
                  sing = sin(gamVde)
                  cosg = cos(gamVde)
                  DoL = para[iaCD, ip] / para[iaCL, ip]
                  Fspec = BW * (sing + cosg * DoL)
                  pare[ieFe, ip] = Fspec / parg[igneng]

                  if (initeng == 0)
                        pare[iembf, ip] = pare[iembf, ip-1]
                        pare[iemblc, ip] = pare[iemblc, ip-1]
                        pare[iembhc, ip] = pare[iembhc, ip-1]
                        pare[iepif, ip] = pare[iepif, ip-1]
                        pare[iepilc, ip] = pare[iepilc, ip-1]
                        pare[iepihc, ip] = pare[iepihc, ip-1]
                        inite = 1


                        # make better estimate for new Tt4, adjusted for new ambient T0
                        dTburn = pare[ieTt4, ip-1] - pare[ieTt3, ip-1]
                        OTR = pare[ieTt3, ip-1] / pare[ieTt2, ip-1]
                        Tt3 = pare[ieT0, ip] * OTR
                        pare[ieTt4, ip] = Tt3 + dTburn + 50.0

                        # make better estimate for new pt5, adjusted for new ambient p0
                        pare[iept5, ip] = pare[iept5, ip-1] * pare[iep0, ip] / pare[iep0, ip-1]

                  else
                        inite = initeng
                  end

                  # use fixed engine geometry, specified Fe
                  icall = 2
                  # use previously-set turbine cooling mass flow
                  icool = 1

                  ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, inite)

                  # store effective thrust, effective TSFC
                  F = pare[ieFe, ip] * parg[igneng]
                  TSFC = pare[ieTSFC, ip]

                  # store integrands for Range and Weight integration
                  FoW[ip] = F / (BW * cosg) - DoL
                  FFC[ip] = F / (W * V * cosg) * TSFC
                  Vgi[ip] = 1.0 / (V * cosg)

                  # if F < 0, then TSFC is not valid, so calculate mdot_fuel directly
                  mfuel = pare[ieff, ip] * pare[iemcore, ip] * parg[igneng]
                  FFC[ip] = gee * mfuel / (W * cosg * V)


                  if (ip > ipdescent1)
                        #  corrector integration step, approximate trapezoidal
                        dh = para[iaalt, ip] - para[iaalt, ip-1]
                        dVsq = pare[ieu0, ip]^2 - pare[ieu0, ip-1]^2

                        FoWavg = 0.5 * (FoW[ip] + FoW[ip-1])
                        FFCavg = 0.5 * (FFC[ip] + FFC[ip-1])
                        Vgiavg = 0.5 * (Vgi[ip] + Vgi[ip-1])

                        dR = para[iaRange, ip] - para[iaRange, ip-1]
                        dt = dR * Vgiavg
                        rW = exp(-dR * FFCavg)

                        para[iatime, ip] = para[iatime, ip-1] + dt
                        para[iafracW, ip] = para[iafracW, ip-1] * rW
                  end

                  if (ip < ipdescentn)

                        # predictor integration step, forward Euler
                        gamVde = para[iagamV, ip+1]
                        cosg = cos(gamVde)
                        W = para[iafracW, ip+1] * WMTO

                        BW = W + para[iaWbuoy, ip+1]

                        CL = para[iaCL, ip+1]
                        rho = pare[ierho0, ip+1]
                        V = sqrt(2 * BW * cosg / (rho * S * CL))
                        pare[ieu0, ip+1] = V

                        dh = para[iaalt, ip+1] - para[iaalt, ip]
                        dVsq = pare[ieu0, ip+1]^2 - pare[ieu0, ip]^2

                        dR = para[iaRange, ip] - para[iaRange, ip-1]
                        dt = dR * Vgi[ip]
                        rW = exp(-dR * FFC[ip])

                        para[iatime, ip+1] = para[iatime, ip] + dt
                        para[iafracW, ip+1] = para[iafracW, ip] * rW

                  end
            end

            # mission fuel fractions and weights
            fracWa = para[iafracW, ipclimb1]
            fracWe = para[iafracW, ipdescentn]
            freserve = parg[igfreserve]
            fburn = fracWa - fracWe
            ffuel = fburn * (1.0 + freserve)
            Wfuel = WMTO * ffuel
            WTO = Wzero + Wfuel

            parm[imWTO] = WTO
            parm[imWfuel] = Wfuel

            # mission PFEI
            Wburn = WMTO * fburn
            parm[imPFEI] = Wburn/gee * pare[iehfuel, ipcruise1] / (parm[imWpay] * parm[imRange])

      end

      return t_prop
end
