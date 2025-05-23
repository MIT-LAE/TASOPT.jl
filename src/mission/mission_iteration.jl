"""
    _mission_iteration!(ac, imission, Ldebug; calculate_cruise = false)

Runs aircraft through mission, calculating fuel burn
and other mission variables. Formerly, `mission!()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: aircraft data storage object
    - `imission::Int64`: mission index
    - `Ldebug::Bool`: debugging flag

NOTE: 
This routine assumes that estimates of the climb-leg flight path 
gamma angles are passed in via para[iagamV,ipclimb1:ipclimbn].
These appear as cos(gamma) factors in the climb equations,
and can be passed in as zero with only a minor error.
They are updated and returned in the same para[iagamV,ip] array.

"""
function _mission_iteration!(ac, imission, Ldebug; calculate_cruise = false)
      #Unpack aircraft
      parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear = unpack_ac(ac, imission) 

      ifirst = true

      # HACK TODO add the para back
      # iairf
      initializes_engine = true

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
      ΔTatmos = parm[imDeltaTatm]
      ip = ipstatic
      altkm = para[iaalt, ip] / 1000.0
      T0, p0, ρ0, a0, μ0 = atmos(altkm, ΔTatmos)

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

      T0, p0, rho0, a0, mu0 = atmos(altkm, ΔTatmos)
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
      T0, p0, rho0, a0, mu0 = atmos(altkm, ΔTatmos)

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
      cosL = cosd(wing.layout.sweep)
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
      sweep = wing.layout.sweep
      cosL = cosd(sweep)
      CLmax = clpmax * cosL^2

      #---- Vs stall speed (takeoff condition)
      rho0 = pare[ierho0, ip]
      a0 = pare[iea0, ip]
      S = wing.layout.S
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
      opt_trim_var = "CL_htail"
      balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

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

            T0, p0, rho0, a0, mu0 = atmos(altkm, ΔTatmos)
            pare[iep0, ip] = p0
            pare[ieT0, ip] = T0
            pare[iea0, ip] = a0
            pare[ierho0, ip] = rho0
            pare[iemu0, ip] = mu0

            rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
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

            V = sqrt(2.0 * BW * cosg / (ρ * S * CL))
            Mach = V / Vsound

            para[iaMach, ip] = Mach
            para[iaReunit, ip] = V * ρ / μ

            pare[ieu0, ip] = V
            pare[ieM0, ip] = Mach

                  # Set pitch trim by adjusting CLh
                  Wf = W - Wzero
                  rfuel = Wf / parg[igWfuel]
                  opt_trim_var = "CL_htail"
                  balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

                  if (ip == ipclimb1)
                        computes_wing_direct = false #use explicitly specified wing cdf, cdp
                  else
                        computes_wing_direct = true #use airfoil database
                  end
                  aircraft_drag!(ac, imission, ip, computes_wing_direct)

                  engine.enginecalc!(ac, "off_design", imission, ip, initializes_engine)

            Ftotal = pare[ieFe, ip] * parg[igneng]
            TSFC = pare[ieTSFC, ip]
            DoL = para[iaCD, ip] / para[iaCL, ip]

            # Calculate improved flight angle
            ϕ = Ftotal / BW
            sing = (ϕ - DoL * sqrt(1.0 - ϕ^2 + DoL^2)) / (1.0 + DoL^2)
            gamV = asin(sing)
            cosg = sqrt(1.0 - sing^2)
            # gamV = atan(sing, cosg)

            para[iagamV, ip] = gamV
            para[iaROC, ip] = sing * V * 60 / ft_to_m #ft per min

            # Store integrands for range and weight integration using a predictor-corrector scheme
            FoW[ip] = Ftotal / (BW * cosg) - DoL

            mfuel = pare[iemfuel, ip]
            FFC[ip] = mfuel * gee / (W * V * cosg)

            Vgi[ip] = 1.0 / (V * cosg)

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

                  dR = (dh + 0.5 * dVsq / gee) / FoWavg
                  dt = dR * Vgiavg
                  rW = exp(-dR * FFCavg)    #Ratio of weights Wᵢ₊₁/Wᵢ

                  para[iaRange, ip] = para[iaRange, ip-1] + dR
                  para[iatime, ip] = para[iatime, ip-1] + dt
                  para[iafracW, ip] = para[iafracW, ip-1] * rW

                  ifirst = false

            end
            if (ip < ipclimbn)
                  # Predictor integration step, forward Euler
                  if (para[iagamV, ip+1] ≤ 0.0)
                        # if gamV guess is not passed in, use previous point as the guess
                        para[iagamV, ip+1] = para[iagamV, ip]
                  end

                  W = para[iafracW, ip+1] * WMTO   # Initial weight fractions have been set in _size_aircraft!
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


                  para[iaRange, ip+1] = para[iaRange, ip] + dR
                  para[iatime, ip+1] = para[iatime, ip] + dt
                  para[iafracW, ip+1] = para[iafracW, ip] * rW

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
      opt_trim_var = "CL_htail"
      balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

      if calculate_cruise #If start of cruise has to be calculated (e.g., in off-design)
            # println("Calculating cruise point")
            # Calculate only if requested since for design mission start of cruise is the des point and ∴ already calcualted 
            # Calculate drag
            computes_wing_direct = true
            aircraft_drag!(ac, imission, ip, computes_wing_direct)
            DoL = para[iaCD, ip] / para[iaCL, ip]
            W = para[iafracW, ip] * WMTO
            BW = W + para[iaWbuoy, ip]
            F = BW * (DoL + para[iagamV, ip])
            pare[ieFe, ip] = F / parg[igneng] #Store required thrust for engine calcs
            Wpay = parg[igWpay]
            
            engine.enginecalc!(ac, "off_design", imission, ip, initializes_engine)

      end

      # set cruise-climb climb angle, from fuel burn rate and atmospheric dp/dz
      TSFC = pare[ieTSFC, ip]
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
      gamVcr1 = DoL * p0 * TSFC / (ρ0 * gee * V - p0 * TSFC)

      para[iagamV, ip] = gamVcr1

      Ftotal = pare[ieFe, ip] * parg[igneng]

      cosg = cos(gamVcr1)

      FoW[ip] = Ftotal / (BW * cosg) - DoL
      mfuel = pare[iemfuel, ip]
      FFC[ip] = mfuel * gee / (W * V * cosg)
      Vgi[ip] = 1.0 / (V * cosg)

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

      T0, p0, ρ0, a0, μ0 = atmos(altkm, ΔTatmos)
      pare[iep0, ip] = p0
      pare[ieT0, ip] = T0
      pare[iea0, ip] = a0
      pare[ierho0, ip] = ρ0
      pare[iemu0, ip] = μ0
      pare[ieM0, ip] = Mach
      pare[ieu0, ip] = Mach * a0
      para[iaReunit, ip] = Mach * a0 * ρ0 / μ0
      para[iaalt, ip] = altd

      ρcab = max(parg[igpcabin], p0) / (RSL * Tref)
      para[iaWbuoy, ip] = (ρcab - ρ0) * gee * parg[igcabVol]

      # Set pitch trim by adjusting CLh
      Wf = para[iafracW, ip] * WMTO - Wzero
      rfuel = Wf / parg[igWfuel]
      opt_trim_var = "CL_htail"
      balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

      # Calc Drag
      computes_wing_direct = true
      aircraft_drag!(ac, imission, ip, computes_wing_direct)
      DoL = para[iaCD, ip] / para[iaCL, ip]
      W = para[iafracW, ip] * WMTO
      BW = W + para[iaWbuoy, ip]
      Ftotal = BW * (DoL + para[iagamV, ip])
      pare[ieFe, ip] = Ftotal / parg[igneng]

      engine.enginecalc!(ac, "off_design", imission, ip, initializes_engine)
      TSFC = pare[ieTSFC, ip]

      V = pare[ieu0, ip]
      p0 = pare[iep0, ip]
      ρ0 = pare[ierho0, ip]
      DoL = para[iaCD, ip] / para[iaCL, ip]

      gamVcr2 = DoL * p0 * TSFC / (ρ0 * gee * V - p0 * TSFC)
      para[iagamV, ip] = gamVcr2

      cosg = cos(gamVcr1)

      FoW[ip] = Ftotal / (BW * cosg) - DoL

      mfuel = pare[iemfuel, ip]
      FFC[ip] = mfuel * gee / (W * V * cosg) 

      Vgi[ip] = 1.0 / (V * cosg)

      ip1 = ipcruise1
      ipn = ipcruisen

      FoWavg = 0.5 * (FoW[ipn] + FoW[ip1])
      FFCavg = 0.5 * (FFC[ipn] + FFC[ip1])
      Vgiavg = 0.5 * (Vgi[ipn] + Vgi[ip1])


      dtcruise = dRcruise * Vgiavg
      rWcruise = exp(-dRcruise * FFCavg)

      para[iaRange, ipn] = para[iaRange, ip1] + dRcruise
      para[iatime, ipn] = para[iatime, ip1] + dtcruise
      para[iafracW, ipn] = para[iafracW, ip1] * rWcruise


      # set intermediate points over cruise, if any, just by interpolating
      for ip = ipcruise1+1:ipcruisen-1
            frac = float(ip - ipcruise1) / float(ipcruisen - ipcruise1)
            Mach = para[iaMach, ip]
            para[iaalt, ip] = altc * (1.0 - frac) + altd * frac
            altkm = para[iaalt, ip] / 1000.0
            T0, p0, rho0, a0, mu0 = atmos(altkm, ΔTatmos)
            pare[iep0, ip] = p0
            pare[ieT0, ip] = T0
            pare[iea0, ip] = a0
            pare[ierho0, ip] = rho0
            pare[iemu0, ip] = mu0
            pare[ieM0, ip] = Mach
            pare[ieu0, ip] = Mach * a0
            para[iaReunit, ip] = Mach * a0 * rho0 / mu0

            rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
            para[iaWbuoy, ip] = (rhocab - rho0) * gee * parg[igcabVol]

            para[iaRange, ip] = para[iaRange, ipcruise1] + dRcruise * frac
            para[iatime, ip] = para[iatime, ipcruise1] + dtcruise * frac
            para[iafracW, ip] = para[iafracW, ipcruise1] * rWcruise^frac

      end

      # Descent
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
            T0, p0, ρ0, a0, μ0 = atmos(altkm, ΔTatmos)
            pare[iep0, ip] = p0
            pare[ieT0, ip] = T0
            pare[iea0, ip] = a0
            pare[ierho0, ip] = ρ0
            pare[iemu0, ip] = μ0

            para[iaRange, ip] = R
            para[iaalt, ip] = alt

            rhocab = max(parg[igpcabin], p0) / (RSL * Tref)
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
            opt_trim_var = "CL_htail"
            balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

            if (ip == ipdescentn)
                  # use explicitly specified wing cdf,cdp
                  computes_wing_direct = false
            else
                  # use airfoil database for wing cdf,cdp
                  computes_wing_direct = true
            end
            aircraft_drag!(ac, imission, ip, computes_wing_direct)

            # set up for engine calculation
            sing = sin(gamVde)
            cosg = cos(gamVde)
            DoL = para[iaCD, ip] / para[iaCL, ip]
            Fspec = BW * (sing + cosg * DoL)
            pare[ieFe, ip] = Fspec / parg[igneng]

            if initializes_engine
                  pare[iembf, ip] = pare[iembf, ip-1]
                  pare[iemblc, ip] = pare[iemblc, ip-1]
                  pare[iembhc, ip] = pare[iembhc, ip-1]
                  pare[iepif, ip] = pare[iepif, ip-1]
                  pare[iepilc, ip] = pare[iepilc, ip-1]
                  pare[iepihc, ip] = pare[iepihc, ip-1]
                  initializes_engine = false #Apparently, this helps convergence
                                             #on descent, where the engine is at lower throttle

                  # make better estimate for new Tt4, adjusted for new ambient T0
                  dTburn = pare[ieTt4, ip-1] - pare[ieTt3, ip-1]
                  OTR = pare[ieTt3, ip-1] / pare[ieTt2, ip-1]
                  Tt3 = pare[ieT0, ip] * OTR
                  pare[ieTt4, ip] = Tt3 + dTburn + 50.0

                  # make better estimate for new pt5, adjusted for new ambient p0
                  pare[iept5, ip] = pare[iept5, ip-1] * pare[iep0, ip] / pare[iep0, ip-1]

            end

            engine.enginecalc!(ac, "off_design", imission, ip, initializes_engine)

            # store effective thrust, effective TSFC
            F = pare[ieFe, ip] * parg[igneng]
            TSFC = pare[ieTSFC, ip]

            # store integrands for Range and Weight integration
            FoW[ip] = F / (BW * cosg) - DoL

            Vgi[ip] = 1.0 / (V * cosg)

            # if F < 0, then TSFC is not valid, so calculate mdot_fuel directly
            mfuel = pare[iemfuel, ip]
            FFC[ip] = mfuel * gee / (W * V * cosg)

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
      Wfvent = parm[imWfvent] #Weight of fuel that is vented from tank

      ffvent = Wfvent/WMTO #weight fraction of vented fuel
      
      fracWa = para[iafracW, ipclimb1]
      fracWe = para[iafracW, ipdescentn]
      freserve = parg[igfreserve]
      fburn = fracWa - fracWe + ffvent #include vented fuel, if any
      ffuel = fburn * (1.0 + freserve)
      Wfuel = WMTO * ffuel
      WTO = Wzero + Wfuel

      parm[imWTO] = WTO
      parm[imWfuel] = Wfuel
      #TODO the above calculation does not account for the effect of venting on the flight profile

      # mission PFEI
      Wburn = WMTO * fburn
      parm[imPFEI] = Wburn/gee * parg[igLHVfuel] / (parm[imWpay] * parm[imRange])

end
