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
function mission!(pari, parg, parm, para, pare,
       NPSS_TS::Base.Process, NPSS_Fan::Base.Process, Ldebug)#, iairf, initeng, ipc1)
    t_prop = 0.0
    calc_ipc1 = true
    ifirst = true

    itergmax::Int64 = 15
    gamVtol  = 1.0e-12
#     gamVtol  = 1.0e-10

    # unpack flags
        iengloc = pari[iiengloc]
        ifclose = pari[iifclose]
        ifuel   = pari[iifuel  ]

    # Mission range
        Rangetot = parm[imRange]
        para[iaRange, ipdescentn] = Rangetot

    # MTOW
        WMTO = parg[igWMTO]

    # Payload fraction for this mission
        rpay = parm[imWpay]/parg[igWpay]
        ξpay = 0.
    # Zero-fuel weight for this mission
        Wzero = WMTO - parg[igWfuel] - parg[igWpay] + parm[imWpay]
    # mission TO weight
        WTO = parm[imWTO]

    # mission TO fuel weight
        WfTO = WTO - Wzero

# set known operating conditions

    # takeoff altitude conditions
           ip = ipstatic
           altkm = para[iaalt,ip]/1000.0
           T_std, p_std, ρ_std, a_std, μ_std = atmos(altkm)
           T0 = parm[imT0TO]
           p0   = p_std
           ρ0   = ρ_std*(T_std/T0)
           a0   = a_std*sqrt(T0/T_std)
           μ0   = μ_std*(T0/T_std) ^ 0.8
     
           pare[iep0  ,ip] = p0
           pare[ieT0  ,ip] = T0
           pare[iea0  ,ip] = a0
           pare[ierho0,ip] = ρ0
           pare[iemu0 ,ip] = μ0
           pare[ieM0  ,ip] = 0.0
           pare[ieu0  ,ip] = 0.0
     
           para[iaMach,ip] = 0.0
           para[iaCL  ,ip] = 0.0
           para[iaReunit,ip] = 0.0
           para[iaWbuoy,ip] = 0.
     
           # ip ∈ iprotate:ipclimb1
           pare[iep0   , iprotate:ipclimb1] .= p0
           pare[ieT0   , iprotate:ipclimb1] .= T0
           pare[iea0   , iprotate:ipclimb1] .= a0
           pare[ierho0 , iprotate:ipclimb1] .= ρ0
           pare[iemu0  , iprotate:ipclimb1] .= μ0
           para[iaWbuoy, iprotate:ipclimb1] .= 0.
      
     
     # Start-of-cruise altitude conditions
           ip = ipcruise1
           Mach  = para[iaMach,ip]
           altkm = para[iaalt,ip]/1000.0
           
           T0, p0, rho0, a0, mu0 = atmos(altkm)
           pare[iep0  ,ip] = p0
           pare[ieT0  ,ip] = T0
           pare[iea0  ,ip] = a0
           pare[ierho0,ip] = rho0
           pare[iemu0 ,ip] = mu0

           pare[ieM0  ,ip]   = Mach
           pare[ieu0  ,ip]   = Mach*a0
           para[iaReunit,ip] = Mach*a0 * rho0/mu0
     
     # End-of-descent altitude conditions
           ip = ipdescentn
           altkm = para[iaalt,ip]/1000.0

           T_std, p_std, rho_std, a_std, mu_std = atmos(altkm)

           T0   = parm[imT0TO]
           p0   = p_std
           rho0 = rho_std*(T_std/T0)
           a0   = a_std*sqrt(T0/T_std)
           mu0  = mu_std*(T0/T_std) ^ 0.8

           pare[iep0  ,ip] = p0
           pare[ieT0  ,ip] = T0
           pare[iea0  ,ip] = a0
           pare[ierho0,ip] = rho0
           pare[iemu0 ,ip] = mu0
           para[iaWbuoy,ip] = 0.
     
     # interpolate CL over climb points, 
     #  between specified ipclimb1+1, ipclimbn values
           CLa = para[iaCL,ipclimb1+1]
           CLb = para[iaCL,ipcruise1 ]
     
           @inbounds for  ip = ipclimb1+1 : ipclimbn
             frac = float(ip       - (ipclimb1+1)) /
                    float(ipclimbn - (ipclimb1+1))
             para[iaCL,ip] = CLa*(1.0-frac^2) +
	                       CLb*     frac^2
           end
     
     #---- interpolate CL over descent points, 
     #-     between specified ipdescent1, ipdescentn-1 values
     #      CLd = para[iaCL,ipcruisen]
     #      CLe = para[iaCL,ipcruisen]
           CLd = para[iaCL,ipcruisen] * 0.96
           CLe = para[iaCL,ipcruisen] * 0.50
     
           @inbounds for  ip = ipdescent1: ipdescentn-1
             frac = float( ip           -ipdescent1) /
	              float((ipdescentn-1)-ipdescent1)
             fb = 1.0-frac
             para[iaCL,ip] = CLd*     fb^2 +
	                       CLe*(1.0-fb^2)
           end
     
#      #---- interpolate altitudes over climb
#            altb = para[iaalt,iptakeoff]
#            altc = para[iaalt,ipcruise1]
#            for ip = ipclimb1+1: ipclimbn
#              frac = float(ip-ipclimb1) / float(ipclimbn-ipclimb1)
#              para[iaalt,ip] = altb*(1.0-frac) + altc*frac
#            end
     
     #---- estimate takeoff speed and set V,Re over climb and descent
           cosL = cos(parg[igsweep]*π/180.0)
           CLTO = para[iaclpmax,iptakeoff]*cosL^2
           # [prash] I think the assumption here is that Wcruise/WTO~1 and
           # just scaling it with rho and CL to get an initial estimate for V
           VTO = pare[ieu0,ipcruise1] *
                  sqrt(pare[ierho0,ipcruise1] / pare[ierho0,iptakeoff]) *
                  sqrt(para[iaCL  ,ipcruise1] / CLTO )
           ReTO = VTO*pare[ierho0,iptakeoff]/pare[iemu0 ,iptakeoff]

            pare[ieu0    , iprotate: ipclimb1] .= VTO
            para[iaReunit, iprotate: ipclimb1] .= ReTO

           @inbounds for  ip = ipclimb1+1: ipclimbn
             frac = float(ip-ipclimb1) / float(ipclimbn-ipclimb1)
             V  =  VTO*(1.0-frac) + pare[ieu0,ipcruise1]*frac
             Re = ReTO*(1.0-frac) + para[iaReunit,ipcruise1]*frac
             pare[ieu0,ip] = V
             para[iaReunit,ip] = Re
           end
           @inbounds for  ip = ipdescent1: ipdescentn
             frac = float(ip-ipdescent1) / float(ipdescentn-ipdescent1)
             V  =  VTO*frac + pare[ieu0,ipcruisen]*(1.0-frac)
             Re = ReTO*frac + para[iaReunit,ipcruisen]*(1.0-frac)
             pare[ieu0,ip] = V
             para[iaReunit,ip] = Re
           end
     
     
    
     #---- takeoff CLmax via section clmax and sweep correction
           ip = iprotate
           clpmax = para[iaclpmax,ip]
           sweep  = parg[igsweep]
           cosL   = cos(sweep*π/180.0)
           CLmax  = clpmax * cosL^2
     
     #---- Vs stall speed (takeoff condition)
           rho0 = pare[ierho0,ip]
           a0   = pare[iea0  ,ip]
           S    = parg[igS]
           Vstall = sqrt(2.0*WTO/(rho0*S*CLmax))
           Mstall = Vstall/a0
           pare[ieu0,ip] = Vstall
           pare[ieM0,ip] = Mstall
           para[iaMach,ip] = Mstall
           para[iaReunit,ip] = Vstall*pare[ierho0,ip] / pare[iemu0 ,ip]
     
     #---- V2 speed per FAR-25  (takeoff,cutback,climb1 condition)
           V2 = Vstall*1.2
           M2 = Mstall*1.2
           CL2 = CLmax/1.2^2
     
           ip = iptakeoff
           pare[ieu0    ,ip] = V2
           pare[ieM0    ,ip] = M2
           para[iaMach  ,ip] = M2
           para[iaReunit,ip] = V2*pare[ierho0,ip] / pare[iemu0 ,ip]
           para[iaCL,ip] = CL2
        
     #---- set pitch trim by adjusting CLh
           Wf = WTO - Wzero
           rfuel = Wf/parg[igWfuel]
           itrim = 1
           balance(pari, parg, view(para, :,ip), rfuel, rpay, ξpay, itrim)

           CLh2 = para[iaCLh,ip]
           xCG2 = para[iaxCG,ip]
           xCP2 = para[iaxCP,ip]
           xNP2 = para[iaxNP,ip]

            para[iaxCG, 1: iprotate] .= xCG2
            para[iaxCP, 1: iprotate] .= xCP2
            para[iaxNP, 1: iprotate] .= xNP2
            para[iaCLh, 1: iprotate] .= 0.

      #---- initial guesses for climb points
           ip = ipcutback
           pare[ieu0    ,ip] = V2
           pare[ieM0    ,ip] = M2
           para[iaMach  ,ip] = M2
           para[iaReunit,ip] = V2*pare[ierho0,ip] / pare[iemu0 ,ip]
           para[iaCL,ip] = CL2
           para[iaxCG,ip] = xCG2
           para[iaxCP,ip] = xCP2
           para[iaxNP,ip] = xNP2
           para[iaCLh,ip] = CLh2

           ip = ipclimb1
           pare[ieu0    ,ip] = V2
           pare[ieM0    ,ip] = M2
           para[iaMach  ,ip] = M2
           para[iaReunit,ip] = V2*pare[ierho0,ip] / pare[iemu0 ,ip]
           para[iaCL,ip] = CL2
           para[iaxCG,ip] = xCG2
           para[iaxCP,ip] = xCP2
           para[iaxNP,ip] = xNP2
           para[iaCLh,ip] = CLh2

     #---- also use V2 speed for end-of-descent condition, with weight correction
           ip = ipdescentn
           Vrat = sqrt(para[iafracW,ip]/para[iafracW,ipclimb1])
           pare[ieu0    ,ip] = V2*Vrat
           pare[ieM0    ,ip] = M2*Vrat
           para[iaMach  ,ip] = M2*Vrat
           para[iaReunit,ip] = V2*Vrat*pare[ierho0,ip] / pare[iemu0 ,ip]
           para[iaCL,ip] = CL2
     #
     # ============================================================================
     #---- set up climb points at equal altitude intervals, from altb to altc
     #-    (takeoff ground temperature is neglected here -- std atmosphere is used)
           altb = para[iaalt,iptakeoff]
           altc = para[iaalt,ipcruise1]
           @inbounds for   ip = ipclimb1+1: ipclimbn
             frac = float(ip-ipclimb1) / float(ipclimbn-ipclimb1)
             para[iaalt,ip] = altb*(1.0-frac) + altc*frac

             altkm = para[iaalt,ip]/1000.0
             
             T0,p0,rho0,a0,mu0 = atmos(altkm)
             pare[iep0  ,ip] = p0
             pare[ieT0  ,ip] = T0
             pare[iea0  ,ip] = a0
             pare[ierho0,ip] = rho0
             pare[iemu0 ,ip] = mu0
     
             rhocab = max( parg[igpcabin] , p0 ) / (RSL*TSL)
             para[iaWbuoy,ip] = (rhocab-rho0)*gee*parg[igcabVol]
           end
     
     #---- set climb Tt4's from fractions
           fT1 = parg[igfTt4CL1]
           fTn = parg[igfTt4CLn]
           Tt4TO = pare[ieTt4,iptakeoff]
           Tt4CR = pare[ieTt4,ipcruise1]
           @inbounds for   ip = ipclimb1: ipclimbn
             frac = float(ip      -ipclimb1) /
	              float(ipclimbn-ipclimb1)
             Tfrac = fT1*(1.0-frac) + fTn*frac
             pare[ieTt4,ip] = Tt4TO*(1.0-Tfrac) + Tt4CR*Tfrac
           end
     #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     #---- initial values for range, time, weight fraction
           para[iaRange, iprotate:ipclimb1] .= 0.0
           para[iatime , iprotate:ipclimb1] .= 0.0
           para[iafracW, iprotate:ipclimb1] .= WTO/WMTO
           para[iaWbuoy, iprotate:ipclimb1] .= 0.

      # Initialize climb integrands
      FoW = zeros(Float64, iptotal)
      FFC = zeros(Float64, iptotal)
      Vgi = zeros(Float64, iptotal)
      dPdt = zeros(Float64, iptotal)
      dygdt = zeros(Float64, iptotal)

      pare[iePLH2, ipclimb1] = 1.5
      pare[ieyg, ipclimb1] = 0.1

      

      # integrate trajectory over climb
      @inbounds for   ip = ipclimb1:ipclimbn
            if(Ldebug)
            printstyled("Climb angle integration - ip = ", ip-ipclimb1+1, "\n"; color=:light_green)
            end
            
            # velocity calculation from CL, Weight, altitude
            W  = para[iafracW, ip] * WMTO
            CL = para[iaCL, ip]
            ρ  = pare[ierho0, ip]
            μ  = pare[iemu0 , ip]
            Vsound = pare[iea0, ip]
            cosg = cos(para[iagamV, ip])
            BW   = W + para[iaWbuoy, ip]
            

            dgamV  = 1.0
            Ftotal = 0.0
            DoL    = 0.0
            mdotf  = 0.0
            V      = 0.0
            if(Ldebug)
                  printstyled(@sprintf("\t%5s  %10s  %10s  %10s  %10s  %10s  %10s \n",
                        "iterg", "dgamV", "gamV", "BW", "Ftotal", "DoL", "V"); color = :light_green)
            end
            @inbounds for  iterg = 1:itergmax
                  V = sqrt(2.0*BW*cosg/(ρ*S*CL))
                  Mach = V/Vsound

                  para[iaMach, ip] = Mach
                  para[iaReunit, ip] = V*ρ/μ

                  pare[ieu0, ip] = V
                  pare[ieM0, ip] = Mach

                  # Set pitch trim by adjusting CLh
                  Wf = W - Wzero
                  rfuel = Wf/parg[igWfuel]
                  itrim = 1
                  balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)

                  if(ip == ipclimb1)
                        icdfun = 0
                  else 
                        icdfun = 1
                  end
                  cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)

                  t_prop += @elapsed Ftotal, η, P, Hrej,
                  mdotf, BSFC,
                  deNOx = PowerTrainOD(NPSS_TS, NPSS_Fan, para[iaalt, ip], Mach, pare[ieTt4, ip],
                                                0.0, 0.0, parpt, parmot, pargen, ifirst, Ldebug)
                  ifirst = false
                  pare[iedeNOx, ip] = deNOx
                  pare[iemdotf, ip] = mdotf
                  pare[ieemot:ieethermal, ip] .= η

                  DoL = para[iaCD, ip]/ para[iaCL, ip]

                  # Calculate improved flight angle
                  ϕ = Ftotal/BW
                  sing = (ϕ - DoL*sqrt(1.0 - ϕ^2 + DoL^2))/(1.0 + DoL^2)
                  cosg = sqrt(1.0 - sing^2)
                  gamV = atan(sing, cosg)

                  dgamV = gamV - para[iagamV, ip]
                  
                  para[iagamV, ip] = gamV
                  if(Ldebug)
                        printstyled(@sprintf("\t%5d  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e\n",
                                 iterg, abs(dgamV), gamV*180/π, BW, Ftotal, DoL, V); color =:light_green)
                  end
                  
                  if(abs(dgamV) < gamVtol) 
                        break
                  end
            end

            if(abs(dgamV) > gamVtol) 
                  println("Climb gamV not converged")
            end
            pare[ieFe, ip] = Ftotal
            # Store integrands for range and weight integration using a predictor-corrector scheme
            FoW[ip] = Ftotal/(BW*cosg) - DoL
            # FFC[ip] = Ftotal*TSFC/(W*V*cosg)
            FFC[ip] = mdotf*gee/(W*V*cosg)
            Vgi[ip] = 1.0/(V*cosg)

            ΔT = pare[ieT0, ip] - 20.0
            P = pare[iePLH2, ip]
            yg = pare[ieyg, ip]
            dygdt[ip] =  mdotf/ρmix(yg, P)
            dPdt[ip] = dPdt_LH2(ΔT, mdotf, P, yg)

            Mach = para[iaMach, ip]
            CL   = para[iaCL  , ip]
            CD   = para[iaCD  , ip]
            gamV = para[iagamV, ip]
            
            if(ip > ipclimb1)
                  # Corrector step
                  dh   = para[iaalt, ip] - para[iaalt, ip-1]
                  dVsq = pare[ieu0, ip]^2 - pare[ieu0, ip-1]^2

                  FoWavg = 0.5*(FoW[ip] + FoW[ip-1])
                  FFCavg = 0.5*(FFC[ip] + FFC[ip-1])
                  Vgiavg = 0.5*(Vgi[ip] + Vgi[ip-1])

                  dPdtavg = 0.5*(dPdt[ip] + dPdt[ip-1])
                  dygdtavg = 0.5*(dygdt[ip] + dygdt[ip-1])


                  dR = (dh + 0.5*dVsq/gee) / FoWavg
                  dt = dR*Vgiavg
                  rW = exp(-dR*FFCavg)    #Ratio of weights Wᵢ₊₁/Wᵢ
                 
                  dP = dPdtavg*Vgiavg *dR
                  dyg = dygdtavg*Vgiavg *dR
               
                  para[iaRange, ip] = para[iaRange, ip-1] + dR
                  para[iatime , ip] = para[iatime , ip-1] + dt
                  para[iafracW, ip] = para[iafracW, ip-1]*rW

                  pare[iePLH2, ip] = pare[iePLH2, ip-1] + dP
                  pare[ieyg  , ip] = pare[ieyg  , ip-1] + dyg

                  ifirst = false
                
            end
            if(ip < ipclimbn)
            # Predictor integration step, forward Euler
                  if(para[iagamV,ip+1] ≤ 0.0)
                  # if gamV guess is not passed in, use previous point as the guess
                        para[iagamV,ip+1] = para[iagamV,ip]
                  end
                  
                  W  = para[iafracW,ip+1]*WMTO   # Initial weight fractions have been set in wsize.jl
                  CL = para[iaCL,ip+1]
                  ρ  = pare[ierho0,ip+1]
                  cosg = cos(para[iagamV,ip+1])
      
                  BW = W + para[iaWbuoy,ip+1]
      
                  V = sqrt(2*BW*cosg/(ρ*S*CL))
                  pare[ieu0,ip+1] = V
                  dh   = para[iaalt,ip+1]   - para[iaalt,ip]
                  dVsq = pare[ieu0 ,ip+1]^2 - pare[ieu0 ,ip]^2 
      
                  dR = (dh + 0.5*dVsq/gee) / FoW[ip]
                  dt = dR*Vgi[ip]
                  rW = exp(-dR*FFC[ip])

                  dP = dPdt[ip]*Vgi[ip] * dR
                  dyg = dygdt[ip]*Vgi[ip] *dR
      
                  para[iaRange,ip+1] = para[iaRange,ip] + dR
                  para[iatime ,ip+1] = para[iatime ,ip] + dt
                  para[iafracW,ip+1] = para[iafracW,ip]*rW

                  pare[iePLH2, ip+1] = pare[iePLH2, ip] + dP
                  pare[ieyg  , ip+1] = pare[ieyg  , ip] + dyg
            end 

      end # done integrating climb

      # First cruise point is last climb point
      para[iaRange, ipcruise1] = para[iaRange, ipclimbn]
      para[iatime , ipcruise1] = para[iatime , ipclimbn]
      para[iafracW, ipcruise1] = para[iafracW, ipclimbn]
      para[iaWbuoy, ipcruise1] = para[iaWbuoy, ipclimbn]

      pare[iePLH2, ipcruise1] = pare[iePLH2, ipclimbn] 
      pare[ieyg  , ipcruise1] = pare[ieyg  , ipclimbn] 

      # Cruise start
      ip = ipcruise1
      # Set pitch trim by adjusting CLh
      Wf = para[iafracW, ip]*WMTO - Wzero
      rfuel = Wf/parg[igWfuel]
      itrim = 1
      balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)

      if(calc_ipc1)
            # Calculate only if requested since for design mission TOC is the des point and ∴ already calcualted 
            # Calculate drag
            icdfun = 1
            cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)
            DoL = para[iaCD, ip]/ para[iaCL, ip]
            W  = para[iafracW, ip]*WMTO
            BW = W + para[iaWbuoy, ip]
            F  = BW*(DoL + para[iagamV, ip])

            Mach = pare[ieM0, ip]
            # Run powertrain [TODO] actually should run this to a required thrust not Tt4
            t_prop += @elapsed Ftotal, η, P, Hrej,
            mdotf, BSFC,
            deNOx = PowerTrainOD(NPSS_TS, NPSS_Fan, para[iaalt, ip], Mach, pare[ieTt4, ip],
                                          0.0, 0.0, parpt, parmot, pargen, ifirst, Ldebug)
            
            pare[iedeNOx, ip] = deNOx
            pare[iemdotf, ip] = mdotf
            pare[ieemot:ieethermal, ip] .= η
            
            # println("Cruise Ftotal, offdes = ", Ftotal)
            # println("Thrust req = ", F)
            # printstyled("Tt4 = $(pare[ieTt4, ip])\n"; color=:red)
            # println("Ftotal diff = ", Ftotal - pare[ieFe, ip])
      end
      V   = pare[ieu0, ip]
      p0  = pare[iep0   ,ip]
      ρ0  = pare[ierho0 ,ip]
      DoL = para[iaCD, ip]/para[iaCL, ip]
      W   = para[iafracW, ip]*WMTO
      BW  = W + para[iaWbuoy, ip]

      # Calculate Cruise-climb angle:
      # ------
      # γ ≈ dh/dR = dh/dp × dp/dW × dW/dR
      #           = -1/(ρg) × p/W × -̇mf g/(Vcosγ)
      #        ∴γ ≈ ̇mf*p/(ρWV)
      # For conventional aircraft can represent this as function of TSFC: (see TASOPT docs Eq 438)
      # gamVcr1 = DoL*p0*TSFC/(ρ0*gee*V - p0*TSFC)
      gamVcr1 = mdotf*p0/(ρ0*W*V) # TSFC not the most useful for turbo-electric systems. Buoyancy weight has been neglected.
      para[iagamV, ip] = gamVcr1
      
      Ftotal = pare[ieFe, ip]
      cosg = cos(gamVcr1)  

      FoW[ip] = Ftotal/(BW*cosg) - DoL
      # FFC[ip] = Ftotal*TSFC/(W*V*cosg)
      FFC[ip] = mdotf*gee/(W*V*cosg)
      Vgi[ip] = 1.0/(V*cosg)

      ΔT = pare[ieT0, ip] - 20.0
      P = pare[iePLH2, ip]
      yg = pare[ieyg, ip]
      dygdt[ip] =  mdotf/ρmix(yg, P)
      dPdt[ip] = dPdt_LH2(ΔT, mdotf, P, yg)

      #---- set end-of-cruise point "d" using cruise and descent angles, 
      #-     and remaining range legs
            gamVde1 = parm[imgamVDE1]
            gamVden = parm[imgamVDEn]
            gamVdeb = 0.5*(gamVde1+gamVden)

            alte = para[iaalt,ipdescentn]
            dRclimb  = para[iaRange,ipclimbn] - para[iaRange,ipclimb1]
            dRcruise = (alte - altc - gamVdeb*(Rangetot-dRclimb)) / (gamVcr1 - gamVdeb)

            altd = altc + gamVcr1*dRcruise

      # Final cruise point
            ip = ipcruisen
            Mach = para[iaMach, ip]
            altkm = altd/1000.0

            T0, p0, ρ0, a0, μ0 = atmos(altkm)
            pare[iep0  ,ip] = p0
            pare[ieT0  ,ip] = T0
            pare[iea0  ,ip] = a0
            pare[ierho0,ip] = ρ0
            pare[iemu0 ,ip] = μ0
 
            pare[ieM0  ,ip]   = Mach
            pare[ieu0  ,ip]   = Mach*a0
            para[iaReunit,ip] = Mach*a0 * ρ0/μ0

            para[iaalt, ip]  = altd
            
            ρcab = max(parg[igpcabin], p0)/ (RSL*TSL)
            para[iaWbuoy, ip] = (ρcab - ρ0)*gee*parg[igcabVol]
            
            # Set pitch trim by adjusting CLh
            Wf = para[iafracW, ip]*WMTO - Wzero
            rfuel = Wf/parg[igWfuel]
            itrim = 1
            balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)
            # Calc Drag
            icdfun = 1
            cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)
            DoL = para[iaCD, ip]/ para[iaCL, ip]
            W  = para[iafracW, ip]*WMTO
            BW = W + para[iaWbuoy, ip]
            F  = BW*(DoL + para[iagamV, ip])

            Mach = pare[ieM0, ip]
            # Run powertrain [TODO] actually should run this to a required thrust not Tt4
            t_prop += @elapsed Ftotal, η, P, Hrej,
            mdotf, BSFC,
            deNOx = PowerTrainOD(NPSS_TS, NPSS_Fan, para[iaalt, ip], Mach, pare[ieTt4, ip],
                                          0.0, 0.0, parpt, parmot, pargen, ifirst, Ldebug)

            pare[iedeNOx, ip] = deNOx
            pare[iemdotf, ip] = mdotf
            pare[ieemot:ieethermal, ip] .= η

            V   = pare[ieu0, ip]
            p0  = pare[iep0   ,ip]
            ρ0  = pare[ierho0 ,ip]
            DoL = para[iaCD, ip]/para[iaCL, ip] 

            gamVcr2 = mdotf*p0/(ρ0*W*V)
            para[iagamV, ip] = gamVcr2

            cosg = cos(gamVcr1)  

            FoW[ip] = Ftotal/(BW*cosg) - DoL
            # FFC[ip] = Ftotal*TSFC/(W*V*cosg)
            FFC[ip] = mdotf*gee/(W*V*cosg)
            Vgi[ip] = 1.0/(V*cosg)

            ΔT = pare[ieT0, ip] - 20.0
            P = pare[iePLH2, ip]
            dygdt[ip] =  mdotf/ρmix(yg, P)
            dPdt[ip] = dPdt_LH2(ΔT, mdotf, P, yg)

            ip1 = ipcruise1
            ipn = ipcruisen

            FoWavg = 0.5*(FoW[ipn] + FoW[ip1])
            FFCavg = 0.5*(FFC[ipn] + FFC[ip1])
            Vgiavg = 0.5*(Vgi[ipn] + Vgi[ip1])

            dPdtavg = 0.5*(dPdt[ipn] + dPdt[ip1])
            dygdtavg = 0.5*(dygdt[ipn] + dygdt[ip1])
            dP = dPdtavg*Vgiavg *dRcruise
            dyg = dygdtavg*Vgiavg *dRcruise
      
            dtcruise = dRcruise*Vgiavg
            rWcruise = exp(-dRcruise*FFCavg)
      
            para[iaRange,ipn] = para[iaRange,ip1] + dRcruise
            para[iatime ,ipn] = para[iatime ,ip1] + dtcruise
            para[iafracW,ipn] = para[iafracW,ip1]*rWcruise

            pare[iePLH2, ipn] = pare[iePLH2, ip1] + dP
            pare[ieyg, ipn] = pare[ieyg, ip1] + dyg

      # Descent
            ip = ipdescent1
            pare[iep0  ,ip]  = pare[iep0  ,ipcruisen]  
            pare[ieT0  ,ip]  = pare[ieT0  ,ipcruisen]  
            pare[iea0  ,ip]  = pare[iea0  ,ipcruisen]  
            pare[ierho0,ip]  = pare[ierho0,ipcruisen]  
            pare[iemu0 ,ip]  = pare[iemu0 ,ipcruisen]  
            pare[ieM0  ,ip]  = pare[ieM0  ,ipcruisen]  
            pare[ieu0  ,ip]  = pare[ieu0  ,ipcruisen]  
      
            para[iaMach  ,ip] = para[iaMach  ,ipcruisen]
            para[iaReunit,ip] = para[iaReunit,ipcruisen]
            para[iaalt   ,ip] = para[iaalt   ,ipcruisen]
      
            para[iaRange,ip] = para[iaRange,ipcruisen]
            para[iatime ,ip] = para[iatime ,ipcruisen]
            para[iafracW,ip] = para[iafracW,ipcruisen]
            para[iaWbuoy,ip] = para[iaWbuoy,ipcruisen]


      # Mission fuel fractions and weights
            fracWa = para[iafracW, ipclimb1  ]
            # fracWe = para[iafracW, ipdescentn] 
            fracWe = para[iafracW, ipdescent1]*(1-(1-0.993)*((para[iaalt, ipdescent1]/ft_to_m)/39858.0))# Pp temp set to 95% of ipdescent1 instead of ipdescentn since descent has not been calculated yet 
            freserve = parg[igfreserve]
            fburn = fracWa - fracWe            # Burnt fuel fraction
            ffuel = fburn*(1.0+freserve)
            Wfuel  = WMTO*ffuel
            WTO = Wzero + Wfuel

            parm[imWTO  ] = WTO
            parm[imWfuel] = Wfuel
            # printstyled("Wfuel = $Wfuel \n fburn = $fburn \n", color=:red)

            Wburn = WMTO*fburn
            parm[imPFEI] = Wburn/gee*120.0e6 / (parm[imWpay]*parm[imRange])

      return t_prop
end
