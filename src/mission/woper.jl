"""
    woper(pari, parg, parm, para, pare, 
          parad, pared, itermax, initeng, NPSS_PT, NPSS)

`woper` runs the aircraft through input off-design missions

!!! compat "Future Changes"
      In an upcoming revision, an `aircraft` struct and auxiliary indices will be passed in lieu of pre-sliced `par` arrays.

"""
function woper(pari, parg, parm, para, pare, 
              parad, pared, itermax, initeng)

    # # Initialze some variables
    # ifirst = true
    # NPSS = Base.Process
    # NPSS_TS = Base.Process
    # NPSS_Fan = Base.Process
    # NPSS_AftFan = Base.Process
    # NPSS_PT = true

    time_propsys = 0.0

    tolerW = 1.0e-8
    errw   = 1.0

    #iterate through mission points
    for ip = 1: iptotal
        #iterate through aero parameters
          for ia = 1: iatotal
                #Setting altitude and mach speeds
                para[ia,ip] = parad[ia,ip]
          end
          #iterate through engine parameters
          for ie = 1: ietotal
                pare[ie,ip] = pared[ie,ip]
          end
    end
    
    para[iaalt, ipcruise1] = 10668
    para[iaalt, ipclimbn] = 10668
    para[iaMach, ipclimbn:ipdescent1] .= 0.78

#------ mission-varying excrescence factors disabled in this version
#-      ( also commented out in getparm.f )
#        para(iafexcdw,ip) = parm[imfexcdw]
#        para(iafexcdt,ip) = parm[imfexcdt]
#        para(iafexcdf,ip) = parm[imfexcdf]

    # Calculates surface velocities, boundary layer, wake 
    fusebl!(pari, parg, para, ipcruise1)

#---- assume K.E., dissipation, drag areas will be the same for all points
    KAfTE   = para[iaKAfTE  , ipcruise1] # Kinetic energy area at T.E.
    DAfsurf = para[iaDAfsurf, ipcruise1] # Surface dissapation area 
    DAfwake = para[iaDAfwake, ipcruise1] # Wake dissapation area
    PAfinf  = para[iaPAfinf , ipcruise1] # Momentum area at ∞

    # Assume K.E., Disspation and momentum areas are const. for all mission points:
    para[iaKAfTE  , :] .= KAfTE
    para[iaDAfsurf, :] .= DAfsurf
    para[iaDAfwake, :] .= DAfwake
    para[iaPAfinf , :] .= PAfinf

# ===================================================================
# ---- max range and this mission range
    Rangemax = parg[igRange]
    Rangetot = parm[imRange]

#---- max TO weight
    WMTO = parg[igWMTO]

# ---- zero-fuel weight for this mission
    Wzero = WMTO-
          parg[igWfuel]-
          parg[igWpay]+
          parm[imWpay]

# ===================================================================
# ---- initial fuel and gross takeoff weight estimates from Breguet, R ~ ln(1+f)
    gmax = log(1.0 + parg[igWfuel]/Wzero)
    gmaxp = gmax * Rangetot/Rangemax
    Wfuel = (exp(gmaxp) - 1.0) * Wzero
    WTO = Wzero + Wfuel

    parm[imWfuel] = Wfuel
    parm[imWTO]   = WTO

#---- scale initial weight fractions by takeoff and descent weight ratios
    rTO = WTO/WMTO
    rDE = Wzero/(WMTO-parg[igWfuel])

    para[iafracW, ipstatic ] = parad[iafracW,ipstatic ]*rTO
    para[iafracW, iprotate ] = parad[iafracW,iprotate ]*rTO
    para[iafracW, iptakeoff] = parad[iafracW,iptakeoff]*rTO
    para[iafracW, ipcutback] = parad[iafracW,ipcutback]*rTO

    # Climb
    @inbounds for ip = ipclimb1:ipclimbn
          para[iafracW,ip] = parad[iafracW,ip] * rTO
    end
    # Cruise
    @inbounds for ip = ipcruise1:ipcruisen
          frac = float(ip       -ipcruise1)/
                float(ipcruisen-ipcruise1)
          rCR = rTO*(1.0-frac) + rDE*frac
          para[iafracW,ip] = parad[iafracW,ip] * rCR
    end
    # Descent
    para[iafracW,ipdescent1:ipdescentn] .= parad[iafracW,ipdescent1:ipdescentn] .* rDE
    para[iagamV,:] .= parad[iagamV,:]

#---- estimate takeoff speed and set V,Re over climb and descent
#-    (needed to start trajectory integration)
    ip = iptakeoff
    VTO = pared[ieu0,ip] * sqrt(pared[ierho0,ip]/pare[ierho0,ip])
    ReTO = VTO*pare[ierho0,ip]/pare[iemu0,ip]

    ip = ipcruise1
    VCR = pared[ieu0,ip]
    ReCR = parad[iaReunit,ip]

    for ip = iprotate: ipclimb1
      pare[ieu0,ip] = VTO
      para[iaReunit,ip] = ReTO
    end
    for ip = ipclimb1+1 : ipclimbn
      frac = float(ip-ipclimb1) / float(ipclimbn-ipclimb1)
      V  =  VTO*(1.0-frac) +  VCR*frac
      Re = ReTO*(1.0-frac) + ReCR*frac
      pare[ieu0,ip] = V
      para[iaReunit,ip] = Re
    end
    for ip = ipdescent1: ipdescentn
      frac = float(ip-ipdescent1) / float(ipdescentn-ipdescent1)
      V  =  VTO*frac +  VCR*(1.0-frac)
      Re = ReTO*frac + ReCR*(1.0-frac)
      pare[ieu0,ip] = V
      para[iaReunit,ip] = Re
    end

    if (initeng)
#----- use design case as initial guess for engine state
          for ip = 1: iptotal
                for ie = 1: ietotal
                      pare[ie,ip] = pared[ie,ip]
                end
          end
    end
  
    for ip = ipstatic: ipdescentn
      para[iaCfnace,ip] = parad[iaCfnace,ip]
    end

#--------------------------------------------------------------------------
#---- set wing pitching moment constants
    b  = parg[igb]
    bs = parg[igbs]
    bo = parg[igbo]
    sweep = parg[igsweep]
    Xaxis = parg[igXaxis]
    λs = parg[iglambdas]
    λt = parg[iglambdat]
    AR = parg[igAR]
    fLo = parg[igfLo]
    fLt = parg[igfLt]

    ip = iptakeoff
    cmpo = para[iacmpo,ip]
    cmps = para[iacmps,ip]
    cmpt = para[iacmpt,ip]
    γt = parg[iglambdat]*para[iarclt,ip]
    γs = parg[iglambdas]*para[iarcls,ip]

    CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                            λt,λs,γt,γs, 
                            AR,fLo,fLt,cmpo,cmps,cmpt)

    para[iaCMw0, ipstatic:ipclimb1] .= CMw0
    para[iaCMw1, ipstatic:ipclimb1] .= CMw1

    ip = ipcruise1
    cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]

    γt = parg[iglambdat]*para[iarclt, ip]
    γs = parg[iglambdas]*para[iarcls, ip]
    
    CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                      λt,λs,γt,γs, 
                      AR,fLo,fLt,cmpo,cmps,cmpt)
    
    para[iaCMw0, ipclimb1+1:ipdescentn-1] .= CMw0
    para[iaCMw1, ipclimb1+1:ipdescentn-1] .= CMw1
    
    ip = ipdescentn
    cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
    γt = parg[iglambdat]*para[iarclt, ip]
    γs = parg[iglambdas]*para[iarcls, ip]

    CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                      λt,λs,γt,γs, 
                      AR,fLo,fLt,cmpo,cmps,cmpt)

    para[iaCMw0, ipdescentn] = CMw0
    para[iaCMw1, ipdescentn] = CMw1

#---- tail pitching moment constants
    bh      = parg[igbh]
    boh     = parg[igboh]
    sweeph  = parg[igsweeph]
    λh      = parg[iglambdah]
    ARh     = parg[igARh]
    fLoh = 0.
    fLth = fLt
    cmph = 0.

    CMh0, CMh1 = surfcm(bh, boh, boh, sweeph, Xaxis, λh, 1.0, λh, 1.0,
    ARh, fLoh, fLth, 0.0, 0.0, 0.0)

    para[iaCMh0, :] .= CMh0
    para[iaCMh1, :] .= CMh1

    # Initialize previous weight iterations
    WMTO1, WMTO2, WMTO3 = zeros(Float64, 3) #1st-previous to 3rd previous iteration weight for convergence criterion

#---- no convergence yet
    Lconv = false

# -------------------------------------------------------    
#                   Weight loop
# -------------------------------------------------------    
  @inbounds for  iterw = 1:itermax

    if iterw == itermax
        println("Reached max iterations in weight sizing loop!")
    end

    rlx = 1.0
    if (iterw > itermax-5)
          rlx = 0.5
    end

    ip = ipcruise1
    altkm = para[iaalt, ip]/1000.0
    T0, p0, ρ0, a0, μ0 = atmos(altkm)
    Mach  = para[iaMach, ip]
    pare[iep0  , ip] = p0 
    pare[ieT0  , ip] = T0
    pare[ierho0, ip] = ρ0
    pare[iemu0 , ip] = μ0
    pare[iea0  , ip] = a0

    pare[ieM0, ip] = Mach
    pare[ieu0, ip] = Mach*a0
    para[iaReunit, ip] = Mach*a0 *ρ0/μ0

    # Calling mission
    time_propsys += mission!(pari, parg, parm, para, pare, false)
    # println(parm[imWfuel,:])
    
#-------------------------------------------------------------------------

# Convergence tests
    
    WMTO = parg[igWMTO]
    errw1 = (WMTO - WMTO1)/WMTO
    errw2 = (WMTO - WMTO2)/WMTO
    errw3 = (WMTO - WMTO3)/WMTO

    errw = max(abs(errw1), abs(errw2), abs(errw3))

    if (errw <= tolerW) 
          Lconv = true
          printstyled("Converged!", "\n"; color=:green)
        
          break
    end

#-----------------------------------------------------------------
#---- set previous-iteration weights for next iteration
    WMTO3 = WMTO2
    WMTO2 = WMTO1
    WMTO1 = parg[igWMTO]

    end

return
end
