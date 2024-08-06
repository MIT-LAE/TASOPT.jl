"""
    woper(ac, mi, itermax, initeng, saveOffDesign)

`woper` runs the aircraft through input off-design missions

!!! details "🔃 Inputs and Outputs"
**Inputs:**
- `ac::aircraft`: Aircraft with first mission being the design mission
- `mi::Int4`: Off design mission to run (Default: 1)
- `itermax::Int64`: Maximum iterations for sizing loop
- `initeng::Boolean`: Use design case as initial guess for engine state if true
- `saveOffDesign::Boolean`: Set true if you want computed quanties to be saved in the selected off design par arrays of the aircraft model

**Outputs:**
- No explicit outputs. Computed quantities are saved to `par` arrays of `aircraft` model for the off design mission selected

"""
function woper(ac, mi = 1; itermax = 35, initeng = true, saveOffDesign = false)

    pari = ac.pari
    parg = ac.parg
    parm = ac.parm[:,mi:mi]
    para = ac.para[:,:,mi:mi]
    pare = ac.pare[:,:,mi:mi]
    parad = ac.parad
    pared = ac.pared

    fuse = ac.fuse
    wing = ac.wing

    time_propsys = 0.0

    tolerW = 1.0e-8
    errw   = 1.0

    para .= parad
    pare .= pared
    
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

    if initeng
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
    sweep = wing.layout.sweep
    Xaxis = parg[igXaxis]
    λs = wing.inboard.layout.λ
    λt = wing.outboard.layout.λ
    AR = wing.layout.AR
    fLo = parg[igfLo]
    fLt = parg[igfLt]

    ip = iptakeoff
    cmpo = para[iacmpo,ip]
    cmps = para[iacmps,ip]
    cmpt = para[iacmpt,ip]
    γt = wing.outboard.layout.λ*para[iarclt,ip]
    γs = wing.inboard.layout.λ*para[iarcls,ip]

    CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                            λt,λs,γt,γs, 
                            AR,fLo,fLt,cmpo,cmps,cmpt)

    para[iaCMw0, ipstatic:ipclimb1] .= CMw0
    para[iaCMw1, ipstatic:ipclimb1] .= CMw1

    ip = ipcruise1
    cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]

    γt = wing.outboard.layout.λ*para[iarclt, ip]
    γs = wing.inboard.layout.λ*para[iarcls, ip]
    
    CMw0, CMw1 = surfcm(b, bs, bo, sweep, Xaxis,
                      λt,λs,γt,γs, 
                      AR,fLo,fLt,cmpo,cmps,cmpt)
    
    para[iaCMw0, ipclimb1+1:ipdescentn-1] .= CMw0
    para[iaCMw1, ipclimb1+1:ipdescentn-1] .= CMw1
    
    ip = ipdescentn
    cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
    γt = wing.outboard.layout.λ*para[iarclt, ip]
    γs = wing.inboard.layout.λ*para[iarcls, ip]

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

    set_ambient_conditions!(ac, ipcruise1)

    # Calling mission
    time_propsys += mission!(pari, parg, parm, para, pare,fuse, wing, false)
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
          if saveOffDesign
            ac.parm[:,mi:mi] = parm
            ac.para[:,:,mi:mi] = para
            ac.pare[:,:,mi:mi] = pare
          end
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