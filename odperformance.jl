function odperf!(pari, parg, parm, para, pare, Wfrac0,
    NPSS_TS::Base.Process, 
    NPSS_Fan::Base.Process, 
    NPSS_AftFan::Base.Process, Ldebug, ifirst)

calc_ipc1 = true
# ifirst = true

itergmax::Int64 = 15
gamVtol  = 1.0e-12

# BLI
    fBLIf = parg[igfBLIf]
    DAfsurf = para[iaDAfsurf, ipcruise1]
    KAfTE = para[iaKAfTE, ipcruise1]

# Zero-fuel weight for this mission
    Wzero = parg[igWMTO] - parg[igWfuel] - parg[igWpay] + parm[imWpay]
# Payload fraction for this mission
    rpay = parm[imWpay]/parg[igWpay]
    ξpay = 0.

S    = parg[igS]
# Flight levels to output in BADA file:
FL = [  0 ,    5 ,   10 ,   15 ,   20 ,   
       30 ,   40 ,   60 ,   80 ,  100 , 
      120 ,  140 ,  160 ,  180 ,  200 ,
      220 ,  240 ,  260 ,  280 ,  290 , 
      310 ,  330 ,  350 ,  370 ,  390 ,
      410 ,  430 ,  431 ]
# FL = LinRange(0,430, 44)

alts = FL*100*ft_to_m
N = length(FL)
T0s = zeros(Float64, N)
p0s = zeros(Float64, N)
ρ0s = zeros(Float64, N)
a0s = zeros(Float64, N)
μ0s = zeros(Float64, N)
V0s = zeros(Float64, N)
Reunits = zeros(Float64, N)
M0s = zeros(Float64, N)
Wbouys = zeros(Float64, N)
Ws = zeros(Float64, N)
γs = zeros(Float64, N)
Tt4s = zeros(Float64, N)

deNOx = zeros(Float64, N)
mdotf = zeros(Float64, N)
ROC   = zeros(Float64, N)
ROCcrz = zeros(Float64, N)
FFcrz  = zeros(Float64, N)
EGTcrz = zeros(Float64, N)

iceil = N
crzmdotf = zeros(Float64, N)
crzTAS   = zeros(Float64, N)

# Initialize climb integrands
FoW = zeros(Float64, N)
FFC = zeros(Float64, N)
Vgi = zeros(Float64, N)

@inbounds for i =1:N
    T0s[i], p0s[i], ρ0s[i], a0s[i], μ0s[i] = atmos(alts[i]/1000)
    rhocab = max( parg[igpcabin] , p0s[i] ) / (RSL*TSL)
    Wbouys[i] = (rhocab-ρ0s[i])*gee*parg[igcabVol]
end
M0s[end] = 0.8
V0s[end] = 0.8 * a0s[end]
Reunits[end] = V0s[end] *ρ0s[end]/μ0s[end]

#---- estimate takeoff speed and set V,Re over climb and descent
cosL = cos(parg[igsweep]*π/180.0)
CLTO = para[iaclpmax,iptakeoff]*cosL^2
# [prash] The assumption here is that Wcruise/WTO~1 and
# just scaling it with rho and CL to get an initial estimate for V
VTO = pare[ieu0,ipcruise1] *
        sqrt(pare[ierho0,ipcruise1] / pare[ierho0,iptakeoff]) *
        sqrt(para[iaCL  ,ipcruise1] / CLTO )
ReTO = VTO*pare[ierho0,iptakeoff]/pare[iemu0 ,iptakeoff]

@inbounds for  ip = 1:N
    frac = (alts[ip] - 0.0)/(alts[N] - 0.0)
    V  =  VTO*(1.0-frac) + V0s[end]*frac
    Re = ReTO*(1.0-frac) + Reunits[end]*frac
    W  = Wfrac0*parg[igWMTO]*(1.0 - frac) + parg[igWMTO]*0.9*frac
    V0s[ip] = V
    Ws[ip] = W
    Reunits[ip] = Re
end
# println(Ws)
#---- set climb Tt4's from fractions
fT1 = 0.2#parg[igfTt4CL1]
fTn = 0.2#parg[igfTt4CLn]
Tt4TO = pare[ieTt4,iptakeoff]
Tt4CR = pare[ieTt4,ipcruise1]
for  ip = 1:N
    frac = (alts[ip] - 0.0)/(alts[N] - 0.0)
    Tfrac = fT1*(1.0-frac) + fTn*frac
    Tt4s[ip] = Tt4TO*(1.0-Tfrac) + Tt4CR*Tfrac
end
# println(Tt4s)
# integrate trajectory over climb
for   i = 1:N
    
    # velocity calculation from CL, Weight, altitude
    W  = Ws[i]
    if i == 1
        ip = ipclimbn # dummy location to store values
        CL = para[iaCL, ipclimb1]
    else
        ip = ipclimbn # dummy location to store values
        CL = para[iaCL, ipcruise1]
    end
    ρ  = ρ0s[i]
    μ  = μ0s[i]
    Vsound = a0s[i]
    cosg = cos(γs[i])
    BW   = W + Wbouys[i]

    dgamV  = 1.0
    Ftotal = 0.0
    DoL    = 0.0
    V      = 0.0
    if (Ldebug)
        printstyled(@sprintf("\t%5s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                        "iterg", "dgamV", "gamV", "BW", "Ftotal", "DoL", "V", "Alt", "M", "Tt4"); color = :light_green)
    end
    for  iterg = 1:itergmax
        V = sqrt(2.0*BW*cosg/(ρ*S*CL))
        Mach = V/Vsound

        V0s[i] = V
        M0s[i] = Mach
        Reunits[i] = V*ρ/μ

        # Set pitch trim by adjusting CLh
        Wf = W - Wzero
        rfuel = Wf/parg[igWfuel]
        itrim = 1
        balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)
        pare[ieM2, ip] = 0.0

        if (i == 1)
            icdfun = 0
        else 
            icdfun = 1
        end
        cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)

        ρ0 = pare[ierho0, ipcruise1]
        u0 = pare[ieu0  , ipcruise1] 
        Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
        Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines

        # ifirst = false
        Ftotal, η, P, Hrej, heatexcess,
        mdotf[i], BSFC,
        deNOx[i], EGT = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt4s[i],
                                    Kinl, Φinl, parpt, parmot, pargen, ifirst, false)
        ifirst = false
        # println(Tt4s[i]*10/18/T0s[i])
        DoL = para[iaCD, ip]/ para[iaCL, ip]
        # println("LoD = $(1/DoL)")
        # Calculate improved flight angle
        ϕ = Ftotal/BW
        sing = (ϕ - DoL*sqrt(1.0 - ϕ^2 + DoL^2))/(1.0 + DoL^2)
        gamV = asin(sing)
        ROC[i] = sing*V/ft_to_m*60 # Rate of climb in m/s
        
        dgamV = gamV - γs[i]
        γs[i] = gamV
        
        # para[iagamV, ip] = gamV
        if (Ldebug)
            printstyled(@sprintf("\t%5d  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4f  %9.4f  %9.4f\n",
                        iterg, abs(dgamV), gamV*180/π, BW, Ftotal, DoL, V, alts[i], Mach, Tt4s[i]); color =:light_green)
        end
        
        if(abs(dgamV) < gamVtol) 
            break
        end
    end
    if ROC[i] < 100
        iceil = i
        ROC[i] = 0.0
        println("Climb ceiling reached at FL",FL[i])
        break
    end
    if (abs(dgamV) > gamVtol) 
            println("Climb gamV not converged")
    end
    pare[ieFe, ip] = Ftotal
    # Store integrands for range and weight integration using a predictor-corrector scheme
    FoW[i] = Ftotal/(BW*cosg) - DoL
    # FFC[ip] = Ftotal*TSFC/(W*V*cosg)
    FFC[i] = mdotf[i]*gee/(W*V*cosg)
    Vgi[i] = 1.0/(V*cosg)

    Mach = M0s[i]
    CL   = para[iaCL  , ip]
    CD   = para[iaCD  , ip]
    gamV = γs[i]

    if (i > 1)
            # Corrector step
            dh   = alts[i]- alts[i-1]
            dVsq = V0s[i]^2 - V0s[i-1]^2

            FoWavg = 0.5*(FoW[i] + FoW[i-1])
            FFCavg = 0.5*(FFC[i] + FFC[i-1])
            Vgiavg = 0.5*(Vgi[i] + Vgi[i-1])

            dR = (dh + 0.5*dVsq/gee) / FoWavg
            dt = dR*Vgiavg
            rW = exp(-dR*FFCavg)    #Ratio of weights Wᵢ₊₁/Wᵢ
    
            para[iaRange, ip] = para[iaRange, ip-1] + dR
            para[iatime , ip] = para[iatime , ip-1] + dt
            Ws[i] = Ws[i-1]*rW

            ifirst = false
        
    end
    if (i < N)
    # Predictor integration step, forward Euler
        # if gamV guess is not passed in, use previous point as the guess
        if (γs[i+1] ≤ 0.0)
            γs[i+1] = γs[i]
        end
        
        W  = Ws[i+1]  # Initial weight fractions have been set in wsize.jl
        CL = para[iaCL,ipcruise1]
        ρ  = ρ0s[i+1]
        cosg = cos(γs[i+1])

        BW = W + Wbouys[i+1]

        V = sqrt(2*BW*cosg/(ρ*S*CL))
        V0s[i+1] = V
        dh   = alts[i+1] - alts[i]
        dVsq = V0s[i+1]^2 - V0s[i]^2

        dR = (dh + 0.5*dVsq/gee) / FoW[i]
        dt = dR*Vgi[i]
        rW = exp(-dR*FFC[i])

        para[iaRange,ip+1] = para[iaRange,ip] + dR
        para[iatime ,ip+1] = para[iatime ,ip] + dt
        Ws[i+1] = Ws[i]*rW


    end
    # Do cruise too 
    if FL[i]≥ 270 && FL[i]≤430
        ip = ipcruise1
        
        Wf = Ws[i] - Wzero
        rfuel = Wf/parg[igWfuel]
        itrim = 1
        balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)
        icdfun = 1
        cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)
        DoL = para[iaCD, ip]/ para[iaCL, ip]
        W  = Ws[i]
        BW = W + Wbouys[i]
        F  = BW*(DoL)

        Mach = 0.8
        V = Mach*a0s[i]
        # Run powertrain [TODO] actually should run this to a required thrust not Tt4
        ρ0 = ρ0s[i]
        u0 = Mach*a0s[i]
        Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
        Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines
        TR = 8.5 # Tt4/ Tamb
        Tt4 = max(2850, min(3550, T0s[i]*TR*18/10)) # convert to [R]
        # println(alts[i], " ", Tt4)
        Ftotal, η, P, Hrej, heatexcess,
        FFcrz[i], BSFC,
        deNOxcrz, EGT = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt4 ,
                                        Kinl, Φinl, parpt, parmot, pargen, ifirst, Ldebug)

        gam = Ftotal/BW - DoL
        ROCcrz[i] = sin(gam)*V*60/0.3048
        while gam>0
            TR = TR*0.99
            Tt4 = max(2850, min(3550, T0s[i]*TR*18/10)) # convert to [R]
            # println(alts[i], " ", Tt4)
            Ftotal, η, P, Hrej, heatexcess,
            crzmdotf[i], BSFC,
            deNOxcrz, EGT = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt4 ,
            Kinl, Φinl, parpt, parmot, pargen, ifirst, Ldebug)
            gam = Ftotal/BW - DoL
            if Tt4 == 2850.0 && gam>0
                gam = 0.0
                println("Min Temp reached")
            end
        end
        EGTcrz[i] = EGT


        crzTAS[i] = V0s[i]/kts_to_mps
    end

end # done integrating climb

return Ws[1], alts[iceil], V0s, ROC, mdotf, crzmdotf, crzTAS, EGTcrz, FFcrz, ROCcrz
end