function odperf!(pari, parg, parm, para, pare, Wfrac, FL, 
    NPSS_TS::Base.Process, 
    NPSS_Fan::Base.Process, 
    NPSS_AftFan::Base.Process, Ldebug, ifirst, NPSS_PT, NPSS::Base.Process)

calc_ipc1 = true
# ifirst = true

itergmax::Int64 = 15
gamVtol  = 1.0e-12

# BLI
    fBLIf   = parg[igfBLIf]
    DAfsurf = para[iaDAfsurf, ipcruise1]
    KAfTE   = para[iaKAfTE, ipcruise1]

# Zero-fuel weight for this mission
    Wzero   = parg[igWMTO] - parg[igWfuel] - parg[igWpay] + parm[imWpay]  #This ensures that this will work for multi-mission in the future
    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
# Payload fraction for this mission
    rpay = parm[imWpay]/parg[igWpay]
    ξpay = 0.

S    = parg[igS]
# Convert FLs to altitude in m
alts = FL*100*ft_to_m

# Initialize arrays
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
ROCmaxcrz = zeros(Float64, N)
FFmaxcrz  = zeros(Float64, N)
EGTcrz = zeros(Float64, N)

iceil = N
crzmdotf = zeros(Float64, N)
crzEINOx = zeros(Float64, N)
crzFAR   = zeros(Float64, N)
clmbEINOx = zeros(Float64, N)
crzTAS   = zeros(Float64, N)
Tt4crz   = zeros(Float64, N)
Tt4crzmax   = zeros(Float64, N)

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

# Determine guesses for start and end of climb weights from fracW
# at end of climb of design mission, but fuel weights adjusted 
# by Wfrac
Wclimb1 = Wfrac*parg[igWMTO]     #Wzero + Wfrac*parg[igWfuel]
Wclimbn = Wfrac*parg[igWMTO]*0.9 #Wzero + Wfrac*parg[igWfuel]*para[iafracW, ipclimbn]
@inbounds for  ip = 1:N
    frac = (alts[ip] - 0.0)/(alts[N] - 0.0)
    V  =  VTO*(1.0-frac) + V0s[end]*frac
    Re = ReTO*(1.0-frac) + Reunits[end]*frac
    W  = Wclimb1*(1.0 - frac) + Wclimbn*frac
    V0s[ip] = V
    Ws[ip] = W
    Reunits[ip] = Re
end
# println(Ws/9.81/1000)
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


println(@sprintf("%12s %12s %12s %12s %12s %12s %12s", "FL", "Tt4max", "Tmetmax", "FFmax", "Tt4cruise", "Tmetcruise", "FFcruise"))
# integrate trajectory over climb
for   i = 1:N
    
    # velocity calculation from CL, Weight, altitude
    W  = Wfrac
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

        ρ0 = pare[ierho0, ip]
        u0 = pare[ieu0  , ip] 
        Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
        Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines

        # ifirst = false
        if NPSS_PT
            NPSS_success, Ftotal, η, P, Hrej, heatexcess, 
            mdotf[i], deNOx_, EINOx1, EINOx2, FAR, Tt3, OPR,
            Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, alts[i], Mach, 0.0, Tt4s[i], 
                Kinl, Φinl, 0.0, 0.0, ifirst, parg, parpt, pare, iptest)
        else
        Ftotal, η, P, Hrej, heatexcess,
        mdotf[i], BSFC,
        deNOx[i], EGT, Tt3, W3, EINOx1, EINOx2, FAR = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt4s[i],
                                    Kinl, Φinl, parpt, parmot, pargen, ifirst, false)
        end
        ifirst = false
        clmbEINOx[i] = EINOx2
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
    if ROC[i] <= 100
        iceil = i
        println("Climb ceiling reached at FL",FL[i])
        println("ROC = $(ROC[i]), Tt4 = $(Tt4s[i]), M = $Mach, W = $BW , CL = $CL, Ftotal = $Ftotal")
        ROC[i:end] .= 0.0
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

    if FL[i]≥ 270 && FL[i]≤431
        ip = ipcruise1

        W = Wfrac*parg[igWMTO]
        Wf = W - Wzero


        rfuel = Wf/parg[igWfuel]
        itrim = 1
        balance(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, itrim)
        icdfun = 1
        cdsum!(pari, parg, view(para, :, ip), view(pare, :, ip), icdfun)
        DoL = para[iaCD, ip]/ para[iaCL, ip]
        # println(1/DoL)
        BW = W + Wbouys[i]
        F  = BW*(DoL) #zero climb angle

        # println("F = $F, BW = $BW, W = $W, Wb = $(Wbouys[i])")

        Mach = 0.8
        V = Mach*a0s[i]
        ρ0 = ρ0s[i]
        u0 = Mach*a0s[i]
        Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
        Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines
        TR = 9.0 # Tt4/ Tamb
        Tt4 = max(2850, min(Tt4s[i], T0s[i]*TR*18/10)) # convert to [R]

        # Calculate max possible thrust at this Tt4, alt and Mach
        # Tt4 = Tt4s[i]
        Tt4crzmax[i] = Tt4
        if NPSS_PT
            NPSS_success, Ftotal, η, P, Hrej, heatexcess, 
            FFmaxcrz[i], deNOx, EINOx1, EINOx2, FAR, Tt3, OPR,
            Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, alts[i], Mach, 0.0, Tt4, 
                Kinl, Φinl, 0.0, 0.0, ifirst, parg, parpt, pare, iptest)
            Tmetmax = pare[ieTmet1, iptest]
        else
        Ftotal, η, P, Hrej, heatexcess,
        FFmaxcrz[i], BSFC,
        deNOxcrz, EGT, Tt3, W3, EINOx1, EINOx2, FAR = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt4 ,
                                        Kinl, Φinl, parpt, parmot, pargen, ifirst, Ldebug)
        end

        gam = Ftotal/BW - DoL
        Tt4crzmax[i] = Tt4
        ROCmaxcrz[i] = sin(gam)*V*60/0.3048 # this is the max possible ROC at this weight, alt etc.
        if Ftotal <F
            println("Ehhhhh Ftotal = ", Ftotal, "Freq = ",F)
        end

        if NPSS_PT
            NPSS_success, Ftotal, η, P, Hrej, heatexcess, 
            crzmdotf[i], deNOx, EINOx1, crzEINOx[i], crzFAR[i], Tt3, OPR,
            Wc3, Tt4crz[i], EGT = NPSS_TEsysOD(NPSS, alts[i], Mach, F, 0.0, 
                Kinl, Φinl, 0.0, 0.0, ifirst, parg, parpt, pare, iptest)
            Tmetcrz= pare[ieTmet1, iptest]
        else
            iter = 1
            itermax = 20
            for iter = 1:itermax
                if abs(Ftotal - F)<1
                    break
                end
                ΔF = F - Ftotal
                Tt4 = Tt4*(1 + ΔF/Ftotal/5) # 5 is just a scale factor so you don't get random oscillations
                # println(Tt4)
                Tt4 = max(2850, min(Tt4, Tt4s[i]))
                Tt4crz[i] = Tt4
                # println("Adjusted --> ",Tt4)
                # TR = TR*0.99
                # Tt4 = max(2850, min(Tt4s[i], T0s[i]*TR*18/10)) # convert to [R]
                # println(alts[i], " ", Tt4)
                Ftotal, η, P, Hrej, heatexcess,
                crzmdotf[i], BSFC,
                deNOxcrz, EGT, Tt3, W3, EINOx1, EINOx2, FAR = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt4 ,
                Kinl, Φinl, parpt, parmot, pargen, ifirst, Ldebug)

                crzEINOx[i] = EINOx2
                crzFAR[i]   = FAR

                gam = Ftotal/BW - DoL
                # println(Wfrac, "--> gam = ", gam, " ROC = ", sin(gam)*V*60/0.3048, " FFsteadylevel = ", crzmdotf[i])
                if Tt4 ≤ 2850 && Ftotal>F    
                    println("gam = ", gam, "ROC = ", sin(gam)*V*60/0.3048)
                    gam = 0.0
                    println("Min Temp reached")
                end
            end
            if iter == itermax && abs(Ftotal - F)>1
                println("Steady level cruise not converged - increase itermax?")
            end
        end

        EGTcrz[i] = EGT

        if ROCmaxcrz[i]<=0
            ROCmaxcrz[i] = 0.0
        end
        crzTAS[i] = V0s[i]/kts_to_mps

        println(@sprintf("%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f", 
        FL[i], Tt4crzmax[i], Tmetmax, FFmaxcrz[i], Tt4crz[i], Tmetcrz, crzmdotf[i]))
    end

end # done integrating climb


return Ws[1], alts[iceil], V0s, ROC, mdotf, crzmdotf, crzTAS, EGTcrz, FFmaxcrz, ROCmaxcrz, Tt4crz, Tt4crzmax, crzEINOx, clmbEINOx, crzFAR
end


function get_cruisespeed(h, Mcr)
    h_ft = h/0.3048
    htrans_ft = Hptrans(280/1.944, 0.8)
    T, P, ρ,  a = atmos(h/1000)

    if 6000<=h_ft<=13999
        CAS = 250/1.944 #Assume 250 kts CAS cruise consistent with BADA
        TAS =  CAS_TAS(CAS, h)
    elseif 14000<=h_ft< htrans_ft
        CAS = 280/1.944  #Assume 280 kts CAS cruise similar to BADA
        TAS =  CAS_TAS(CAS, h)
    elseif h_ft>= htrans_ft
        TAS =  MNcr * a
        CAS =  TAS_CAS(TAS, h)
    end
    
     MN =  TAS/ a
     return TAS, MN
end

function Hptrans(CAS, MNcr)

    gam = gamSL
    gmi = gam - 1
    βT = -0.0065
    R = 287.05287

    δtrans = ((1+0.5*gmi*(CAS/aSL)^2)^(gam/gmi) - 1)/((1+0.5*gmi*MNcr^2)^(gam/gmi) - 1)
    θtrans = δtrans^(-βT*R/gee)

    htrans_ft = 1000/0.3048/6.5 * TSL*(1-θtrans)

    return htrans_ft
end

function CAS_TAS(CAS, h)
    T, P, ρ,  a = atmos(h/1000) 
    gam = gamSL
    gmi = gam - 1
    k = gmi/gam

    temp1 = (1 + k/2*ρ0/PSL*CAS^2)^(1/k)
    temp2 = (1 + PSL/P * (temp1 - 1))^k
    TAS   = (2/k *P/ρ * (temp2 - 1))^0.5
    return TAS
end

function TAS_CAS(TAS, h)
    T, P, ρ,  a = atmos(h/1000) 
    gam = gamSL
    gmi = gam - 1
    k = gmi/gam

    temp1 = (1 + k/2*ρ/P*TAS^2)^(1/k)
    temp2 = (1 + P/PSL * (temp1 - 1))^k
    CAS   = (2/k *PSL/ρ0 * (temp2 - 1))^0.5
    return CAS
end
