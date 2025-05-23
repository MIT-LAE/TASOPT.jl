"""
`odperf!` runs the aircraft in off-design mode to generate a `BADA`-like 
`PTF` file for use in `AEIC`.

!!! compat "Future Changes"
    This function will be overhauled and renamed in an upcoming revision. Neither NPSS nor turboelectric compatibility are currently in the scope.

"""
function odperf!(pari, parg, parm, para, pare, Wfrac, FL, 
    NPSS_TS::Base.Process, 
    NPSS_Fan::Base.Process, 
    NPSS_AftFan::Base.Process, Ldebug, ifirst, NPSS_PT, NPSS::Base.Process)

@warn "The function `odperf!` will be overhauled and renamed in an upcoming revision. Neither NPSS nor turboelectric compatibility are currently in the scope."

calc_ipc1 = true
# ifirst = true
# Ldebug = true

itergmax::Int64 = 20
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
γs = ones(Float64, N)*0.015
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
desTAS   = zeros(Float64, N)
Tt4crz   = zeros(Float64, N)
Tt4crzmax   = zeros(Float64, N)

@inbounds for i =1:N
    T0s[i], p0s[i], ρ0s[i], a0s[i], μ0s[i] = atmos(alts[i]/1000)
    rhocab = max( parg[igpcabin] , p0s[i] ) / (RSL*TSL) #Should be T0s to be more accurate?
    Wbouys[i] = (rhocab-ρ0s[i])*gee*parg[igcabVol]
end

MNcr = para[iaMach, ipcruise1]
Tt4max = maximum(pare[ieTt4, :])
println(@sprintf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s", 
"FL", "TAS", "CAS", "Mach", "Fn", "L/D", "Tt4max", "Tmetmax", "FFmax", "Tt4cruise", "Tmetcruise", "FFcruise", "CL", "CLh"))
# integrate trajectory over climb
for   i = 1:N
    
    # velocity calculation from CL, Weight, altitude
    W  = Wfrac*parg[igWMTO]

    #Get flight speed from climb schedule
    CAScl, TAScl, Mcl = get_climbspeed(alts[i], MNcr)
    # println(@sprintf("%.2f, %.2f, %.2f",CAScl, TAScl, Mcl))
    CASdes, TASdes, Mdes = get_descentspeed(alts[i], MNcr)

    if i == 1 #assume climb CL below 
        ip = iptest # dummy location to store values
        CL = para[iaCL, ipclimb1]
        para[:, ip] .= para[:, ipclimb1] #initialize iptest location
        pare[:, ip] .= pare[:, ipclimb1] #initialize iptest location

    else
        ip = iptest # dummy location to store values
        CL = para[iaCL, ipcruise1]
        para[:, ip] .= para[:, ipcruise1]
        pare[:, ip] .= pare[:, ipcruise1]
    end

    ρ  = ρ0s[i]
    μ  = μ0s[i]
    Vsound = a0s[i]
    gamV = 0.0
    cosg = cos(gamV)
    γs[i] = gamV
    BW   = W + Wbouys[i]

    dgamV  = 1.0
    Ftotal = 0.0
    DoL    = 0.0
    V      = 0.0

    if (Ldebug)
        printstyled(@sprintf("\t%5s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                        "iterg", "dgamV", "gamV", "cosg", "CL", "BW", "Tt4", "Ftotal", "DoL", "V", "Alt", "M", "Tt4"); color = :light_green)
    end
    for  iterg = 1:itergmax

        q = 0.5*ρ*TAScl^2
        CL = BW*cosg/(q*S)
        para[iaCL, ip] = CL

        V = TAScl
        Mach = V/Vsound
        para[iaMach, ip] = Mach

        V0s[i] = V # Save climb TAS (m/s) to V0s
        desTAS[i] = TASdes/kts_to_mps # Save descent TAS separately, convert to knots
        M0s[i] = Mach
        Reunits[i] = V*ρ/μ
        
        para[iaMach, ip] = Mach
        para[iaReunit, ip] = V*ρ/μ
        pare[ieu0, ip] = V
        pare[ieM0, ip] = Mach

        # Set pitch trim by adjusting CLh
        Wf = W - Wzero
        rfuel = Wf/parg[igWfuel]*0
        opt_trim_var = "CL_htail"
        balance_aircraft!(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

        
        # Calculate Drag
        if (i == 1)
            computes_wing_direct = false
        else 
            computes_wing_direct = true
        end
        if CL > 0.9
            computes_wing_direct = false
        end
        aircraft_drag!(pari, parg, view(para, :, ip), view(pare, :, ip), computes_wing_direct)

        #BLI parameters
        ρ0 = pare[ierho0, ip]
        u0 = pare[ieu0  , ip] 
        Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
        Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines

        # ifirst = false
        #Run propsys to MAX THRUST
        TR = 11.5
        Tt4s[i] =  max(2850, min(Tt4max, T0s[i]*TR*18/10)) # convert to [R]
        if NPSS_PT
            NPSS_success, Ftotal, η, P, Hrej, heatexcess, 
            mdotf[i], deNOx_, EINOx1, EINOx2, FAR, Tt3, OPR,
            Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, alts[i], Mach, 0.0, Tt4max, 
                Kinl, Φinl, 0.0, 0.0, ifirst, parg, parpt, pare, iptest)
        else
            Ftotal, η, P, Hrej, heatexcess,
            mdotf[i], BSFC,
            deNOx[i], EGT, Tt3, W3, EINOx1, EINOx2, FAR = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt4s[i],
                                        Kinl, Φinl, parpt, parmot, pargen, ifirst, false)
        end
        ifirst = false
        clmbEINOx[i] = EINOx2
        DoL = para[iaCD, ip]/ para[iaCL, ip]

        # Calculate improved flight angle
        ϕ = Ftotal/BW
        sing = (ϕ - DoL*sqrt(1.0 - ϕ^2 + DoL^2))/(1.0 + DoL^2)
        gamV = asin(sing)
        cosg = sqrt(1 - sing)
        ROC[i] = sing*V/ft_to_m*60 # Rate of climb in m/s
        
        dgamV = gamV - γs[i]
        γs[i] = gamV
        
        para[iagamV, ip] = gamV

        if (Ldebug)
            printstyled(@sprintf("\t%5d  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4f  %9.4f  %9.4f\n",
                        iterg, abs(dgamV), gamV*180/π, cosg, CL, BW, Tt4max, Ftotal, DoL, V, alts[i], Mach, Tt4s[i]); color =:light_green)
        end
        
        if(abs(dgamV) < gamVtol) 
            break
        end
    end

    if ROC[i] <= 00
        iceil = i
        println("Climb ceiling reached at FL",FL[i])
        println("ROC = $(ROC[i]), Tt4 = $(Tt4s[i]), M = $Mach, W = $BW , CL = $CL, Ftotal = $Ftotal")
        ROC[i:end] .= 0.0
        break
    end
    if (abs(dgamV) > gamVtol) 
            println("Climb gamV not converged. dgamV = $dgamV")
    end
    pare[ieFe, ip] = Ftotal
  
    # Cruise Section
    if FL[i]≥ 30 && FL[i]≤431
        ip = iptest
        #Get flight speed from climb schedule
        CAScr, TAScr, Mcr = get_cruisespeed(alts[i], MNcr)
        # Mcr = para[iaMach, ipcruise1]
        # TAScr = Mcr*a0s[i]
        # CAScr = TAS_CAS(TAScr, alts[i])
        crzTAS[i] = TAScr/kts_to_mps
        
        BW = W + Wbouys[i]
        q = 0.5*ρ0s[i]*TAScr^2
        CL = BW/(q*S)
        para[iaCL, ip] = CL

        V = TAScr
        # Commented out: 
        # Because climb/descent TAS has been overwritten by cruise TAS if run for same altitude. "crzTAS" is separately saved and sent to printBADA
        # V0s[i] = TAScr 
        Mach = Mcr
        para[iaMach, ip] = Mach
        para[iaReunit, ip] = V*ρ/μ
        pare[ieu0, ip] = V
        pare[ieM0, ip] = Mach

        Wf = W - Wzero
        rfuel = Wf/parg[igWfuel]*0

        #Trim aircraft
        opt_trim_var = "CL_htail"
        balance_aircraft!(pari, parg, view(para, :, ip), rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)
        computes_wing_direct = true
        if CL > 1.0
            println("CL during cruise is $CL")
            computes_wing_direct = false
        end
        #Get Drag
        aircraft_drag!(pari, parg, view(para, :, ip), view(pare, :, ip), computes_wing_direct)
        DoL = para[iaCD, ip]/ para[iaCL, ip]

        F  = BW*(DoL) #zero climb angle for cruise
        # println("F = $F, BW = $BW, W = $W, Wb = $(Wbouys[i]), Wf = $(Wf)")


        ρ0 = ρ0s[i]
        u0 = Mach*a0s[i]
        Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
        Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines
        TR = 11.0 # Tt4/ Tamb
        Tt4 = max(2850, min(Tt4max, T0s[i]*TR*18/10)) # convert to [R]

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
        ROCmaxcrz[i] = sin(gam)*V*60/0.3048 # this is the max possible ROC at this weight, alt and cruise speed etc.
        if Ftotal <F
            println("Ehhhhh Ftotal = ", Ftotal, "Freq = ",F)
            break
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

        println(@sprintf("%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f", 
        FL[i], TAScr/kts_to_mps, CAScr/kts_to_mps, Mcr, F, 1/DoL, Tt4crzmax[i], Tmetmax, FFmaxcrz[i], Tt4crz[i], Tmetcrz, crzmdotf[i], CL, para[iaCLh, ip]))
    end #cruise section done

end #outer loop

return Wfrac*parg[igWMTO], alts[iceil], V0s, desTAS, ROC, mdotf, crzmdotf, crzTAS, EGTcrz, FFmaxcrz, ROCmaxcrz, Tt4crz, Tt4crzmax, crzEINOx, clmbEINOx, crzFAR
end

function get_cruisespeed(h, MNcr)
    h_ft = h/0.3048
    Vcr2 = 280/1.944 # From BADA 737__.APF
    htrans_ft = Hptrans(Vcr2, 0.8)
    T, P, ρ,  a = atmos(h/1000)
    TAS = MNcr*a
    if h_ft<3000
        CAS = 170/1.944
        TAS = CAS_TAS(CAS, h)
    elseif 3000<=h_ft<6000
        CAS = 220/1.944
        TAS = CAS_TAS(CAS, h)
    elseif 6000<=h_ft<14000
        CAS = 250/1.944 #Assume 250 kts CAS cruise consistent with BADA
        TAS =  CAS_TAS(CAS, h)
    elseif 14000<=h_ft< htrans_ft
        CAS = Vcr2  #Assume CAS cruise similar to BADA
        TAS =  CAS_TAS(CAS, h)
    elseif h_ft>= htrans_ft
        TAS =  MNcr * a
        CAS =  TAS_CAS(TAS, h)
    end
    
    MN =  TAS/ a
    return CAS, TAS, MN
end

function get_climbspeed(h,MNcr)
    h_ft = h/0.3048
    T, P, ρ,  a = atmos(h/1000)
    TAS = a*MNcr
    Vcl1 = 300/1.944 # From BADA 737__.APF
    Vcl2 = 300/1.944 # From BADA 737__.APF
    htrans_ft = Hptrans(Vcl2/1.944, 0.8)

    VstallTO = 117/1.944 # From 738__.PTF file BADA, flap 5 (Standard Takeoff Setting for B738)
    
    CVmin = 1.3 # From BADA manual
    VdCL1, VdCL2, VdCL3, VdCL4, VdCL5 = [5, 10, 30, 60, 80]/1.944

    if 0<=h_ft<1500
        CAS = CVmin*VstallTO + VdCL1
        TAS = CAS_TAS(CAS, h)
    elseif 1500<=h_ft<3000
        CAS = CVmin*VstallTO + VdCL2
        TAS = CAS_TAS(CAS, h)
    elseif 3000<=h_ft<4000
        CAS = CVmin*VstallTO + VdCL3
        TAS = CAS_TAS(CAS, h)
    elseif 4000<=h_ft<5000
        CAS = CVmin*VstallTO + VdCL4
        TAS = CAS_TAS(CAS, h)
    elseif 5000<=h_ft<6000
        CAS = CVmin*VstallTO + VdCL5
        TAS = CAS_TAS(CAS, h)
    elseif 6000<=h_ft<10000 
        CAS = min(Vcl1, 250/1.944)
        TAS = CAS_TAS(CAS, h)
    elseif 10000<=h_ft<htrans_ft
        CAS = Vcl2
        TAS = CAS_TAS(CAS, h)
    elseif h_ft>=htrans_ft
        TAS = MNcr*a
        CAS = TAS_CAS(TAS, h)
    end

    MN = TAS/a
    #Limit MN to be less than MNcr at design
    if MN>MNcr
        MN = MNcr
        TAS = MN*a
    end

    return CAS, TAS, MN
end

function get_descentspeed(h,MNcr)
    h_ft = h/0.3048
    T, P, ρ,  a = atmos(h/1000)
    TAS = a*MNcr
    Vdes1 = 290/1.944 # From BADA 737__.APF
    Vdes2 = 290/1.944 # From BADA 737__.APF
    htrans_ft = Hptrans(Vdes2/1.944, 0.8)

    VstallLD = 107/1.944 # From 738__.PTF file BADA, flap 30 (Standard Landing Setting for B738)
    
    CVmin = 1.3 # From BADA manual
    VdCL1, VdCL2, VdCL3, VdCL4 = [5, 10, 20, 50]/1.944

    if 0<=h_ft<1000
        CAS = CVmin*VstallLD + VdCL1
        TAS = CAS_TAS(CAS, h)
    elseif 1000<=h_ft<1500
        CAS = CVmin*VstallLD + VdCL2
        TAS = CAS_TAS(CAS, h)
    elseif 1500<=h_ft<2000
        CAS = CVmin*VstallLD + VdCL3
        TAS = CAS_TAS(CAS, h)
    elseif 2000<=h_ft<3000
        CAS = CVmin*VstallLD + VdCL4
        TAS = CAS_TAS(CAS, h)
    elseif 3000<=h_ft<6000
        CAS = min(Vdes1, 220/1.944)
        TAS = CAS_TAS(CAS, h)
    elseif 6000<=h_ft<10000 
        CAS = min(Vdes1, 250/1.944)
        TAS = CAS_TAS(CAS, h)
    elseif 10000<=h_ft<htrans_ft
        CAS = Vdes2
        TAS = CAS_TAS(CAS, h)
    elseif h_ft>=htrans_ft
        TAS = MNcr*a
        CAS = TAS_CAS(TAS, h)
    end

    MN = TAS/a
    #Limit MN to be less than MNcr at design
    if MN>MNcr
        MN = MNcr
        TAS = MN*a
    end

    return CAS, TAS, MN
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

    temp1 = (1 + k/2*ρSL/pSL*CAS^2)^(1/k)
    temp2 = (1 + pSL/P * (temp1 - 1))^k
    TAS   = (2/k *P/ρ * (temp2 - 1))^0.5
    return TAS
end

function TAS_CAS(TAS, h)
    T, P, ρ,  a = atmos(h/1000) 
    gam = gamSL
    gmi = gam - 1
    k = gmi/gam

    temp1 = (1 + k/2*ρ/P*TAS^2)^(1/k)
    temp2 = (1 + P/pSL * (temp1 - 1))^k
    CAS   = (2/k *pSL/ρSL * (temp2 - 1))^0.5
    return CAS
end

using Plots 
function show_ff_sens(ff, H)
    # Normalize `ff` by dividing each column by its minimum
    ffnorm = ff' ./ minimum.(eachcol(ff))
    
    # Create the plot
    plot(
        ffnorm', H,
        xlabel = "Normalized ff",
        ylabel = "H",
        legend = false,
        grid = true
    )
end