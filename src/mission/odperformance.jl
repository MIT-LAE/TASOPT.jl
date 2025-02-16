function odperf!(ac, W, FL, Ldebug)
    calc_ipc1 = true
    # ifirst = true
    # Ldebug = true

    itergmax::Int64 = 50
    gamVtol  = 1.0e-12
    method = "cubic"
    imission = 1

    # BLI
        fBLIf   = ac.parg[igfBLIf]
        DAfsurf = ac.parad[iaDAfsurf, ipcruise1]
        KAfTE   = ac.parad[iaKAfTE, ipcruise1]

    # Zero-fuel weight for this mission
        Wzero   = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay] + ac.parmd[imWpay]  #This ensures that this will work for multi-mission in the future
        Wempty  = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay]
    # Payload fraction for this mission
        rpay = ac.parmd[imWpay]/ac.parg[igWpay]
        ξpay = 0.

    S    = ac.wing.layout.S
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
    Tt41s = zeros(Float64, N)

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
    Tt41crz   = zeros(Float64, N)
    Tt41crzmax   = zeros(Float64, N)

    @inbounds for i =1:N
        T0s[i], p0s[i], ρ0s[i], a0s[i], μ0s[i] = atmos(alts[i]/1000)
        rhocab = max( ac.parg[igpcabin] , p0s[i] ) / (RSL*TSL) #Should be T0s to be more accurate?
        Wbouys[i] = (rhocab-ρ0s[i])*gee*ac.parg[igcabVol]
    end

    MNcr = ac.parad[iaMach, ipcruise1]
    Tt41max = maximum(ac.pared[ieTt41, :])  #TODO set to max temp
    Tt41CR = ac.pared[ieTt41,ipcruise1]
    println(@sprintf("%12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s", 
    "FL", "TAS", "CAS", "Mach", "Fn", "L/D", "Tt4max", "Tmetmax", "FFmax", "Tt4cruise", "Tmetcruise", "FFcruise", "CL", "CLh"))
    # integrate trajectory over climb
    for   i = 1:N
        
        # velocity calculation from CL, Weight, altitude
        #Get flight speed from climb schedule
        if ac.aircraft_type == "2" || lowercase(string(ac.aircraft_type)) == "wide body aircraft"

            CAScl, TAScl, Mcl = get_climbspeed_77W(alts[i], MNcr)
            #Get flight speed from descent schedule
            CASdes, TASdes, Mdes = get_descentspeed_77W(alts[i], MNcr)
        else
            CAScl, TAScl, Mcl = get_climbspeed(alts[i], MNcr)
            #Get flight speed from descent schedule
            CASdes, TASdes, Mdes = get_descentspeed(alts[i], MNcr)
        end
        # println(@sprintf("%.2f, %.2f, %.2f",CAScl, TAScl, Mcl))
        # println(@sprintf("%.2f, %.2f, %.2f",CASdes, TASdes, Mdes))

        if alts[i] <= 2000/3.28084 #assume climbout CL as below 
            ip = iptest # dummy location to store values
            CL = ac.parad[iaCL, ipclimb1]
            ac.pared[:, ip] .= ac.pared[:, ipclimb1]
            ac.parad[:, ip] .= ac.parad[:, ipclimb1]
            # ac.parad[:, ip] .= ac.parad[:, ipclimb1] #initialize iptest location
            # ac.pared[:, ip] .= ac.pared[:, ipclimb1] #initialize iptest location

        else
            ip = iptest # dummy location to store values
            CL = ac.parad[iaCL, ipcruise1] # To be used by cruise values (Climb CL will be recalc. and override this)
            # ac.parad[:, ip] .= ac.parad[:, ipcruise1] #initialize iptest location
            # ac.pared[:, ip] .= ac.pared[:, ipclimb1] #initialize iptest location
            ac.pared[:, ip] .= ac.pared[:, ipclimb1]
            ac.parad[:, ip] .= ac.parad[:, ipcruise1]
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
            printstyled(@sprintf("\t%5s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n",
                            "iterg", "dgamV", "gamV", "cosg", "CL", "BW", "Ftotal", "DoL", "V", "Alt", "M", "Tt41"); color = :light_green)
        end
        for  iterg = 1:itergmax
            
            q = 0.5*ρ*TAScl^2
            CL = BW*cosg/(q*S)

            ac.parad[iaCL, ip] = CL

            V = TAScl
            Mach = V/Vsound
            ac.parad[iaMach, ip] = Mach

            V0s[i] = V
            desTAS[i] = TASdes/kts_to_mps
            M0s[i] = Mach
            Reunits[i] = V*ρ/μ
            
            ac.parad[iaMach, ip] = Mach
            ac.parad[iaReunit, ip] = V*ρ/μ
            ac.pared[ieu0, ip] = V
            ac.pared[ieM0, ip] = Mach

            # Set pitch trim by adjusting CLh
            Wf = W - Wzero
            rfuel = Wf/ac.parg[igWfuel]*0
            itrim = 1
            balance(ac, imission, ip, rfuel, rpay, ξpay, itrim)
            
            # Calculate Drag
            if (i == 1)
                icdfun = 0
            else 
                icdfun = 1
            end
            if alts[i] < 5000/3.28084 # Under 5000ft climb case
                # Transient/unsteady with flap etc roattion phase (takeoff/ rotation/ don't try to calulate CD and just use what you have)
                # Maybe in the future we want to incorporate different CL depending on load and flap schedule, altitude etc.
                # HILOAD case obviously flap is extended until 4~5000 ft climb
                icdfun = 0
                if  CL < 0.875 # if it is a low enough load case with low enough CL, if only impose this, it oscillates if one iter gives CL = 0.9 and other give 0.85 for example. 
                    icdfun = 1 
                end
            end
            cdsum!(ac, imission, ip, icdfun)

            #BLI ac.parameters
            ρ0 = ac.pared[ierho0, ip]
            u0 = ac.pared[ieu0  , ip] 
            Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
            Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines

            # ifirst = false
            #Run propsys to MAX THRUST
            TR = 11.5
            OEW = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay]

            # Climb Tt41 schedule
            if ac.aircraft_type == "2" || lowercase(string(ac.aircraft_type)) == "wide body aircraft"
                if W > 0.95 * ac.parg[igWMTO]
                    Tt41frac = LinRange(0.965, 0.97, length(alts)-5) # length(alts)-"5" is to account for alt at or under 2000ft
                elseif W < 1.05*(1.2*OEW) # Low weight
                    Tt41frac = LinRange(0.86, 0.955, length(alts)-5) # length(alts)-"5" is to account for alt at or under 2000ft        
                else # mid-weight (reference mass)
                    Tt41frac = LinRange(0.925, 0.965, length(alts)-5) # length(alts)-"5" is to account for alt at or under 2000ft
                end
            else
                if W > 0.95 * ac.parg[igWMTO]
                    Tt41frac = LinRange(0.976,0.976, length(alts)-5) # length(alts)-"5" is to account for alt at or under 2000ft
                elseif W < 1.05*(1.2*OEW) # Low weight
                    Tt41frac = LinRange(0.92,0.963, length(alts)-5) # length(alts)-"5" is to account for alt at or under 2000ft        
                else # mid-weight (reference mass)
                    Tt41frac = LinRange(0.94,0.975, length(alts)-5) # length(alts)-"5" is to account for alt at or under 2000ft
                end
            end

            if alts[i] <= (2000/3.28084) * 1.1 # At or under 2000 [ft]
                Tt41s[i] = Tt41frac[1]*Tt41max # convert to [R]

            else # climb thrust
                Tt41s[i] = Tt41frac[i-5]*Tt41max # convert to [R], Tt41frac[i-5] for same reason as above
            end

            # if NPSS_PT
                # NPSS_success, Ftotal, heatexcess, 
                # mdotf[i], EINOx1, FAR, Mtip, Tblade, Tt3, OPR, BPR,
                # Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, alts[i], Mach, 0.0, Tt41s[i], mofft, Pofft, 
                #     ifirst, ac.parg, parpt, ac.pare, iptest)
            # else
            #     Ftotal, η, P, Hrej, heatexcess,
            #     mdotf[i], BSFC,
            #     deNOx[i], EGT, Tt3, W3, EINOx1, EINOx2, FAR = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt41s[i],
            #                                 Kinl, Φinl, parpt, ac.parmot, ac.pargen, ifirst, false)
            # end

            # NON NPSS
            ip = iptest

            # ac.pared[:, ip] = ac.pared[:, ipstatic]
            # ac.parad[:, ip] = ac.parad[:, ipstatic]
            # alts[i]
            # Tt41s[i]

            icall = 1
            icool = 1
            initeng = false
            ac.pared[ieTt41,ip] = Tt41s[i]
            ac.parad[iaalt,ip] = alts[i]

            TASOPT.tfcalc!(ac.pari, ac.parg, view(ac.parad, :, ip), 
                                        view(ac.pared, :, ip), ac.wing, ip, icall, icool, initeng)

            mdotf[i] = ac.pared[ieff, iptest] * ac.pared[iemcore]     
            EINOx1 = EINOx(ac, iptest; method=method)
            Ftotal =  ac.pared[ieFe, iptest] * ac.parg[igneng]
            ifirst = false
            clmbEINOx[i] = EINOx1
            DoL = ac.parad[iaCD, ip]/ ac.parad[iaCL, ip]

            # Calculate improved flight angle
            ϕ = Ftotal/BW
            sing = (ϕ - DoL*sqrt(1.0 - ϕ^2 + DoL^2))/(1.0 + DoL^2)
            gamV = asin(sing)
            cosg = sqrt(1 - sing)
            ROC[i] = sing*V/ft_to_m*60 # Rate of climb in m/s
            
            dgamV = gamV - γs[i]
            γs[i] = gamV
            
            ac.parad[iagamV, ip] = gamV

            if (Ldebug)
                printstyled(@sprintf("\t%5d  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4e  %9.4f  %9.4f  %9.4f\n",
                            iterg, abs(dgamV), gamV*180/π, cosg, CL, BW, Tt41max, Ftotal, DoL, V, alts[i], Mach, Tt41s[i]); color =:light_green)
            end
            
            if(abs(dgamV) < gamVtol) 
                break
            end
        end

        if ROC[i] <= 00
            iceil = i
            println("Climb ceiling reached at FL",FL[i])
            println("ROC = $(ROC[i]), Tt41 = $(Tt41s[i]), M = $Mach, W = $BW , CL = $CL, Ftotal = $Ftotal") #CD, total drag
            ROC[i:end] .= 0.0
            break
        end
        if (abs(dgamV) > gamVtol) 
            println("Climb gamV not converged. dgamV = $dgamV")
            # break
        end
        ac.pared[ieFe, ip] = Ftotal
    
        # Cruise Section
        FLcrzmax = 410
        if ac.aircraft_type == "2" || lowercase(string(ac.aircraft_type)) == "wide body aircraft"
            FLcrzmax = 431
        end
        if FL[i]≥ 60 && FL[i]≤FLcrzmax
            ip = iptest
            #Get flight speed from climb schedule
            if ac.aircraft_type == "2" || lowercase(string(ac.aircraft_type)) == "wide body aircraft"
                CAScr, TAScr, Mcr = get_cruisespeed_77W(alts[i], MNcr)
            else
                CAScr, TAScr, Mcr = get_cruisespeed(alts[i], MNcr)
            end
            # Mcr = ac.parad[iaMach, ipcruise1]
            # TAScr = Mcr*a0s[i]
            # CAScr = TAS_CAS(TAScr, alts[i])
            crzTAS[i] = TAScr/kts_to_mps
            
            BW = W + Wbouys[i]
            q = 0.5*ρ0s[i]*TAScr^2
            CL = BW/(q*S)
            ac.parad[iaCL, ip] = CL

            V = TAScr
            # Commented out: 
            # Because climb/descent TAS has been overwritten by cruise TAS if run for same altitude. "crzTAS" is separately saved and sent to printBADA
            # V0s[i] = TAScr 
            Mach = Mcr
            ac.parad[iaMach, ip] = Mach
            ac.parad[iaReunit, ip] = V*ρ/μ
            ac.pared[ieu0, ip] = V
            ac.pared[ieM0, ip] = Mach

            Wf = W - Wzero
            rfuel = Wf/ac.parg[igWfuel]*0

            #Trim aircraft
            itrim = 1
            balance(ac, imission, ip, rfuel, rpay, ξpay, itrim)
            icdfun = 1

            # printstyled(@sprintf("%9.4e  %9.4e n", FL[i], crzTAS[i]); color =:light_red)

            if (FL[i] < 50 && CL > 0.8)
                # NEED SIMILAR THING AS CLIMB FOR HIGH LOAD CASE CRUISE, ex 3000ft CRZ at 230 kts cannot be flown clean wing
                # FOR NOW, REVERT TO WHAT YOU HAVE

                printstyled(@sprintf("%9.4e  %9.4e \n", CL, BW); color =:light_red)
                println("CL during cruise is $CL")

                icdfun = 0
            end
            #Get Drag
            cdsum!(ac, imission, ip, icdfun)
            
            DoL = ac.parad[iaCD, ip]/ ac.parad[iaCL, ip]

            F  = BW*(DoL) #zero climb angle for cruise
            # println("F = $F, BW = $BW, W = $W, Wb = $(Wbouys[i]), Wf = $(Wf)")


            ρ0 = ρ0s[i]
            u0 = Mach*a0s[i]
            Φinl = 0.5*ρ0*u0^3 * (DAfsurf*fBLIf)/2.0 
            Kinl = 0.5*ρ0*u0^3 * (KAfTE  *fBLIf)/2.0 # Assume 2 engines

            TR = 11.0 # Tt41/ Tamb
            Tt41 = max(Tt41CR, min(Tt41max, T0s[i]*TR)) 

            # Calculate max possible thrust at this Tt4, alt and Mach
            # Tt4 = Tt4s[i]
            Tt41crzmax[i] = Tt41
            # if NPSS_PT
            #     NPSS_success, Ftotal, heatexcess, 
            #     FFmaxcrz[i], EINOx1, FAR, Mtip, Tblade, Tt3, OPR, BPR,
            #     Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, alts[i], Mach, 0.0, Tt41s[i], mofft, Pofft, 
            #         ifirst, ac.parg, parpt, ac.pare, iptest)
            #     Tmetmax = ac.pared[ieTmet1, iptest]

            # else
            #     Ftotal, η, P, Hrej, heatexcess,
            #     FFmaxcrz[i], BSFC,
            #     deNOxcrz, EGT, Tt3, W3, EINOx1, EINOx2, FAR = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt41 ,
            #                                     Kinl, Φinl, parpt, ac.parmot, ac.pargen, ifirst, Ldebug)
            # end

            icall = 1
            icool = 1
            initeng = false
            ac.pared[ieTt41,ip] = Tt41s[i]
            ac.parad[iaalt,ip] = alts[i]

            TASOPT.tfcalc!(ac.pari, ac.parg, view(ac.parad, :, ip), 
                                        view(ac.pared, :, ip), ac.wing, ip, icall, icool, initeng)

            Ftotal = ac.pared[ieFe, iptest] * ac.parg[igneng]
            FFmaxcrz[i] = ac.pared[ieff, iptest] * ac.pared[iemcore]     
            EINOx1 = EINOx(ac, iptest; method=method)
            Tmetmax = ac.pared[ieTmet1, iptest]
            Tt41 = ac.pared[ieTt41, iptest]

            gam = Ftotal/BW - DoL
            Tt41crzmax[i] = Tt41
            ROCmaxcrz[i] = sin(gam)*V*60/0.3048 # this is the max possible ROC at this weight, alt and cruise speed etc.
            if Ftotal <F
                println("Ehhhhh Ftotal = ", Ftotal, "Freq = ",F)
                break
            end

            

            icall = 2
            icool = 1
            initeng = false
            ac.pared[ieFe, ip] = F
            TASOPT.tfcalc!(ac.pari, ac.parg, view(ac.parad, :, ip), 
                                        view(ac.pared, :, ip), ac.wing, ip, icall, icool, initeng)
            
            Ftotal = ac.pared[ieFe, iptest] * ac.parg[igneng]
            crzmdotf[i] = ac.pared[ieff, iptest] * ac.pared[iemcore]     
            EINOx1 = EINOx(ac, iptest; method=method)
            FAR = ac.pared[ieff, iptest]
            Tt3   = ac.pared[ieTt3, iptest]
            Tt41 = ac.pared[ieTt41, iptest]
            Tmetcrz= ac.pared[ieTmet1, iptest]
            crzEINOx[i] = EINOx1
            crzFAR[i]   = FAR
            # Velocity to mach
            _, _, _, a = atmos(alts[i]/1000)
            mach_num = ac.pared[ieu9, iptest]/a
            # Tt9 to T9 
            Tt9oT9 = 1 + (((1.4-1)/2)*mach_num^2)
            T9 = Tt9oT9/ac.pared[ieTt9, iptest]
            EGT = T9
            EGTcrz[i] = EGT
            # if NPSS_PT
            #     NPSS_success, Ftotal, heatexcess, 
            #     crzmdotf[i], EINOx1, FAR, Mtip, Tblade, Tt3, OPR, BPR,
            #     Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, alts[i], Mach, F, Tt41s[i], mofft, Pofft, 
            #         ifirst, ac.parg, parpt, ac.pare, iptest)
            #     Tmetcrz= ac.pared[ieTmet1, iptest]
            #     crzEINOx[i] = EINOx1
            #     crzFAR[i]   = FAR
            # else
            #     iter = 1
            #     itermax = 100
            #     for iter = 1:itermax
            #         if abs(Ftotal - F)<1
            #             break
            #         end
            #         ΔF = F - Ftotal
            #         Tt41 = Tt41*(1 + ΔF/Ftotal/5) # 5 is just a scale factor so you don't get random oscillations
            #         # println(Tt4)
            #         Tt41 = max(Tt41CR, min(Tt41, Tt41s[i]))
            #         Tt41crz[i] = Tt41
            #         # println("Adjusted --> ",Tt41)
            #         # TR = TR*0.99
            #         # Tt41 = max(Tt41CR, min(Tt41s[i], T0s[i]*TR*18/10)) # convert to [R]
            #         # println(alts[i], " ", Tt41)
            #         Ftotal, η, P, Hrej, heatexcess,
            #         crzmdotf[i], BSFC,
            #         deNOxcrz, EGT, Tt3, W3, EINOx1, EINOx2, FAR = PowerTrainOD(NPSS_TS, NPSS_Fan, NPSS_AftFan, alts[i], Mach, Tt41 ,
            #         Kinl, Φinl, parpt, ac.parmot, ac.pargen, ifirst, Ldebug)

            #         crzEINOx[i] = EINOx2
            #         crzFAR[i]   = FAR

            #         gam = Ftotal/BW - DoL
            #         # println(Wfrac, "--> gam = ", gam, " ROC = ", sin(gam)*V*60/0.3048, " FFsteadylevel = ", crzmdotf[i])
            #         if Tt41 ≤ Tt41CR && Ftotal>F    
            #             println("gam = ", gam, "ROC = ", sin(gam)*V*60/0.3048)
            #             gam = 0.0
            #             println("Min Temp reached")
            #         end
            #     end
            #     if iter == itermax && abs(Ftotal - F)>1
            #         println("Steady level cruise not converged - increase itermax?")
            #     end
            # end

            

            if ROCmaxcrz[i]<=0
                ROCmaxcrz[i] = 0.0
                println("Setting ROC maxcrz to 0.0")
            end
            println(@sprintf("%12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f", 
            FL[i], TAScr/kts_to_mps, CAScr/kts_to_mps, Mcr, F, 1/DoL, Tt41crzmax[i], Tmetmax, FFmaxcrz[i], Tt41crz[i], Tmetcrz, crzmdotf[i], CL, ac.parad[iaCLh, ip]))
        end #cruise section done

    end #outer loop

    # println("crzmdotf = ",crzmdotf)
    maxcruisealt = alts[maximum(findall(crzmdotf -> crzmdotf != 0, crzmdotf))]

    return W, alts[iceil], V0s, desTAS, ROC, mdotf, crzmdotf, crzTAS, EGTcrz, FFmaxcrz, ROCmaxcrz, Tt41crz, Tt41crzmax, crzEINOx, clmbEINOx, crzFAR, maxcruisealt
end

function get_cruisespeed(h, MNcr)
    h_ft = round(h/0.3048/1000)*1000
    Vcr2 = 280/1.944
    htrans_ft = Hptrans(Vcr2, 0.78)
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
    htrans_ft = Hptrans(Vcl2/1.944, 0.78)
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
    htrans_ft = Hptrans(Vdes2/1.944, 0.78)

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

function show_ff_sens(ff, H; ax = nothing)

    if ax === nothing
        fig, ax = plt.subplots()
    end
    ffnorm = ff'./minimum.(eachcol(ff))

    ax.plot(ffnorm', H)

end


## FOR 77W__
function get_cruisespeed_77W(h, MNcr)
    h_ft = round(h/0.3048/1000)*1000
    Vcr2 = 310/1.944
    htrans_ft = Hptrans(Vcr2, 0.84)
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

function get_climbspeed_77W(h,MNcr)
    h_ft = h/0.3048
    T, P, ρ,  a = atmos(h/1000)
    TAS = a*MNcr
    Vcl1 = 310/1.944 # From BADA 77W__.APF
    Vcl2 = 310/1.944 # From BADA 77W__.APF
    htrans_ft = Hptrans(Vcl2/1.944, 0.84)
    VstallTO = 139/1.944 # From 77W__.PTF file BADA, flap 5 (Standard Takeoff Setting for B77W)
    
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
    
function get_descentspeed_77W(h,MNcr)
    h_ft = h/0.3048
    T, P, ρ,  a = atmos(h/1000)
    TAS = a*MNcr
    Vdes1 = 310/1.944 # From BADA 77W__.APF
    Vdes2 = 310/1.944 # From BADA 77W__.APF
    htrans_ft = Hptrans(Vdes2/1.944, 0.84)

    VstallLD = 122/1.944 # From 77W__.PTF file BADA, flap 30 (Standard Landing Setting for B77W)
    
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