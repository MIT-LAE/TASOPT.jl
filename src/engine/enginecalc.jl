function enginecalc!(ac, case, engine_type, ip, initeng, iterw = 0)
    #Unpack data storage arrays
    pari = ac.pari
    parg = ac.parg
    para = ac.parad
    pare = ac.pared  
    if engine_type == "turbofan"
        if case == "design"
            icall = 0
            icool = 1
            if (iterw == 1 || initeng == 0)
                # initialize engine state
                inite1 = 0
            else
                # start with current engine state
                inite1 = 1
            end

            ichoke5, ichoke7 = tfcalc!(pari,parg,view(para, :, ip), view(pare, :, ip), ip, icall, icool, inite1)

            # store engine design-point parameters for all operating points
            parg[igA5] = pare[ieA5, ip] / pare[ieA5fac, ip]
            parg[igA7] = pare[ieA7, ip] / pare[ieA7fac, ip]
            
            pare[ieA2, :] .= pare[ieA2, ip]
            pare[ieA25, :] .= pare[ieA25, ip]
            pare[ieA5, :] .= parg[igA5] .* pare[ieA5fac, :]
            pare[ieA7, :] .= parg[igA7] .* pare[ieA7fac, :]

            pare[ieNbfD, :] .= pare[ieNbfD, ip]
            pare[ieNblcD, :] .= pare[ieNblcD, ip]
            pare[ieNbhcD, :] .= pare[ieNbhcD, ip]
            pare[ieNbhtD, :] .= pare[ieNbhtD, ip]
            pare[ieNbltD, :] .= pare[ieNbltD, ip]

            pare[iembfD, :] .= pare[iembfD, ip]
            pare[iemblcD, :] .= pare[iemblcD, ip]
            pare[iembhcD, :] .= pare[iembhcD, ip]
            pare[iembhtD, :] .= pare[iembhtD, ip]
            pare[iembltD, :] .= pare[iembltD, ip]

            pare[iepifD, :] .= pare[iepifD, ip]
            pare[iepilcD, :] .= pare[iepilcD, ip]
            pare[iepihcD, :] .= pare[iepihcD, ip]
            pare[iepihtD, :] .= pare[iepihtD, ip]
            pare[iepiltD, :] .= pare[iepiltD, ip]
            
        elseif case == "off_design"
            if ip in range(ipstatic, ipclimbn)
                icall = 1
                icool = 1
            else
                icall = 2
                icool = 1
            end
            ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, initeng)
        elseif case == "cooling_sizing"
            icall = 1
            icool = 2
            ichoke5, ichoke7 = tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, initeng)
        end

    elseif engine_type == "ducted_fan"
        if case == "design"
            #Design ducted fan for start of cruise
            icall = 0

            ductedfancalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, initeng)

            parg[igA7] = pare[ieA7, ip] / pare[ieA7fac, ip]
            
            pare[ieA2, :] .= pare[ieA2, ip]
            pare[ieA7, :] .= parg[igA7] .* pare[ieA7fac, :]

            pare[iembfD, :] .= pare[iembfD, ip]
            pare[iepifD, :] .= pare[iepifD, ip]
            pare[ieNbfD, :] .= pare[ieNbfD, ip]

            #Design fuel cell for static conditions
            Pfanmax = pare[iePfanmax, ipstatic]

            ## Model of electric machine to deliver Pfanmax
            Pmotormax = Pfanmax #100% efficiency for now
            ##
            pare[iePfcdes, :] .= Pmotormax

            powersizing!(ac, engine_type, ipstatic)

            #Evaluate state at start of cruise
            Pfan = pare[iePfan, ip]
            ## Model of electric machine to deliver Pfan
            Pmotor = Pfan #100% efficiency for now
            ##
            pare[iePfc, ip] = Pmotor
            poweroper!(ac, engine_type, ip)

        else
            if ip in range(ipstatic, ipcutback)
                icall = 1

            elseif ip in range(ipclimb1, ipclimbn)
                icall = 2

                neng = parg[igneng]
                WMTO = parg[igWMTO]
                S = parg[igS]
                DoL = para[iaCD, ip] / para[iaCL, ip]
                W = para[iafracW, ip] * WMTO
                BW = W + para[iaWbuoy, ip]
                CL = para[iaCL, ip]
                ρ = pare[ierho0, ip]
                #Ftotal = BW * (DoL + para[iagamVdes, ip])
                ROC = para[iaROCdes, ip]

                #Solve with numerical & analytical solution
                f(γ) = ROC - sin(γ) * sqrt(2*BW*cos(γ) / (ρ*S*CL))
                γ = find_zero(f, para[iagamV, ip])
                A = sin(γ)
                B = DoL
                ϕ = (sqrt(-A^2*B^6 - 2*A^2*B^4 - A^2*B^2 + B^6 + 2*B^4 + B^2) + A*B^2 + A)/(B^2 + 1)
                Ftotal = BW * ϕ
                pare[ieFe, ip] = Ftotal / neng

            else
                icall = 2

            end

            ductedfancalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, initeng)

            Pfan = pare[iePfan, ip]
            ## Model of electric machine to deliver Pfan
            Pmotor = Pfan #100% efficiency for now
            ##
            pare[iePfc, ip] = Pmotor

            poweroper!(ac, engine_type, ip)
        end
    end
end