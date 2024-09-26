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
            powersizing!(ac, engine_type, ipstatic)
            icall = 0

            ductedfancalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, initeng)

            parg[igA7] = pare[ieA7, ip] / pare[ieA7fac, ip]
            
            pare[ieA2, :] .= pare[ieA2, ip]
            pare[ieA7, :] .= parg[igA7] .* pare[ieA7fac, :]

            pare[iembfD, :] .= pare[iembfD, ip]
            pare[iepifD, :] .= pare[iepifD, ip]

            poweroper!(ac, engine_type, ip)

        else
            if ip in range(ipstatic, ipclimbn)
                icall = 1
            else
                icall = 2
            end

            ductedfancalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, initeng)

            poweroper!(ac, engine_type, ip)
        end
    end
end