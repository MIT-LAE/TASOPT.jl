function enginecalc!(ac, case, engine_type, initeng, iterw)
    #Unpack data storage arrays
    pari = ac.pari
    parg = ac.parg
    para = ac.parad
    pare = ac.pared  
    if engine_type == "turbofan"
        if case == "sizing"
            ip = ipcruise1
            icall = 0
            icool = 1
            if (iterw == 1 || initeng == 0)
                # initialize engine state
                inite1 = 0
            else
                # start with current engine state
                inite1 = 1
            end

            ichoke5, ichoke7 = tfcalc(pari,parg,view(para, :, ip), view(pare, :, ip), ip, icall, icool, inite1)

            # store engine design-point parameters for all operating points
            parg[igA5] = pare[ieA5, ip] / pare[ieA5fac, ip]
            parg[igA7] = pare[ieA7, ip] / pare[ieA7fac, ip]
            for jp = 1:iptotal
                pare[ieA2, jp] = pare[ieA2, ip]
                pare[ieA25, jp] = pare[ieA25, ip]
                pare[ieA5, jp] = parg[igA5] * pare[ieA5fac, jp]
                pare[ieA7, jp] = parg[igA7] * pare[ieA7fac, jp]

                pare[ieNbfD, jp] = pare[ieNbfD, ip]
                pare[ieNblcD, jp] = pare[ieNblcD, ip]
                pare[ieNbhcD, jp] = pare[ieNbhcD, ip]
                pare[ieNbhtD, jp] = pare[ieNbhtD, ip]
                pare[ieNbltD, jp] = pare[ieNbltD, ip]

                pare[iembfD, jp] = pare[iembfD, ip]
                pare[iemblcD, jp] = pare[iemblcD, ip]
                pare[iembhcD, jp] = pare[iembhcD, ip]
                pare[iembhtD, jp] = pare[iembhtD, ip]
                pare[iembltD, jp] = pare[iembltD, ip]

                pare[iepifD, jp] = pare[iepifD, ip]
                pare[iepilcD, jp] = pare[iepilcD, ip]
                pare[iepihcD, jp] = pare[iepihcD, ip]
                pare[iepihtD, jp] = pare[iepihtD, ip]
                pare[iepiltD, jp] = pare[iepiltD, ip]
            end
        elseif case == "offdesign"
        end
    end
end