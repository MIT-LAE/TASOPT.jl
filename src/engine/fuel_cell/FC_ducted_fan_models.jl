function calculate_fuel_cell_with_ducted_fan(ac, case, imission, ip, initializes_engine, iterw = 0) 
    #Unpack aircraft data
    pari, parg, _, para, pare, _, _, _, _, _, _ = unpack_ac(ac, imission)
    fcdata = ac.engine.data #Extract fuel cell data

    if case == "design"
        #Design ducted fan for start of cruise
        icall = 0

        ductedfancalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, initializes_engine)

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
        fcdata.design_power = Pmotormax

        ipdes = ipstatic
        size_fuel_cell!(ac, ipdes, imission)

        #Evaluate state at start of cruise
        ip = ipcruise1

        Pfan = pare[iePfan, ip]
        ## Model of electric machine to deliver Pfan
        Pmotor = Pfan #100% efficiency for now
        ##
        fcdata.fuel_cell_power[ip, imission] = Pmotor
        operate_fuel_cell!(ac, ip, imission)

    else
        if ip in range(ipstatic, ipcutback)
            icall = 1

        elseif ip in range(ipclimb1, ipclimbn)
            icall = 2

            calculate_thrust_from_ROC!(ac, ip, imission) #Calculate thrust required for climb

        else
            icall = 2

        end

        ductedfancalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, initeng)

        Pfan = pare[iePfan, ip]
        ## Model of electric machine to deliver Pfan
        Pmotor = Pfan #100% efficiency for now
        ##
        fcdata.fuel_cell_power[ip, imission] = Pmotor

        operate_fuel_cell!(ac, ip, imission)
    end
end