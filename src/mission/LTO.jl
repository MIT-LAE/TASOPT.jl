function LTO(name, pari,parg,para,pare; fileout = stdout)

    #LTO values
    LTOpoints = [1.0, 0.85, 0.3, 0.07]
    # LTOEINOx = zero(4)
    # LTOmdotf = zero(4)
    # LHV = parg[igLHVfuel]
    # Get Foo based on Tt41 or T/O thrust
    # Tt41 = pare[ieTt41, ipstatic]
    Foo = pare[ieFe , ipstatic] # This is for aircraft!
    ip = iptest
    
    if Foo == 0.0
        icall = 1 #Tt41 specified
        # println("Tt41 specified mode")
    else
        icall = 2 #Fn specified
        # println("Fn specified mode")
    end

    
    icool = 1
    inite1 = 0

    _, _ = TASOPT.tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, inite1)
    Ftotal = pare[ieFe, ip]
    mdotf = pare[iemdotf, ip]
    EINOx1  = pare[ieEINOx1, ip]

    # NPSS_success, Foo, heatexcess, 
    # mdotf, EINOx1, FAR, Mtip, Tblade, Tt3, OPR, BPR,
    # Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, 0.0, 0.0,
    # Foo, 0.0, mofft, Pofft, ifirst, parg, parpt, pare, ip)
    
    println(fileout, "Aircraft LTO characteristics \t\t"*Dates.format(now(), DateFormat("u dd yyyy")))
    println(fileout, name)
    println(fileout, @sprintf("Foo = %.3f [kN]", Foo/1000))
    println(fileout, @sprintf("%2s) %4s %15s %15s %15s %15s %8s %16s", 
    "#", "TS", "Fn[kN]", "̇mf[kg/s]", "EI(NOₓ)[g/kg]", "NOₓ[g/s]", "deNOₓ[%]", "JetAeqEI(NOₓ)" ))
    
        for (i,TS) in enumerate(LTOpoints)
            # NPSS_success, Ftotal, heatexcess, 
            #mdotf, EINOx1, FAR, Mtip, Tblade, Tt3, OPR, BPR, Wc3, Tt41, EGT = NPSS_TFsysOD(NPSS, 0.0, 0.0, Foo*TS, 0.0, mofft, Pofft, ifirst, parg, parpt, pare, ip)
            _, _ = TASOPT.tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, inite1)
            
            Ftotal = pare[ieFe, ip]
            mdotf = pare[iemdotf, ip]
            EINOx1  = pare[ieEINOx1, ip]

            ifirst  = false
    
            println(fileout, 
            @sprintf("%2d) %4.2f %15.5f %15.5f %15.5f %15.5f",
            i, TS, Ftotal/1000, mdotf, EINOx1, EINOx1*mdotf))
        end
    
    
    
    end