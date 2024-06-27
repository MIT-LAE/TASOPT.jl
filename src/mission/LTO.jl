function LTO(name, pari,parg,para,pare; fileout = stdout)

    #LTO values
    LTOpoints = [1.0, 0.85, 0.3, 0.07]
    # LTOEINOx = zero(4)
    # LTOmdotf = zero(4)
    # LHV = parg[igLHVfuel]
    # Get Foo based on Tt41 or T/O thrust
    # Tt41 = pare[ieTt41, ipstatic]
    Foo = pare[ieFe , ipstatic] # This is for aircraft!
    ip = 1#iptest
    
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
    Foo = pare[ieFe, ip]
    mdotf = pare[iemcore, ip]#/1e2
    Pt3 = pare[iept3]
    Tt3 = pare[ieTt3]
    EINOx1  = get_EINOx_4(Pt3/1e3, Tt3, 0)

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
            pare[ieFe, ip] = Foo*TS
            _, _ = TASOPT.tfcalc!(pari, parg, view(para, :, ip), view(pare, :, ip), ip, icall, icool, inite1)
            
            Ftotal = pare[ieFe, ip]
            mdotf = pare[iemcore, ip]#/1e2
            Pt3 = pare[iept3]
            Tt3 = pare[ieTt3]
            EINOx1  = get_EINOx_4(Pt3/1e3, Tt3, 0)
            # pare[ieEINOx1, ip]

            ifirst  = false
    
            println(fileout, 
            @sprintf("%2d) %4.2f %15.5f %15.5f %15.5f %15.5f",
            i, TS, Ftotal/1000, mdotf, EINOx1, EINOx1*mdotf))
        end
    
    
    
    end

function get_EINOx_4(Pt3, Tt3, sp_humidity)
    # Pt3 in kPa and Tt3 in K
    # Coefficients
    a4 = 4.85354237e-11
    b4 = -6.51089333e-08
    c4 =  7.19366066e-06
    d4 = 2.06850617e-02
    e4 =  -6.69412110e+00
    H4 = -19.0 * (sp_humidity - 0.00634)

    EINOx = exp(H4) * Pt3^0.4 * (a4 * Tt3^4 + b4 * Tt3^3 + c4 * Tt3^2 + d4 * Tt3 + e4)

    return EINOx
end
    