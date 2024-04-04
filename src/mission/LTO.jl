"""
    LTO(name, ac; fileout = stdout)

Prints out LTO EI(NOₓ) values
"""
function LTO(name, ac; fileout = stdout)

    #LTO values
    LTOpoints = [1.0, 0.85, 0.3, 0.07]

    LHV = ac.parg[igLHVfuel]
    # Get Foo based on Tt4 or T/O thrust
    Tt4 = ac.pared[ieTt4, ipstatic]
    Foo = ac.pared[ieFe , ipstatic]

    ip = iptest

    ac.pared[:, ip] = ac.pared[:, ipstatic]
    ac.parad[:, ip] = ac.parad[:, ipstatic]

    EIs = zeros(4)
    mfs = zeros(4)

    icall = 2
    icool = 1
    initeng = 0

    tfcalc!(ac.pari, ac.parg, view(ac.parad, :, ip), 
                                view(ac.pared, :, ip), ip, icall, icool, initeng)


    println(fileout, "Aircraft LTO characteristics \t\t"*Dates.format(now(), DateFormat("u dd yyyy")))
    println(fileout, name)
    println(fileout, @sprintf("Foo = %.3f [kN]", Foo/1000))
    println(fileout, @sprintf("%2s) %4s %15s %15s %15s %15s", 
    "#", "TS", "Fn[kN]", "̇mf[kg/s]", "EI(NOₓ)[g/kg]",  "NOₓ[g/s]" ))

        for (i,TS) in enumerate(LTOpoints)
            # println(@sprintf("Running LTO point %.1f %%", TS*100))
            ac.pared[ieFe, iptest] = Foo*TS
            tfcalc!(ac.pari, ac.parg, view(ac.parad, :, iptest), 
                                    view(ac.pared, :, iptest), 
                                    ip, icall, icool, initeng)
            
            F_total = ac.pared[ieFe, iptest]
            mdotf = ac.pared[ieff, iptest] * ac.pared[iemcore]
            P3_kPa = ac.pared[iept3, iptest]/1000.0
            T3_K   = ac.pared[ieTt3, iptest]
            EINOx = EINOx3(ac, iptest)
            EIs[i] = EINOx
            mfs[i] = mdotf
        
            println(fileout, 
            @sprintf("%2d) %4.2f %15.5f %15.5f %15.5f %15.5f",
            i, TS, F_total/1000, mdotf, EINOx, EINOx*mdotf))

        end

        println(fileout, "------------------------- AEIC ENG_EI output --------------------------")
        for (i,j) in zip([3,2,1,4], [1,2,3,4])
            println(fileout,@sprintf("%s, %d, 0.00, 0.00, %10.5f, 0.00, 1.00,  %10.5f, %s, 1.0",
                                        name, j, EIs[i], mfs[i], name))
        end

        return EIs, mfs
end

"""
    EINOx3(P3_kPa,T3_K, sp_humidity = 0.00634)


Returns the EI(NOₓ) for a CFM56 level engine.
Uses a third order polynomial in T₃ 
"""
function EINOx3(P3_kPa, T3_K, sp_humidity = 0.00634)
    # Constants derived using a CRN model for a CFM56 tech level engine
    a = 6.25528852e-08
    b = -1.17064467e-04
    c = 7.36953400e-02
    d = -1.50392850e+01

    H = -19.0*(sp_humidity - 0.00634)
    
    EI_NOx = exp(H)*P3_kPa^0.4 * (a*T3_K^3 + b*T3_K^2 + c*T3_K + d)
   
    return EI_NOx
end
"""
    EINOx4(P3_kPa, T3_K, sp_humidity = 0.00634)

Returns the EI(NOₓ) for a CFM56 level engine.
Uses a fourth order polynomial which behaves a little better at low T₃ than the cubic
"""
function EINOx4(P3_kPa, T3_K, sp_humidity = 0.00634)
    a = 4.85354237e-11
    b = -6.51089333e-08
    c =  7.19366066e-06
    d = 2.06850617e-02
    e =  -6.69412110e+00

    H = -19.0*(sp_humidity - 0.00634)

    EI_NOx = exp(H)*P3_kPa^0.4 * (a*T3_K^4 + b*T3_K^4 + c*T3_K^2 + d*T3_K + e)
    return EI_NOx
end  # function EINOx4

function EINOx3(ac::aircraft, ip::Int, sp_humidity = 0.00634)
       P3_kPa = ac.pared[iept3, ip]/1000.0
       T3_K   = ac.pared[ieTt3, ip]
       EINOx3(P3_kPa, T3_K, sp_humidity)
end