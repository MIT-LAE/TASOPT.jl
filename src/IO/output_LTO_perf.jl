"""
    output_LTO_perf(name, ac; fileout = stdout)

Prints out LTO EI(NOₓ) values
"""
function output_LTO_perf(name, ac; fileout = stdout, method = "cubic", extra_points = false)

    #LTO values
    if extra_points
        LTOpoints = [1.0, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.07]
    else
        LTOpoints = [1.0, 0.85, 0.3, 0.07]
    end

    LHV = ac.parg[igLHVfuel]
    # Get Foo based on Tt4 or T/O thrust
    Tt4 = ac.pared[ieTt4, ipstatic]
    Foo = ac.pared[ieFe , ipstatic]# * ac.parg[igneng]

    ip = iptest

    ac.pared[:, ip] = ac.pared[:, ipstatic]
    ac.parad[:, ip] = ac.parad[:, ipstatic]

    ltopts = size(LTOpoints)[1]
    EIs = zeros(ltopts)
    mfs = zeros(ltopts)

    icall = 2
    icool = 1
    initeng = false

    TASOPT.tfcalc!(ac.pari, ac.parg, view(ac.parad, :, ip), 
                                view(ac.pared, :, ip), ac.wing, ip, icall, icool, initeng)


    println(fileout, "Aircraft LTO characteristics \t\t"*Dates.format(now(), DateFormat("u dd yyyy")))
    println(fileout, name)
    println(fileout, @sprintf("Foo = %.3f [kN]", Foo/1000))
    println(fileout, @sprintf("%2s) %4s %15s %15s %15s %15s", 
    "#", "TS", "Fn[kN]", "̇mf[kg/s]", "EI(NOₓ)[g/kg]",  "NOₓ[g/s]" ))

        for (i,TS) in enumerate(LTOpoints)
            # println(@sprintf("Running LTO point %.1f %%", TS*100))
            ac.pared[ieFe, iptest] = Foo*TS
            TASOPT.tfcalc!(ac.pari, ac.parg, view(ac.parad, :, iptest), 
                                    view(ac.pared, :, iptest), ac.wing,
                                    ip, icall, icool, initeng)
            
            F_total = ac.pared[ieFe, iptest]# * ac.parg[igneng]
            mdotf = ac.pared[ieff, iptest] * ac.pared[iemcore]
            P3_kPa = ac.pared[iept3, iptest]/1000.0
            T3_K   = ac.pared[ieTt3, iptest]
            EI = EINOx(ac, iptest; method=method)
            EIs[i] = EI
            mfs[i] = mdotf
        
            println(fileout, 
            @sprintf("%2d) %4.2f %15.5f %15.5f %15.5f %15.5f",
            i, TS, F_total/1000, mdotf, EI, EI*mdotf))

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
Assumes a default specific humidity of 0.00634 kg water/kg dry air per 
ICAO Annex 16 Vol. II (part 2.1.4.1)
"""
function EINOx3(P3_kPa, T3_K, sp_humidity = 0.00634, ac_type = "1")
    # Constants derived using a CRN model for a CFM56 tech level engine
        
    if ac_type == :wide
        a = 1.01084407e-07
        b = -2.12716481e-04
        c = 1.50618950e-01
        d = -3.49737491e+01
    else
        a = 6.25528852e-08
        b = -1.17064467e-04
        c = 7.36953400e-02
        d = -1.50392850e+01
    end

    H = -19.0*(sp_humidity - 0.00634)
    
    EI_NOx = exp(H)*P3_kPa^0.4 * (a*T3_K^3 + b*T3_K^2 + c*T3_K + d)
   
    return EI_NOx
end
"""
    EINOx4(P3_kPa, T3_K, sp_humidity = 0.00634)

Returns the EI(NOₓ) for a CFM56 level engine.
Uses a fourth order polynomial which behaves a little better at low T₃ than the cubic
Assumes a default specific humidity of 0.00634 kg water/kg dry air per 
ICAO Annex 16 Vol. II (part 2.1.4.1)
"""
function EINOx4(P3_kPa, T3_K, sp_humidity = 0.00634, ac_type = "1")
    if ac_type == :wide
        a = 9.15868718e-10
        b = -2.60435361e-06
        c = 2.75257358e-03
        d = -1.27746983e+00
        e = 2.19853429e+02
    else
        a = 4.85354237e-11
        b = -6.51089333e-08
        c =  7.19366066e-06
        d = 2.06850617e-02
        e =  -6.69412110e+00
    end

    H = -19.0*(sp_humidity - 0.00634)

    EI_NOx = exp(H)*P3_kPa^0.4 * (a*T3_K^4 + b*T3_K^3 + c*T3_K^2 + d*T3_K + e)
    return EI_NOx
end  # function EINOx4

function EINOx(ac::aircraft, ip::Int; sp_humidity = 0.00634, method="cubic")
       P3_kPa = ac.pared[iept3, ip]/1000.0
       T3_K   = ac.pared[ieTt3, ip]
       if lowercase(method) == "cubic"
            EI = EINOx3(P3_kPa, T3_K, sp_humidity,ac.aircraft_type)
       elseif lowercase(method) == "quartic"
            EI = EINOx4(P3_kPa, T3_K, sp_humidity,ac.aircraft_type)
       end
       return EI
end