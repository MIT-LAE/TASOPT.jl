function LTO(name, NPSS::Base.Process, ifirst, pare; fileout = stdout)
#LTO values
LTOpoints = [1.0, 0.85, 0.3, 0.07]
LTOEINOx = zero(4)
LTOmdotf = zero(4)
LHV = parg[igLHVfuel]
# Get Foo based on Tt4 or T/O thrust
Tt4 = pare[ieTt4, ipstatic]
Foo = pare[ieFe , ipclimb1]*1.5  # Temporary need to update with actual take-off Foo
ip = iptest

NPSS_success, Foo, η, P, Hrej, heatexcess, 
mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR,
Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, 0.0, 0.0,
Foo, 0.0, 0.0, 0.0, 0.0, 0.0, ifirst, parg, parpt, pare, ip)
println(fileout, "Aircraft LTO characteristics \t\t"*Dates.format(now(), DateFormat("u dd yyyy")))
println(fileout, name)
println(fileout, @sprintf("Foo = %.3f [kN]", Foo/1000))
println(fileout, @sprintf("%2s) %4s %15s %15s %15s %15s %8s %16s", 
"#", "TS", "Fn[kN]", "̇mf[kg/s]", "EI(NOₓ)[g/kg]", "NOₓ[g/s]", "deNOₓ[%]", "JetAeqEI(NOₓ)" ))

    for (i,TS) in enumerate(LTOpoints)
        NPSS_success, Ftotal, η, P, Hrej, heatexcess, 
        mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, 0.0, 0.0, Foo*TS, 0.0, 0.0, 0.0, 0.0, 0.0, ifirst, parg, parpt, pare, ip)
        ifirst  = false

        println(fileout, 
        @sprintf("%2d) %4.2f %15.5f %15.5f %15.5f %15.5f %8.2f %16.5f",
        i, TS, Ftotal/1000, mdotf, EINOx2, EINOx2*mdotf, deNOx*100, EINOx2*43/LHV))
    end



end