function LTO(name, NPSS::Base.Process, ifirst, pare; fileout = stdout)
nTshaft = parpt[ipt_nTshaft]
#LTO values
LTOpoints = [1.0, 0.85, 0.3, 0.07]
LTOEINOx = zero(4)
LTOmdotf = zero(4)
LHV = parg[igLHVfuel]
# Get Foo based on Tt4 or T/O thrust
Tt4 = pare[ieTt4, ipstatic]
Foo = pare[ieFe , ipstatic]
ip = iptest
EIs = zeros(4)
mfs = zeros(4)
NPSS_success, Foo, η, P, Hrej, heatexcess, 
mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR,
Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, 0.0, 0.0,
Foo, 0.0, 0.0, 0.0, 0.0, 0.0, ifirst, parg, parpt, pare, ip)
println(fileout, "Aircraft LTO characteristics \t\t"*Dates.format(now(), DateFormat("u dd yyyy")))
println(fileout, name)
println(fileout, @sprintf("Foo = %.3f [kN]", Foo/1000))
println(fileout, @sprintf("%2s) %4s %15s %15s %15s %15s %15s %15s %8s %16s", 
"#", "TS", "Fn[kN]", "̇mf_total[kg/s]", "̇mf[kg/s]", "EI(NOₓ)[g/kg]", "Total_NOₓ[g/s]", "NOₓ[g/s]", "deNOₓ[%]", "JetAeqEI(NOₓ)" ))

    for (i,TS) in enumerate(LTOpoints)
        println(@sprintf("Running LTO point %.1f %%", TS*100))
        NPSS_success, Ftotal, η, P, Hrej, heatexcess, 
        mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, 0.0, 0.0, Foo*TS, 0.0, 0.0, 0.0, 0.0, 0.0, ifirst, parg, parpt, pare, ip)
        ifirst  = false

        println(fileout, 
        @sprintf("%2d) %4.2f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %8.2f %16.5f",
        i, TS, Ftotal/1000, mdotf, mdotf/nTshaft, EINOx2, EINOx2*mdotf, EINOx2*mdotf/nTshaft, deNOx*100, EINOx2*43/LHV))
        EIs[i] = EINOx2
        mfs[i] = mdotf/nTshaft
    end
    # for TS in LTOpoints[end:-1:1]
    #     NPSS_success, Ftotal, η, P, Hrej, heatexcess, 
    #     mdotf, deNOx, EINOx1, EINOx2, FAR, Tt3, OPR, Wc3, Tt41, EGT = NPSS_TEsysOD(NPSS, 0.0, 0.0, Foo*TS, 0.0, 0.0, 0.0, 0.0, 0.0, ifirst, parg, parpt, pare, ip)
    # end
    for (i,j) in zip([3,2,1,4], [1,2,3,4])
        println(@sprintf("ZIAENG, %d, 0.00, 0.00, %10.5f, 0.00, 1.00,  %10.5f, ZIAENG, 1.0", j, EIs[i], mfs[i]))
    end

end