# This file will help print useful aircraft level parameters/ performance 

"""
Cruise characteristics
Aircraft deck for contrail avoidence
"""
function cruisechar(io, name, fracW0, M0, FL, FF, FFcrz, ROC, EGT)

    N = length(FL)
    println(io, "Aircraft cruise characteristics                           "*Dates.format(now(), DateFormat("u dd yyyy")))
    println(io, name)
    println(io, @sprintf("Cruise Mach : %3.2f\tW0 = %10.2f kg", M0, para[iafracW, ipcruise1].*parg[igWMTO]/9.81))
    
    println(io,"")
    println(io, @sprintf("%15s  %-7s%6s%7s", "FFmax [kg/s]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥260
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], FF[1, i], FF[2, i], FF[3, i]))
        end
    end
 
    println(io,"")
    println(io, @sprintf("%15s  %-7s%6s%7s", "FFcrz [kg/s]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥260
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], FFcrz[1, i], FFcrz[2, i], FFcrz[3, i]))
        end
    end

    println(io,"")
    println(io, @sprintf("%15s  %-7s%6s%7s", "ROCmax [ft/min]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥260
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], ROC[1, i], ROC[2, i], ROC[3, i]))
        end
    end

    println(io, "")
    println(io, @sprintf("%15s  %-7s%6s%7s", "Cruise EGT [K]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥260
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], EGT[1, i], EGT[2, i], EGT[3, i]))
        end
    end
end