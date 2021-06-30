# This file will help print useful aircraft level parameters/ performance 

"""
Cruise characteristics
Aircraft deck for contrail avoidence
"""
function cruisechar(io, name, fracW0, M0, FL, FFmax, FFcrz, ROC, EGT, Tt4crz, Tt4crzmax, crzFAR)

    N = length(FL)
    println(io, "Aircraft cruise characteristics              "*Dates.format(now(), DateFormat("u dd yyyy")))
    println(io, name)
    println(io, @sprintf("Cruise Mach = %3.2f\n%-22s = %10.2f kg\n%-22s = %10.2f kg\n%-22s = %10.2f kg\n%-22s = %10.2f kg\n(Fuel reserves = %.1f%%)\n%-22s = %10.2f kg\n%-22s = %10.2f kg",
     M0,"W @TOC", para[iafracW, ipcruise1].*parg[igWMTO]/9.81, "WMTO", parg[igWMTO]/9.81,"Wpay", parg[igWpay]/9.81, "Wfuel (incl. reserves)", parg[igWfuel]/9.81, parg[igfreserve]*100, 
     "Wfmax", parg[igWfmax]/9.81,
     "ZFW = WMTO - Wfuel", (parg[igWMTO] - parg[igWfuel])/9.81 ))
     println(io, "Note: Following values assume fixed 100% payload")
    
    println(io,"")
    println(io, @sprintf("%15s  %-7s%6s%7s", "FFmax [kg/s]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥270
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], FFmax[1, i], FFmax[2, i], FFmax[3, i]))
        end
    end
 
    println(io,"")
    println(io, @sprintf("%15s  %-7s%6s%7s", "FFcrz [kg/s]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥270
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], FFcrz[1, i], FFcrz[2, i], FFcrz[3, i]))
        end
    end

    println(io,"")
    println(io, @sprintf("%15s  %-7s%6s%7s", "ROCmax [ft/min]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥270
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], ROC[1, i], ROC[2, i], ROC[3, i]))
        end
    end

    println(io, "")
    println(io, @sprintf("%15s  %-7s%6s%7s", "Cruise EGT [K]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥270
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], EGT[1, i], EGT[2, i], EGT[3, i]))
        end
    end

    println(io, "")
    println(io, @sprintf("%15s  %-7s%6s%7s", "Cruise FAR", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥270
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], crzFAR[1, i], crzFAR[2, i], crzFAR[3, i]))
        end
    end
    println(io, "")
    println(io, @sprintf("%15s  %-7s%6s%7s", "MaxClmb Tt4 [R]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥270
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], Tt4crzmax[1, i], Tt4crzmax[2, i], Tt4crzmax[3, i]))
        end
    end
    println(io, "")
    println(io, @sprintf("%15s  %-7s%6s%7s", "Cruise Tt4 [R]", "--","Wfracs", "--"))
    println(io, @sprintf("%4s  %9.2f  %9.2f  %9.2f", "FL", fracW0[1], fracW0[2], fracW0[3]))
    for i = 1:N
        if FL[i]≥270
            println(io, @sprintf("%4.0f  %9.4f  %9.4f  %9.4f", FL[i], Tt4crz[1, i], Tt4crz[2, i], Tt4crz[3, i]))
        end
    end
end