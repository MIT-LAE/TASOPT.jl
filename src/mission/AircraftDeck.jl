# This file will help print useful aircraft level parameters/ performance 

"""
    cruisechar(io, name, fracW0, M0, FL, TAS, 
    FFmax, FFcrz, ROC, EGT, Tt4crz, Tt4crzmax, crzFAR)

Cruise characteristics of a sized aircraft. 
Outputs an aircraft deck for contrail avoidence simulations.
"""
function cruisechar(io, name, fracW0, M0, FL, TAS, FFmax, FFcrz, ROC, EGT, Tt4crz, Tt4crzmax, crzFAR)
    FLcutoff = 150
    N = length(FL)
    println(io, "Aircraft cruise characteristics              "*Dates.format(now(), DateFormat("u dd yyyy")))
    println(io, name)
    println(io, @sprintf("Cruise Mach = %3.2f\n%-22s = %10.2f kg\n%-22s = %10.2f kg\n%-22s = %10.2f kg\n%-22s = %10.2f kg\n(Fuel reserves = %.1f%%)\n%-22s = %10.2f kg\n%-22s = %10.2f kg",
     M0,"W @TOC", para[iafracW, ipcruise1].*parg[igWMTO]/9.81, "WMTO", parg[igWMTO]/9.81,"Wpay", parg[igWpay]/9.81, "Wfuel (incl. reserves)", parg[igWfuel]/9.81, parg[igfreserve]*100, 
     "Wfmax", parg[igWfmax]/9.81,
     "ZFW = WMTO - Wfuel", (parg[igWMTO] - parg[igWfuel])/9.81 ))
     println(io, "Note: Following values assume fixed 100% payload")
    
     println(io,"")
     println(io, @sprintf("%15s", "TAS[m/s]",))
     println(io, @sprintf("%4s  %9s", "FL", "TAS"))
     for i = 1:N
         if FL[i]≥FLcutoff
             println(io, @sprintf("%4.0f  %9.4f", FL[i], TAS[i]))
         end
     end

    print_variables(io, fracW0, FL, FLcutoff, FFmax, "FFmax [kg/s]")
    print_variables(io, fracW0, FL, FLcutoff, FFcrz, "FFcrz [kg/s]")
    print_variables(io, fracW0, FL, FLcutoff, ROC, "ROCmax [ft/min]")
    print_variables(io, fracW0, FL, FLcutoff, EGT, "Cruise EGT [K]")
    print_variables(io, fracW0, FL, FLcutoff, crzFAR, "Cruise FAR")
    print_variables(io, fracW0, FL, FLcutoff, Tt4crzmax, "MaxClmb Tt4 [R]")
    print_variables(io, fracW0, FL, FLcutoff, Tt4crz, "Cruise FAR")

end

"""
Helper function to print variables in a given format
"""
function print_variables(io, W_lst, FL, FLcutoff, Var, Varname)
    println(io,"")
    println(io, @sprintf("%15s  %-7s%6s%7s", Varname, "--","Wfracs", "--"))
    print(io, @sprintf("%4s", "FL"))
    for W in W_lst
        print(io, @sprintf("%12.4f", W))
    end
    println(io,"")
    
    for i = 1:length(FL)
        if FL[i]≥FLcutoff
            print(io, @sprintf("%4.0f", FL[i]))
            for var in Var[i,:]
                print(io, @sprintf("%12.4f", var))
            end
            println(io,"")
        end
    end

end