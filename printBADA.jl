function printBADA(io, name, W0, flightceil, TAS, ROC, ffpmin, crzf, crzTAS )

    alt = flightceil
    γdes = 3.0 * π/180.0
    sing = sin(γdes)
# Flight levels to output in BADA file:
    FL = [  0 ,    5 ,   10 ,   15 ,   20 ,   
        30 ,   40 ,   60 ,   80 ,  100 , 
        120 ,  140 ,  160 ,  180 ,  200 ,
        220 ,  240 ,  260 ,  280 ,  290 , 
        310 ,  330 ,  350 ,  370 ,  390 ,
        410 ,  430 ,  431 ]
    maxalt = alt/ft_to_m + 7000 #idk why but AEIC subtracts 7000 ft

    println(io, "BADA PERFORMANCE FILE                                    "*Dates.format(now(), DateFormat("u dd yyyy")) )
    println(io, @sprintf(""))
    println(io, @sprintf("AC/Type: %6s", name))
    println(io, "                              Source OPF File:               "*Dates.format(now(), DateFormat("u dd yyyy")) )
    println(io, "                              Source APF file:               "*Dates.format(now(), DateFormat("u dd yyyy")) )
    println(io, @sprintf(""))
    println(io, @sprintf(" Speeds:   CAS(LO/HI)  Mach   Mass Levels [kg]         Temperature:  ISA"))
    println(io, @sprintf(" climb   - 250/310     0.80   low     -  %6.0f", W0[1]/9.81))
    println(io, @sprintf(" cruise  - 250/310     0.80   nominal -  %6.0f        Max Alt. [ft]:%7.0f", W0[2]/9.81, maxalt))
    println(io, @sprintf(" descent - 250/310     0.80   high    -  %6.0f", W0[3]/9.81))
    println(io, @sprintf("=========================================================================================="))
    println(io, @sprintf(" FL |          CRUISE           |               CLIMB               |       DESCENT       "))
    println(io, @sprintf("    |  TAS          fuel        |  TAS          ROCD         fuel   |  TAS  ROCD    fuel  "))
    println(io, @sprintf("    | [kts]       [kg/min]      | [kts]        [fpm]       [kg/min] | [kts] [fpm] [kg/min]"))
    println(io, @sprintf("    |          lo   nom    hi   |         lo    nom    hi    nom    |        nom    nom   "))
    println(io, @sprintf("=========================================================================================="))
   
    n = length(FL)
    for i = 1:n
        if FL[i]≥ 330 && FL[i]≤430
    println(io, @sprintf("%3.0f |  %3.0f   %4.1f  %5.1f  %5.1f |  %3.0f     %4.0f  %4.0f  %4.0f   %4.1f  |  %3.0f   %4.0f  %4.1f  ", 
            FL[i], crzTAS[i], crzf[1, i], crzf[2, i], crzf[3, i], 
                                                                TAS[i],  ROC[1, i], ROC[2, i], ROC[3, i], ffpmin[i], TAS[i], TAS[i]*sing*60/ft_to_m, ffpmin[i]*0.1))
    println(io, @sprintf("    |                           |                                   | "))
        else
    println(io, @sprintf("%3.0f |                           |  %3.0f     %4.0f  %4.0f  %4.0f   %4.1f  |  %3.0f   %4.0f  %4.1f  ", 
                         FL[i],                               TAS[i],  ROC[1, i], ROC[2, i], ROC[3, i], ffpmin[i], TAS[i], TAS[i]*sing*60/ft_to_m, ffpmin[i]*0.1))
    println(io, @sprintf("    |                           |                                   | "))
        end
    end
    println(io, @sprintf("=========================================================================================="))
    


end