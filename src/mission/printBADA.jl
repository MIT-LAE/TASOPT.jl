function printBADA(io, name, W0, maxalt, TAS, desTAS, ROC, ffpmin, crzf, crzTAS, FLin, Wpaymax; NOx = false, wide = false, crzNoxEI =[], ffpminNoxEI=[] )

    γdes = 3.0 * π/180.0
    sing = sin(γdes)
    kts_to_ft = 1.68781
# Flight levels to output in BADA file:
    if wide
        FL = float([  0 ,    5 ,   10 ,   15 ,   20 ,   
        30 ,   40 ,   60 ,   80 ,  100 , 
        120 ,  140 ,  160 ,  180 ,  200 ,
        220 ,  240 ,  260 ,  280 ,  290 , 
        310 ,  330 ,  350 ,  370 ,  390 ,
        410, 430, 431])
    else
        FL = float([  0 ,    5 ,   10 ,   15 ,   20 ,   
            30 ,   40 ,   60 ,   80 ,  100 , 
            120 ,  140 ,  160 ,  180 ,  200 ,
            220 ,  240 ,  260 ,  280 ,  290 , 
            310 ,  330 ,  350 ,  370 ,  390 ,
            410])
    end

    maxalt = maxalt/ft_to_m  #idk why but AEIC subtracts 7000 ft

    println(io, "TASOPT PERFORMANCE FILE                                     "*Dates.format(now(), DateFormat("u dd yyyy")) )
    println(io, @sprintf(" "))
    println(io, @sprintf("AC/Type: %6s", name))
    println(io, "                              Source OPF File:               "*Dates.format(now(), DateFormat("u dd yyyy")) )
    println(io, "                              Source APF file:               "*Dates.format(now(), DateFormat("u dd yyyy")) )
    println(io, @sprintf(" "))
    println(io, @sprintf(" Speeds:   CAS(LO/HI)  Mach   Mass Levels [kg]         Temperature:  ISA"))
    println(io, @sprintf(" climb   - 250/300     0.80   low     -  %6.0f", W0[1]/gee))
    println(io, @sprintf(" cruise  - 250/280     0.80   nominal -  %6.0f        Max Alt. [ft]:%7.0f", W0[2]/gee, maxalt))
    println(io, @sprintf(" descent - 250/290     0.80   high    -  %6.0f        Max Payload [kg]:%7.0f", W0[3]/gee, Wpaymax/gee))
    
    if NOx
    println(io, @sprintf("==================================================================================================================================="))
    println(io, @sprintf(" FL |                        CRUISE                  |                     CLIMB                   |            DESCENT            "))
    println(io, @sprintf("    |  TAS          fuel                EINOx        |  TAS          ROCD         fuel      NOx    |  TAS  ROCD    fuel      NOx   "))
    println(io, @sprintf("    | [kts]       [kg/min]             [g/kg]        | [kts]        [fpm]       [kg/min]  [g/kg]   | [kts] [fpm] [kg/min]  [g/kg]  "))
    println(io, @sprintf("    |          lo   nom    hi       lo   nom    hi   |         lo    nom    hi    nom       nom    |        nom    nom       nom   "))
    println(io, @sprintf("==================================================================================================================================="))
    else
    println(io, @sprintf("=========================================================================================="))
    println(io, @sprintf(" FL |          CRUISE           |               CLIMB               |       DESCENT       "))
    println(io, @sprintf("    |  TAS          fuel        |  TAS          ROCD         fuel   |  TAS  ROCD    fuel  "))
    println(io, @sprintf("    | [kts]       [kg/min]      | [kts]        [fpm]       [kg/min] | [kts] [fpm] [kg/min]"))
    println(io, @sprintf("    |          lo   nom    hi   |         lo    nom    hi    nom    |        nom    nom   "))
    println(io, @sprintf("=========================================================================================="))
    end
    
   
    n = length(FL)
    for i = 1:n
        iin = findall(x->x ≈ FL[i], FLin)[1]
        if FL[i]≥ 60 && FL[i]≤410

            # For AEIC ingestion
            if crzf[1, iin] == 0.0
                crzf[1, iin] = crzf[1, iin-1]
            end
            if crzf[2, iin] == 0.0
                crzf[2, iin] = crzf[1, iin]
            end
            if crzf[3, iin] == 0.0
                crzf[3, iin] = crzf[2, iin]
            end

            # For AEIC ingestion
            if ROC[2, iin] < 100
                ROC[2, iin] = ROC[1, iin]
            end
            
            # if nominal weight climb fuelburn 0 (not defined)
            if ffpmin[iin] == 0.0
                ffpmin[iin] = ffpmin_backup[iin]
            end
            
            if NOx
                println(io, @sprintf("%3.0f |  %3.0f   %5.1f %5.1f %5.1f     %5.1f %5.1f %5.1f |  %3.0f    %4.0f  %4.0f  %4.0f %6.1f   %6.1f    |  %3.0f   %4.0f  %4.1f    %4.1f   ", 
                        FLin[iin], crzTAS[iin], crzf[1, iin], crzf[2, iin], crzf[3, iin], crzNoxEI[1, iin], crzNoxEI[2, iin], crzNoxEI[3, iin],
                                                                            TAS[iin],  ROC[1, iin], ROC[2, iin], ROC[3, iin], ffpmin[iin], ffpminNoxEI[iin], desTAS[iin], desTAS[iin]*sing*60*kts_to_ft, ffpmin[iin]*0.1, ffpminNoxEI[iin]*0.1))
            else
                println(io, @sprintf("%3.0f |  %3.0f    %5.2f %5.2f %5.2f |  %3.0f    %4.0f  %4.0f  %4.0f  %6.2f  |  %3.0f   %4.0f   %5.2f   ", 
                        FLin[iin], crzTAS[iin], crzf[1, iin], crzf[2, iin], crzf[3, iin], 
                                                                            TAS[iin],  ROC[1, iin], ROC[2, iin], ROC[3, iin], ffpmin[iin], desTAS[iin], desTAS[iin]*sing*60*kts_to_ft, ffpmin[iin]*0.1))
            end
            println(io, @sprintf("    |                                                |                                             | "))
        else
            if NOx
                println(io, @sprintf("%3.0f |                                                |  %3.0f    %4.0f  %4.0f  %4.0f %6.1f   %6.1f    |  %3.0f   %4.0f   %4.1f    %4.1f   ", 
                FLin[iin],                               TAS[iin],  ROC[1, iin], ROC[2, iin], ROC[3, iin], ffpmin[iin], ffpminNoxEI[iin], desTAS[iin], desTAS[iin]*sing*60*kts_to_ft, ffpmin[iin]*0.3, ffpminNoxEI[iin]*0.3)) #landing configuration, N1~50%
            else
                println(io, @sprintf("%3.0f |                           |  %3.0f    %4.0f  %4.0f  %4.0f  %6.2f  |  %3.0f   %4.0f   %5.2f   ", 
                                    FLin[iin],                               TAS[iin],  ROC[1, iin], ROC[2, iin], ROC[3, iin], ffpmin[iin], desTAS[iin], desTAS[iin]*sing*60*kts_to_ft, ffpmin[iin]*0.3)) #landing configuration, N1~50%
            end
            println(io, @sprintf("    |                                                |                                             | "))
        end
    end
    println(io, @sprintf("=========================================================================================="))

end