"""
A trivial wrapper that runs a NPSS batch file `bat_file` at the specifed `dir`

>Note: 
>Remember to use `/` as separators
"""
function NPSS_run(dir, bat_file)
    global time_run_NPSS += @elapsed run(`cmd /c cd $dir '&&' $bat_file `)
    return nothing
end

"""
This function starts up and returns an NPSS process that can then be written to
"""
function startNPSS(dir, bat_file)
    NPSS = open(`cmd /c cd $dir '&&' $bat_file `, "w")
    return NPSS
end

"""
Ends NPSS process that can then be written to
"""
function endNPSS(NPSS)
    write(NPSS, "999 \\n")
    if(NPSS.exitcode == 0)
        close(NPSS)
    end
end

"""
# Design mode

Writes an input file for NPSS Turboshaft model 

    - Abmient altitude and mach number
    - Specifies πᵢ ∀ i ∈ {HPC, LPC}
    - Shaft power demand, Tt41 
    - DeNOx target
    - SCR paramters -> CPSI, w, l

"""
function NPSS_TShaft_input(alt_in, MN_in, 
    SHP_dmd, Tt41, 
    LPC_PR, HPC_PR,
    cpsi, w, lcat, deNOx, first ;
    LHV = 120, file_name = "NPSS_Turboshaft/EngineInputs.inp")

    open(file_name, "w") do io
        if(first)
            println(io, "int first = 1;")
        else 
            println(io, "int first = 0;")
        end
        println(io, "\n// Abmient conditions")
        println(io, "Eng.Amb.alt_in = ", alt_in/ft_to_m, ";" )
        println(io, "Eng.Amb.MN_in  = ", MN_in, ";" )

        println(io, "\n// Targets")
        println(io, "real SHP_dmd = ", SHP_dmd, ";" )
        println(io, "real Tt41    = ", Tt41, ";" )
        println(io, "real deNOx_target = ", deNOx, ";")

        println(io, "\n// LHV based on fuel")
        println(io, "Eng.FusEng.LHV = ", LHV*429.923, ";")

        println(io, "\n// Design parameters")
        println(io, "Eng.CmpL.PRdes = ", LPC_PR, ";")
        println(io, "Eng.CmpH.PRdes = ", HPC_PR, ";")

        println(io, "\n// PCEC parameters")
        println(io, "Eng.PCEC.l = ", lcat, ";")
        println(io, "Eng.PCEC.w = ", w, ";")
        println(io, "Eng.PCEC.cpsi = ", cpsi, ";")
    end

end

"""
OffdesMode

Writes an input file for NPSS Turboshaft model in off-des conditions

    - Abmient altitude and mach number
    - Tt41 

"""
function NPSS_TShaft_input(alt_in, MN_in, 
                            Tt41, Nshaft, first ; 
                            file_name = "NPSS_Turboshaft/OffDesInputs.inp")

    open(file_name, "w") do io
        if(first)
            println(io, "int first = 1;")
        else 
            println(io, "int first = 0;")
        end
        println(io, "\n// Abmient conditions")
        println(io, "Eng.Amb.alt_in = ", alt_in/ft_to_m, ";")
        println(io, "Eng.Amb.MN_in  = ", MN_in , ";")

        println(io, "\n// Targets")
        println(io, "real Tt41   = ", Tt41  , ";")
        println(io, "real N2_dmd = ", Nshaft, ";")
        
    end

end

"""
Fan - Design mode:

Writes the desired conditions to the external NPSS process
"""
function runNPSS_Fan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Fn::Float64,
                    Kinl::Float64, Φinl::Float64, 
                    πfan::Float64, first)
    input_string = "111 "*
                "first = $first;" *
                "DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
                "DuctedFan.Amb.MN_in  = $MN_in ;" *
                "Fn_target = $Fn ;" *
                "DuctedFan.InEng.Kinl = $Kinl;" *
                "DuctedFan.Phiinl     = $Φinl;" *
                "DuctedFan.Fan.PRdes = $πfan;" *
                "\n"

    write(NPSS, input_string)
    # sleep(0.5)

end

"""
Fan - Off-design mode:

Writes the desired conditions to the external NPSS process
"""
function runNPSS_Fan(NPSS::Base.Process, alt_in::Float64, MN_in::Float64, Pin::Float64,
                    Kinl::Float64, Φinl::Float64, first)
    input_string = "222 "*
                "first = $first;" *
                "DuctedFan.Amb.alt_in = $(alt_in/ft_to_m);" *
                "DuctedFan.Amb.MN_in  = $MN_in ;" *
                "ShP_input = $(-Pin/745.7) ;" *
                "DuctedFan.InEng.Kinl = $Kinl;" *
                "DuctedFan.Phiinl     = $Φinl;" *
                "\n"

    write(NPSS, input_string)
    # sleep(0.5)
end
        # println(io, "DuctedFan.setOption(\"switchDes\",\"OFFDESIGN\");")

        println(io, "\n// Abmient conditions")
        println(io, "DuctedFan.Amb.alt_in = ", alt_in/ft_to_m, ";" )
        println(io, "DuctedFan.Amb.MN_in  = ", MN_in, ";" )

        println(io, "\n// Input power")
        println(io, "real ShP_input = ", -Pin/745.7, ";")

        
        println(io, "\n// BLI inputs")
        println(io, "DuctedFan.InEng.Kinl = ", Kinl, ";")
        println(io, "DuctedFan.Phiinl     = ", Φinl, ";")


    end

end