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
    NPSS = open(`cmd /c cd $dir '&&' $bat_file `, "w+")
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

    - Abmient altitude [in m]
    - Flight mach number
    - Tt41 
    - the flag first sets whether NPSS should restart or use previous value

"""
function NPSS_TShaft_run(NPSS_TS, alt_in, MN_in, 
                            Tt41, first)

    write(NPSS_TS, "111 Tt41=$Tt41; alt=$(alt_in/0.3048); M0 = $MN_in;first=$first;\n")
    NPSS_success = parse(Bool, String(readavailable(NPSS_TS.out))) # `readavailable(stream)` blocks until data is available
    return NPSS_success
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
    NPSS_success = parse(Bool, String(readavailable(NPSS.out))) # `readavailable(stream)` is blocking only if no data is available
    return NPSS_success
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
    NPSS_success = parse(Bool, String(readavailable(NPSS.out)))
    return NPSS_success
end
