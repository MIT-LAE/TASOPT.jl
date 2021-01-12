"""
A trivial wrapper that runs a NPSS batch file `bat_file` at the specifed `dir`

>Note: 
>Remember to use `/` as separators
"""
function NPSS_run(dir, bat_file)
    run(`cmd /c cd $dir '&&' $bat_file `)
    return nothing
end

"""
# Design mode

Writes an input file for NPSS Turboshaft model
"""
function NPSS_TShaft_input(alt_in, MN_in, 
    SHP_dmd, Tt41, 
    LPC_PR, HPC_PR,
    lcat, deNOx ; PCEC = true,
    LHV = 120, file_name = "NPSS_Turboshaft/EngineInputs.inp")

    open(file_name, "w") do io
        println(io, "// Design State")
        println(io, "Eng.setOption(\"switchDes\",\"DESIGN\");")
        PCEC ? println(io, "Eng.PCEC.setOption(\"switchMode\",\"ON\");") : 

        println(io, "\n// Abmient conditions")
        println(io, "Eng.Amb.alt_in = ", alt_in, ";" )
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

    end

end

"""
OffdesMode
"""
function NPSS_TShaft_input(alt_in, MN_in, 
                            SHP_dmd ; 
                            file_name = "NPSS_Turboshaft/EngineInputs.inp",
                            desScl_name = "DesScl.int")

    open(file_name, "w") do io
        println(io, "// Design State")
        println(io, "Eng.setOption(\"switchDes\",\"OFFDESIGN\");")
        println(io, "Eng.PCEC.setOption(\"switchMode\",\"ON\");") 
        
        println(io, "\n// Abmient conditions")
        println(io, "Eng.Amb.alt_in = ", alt_in, ";" )
        println(io, "Eng.Amb.MN_in  = ", MN_in, ";" )

        println(io, "\n// Targets")
        println(io, "real SHP_dmd = ", SHP_dmd, ";" )

        println(io, "\n// DesPt scalars")
        println(io, "#include \"", desScl_name, "\"")
        
    end

end

"""
Writes an input file for NPSS ducted fan model
"""
function NPSS_Fan_input(alt_in::Float64, MN_in::Float64, Fn::Float64, πfan::Float64 ; file_name = "NPSS_Turboshaft/FanInputs.inp")

    open(file_name, "w") do io
        println(io, "// Design State")
        println(io, "DuctedFan.setOption(\"switchDes\",\"DESIGN\");")

        println(io, "\n// Abmient conditions")
        println(io, "DuctedFan.Amb.alt_in = ", alt_in, ";" )
        println(io, "DuctedFan.Amb.MN_in  = ", MN_in, ";" )

        println(io, "\n// Thrust Target")
        println(io, "real Fn_target = ", Fn, ";")

        println(io, "\n// Design parameters")
        println(io, "DuctedFan.Fan.PRdes = ", πfan, ";")
    end

end

function NPSS_Fan_input(alt_in::Float64, MN_in::Float64, Fn::Float64,
                        MapScalars::Array{Float64, 1}, NozArea::Float64; 
                        file_name = "NPSS_Turboshaft/FanInputs.inp")

    open(file_name, "w") do io
        println(io, "// Design State")
        println(io, "DuctedFan.setOption(\"switchDes\",\"OFFDESIGN\");")

        println(io, "\n// Abmient conditions")
        println(io, "DuctedFan.Amb.alt_in = ", alt_in, ";" )
        println(io, "DuctedFan.Amb.MN_in  = ", MN_in, ";" )

        println(io, "\n// Thrust Target")
        println(io, "real Fn_target = ", Fn, ";")

        println(io, "\n// Map scalars")
        println(io, "DuctedFan.Fan.S_map.s_effDes = ", MapScalars[1], ";")
        println(io, "DuctedFan.Fan.S_map.s_PRdes  = ", MapScalars[2], ";")
        println(io, "DuctedFan.Fan.S_map.s_WcDes  = ", MapScalars[3], ";")
        println(io, "DuctedFan.Fan.S_map.s_NcDes  = ", MapScalars[4], ";")

        println(io, "\n// Nozzle Area")
        println(io, "DuctedFan.FanNozzle.AthCold  = ", NozArea, ";")


    end

end