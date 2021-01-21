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

    - Abmient altitude and mach number
    - Specifies πᵢ ∀ i ∈ {HPC, LPC}
    - Shaft power demand, Tt41 
    - DeNOx target
    - SCR paramters -> CPSI, w, l

"""
function NPSS_TShaft_input(alt_in, MN_in, 
    SHP_dmd, Tt41, 
    LPC_PR, HPC_PR,
    cpsi, w, lcat, deNOx ;
    LHV = 120, file_name = "NPSS_Turboshaft/EngineInputs.inp")

    open(file_name, "w") do io

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
                            Tt41, Nshaft ; 
                            file_name = "NPSS_Turboshaft/OffDesInputs.inp")

    open(file_name, "w") do io

        println(io, "\n// Abmient conditions")
        println(io, "Eng.Amb.alt_in = ", alt_in, ";")
        println(io, "Eng.Amb.MN_in  = ", MN_in , ";")

        println(io, "\n// Targets")
        println(io, "real Tt41   = ", Tt41  , ";")
        println(io, "real N2_dmd = ", Nshaft, ";")
        
    end

end

"""
Writes an input file for NPSS ducted fan model
"""
function NPSS_Fan_input(alt_in::Float64, MN_in::Float64, Fn::Float64,
                        Kinl::Float64, Φinl::Float64,
                        πfan::Float64 ; file_name = "NPSS_Turboshaft/FanInputs.inp")

    open(file_name, "w") do io
        println(io, "// Design State")
        # println(io, "DuctedFan.setOption(\"switchDes\",\"DESIGN\");")

        println(io, "\n// Abmient conditions")
        println(io, "DuctedFan.Amb.alt_in = ", alt_in, ";" )
        println(io, "DuctedFan.Amb.MN_in  = ", MN_in, ";" )

        println(io, "\n// Thrust Target")
        println(io, "real Fn_target = ", Fn, ";")

        println(io, "\n// BLI inputs")
        println(io, "DuctedFan.InEng.Kinl = ", Kinl, ";")
        println(io, "DuctedFan.Phiinl     = ", Φinl, ";")

        println(io, "\n// Design parameters")
        println(io, "DuctedFan.Fan.PRdes = ", πfan, ";")
    end

end

function NPSS_Fan_input(alt_in::Float64, MN_in::Float64, Pin::Float64,
                        Kinl::Float64, Φinl::Float64;
                        file_name = "NPSS_Turboshaft/FanInputs.inp")

    open(file_name, "w") do io
        println(io, "// Design State")
        # println(io, "DuctedFan.setOption(\"switchDes\",\"OFFDESIGN\");")

        println(io, "\n// Abmient conditions")
        println(io, "DuctedFan.Amb.alt_in = ", alt_in, ";" )
        println(io, "DuctedFan.Amb.MN_in  = ", MN_in, ";" )

        println(io, "\n// Input power")
        println(io, "real ShP_input = ", -Pin/745.7, ";")

        
        println(io, "\n// BLI inputs")
        println(io, "DuctedFan.InEng.Kinl = ", Kinl, ";")
        println(io, "DuctedFan.Phiinl     = ", Φinl, ";")


    end

end