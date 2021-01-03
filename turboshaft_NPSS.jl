"""
A trivial wrapper that runs a NPSS batch file `bat_file` at the specifed `dir`

>Note: 
>Remember to use `/` as separators
"""
function turboshaft_NPSS(dir, bat_file)
    run(`cmd /c cd $dir '&&' $bat_file `)
    return nothing
end

"""
Writes an input file for NPSS
"""
function NPSS_input(SHP_dmd, Tt41, LPC_PR, HPC_PR ; LHV = 120e6, file_name = "NPSS_Turboshaft/EngineInputs.inp")

    open(file_name, "w") do io
        println(io, "// Targets")
        println(io, "real SHP_dmd = ", SHP_dmd, ";" )
        println(io, "real Tt41    = ", Tt41, ";" )

        println(io, "// LHV based on fuel")
        println(io, "Eng.FusEng.LHV = ", LHV/430, ";")

        println(io, "// Design parameters")
        println(io, "Eng.CmpL.PRdes = ", LPC_PR, ";")
        println(io, "Eng.CmpH.PRdes = ", HPC_PR, ";")
    end

end