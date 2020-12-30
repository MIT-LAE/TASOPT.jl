"""
A trivial wrapper that runs a NPSS batch file `bat_file` at the specifed `dir`

>Note: 
>Remember to use `/` as separators
"""
function turboshaft_NPSS(dir, bat_file)
    run(`cmd /c cd $dir '&&' $bat_file `)
    return nothing
end