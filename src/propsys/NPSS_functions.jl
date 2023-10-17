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
Starts up and returns an NPSS process that can then be written to.
"""
function startNPSS(dir, bat_file)
    if Sys.iswindows()
        NPSS = open(`cmd /c cd $dir '&&' $bat_file `, "w+")
    elseif Sys.islinux()
        cd(dir) # quick fix cant figure out how to cd in open()
        NPSS = open(`bash $bat_file `, "w+")
        cd("../")
    end
    return NPSS
end

"""
Ends the NPSS process.
"""
function endNPSS(NPSS)
    write(NPSS, "999 \\n")
    if(NPSS.exitcode == 0)
        close(NPSS)
    end
end

