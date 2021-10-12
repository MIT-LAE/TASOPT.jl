# Runs existing models
function load_aircraft(name, dir, saveOD; run = false)
    include("index.inc")
    include(dir *"/"* name)
    fig = plot_details(parg, pari, para, parm)
    if run
        run_wsize(35, 0, false, true, saveOD)
    end

end