using TASOPT
include(joinpath(TASOPT.__TASOPTroot__, "misc/index.inc"))
export read_mdl

"""
    read_mdl(filepath::String, n_flight_pts::Integer=18, n_missions::Integer = 5)

TODO: write docstring lol
"""
function read_mdl(filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/B77-300ER_edit.mdl"),
    n_flight_pts::Integer=18, n_missions::Integer = 5)
    
    #generate global par_ arrays
    init_par_arrays(n_flight_pts=n_flight_pts, n_missions=n_missions)

    #execute .mdl file to populate arrays
    include(filepath)

    #TODO: fill out empty parts of arrays? i.e., repeat same aero spec for all missions if only one aero given?

    #generate the aircraft struct and return
    name = basename(filepath)
    description = "Imported using read_mdl() from: "*filepath*". Assumed and marked unsized on import."
    sized = [false]
    return TASOPT.aircraft(name, description,
        pari, parg, parm, para, pare, sized)
end