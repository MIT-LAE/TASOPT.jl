include(__TASOPTindices__)
export read_mdl, read_jl

"""
    read_mdl(filepath::String, n_flight_pts::Integer=18, n_missions::Integer = 5)

includes old-fashioned `.mdl` files; generates and returns an `aircraft` struct, assumed unsized
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
    is_sized = [false]
    return TASOPT.aircraft(name, description,
        pari, parg, parm, para, pare, is_sized)
end

"""
    read_jl(filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/777TFinput_edited.jl"))

includes `.jl` files; generates and returns an `aircraft` struct, assumed unsized
"""
function read_jl(filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/777TFinput_edited.jl"))
    #execute .jl file to populate arrays
    include(filepath)

    #generate the aircraft struct and return
    name = basename(filepath)
    description = "Imported using read_jl() from: "*filepath*". Assumed and marked unsized on import."
    is_sized = [false]
    return TASOPT.aircraft(name, description,
        pari, parg, parm, para, pare, is_sized)
end