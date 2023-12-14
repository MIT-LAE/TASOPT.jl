export init_par_arrays, autofill_across_mission
include(joinpath(TASOPT.__TASOPTroot__, "misc/index.inc"))

"""
    init_par_arrays(;n_flight_pts::Int, n_missions::Int)

initializes global arrays pari, parg, etc., to enable certain imports (e.g., via `mdl_read`).

Inputs: 
- `n_flight_pts::Int`: number of flight points modelled
- `n_mision::Int`: number of missions modelled
"""
function init_par_arrays(;n_flight_pts::Int, n_missions::Int)
    global pari = Array{Int64, 1}(undef, iitotal)
    global parg = Array{Float64, 1}(undef, igtotal)
    global parm = Array{Float64, 2}(undef, imtotal, n_missions)
    global pare = Array{Float64, 3}(undef, ietotal, n_flight_pts, n_missions)
    global para = Array{Float64, 3}(undef, iatotal, n_flight_pts, n_missions)

end

"""
    autofill_across_mission()

copies the performance specs in pare, para to missions that have none prescribed,
taking the latest valid performance spec (i.e., performance specs are copied into blanks to their "right").

uses global par arrays
"""
function autofill_across_mission()
    
end