export init_par_arrays, fill_par_entry
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
    fill_par_entry(value_or_values, n_ip, n_im)

repeats given value or vector of values to fill a 2d array of size (n_ip, n_im)
used for assigning constants across flight segments conveniently.

accounts for multi-mission input (read as array in value_or_values) by reducing repeats if present
"""
#TODO: make mission repeat assume something when nmisx != size(value_or_values,1)
function fill_par_entry(value_or_values, n_ip, n_im)
    return repeat(transpose(value_or_values), n_ip, 1 + (n_im-size(value_or_values,1)))
end

# Overloads repeat fxn to handle single elements (not default behavior)
# (edge case for some multi-point/multi-mission handling in read_input)
import Base.repeat
function repeat(single_elem::Real, counts::Integer...)
    repeat([single_elem], counts...)
end