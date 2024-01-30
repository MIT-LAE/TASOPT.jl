export init_par_arrays, fill_par_entry, generate_par_indname
global par_indname

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


"""

generates a vector that maps par_array indices to their index name (declared in `index.inc`)

e.g., the call above would result in pari_indname[1] evaluating to "iifuel"
"""
function generate_par_indname(par_suffix::String)
    #ensure index vars are available
    include(joinpath(TASOPT.__TASOPTroot__, "misc/index.inc"))
    #get all variables in global scope
    all_vars = names(TASOPT,all=true)

    # Filter for variables that are par array indices
    # i.e., they start with "i"*par_suffix and evaluate to Integers
    filtered_variables = filter(var_elem -> length(String(var_elem))>1 
                                    && (String(var_elem)[1:nextind(String(var_elem),1)] == "i"*par_suffix)  #using nextind instead of 2 because Unicode chars
                                    && typeof(getfield(TASOPT,var_elem)) <: Int, all_vars)

    #get the indices and the maximum to initialize the output array
    # while preserving the index position
    filtered_values = [getfield(TASOPT, sym) for sym in filtered_variables]
    max_value = maximum(filtered_values)

    #initialize the Vector of Strings that will be output
    par_indname = fill("",max_value)

    #populate the Vector with variable names of the indices
    #e.g., if parvar is Symbol("iifuel"), par_indname[1] = "iifuel"
    for parvar in filtered_variables
        if String(parvar) == "i"*par_suffix*"total"
            continue
        end
        par_indname[getfield(TASOPT,parvar)] = String(parvar)
    end

    return par_indname

end

function generate_par_indname(par_suffix::AbstractVector{String}=["i","g","m","a","e"])
    # generate_par_indname.(par_suffix) 
    #^ this would be the def'n if we could get the global to work

    #initialize output dictionary
    output = Dict()
    #populate from provided suffixes
    for par_suff in par_suffix
        #for each suffix, get the par()_indname array
        output["par"*par_suff] = generate_par_indname(par_suff)
    end
    return output
end