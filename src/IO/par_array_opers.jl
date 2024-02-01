export init_par_arrays, fill_par_entry, generate_par_indname
global par_indname
export array2nestedvectors, nestedvectors2array

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

"""
    array2nestedvectors(A::Array)

Converts any n-dimensional array to a nested Vector representation and returns it.
Used in .csv and .toml IO.
"""
function array2nestedvectors(A::Array)
    if ndims(A) == 1
        # If the array is 1D, return it as is
        return A
    else
        # For higher-dimensional arrays, recursively convert subarrays
        return [array2nestedvectors(copy(subarr)) for subarr in eachslice(A, dims = 1)]
    end
end

"""
    nestedvectors2array(nested_vectors::AbstractVector)

Converts vectors nested within vectors to a multi-dimensional array.
Used in .toml IO.
"""
function nestedvectors2array(nested_vectors::AbstractVector)
    #identify how many levels deep it goes, obtaining lengths
    first_elem = nested_vectors[1]
    lengths = [length(nested_vectors)]
        #while the vector's first element is /still/ a vector
        while typeof(first_elem) <: Vector
            #save length and
            #go a level deeper
            push!(lengths,length(first_elem))
            first_elem = first_elem[1]
        end #once the first element isn't a vector (i.e., reached bottom level)

    #recursively generate flat_vector
    flat_vector = []
    function flatten_vector(v)
        for el in v
            if el isa Vector
                flatten_vector(el)
            else
                push!(flat_vector, el)
            end
        end
    end
    flatten_vector(nested_vectors)

    #reshape flat vector with length information 
    #s.t. output is an Array{Float64,n} where n = length(lengths)
    # result_array = reshape(flat_vector, reverse(lengths)...)
    result_array = permutedims(reshape(flat_vector, reverse(lengths)...),length(lengths):-1:1)

    return result_array
end