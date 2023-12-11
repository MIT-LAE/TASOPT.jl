using TASOPT
using TOML
export quicksave_aircraft, quickload_aircraft


"""
    quicksave_aircraft(ac::TASOPT.aircraft, datafile)



TODO: write docstring lol
"""
function quicksave_aircraft(ac::TASOPT.aircraft=TASOPT.read_aircraft_model(),
    datafile=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_quicksave_aircraft.toml"))

    #convert structure into dictionary
    acdict = Dict(fieldnames(typeof(ac)) .=> getfield.(Ref(ac), fieldnames(typeof(ac))))

    #check dict values, replace 2d arrays (toml-unsupported)
    # with nested arrays via helper fxn
    acdict_fixed = fix_dict_for_toml(acdict)

    #generate io file
    open(datafile, "w") do io
      #write out aircraft-turned-dictionary
        TOML.print(io,acdict_fixed)
    end
end

"""
    quickload_aircraft()

TODO: write docstring
"""
function quickload_aircraft(datafile=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_quicksave_aircraft.toml"))
    data = TOML.parsefile(datafile)

    # for (key, value) in data
    #     if typeof(value) == String
    #         value = repr(value)
    #     end
    #     eval(Meta.parse("const $key = nestedvector2array($value)"))
    # end

    @unpack name, description, pari, parg, parm, para, pare, sized = data

    parm = nestedvectors2array(parm) 
    para = nestedvectors2array(para) 
    pare = nestedvectors2array(pare)

    return TASOPT.aircraft(name, description,
        pari, parg, parm, para, pare, sized)
end

function print_nested_dict(dict, indent = "    ")
    for (key, value) in dict
        if isa(value, Dict)
            println("$indent$key (Type: Dict)")
            print_nested_dict(value, "$indent  ")
        else
            value_type = typeof(value)
            println("$indent$key (Type: $value_type): $value")
        end
    end
end

#helper fxn: converts matrix into array of arrays
# function slicematrix(A::AbstractMatrix)
#     return [A[i, :] for i in 1:size(A,1)]
# end
"""

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

"""
function nestedvectors2array(nested_vectors)
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


"""

helper fxn: 
inspects dictionary for 2d arrays and replaces them with nested arrays (for toml-compatibility)
called recursively for internal dicts


"""
function fix_dict_for_toml(dict::Dict)
    #deep copy of dictionary
    dict = deepcopy(dict)
    for (key, value) in dict
        if isa(value, Array) && ndims(value) >= 2
            dict[key] = array2nestedvectors(value)
        elseif isa(value, Dict)
            dict[key] = fix_dict_for_toml(value)
        end
    end
    return dict
end
