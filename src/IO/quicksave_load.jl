export quicksave_aircraft, quickload_aircraft

"""
    quicksave_aircraft(ac::TASOPT.aircraft=TASOPT.read_aircraft_model(),
        filepath=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_quicksave_aircraft.toml"))

Converts aircraft `struct` into a dictionary with fields as keys and saves as a toml file.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::TASOPT.aircraft`: TASOPT aircraft `struct` containing model in any state. 
    - `filepath::String`: path and name of .toml file to be written.
    **Outputs:**
    - None.
"""
function quicksave_aircraft(ac::TASOPT.aircraft=TASOPT.read_aircraft_model();
    filepath=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_quicksave_aircraft.toml"))

    #convert structure into dictionary
    acdict = Dict(fieldnames(typeof(ac)) .=> getfield.(Ref(ac), fieldnames(typeof(ac))))

    #check dict values, replace 2d arrays (toml-unsupported)
    # with nested arrays via helper fxn
    acdict_fixed = fix_dict_for_toml(acdict)

    #generate io file
    open(filepath, "w") do io
      #write out aircraft-turned-dictionary
        TOML.print(io,acdict_fixed)
    end
end

"""
    quickload_aircraft(datafile=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_quicksave_aircraft.toml"))

Reads a .toml file generated by `quicksave_aircraft()` and returns the generated `aircraft` structure.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `filepath::String`: path and name of .toml file to be written.
    **Outputs:**
    - `ac::TASOPT.aircraft`: TASOPT aircraft `struct` with fields specified by the given .toml file. 
"""
function quickload_aircraft(filepath=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_quicksave_aircraft.toml"))
    data = TOML.parsefile(filepath)

    name = data["name"]
    description = data["description"]
    pari = data["pari"]
    parg = data["parg"]
    parm = data["parm"]
    para = data["para"]
    pare = data["pare"]
    sized = data["sized"]

    parm = convert(Array{Float64},nestedvectors2array(parm)) #conversion for type-consistency w constructors
    para = convert(Array{Float64},nestedvectors2array(para))
    pare = convert(Array{Float64},nestedvectors2array(pare))

    return TASOPT.aircraft(name, description,
        pari, parg, parm, para, pare, sized)
end


#helper fxn: converts matrix into array of arrays
# function slicematrix(A::AbstractMatrix)
#     return [A[i, :] for i in 1:size(A,1)]
# end


"""
    fix_dict_for_toml(dict::Dict)

Prepares a `dict` for output into .toml by:
    - replacing any multi-dimensional arrays with nested arrays (TOML library compatibility restriction),
    - and recursively applying this for any nested `dict`s.
    - Also replaces structs with dictionaries, then recursively applying this function.
"""
function fix_dict_for_toml(dict::Dict)
    #deep copy of dictionary
    dict = deepcopy(dict)
    for (key, value) in dict
        #if it's an array, convert to nested vectors
        if isa(value, Array) && ndims(value) >= 2
            dict[key] = array2nestedvectors(value)
        #if it's another dict, pass it through this fxn
        elseif isa(value, Dict)
            dict[key] = fix_dict_for_toml(value)
        #if it's a TOML-compatible value, assign directly
        #TODO: convert structs to dictionaries for output 
        elseif (value isa TOML.Internals.Printer.TOMLValue)
            dict[key] = value
        else #if not TOML-compatible
            @warn "TOML cannot output the following `aircraft` field and will skip it: "*String(key)
            # println("Field value: ", value)
            delete!(dict, key)
        end
    end
    return dict
end
