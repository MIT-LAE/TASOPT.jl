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
function quicksave_aircraft(ac::TASOPT.aircraft=TASOPT.read_aircraft_model(), 
    filepath=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_quicksave_aircraft.toml"))

    #convert structure into dictionary
    # acdict = Dict(fieldnames(typeof(ac)) .=> getfield.(Ref(ac), fieldnames(typeof(ac))))
    ## We might want to consider something like:
    acdict = struct2dict(ac)
    ## But this works for the output portion but need to be *very* careful in reading it. 
    # The dictionary stores by sorting the hash or something so the order is going to be different
    # So the only way to effectively load this would be to enforce kwargs explicitly either 
    # by ensuring all the structs are defined with the @kwargs macro or an explicit kw constructor is provided.

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

    #TODO:
    #add fuselage and wing and v/htail
    fuse_tank = dict2struct(data["fuse_tank"], struct_keys["fuse_tank"])
    fuselage = dict2struct(data["fuselage"], struct_keys["fuselage"])
    wing = dict2struct(data["wing"], struct_keys["wing"])
    htail = dict2struct(data["htail"], struct_keys["htail"])
    vtail = dict2struct(data["vtail"], struct_keys["vtail"])

    return TASOPT.aircraft(name, description,
        pari, parg, parm, para, pare, sized,
        fuse_tank, fuselage,
        wing, htail, vtail)
end

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

function struct2dict(obj)
    type = typeof(obj)
    # println("T = $type")
    d = Dict()
    for key ∈ fieldnames(type)
        # println("field = $key")
        if fieldtype(type, key) <: Union{AbstractArray,Number,String}
            # println("pushing...")
            push!(d, key => getfield(obj, key))
        else
            # println("recursing...")
            push!(d, key => struct2dict(getfield(obj, key)))
        end
    end
    return d
end

# Dict to convert the keyword keys in TOMLs to their corresponding struct types
struct_keys = Dict( 
        "fuselage" => TASOPT.Fuselage,
        "layout"    => TASOPT.structures.FuselageLayout,
        "cabin"     => TASOPT.structures.Cabin,
        "material"  => TASOPT.structures.StructuralAlloy,
        "skin"      => TASOPT.structures.StructuralMember,
        "shell"     => TASOPT.structures.StructuralMember,
        "cone"      => TASOPT.structures.StructuralMember,
        "floor"     => TASOPT.structures.StructuralMember,
        "insulation"=> TASOPT.structures.Weight,
        "window"    => TASOPT.structures.Weight,
        "bendingmaterial_h" => TASOPT.structures.StructuralMember,
        "bendingmaterial_v" => TASOPT.structures.StructuralMember,
        "APU"       => TASOPT.structures.Weight,
        "seat"      => TASOPT.structures.Weight,
        "added_payload"     => TASOPT.structures.Weight,
        "HPE_sys"   => TASOPT.structures.Weight,
        "fixed"     => TASOPT.structures.Weight,
        "Weight"    => TASOPT.structures.Weight,
        # "weight"    => TASOPT.structures.Weight,
        "frame"     => TASOPT.structures.Frame,
        "cross_section"=> TASOPT.structures.AbstractCrossSection,

        "fuse_tank" => TASOPT.fuselage_tank,
        "inner_material"    => TASOPT.StructuralAlloy,
        "outer_material"    => TASOPT.StructuralAlloy,
        "material_insul"    => TASOPT.ThermalInsulator, #this is a vector

        "wing"      => TASOPT.structures.Wing,
        "htail"     => TASOPT.structures.Tail,
        "vtail"     => TASOPT.structures.Tail,
        "strut"     => TASOPT.structures.Strut,
        "inboard"   => TASOPT.structures.WingSection,
        "outboard"  => TASOPT.structures.WingSection,
        "airsection"=> TASOPT.aerodynamics.airfoil,
    )

"""
    dict2struct(dict::Dict, struct_type::DataType)

creates a struct from a dictionary and a struct type (recursively). 
references internal dictionary `struct_keys` for keys that should create structs.
Unfortunately, there's a lot of edge cases from inconsistent struct construction and some name repeats.
Recommend using unique names for fields that contain structs, and making all constructors keyword-based (not order-based).
"""
function dict2struct(dict::Dict, struct_type::DataType; 
    struct_context="")
    # Prepare a dictionary to hold arguments for struct construction
    struct_args = Dict{Symbol,Any}()

    #if the dict is a single key-value pair, 
    # check if the user meant to lose this outer dict layer
    # by seeing if the key matches the struct type
    # e.g., calling dict2struct(dict, TASOPT.Fuselage) where dict["fuselage"] = dict
    if length(dict) == 1
        (key, value) = first(dict)
        if struct_keys[key] == struct_type
            #if so, lose the outer layer
            dict = value
        end
    end

    #if we have multiple key-value pairs, we have fields to populate a struct
    for (key, value) in dict
       #1) sort out any structs
        # If the key corresponds to a struct per struct_keys...
        if haskey(struct_keys, key)
            #...and the value is already a struct, assign it directly
            if isa(value, struct_keys[key]) 
                struct_args[Symbol(key)] = value
            #...and the value is a dict (i.e., we have field names and values), recurse
            elseif isa(value, Dict)
                struct_args[Symbol(key)] = dict2struct(value, struct_keys[key], 
                                            struct_context = struct_context * "." * key)
            #...and the value is a list
            elseif isa(value, Array)
                #...and the list is empty, assign it directly
                if isempty(value)
                    struct_args[Symbol(key)] = value
                #...and the list is not empty, recurse
                else
                    struct_args[Symbol(key)] = [dict2struct(v, struct_keys[key],
                                            struct_context = struct_context * "." * key) for v in value]
                end
            #...and the value is anything else, throw an error
            else    
                throw(ArgumentError("Expected a dictionary for field $key, but got a $(typeof(value))"))
            end
        
       #2) sort out any arrays
        # If the value is an array, convert it (if needed) and assign it
        elseif isa(value, Array)
            #prune cases with value = Union{}[] (i.e., oddly untyped arrays)
            if length(value) == 0
                value = zeros(size(value))
            end
            #assign, converting from nested vectors to a mtulti-dim array if needed
            struct_args[Symbol(key)] = nestedvectors2array(value)

        
       #3) just save it
       # assign the value directly
        else
            struct_args[Symbol(key)] = value
        end
    end #end for loop

    # Edge case: frames - constructed without an id
    if struct_type == TASOPT.structures.Frame
        if struct_args[:origin] == [0.0, 0.0, 0.0] #usually, world frame
            return TASOPT.WORLD #avoid constructing a new frame
        else
            return TASOPT.structures.Frame(struct_args[:origin])
        end
    end

    # Edge case: cross_section
    # (Wing) or (Fuselage [Single or MultiBubble])
    if struct_type == TASOPT.structures.AbstractCrossSection
        if haskey(dict, "thickness_to_chord")
            struct_type = TASOPT.structures.WingCrossSection
        elseif haskey(dict, "n_webs")
            struct_type = TASOPT.MultiBubble
        else
            struct_type = TASOPT.SingleBubble
        end
    end

    # Edge case: layout - WingLayout or FuselageLayout since both are ___.layout
    if struct_type == TASOPT.structures.FuselageLayout
        if haskey(dict, "x_nose")
            struct_type = TASOPT.structures.FuselageLayout
        else
            struct_type = TASOPT.structures.WingLayout
        end
    end 

    # Construct the struct using the collected arguments
    return struct_type(; struct_args...)
end #end fxn