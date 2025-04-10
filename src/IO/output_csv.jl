export output_csv, default_output_indices
export output_indices_all, output_indices_wGeom, output_indices_wEngine

"""
    output_csv(ac::TASOPT.aircraft=TASOPT.load_default_model(), 
            filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/default_output.csv");
            overwrite::Bool = false, indices::Dict = default_output_indices,
            includeMissions::Union{AbstractVector,Colon,Bool,Integer} = false, 
            includeFlightPoints::Union{AbstractVector,Colon,Bool,Integer} = false,
            forceMatrices::Bool = false)

writes the values of `ac` to CSV file `filepath` with index variables as headers. 
A typical set of values is output by default for the design mission at the first cruise point.
Appends to extant `filepath` if headers are compatible, appending integer suffixes to filename when not.

Output is customizable by:
  - `indices`: specifies the desired indices of each par array,
  - `includeMissions`: allows output of all missions (i.e., =true) or specifiable indices (e.g., =[1,2,3]),
  - `includeFlightPoints`: allows output of all flight points (as for `includeMissions`).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::TASOPT.aircraft`: TASOPT aircraft `struct` containing model in any state. 
    - `filepath::String`: path and name of .csv file to be written.
    - `overwrite::Bool`: deletes existing file at filepath when true, default is false.
    - `indices::Dict{String => Union{AbstractVector,Colon(), Integer}}`: specifies desired indices of par arrays given as keys. Customizable; built-in options: [`default_output_indices`](@ref), `output_indices_all`, and `output_indices_wEngine`.
    - `includeMissions::Union{AbstractVector,Colon,Bool,Integer}`: saves all mission entries as an array in a CSV cell when true, default is false, inner nested array when flight points are also output. specific indices can also be specified as Vectors of Ints.
    - `includeFlightPoints::Union{AbstractVector,Colon,Bool,Integer}`: saves all flight point entries as an array in a CSV cell when true, default is false, outer nested array when missions are also output. specific indices can also be specified as Vectors of Ints.
    - `forceMatrices::Bool`: forces all entries that vary with mission and flight point to follow nested array structure
    - `struct_excludes::AbstractVector{String}`: names/substrings of fields to exclude from output of ac struct. Default is `[]`, excluding nothing.
    **Outputs:**
    - `newfilepath::String`: actual output filepath; updates in case of header conflicts. same as input filepath if `overwrite = true`.
"""
function output_csv(ac::TASOPT.aircraft=TASOPT.load_default_model(), 
    filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/default_output.csv");
    overwrite::Bool = false, indices::Dict = default_output_indices,
    includeMissions::Union{AbstractVector,Colon,Bool,Integer} = false, 
    includeFlightPoints::Union{AbstractVector,Colon,Bool,Integer} = false,
    forceMatrices::Bool = false,
    struct_excludes::AbstractVector = [])

    #initialize row to write out
    csv_row = []
    #initalize header row to assess compatibility, and maybe write out
    header = ["Model Name","Description","Sized?"]
    par_array_names = ["parg", "parm", "para", "pare"] #used for checking/populating indices Dict, deleting from struct

    #pull name and description, removing line breaks and commas
    append!(csv_row, [ac.name, 
                     replace(ac.description, r"[\r\n,]" => " "),
                     ac.is_sized[1]])

    #if the entire struct is not excluded from output
    if length(intersect(struct_excludes, ["aircraft","ac", "all"])) < 1
        #flatten the struct into a list of names and values
        titles, values = flatten_struct(ac) 
        #remove par arrays from titles and values (the contents will be appended later)
        par_array_indices = findall(title -> any( contains.(title, par_array_names) ), titles)
        titles = deleteat!(titles, par_array_indices)
        values = deleteat!(values, par_array_indices)
        #remove any fields that contain the substrings in struct_excludes
        if !isempty(struct_excludes)
            exclude_indices = findall(title -> any( contains.(title, struct_excludes) ), titles)
            titles = deleteat!(titles, exclude_indices)
            values = deleteat!(values, exclude_indices)
        end
        #append the rest of the struct to the row
        append!(csv_row, values)
        append!(header, titles)
    end

    #if all outputs are desired, map all par arrays to colon() for data pull below
    if indices in [Colon(), Dict()] #if indices is Colon() or an EMPTY dict()
        indices = Dict(name => Colon() for name in par_array_names)
    #if desired output indices are indicated...
    elseif indices isa Dict
        for indexkey in keys(indices)
            # ...if key is not in valid names dict, remove it
            if !(indexkey in par_array_names)
                @warn "indices Dict for csv_output() contains an unparsable key: "*String(indexkey)*". Removing key and continuing."
                delete!(indices,indexkey)
            # ...else if key is not Colon(), ensure there are no repeats
            elseif indices[indexkey] != Colon()
                #remove duplicate indices in array because CSV.jl can't handle them
                oldlen = length(indices[indexkey])
                indices[indexkey] = union(indices[indexkey])
                if length(indices[indexkey]) < oldlen
                    @warn "duplicate indices identified and removed from indices Dict["*String(indexkey)*"]" 
                end
            end
        end
    else #complain if not as expected
        throw(ArgumentError("indices parameter must be Dict or Colon()"))
    end

  #mission and flight point indexing logic
    #default indices (design mission, cruise point)
    default_imiss = 1 
    default_ipts = ipcruise1

    #convert boolean selections into index arrays 
    # append all mission/flight points (i.e., Colon()) if include is true, 
    # else just first mission (design)/cruise point
    if includeMissions isa Bool
        includeMissions = includeMissions ? Colon() : default_imiss
    end
    if includeFlightPoints isa Bool
        includeFlightPoints = includeFlightPoints ? Colon() : default_ipts
    end
    #if desired, force matrix output when value is mission- and flight-point-dependent
    # (prevent flattening into 1d vector if only 1 index is requested; done by creating a range)
    if forceMatrices
        #convert length-1 vectors to integers
        if (includeMissions isa Vector) && (length(includeMissions) == 1); includeMissions = includeMissions[1]; end
        if (includeFlightPoints isa Vector) && (length(includeFlightPoints) == 1); includeFlightPoints = includeFlightPoints[1]; end

        #convert integers to ranges to force matrix output
        if includeMissions isa Integer; includeMissions = (includeMissions:includeMissions); end
        if includeFlightPoints isa Integer; includeFlightPoints = (includeFlightPoints:includeFlightPoints); end
    end

    #append mission- and flight-point-independent arrays
    append!(csv_row, ac.parg[indices["parg"]])

    #append mission-dependent arrays
    append!(csv_row, array2nestedvectors(ac.parm[indices["parm"], includeMissions])) #append to csv_rows

    #append flight-point- and mission-dependent arrays
    if (indices["para"]==Colon()) || !isempty(indices["para"])    #explicit check required, lest an empty [] be appended
        append!(csv_row,array2nestedvectors(ac.para[indices["para"], includeFlightPoints , includeMissions]))
    end
    if (indices["pare"]==Colon()) || !isempty(indices["pare"])    # ( " )
        append!(csv_row,array2nestedvectors(ac.pare[indices["pare"], includeFlightPoints , includeMissions]))
    end

    #jUlIa iS cOluMn MaJoR
    csv_row=reshape(csv_row,1,:) 

    #get lookup dict for column names from variable indices in indices Dict()
    suffixes = ["g","m","a","e"]
    par_indname_dict = generate_par_indname(suffixes) #gets dictionary with all index var names via global scope

    #get column names in order
    for suffix in suffixes
        append!(header,par_indname_dict["par"*suffix][indices["par"*suffix]])
    end
    
    #if overwrite indicated (not default), do so
    if overwrite
        newfilepath = filepath
        bool_append = false
        if isfile(filepath)
            rm(filepath)
            @info "overwriting "*filepath
        end

    # else, assign new filepath, checking header append-compatibility
    # i.e, if can append here, do so. else, add a suffix and warn
    else 
        newfilepath, bool_append = check_file_headers(filepath, header)
    end

    io = open(newfilepath,"a+")
    CSV.write(io, Tables.table(csv_row), header = header, append = bool_append,
                    delim=',', quotechar='"')
    close(io)
    return newfilepath
end



function flatten_struct(s, prefix="")
    T = typeof(s)
    #get all field names and values of current struct
    field_names = String.(collect(fieldnames(T))) 
    field_values = [getfield(s, f) for f in fieldnames(T)]

    flat_names = []
    flat_values = []
    
    #iterate over field-value pairs
    for (name, value) in zip(field_names, field_values)
        full_name = prefix * name #add prefix to name if this struct is a field of another struct
        #if value is a struct, call recursively
        if isstructtype(typeof(value)) && !(value isa AbstractArray)
            sub_names, sub_values = flatten_struct(value, full_name * ".")
            append!(flat_names, sub_names)
            append!(flat_values, sub_values)
        else
            push!(flat_names, full_name)
            push!(flat_values, value)
            # push!(flat_values, string(value))   #ensure value is a string for CSV writing - unclear if necessary
        end
    end

    return flat_names, flat_values
end

"""
    check_file_headers(filepath, header)

checks a filepath for existence. returns suffixed filename if extant .csv contains 
    column headers that conflict with header parameter. recursively called if suffix is needed/added

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `filepath::String`: path and name of .csv file to be checked for header consistency.
    - `header::Vector{String}: column headers intended for output, will be checked with filepath's header`
    **Outputs:**
    - `newfilepath::String`: path and name of .csv file cleared for writing.
    - `bool_append::Bool`: true if newfilepath is appendable; false if newfilepath is new and needs headers written

"""
function check_file_headers(filepath, header)
    #does file already exist?
    #if not
    if !isfile(filepath)    
        @info "output to .csv creating new file: "*basename(filepath)
        return filepath, false  #write to the given filepath, append = false
    else #if file already exists
        #check if headers match, 
         #replacing any default Column names with blanks to match
        header_old = replace.(String.(propertynames(CSV.File(filepath))), 
                               r"Column\d+" .=> "")
        
        if header == header_old #if so, simple append to given filepath, append = true
            @info "output to .csv appending to: "*basename(filepath)
            return filepath, true
        else    #if not, add a suffix to the filename
            #split filepath into components
            dir, file = splitdir(filepath)
            filebase, ext = split(file,".")
   
            #detect if end of filename is already an int suffix
            regexp_pattern = r"_(\d+)$"
            matches = match(regexp_pattern, filebase)
           
            #if not, add suffix "_1"
            if matches===nothing 
               suffix = "_1"
               filebase = filebase*suffix
            #if so, increment suffix by 1
            else
                int_suffix = parse(Int, matches[1])
                new_int_suffix = int_suffix + 1
                filebase = replace(filebase, regexp_pattern => "_$new_int_suffix")
            end

            #collate new path with modified filebase and check THAT file recursively
            newfilepath = joinpath(dir,filebase*"."*ext)
            newfilepath, bool_append = check_file_headers(newfilepath,header)

            return newfilepath, bool_append
        end
    end
end

"""
Indices of quantities of par arrays selected to be output by default by [`output_csv()`](@ref).
Formatted as a Dict() where keys are par array names and values are arrays of indices.
Note that this selection is only along the first dimension of par arrays (i.e., selecting quantities, not missions or flight points).

Custom Dicts() can be passed to [`output_csv()`](@ref) following the format: `Dict(){String => Union{AbstractVector,Colon(), Integer}}`. 

For example:
1 : `default_output_indices["parg"] = [1, 5, igtotal]`
        > when submitted to `output_csv()`, will output as columns the first, fifth, and last entry of `parg`
    
2 : `default_output_indices["para"] = Colon() #NOT [Colon()]`
        > when submitted to `output_csv()`, will output all indices of `para` (whether all flights or missions are included is controlled by a separate parameter)
"""
default_output_indices = 
    Dict("parg" => [#weights
                    igWMTO, igWfuel, igWpay, igWpaymax,
                    igWtesys, igWftank, 
                    #other
                    igdfan, igGearf,
                    #performance, different from im?
                    igRange
                    ],
        "parm" => [imRange, imWpay, imWTO, imWfuel, imPFEI, 
                    ],
        "para" => [iaalt, iaMach, iaReunit,iagamV,      #flight conditions
                    iaCL, iaCD, iaCDi, iaCLh,iaspaneff, #performance
                    iaxCG, iaxCP,                       #balance
                    ],
        "pare" => [iehfuel, ieTfuel, ieff, ieTSFC, ieBPR, ieOPR, iepif, iepilc, iepihc
        ])

output_indices_wEngine = deepcopy(default_output_indices)
#append engine params
pareinds = output_indices_wEngine["pare"]
pare_toadd = [
            ieFe, ieFsp,
            ieN1, ieN2,#spool speeds
            #core flow
            ieTt0,ieTt19,ieTt3,ieTt4,ieTt45, ieTt49, ieTt5, 
            iept0,iept19,iept3,iept4,iept45, iept49, iept5, 
            #bypass flow
            ieTt2, ieTt21, ieTt9,
            iept2, iept21, iept9,
            
            iedeNOx, iemdotf, ieEINOx1, ieEINOx2,
            ieFAR, ieOPR
            ]
append!(pareinds, pare_toadd)
output_indices_wEngine["pare"] = pareinds

output_indices_all =  Dict("parg" => Colon(),
                            "parm" => Colon(),
                            "para" => Colon(),
                            "pare" => Colon())