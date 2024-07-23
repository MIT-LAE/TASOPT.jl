export output_csv, default_output_indices
export output_indices_all, output_indices_wGeom, output_indices_wEngine

"""
    output_csv(ac::TASOPT.aircraft=TASOPT.load_default_model(), 
            filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_output.csv");
            overwrite::Bool = false, indices::Dict = default_output_indices,
            includeMissions::Union{AbstractVector,Colon,Bool,Integer} = false, 
            includeFlightPoints::Union{AbstractVector,Colon,Bool,Integer} = false,
            forceMatrices::Bool = false)

writes the values of `ac` to CSV file `filepath` with index variables as headers. 
A typical set of values is output by default for the design mission at the first cruise point.
Appends to extant `filepath` if headers are compatible, appending integer suffixes to filename when not.

Output is customizable by:
    - `indices` dictionary mapping par arrays to desired indices,
    - `includeMissions` allows output of all missions (i.e., =true) or specifiable indices (e.g., =[1,2,3]),
    - `includeFlightPoints` allows output of all flight points (as for `includeMissions`).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::TASOPT.aircraft`: TASOPT aircraft `struct` containing model in any state. 
    - `filepath::String`: path and name of .csv file to be written.
    - `overwrite::Bool`: deletes existing file at filepath when true, default is false.
    - `indices::Dict{String => Union{AbstractVector,Colon(), Integer}}`: See [`default_output_indices`](@ref)
    - `includeMissions::Union{AbstractVector,Colon,Bool,Integer}`: saves all mission entries as an array in a CSV cell when true, default is false, inner nested array when flight points are also output. specific indices can also be specified as Vectors of Ints.
    - `includeFlightPoints::Union{AbstractVector,Colon,Bool,Integer}`: saves all flight point entries as an array in a CSV cell when true, default is false, outer nested array when missions are also output. specific indices can also be specified as Vectors of Ints.
    - `forceMatrices::Bool`: forces all entries that vary with mission and flight point to follow nested array structure
    **Outputs:**
    - `newfilepath::String`: actual output filepath; updates in case of header conflicts. same as input filepath if `overwrite = true`.

TODO:
- human-readable column names - lookup with /another/ dict
    currently column headers are index vars - e.g., iifuel, ieA5fac
    desire human readable. can generate a lookup dict: index2string["pari"][iifuel] = "fuel_type"
- bug: bonus rows of last #? (comes from mismatch of headers and cells)
"""
function output_csv(ac::TASOPT.aircraft=TASOPT.load_default_model(), 
    filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_output.csv");
    overwrite::Bool = false, indices::Dict = default_output_indices,
    includeMissions::Union{AbstractVector,Colon,Bool,Integer} = false, 
    includeFlightPoints::Union{AbstractVector,Colon,Bool,Integer} = false,
    forceMatrices::Bool = false)

    #initialize row to write out
    csv_row = []

    #pull name and description, removing line breaks and commas
    append!(csv_row, [ac.name, 
                     replace(ac.description, r"[\r\n,]" => ""),
                     ac.sized[1]])

    par_array_names = ["pari", "parg", "parm", "para", "pare"] #used for checking/populating indices Dict
    #if all outputs are desired, map all par arrays to colon() for data pull below
    if indices in [Colon(), Dict()]
        indices = Dict(name => Colon() for name in par_array_names)
    #if desired output indices are indicated, check inputs
    elseif indices isa Dict
        for indexkey in keys(indices)
            if !(indexkey in par_array_names)
                @warn "indices Dict for csv_output() contains an unparsable key: "*String(indexkey)*". Removing key and continuing."
                delete!(indices,indexkey)
            elseif indices[indexkey] != Colon() #if indexkey is okay and not Colon()
                #remove duplicate indices in array because CSV.jl can't handle them
                oldlen = length(indices[indexkey])
                indices[indexkey] = union(indices[indexkey])
                if length(indices[indexkey]) < oldlen
                    @warn "duplicate indices identified and removed from indices Dict["*String(indexkey)*"]" 
                end
            end
        end
    else
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
    append!(csv_row,ac.pari[indices["pari"]],
                    ac.parg[indices["parg"]])

    #append mission-dependent arrays
    append!(csv_row,array2nestedvectors(ac.parm[indices["parm"], includeMissions]))

    if (indices["para"]==Colon()) || !isempty(indices["para"])    #explicit check required, lest an empty [] be appended
        append!(csv_row,array2nestedvectors(ac.para[indices["para"], includeFlightPoints , includeMissions]))
    end
    if (indices["pare"]==Colon()) || !isempty(indices["pare"])    # ( " )
        append!(csv_row,array2nestedvectors(ac.pare[indices["pare"], includeFlightPoints , includeMissions]))
    end

    #jUlIa iS cOluMn MaJoR
    csv_row=reshape(csv_row,1,:) 

    #get lookup dict for column names from variable indices in indices Dict()
    suffixes = ["i","g","m","a","e"]
    par_indname_dict = generate_par_indname(suffixes) #gets dictionary with all index var names via global scope

    #get column names in order
    header = ["Model Name","Description","Sized?"]
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
    CSV.write(io, Tables.table(csv_row), header = header, append = bool_append)
    close(io)
    return newfilepath
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

#TODO: include this via ```@docs default_output_indices```
"""
par arrays and indices output by default by `output_csv()`.

Format: Dict(){String => Union{AbstractVector,Colon(), Integer}}

    use ex.1 : default_output_indices["pari"] = [1, 5, iitotal]
    
    use ex.2 : default_output_indices["para"] = Colon() #NOT [Colon()]
"""
default_output_indices = 
    Dict("pari" => Colon(),
        "parg" => [#weights
                    igWMTO, igWfuel, igWpay, igWpaymax, igflgnose,igflgmain,
                    igWweb,igWcap, igfflap,igfslat, igfaile, igflete, igfribs, igfspoi, igfwatt,
                 igWwing, igWvtail, igWhtail, igWtesys, igWftank, 
                    
                    #other
                    igdfan, igGearf,

                    #performance, different from im?
                    igPFEI, igRange
                    ],
        "parm" => [imRange, imWpay, imWTO, imWfuel, imPFEI, 
                    ],
        "para" => [iaalt, iaMach, iaReunit,iagamV,      #flight conditions
                    iaCL, iaCD, iaCDi, iaCLh,iaspaneff, #performance
                    iaxCG, iaxCP,                       #balance
                    ],
        "pare" => [iehfuel, ieTfuel, ieff, ieTSFC, ieBPR, ieOPR, iepif, iepilc, iepihc
        ])
                     
output_indices_wGeom = deepcopy(default_output_indices)
#append geometry params
parginds = output_indices_wGeom["parg"]
parg_toadd = [igAR, igS, igb, igbo, igbs, igco, iglambdat, iglambdas, igsweep, #wing geom
                #stabilizers geom
                igVh, igARh, igSh, igbh, igboh, iglambdah, igcoh, igsweeph, 
                igVv, igARv, igSv, igbv, igbov, iglambdav, igcov, igsweepv,]
append!(parginds, parg_toadd)
output_indices_wGeom["parg"] = parginds

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

output_indices_all =  Dict("pari" => Colon(),
                            "parg" => Colon(),
                            "parm" => Colon(),
                            "para" => Colon(),
                            "pare" => Colon())