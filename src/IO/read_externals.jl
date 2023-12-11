using TASOPT
using TOML
include(joinpath(TASOPT.__TASOPTroot__, "misc/index.inc"))

"""
    read_mdl(filepath::String)



TODO: write docstring lol
"""
function read_mdl(filepath::String=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/A320_L1B.mdl"),
        max_arraydim_capacity::Int=1000)

    mdl_file_lines = readlines(filepath)

    pari = Array{Int64, 1}(undef, iitotal)
    parg = Array{Float64, 1}(undef, igtotal)
    parm = Array{Float64, 2}(undef, imtotal, max_arraydim_capacity)
    pare = Array{Float64, 3}(undef, iatotal, max_arraydim_capacity, max_arraydim_capacity)
    para = Array{Float64, 3}(undef, ietotal, max_arraydim_capacity, max_arraydim_capacity)

    for line in mdl_file_lines
        # Extract matches from the string using the regular expression
        # pattern = r"([a-zA-Z_][a-zA-Z0-9_]*)\[(.*?)\].*?\.=(.*?)\s*"
        # pattern = r"([a-zA-Z_][a-zA-Z0-9_]*)\[(.*?)\]\s*\.=\s*(\[.*?\])\s*"
        # match_ = match(pattern, line)
        # if match_ !== nothing
        #     variable_name = match_.captures[1]
        #     index = parse.(Int, split(match_.captures[2], ','))
        #     value = parse.(Float64, split(match_.captures[3], ','))

        #     # update the variable at index with value
        #     variable = get!(variables, variable_name, zeros(length(index)...))
        #     setindex!(variable, value, index...)

        if line[1:5] in ["para[","pare[","parg[","parm[","pari["]
            eval(Meta.parse(line))


        end
    end

    name = basename(filepath)
    description = "Imported by read_mdl() from: "*filepath
    sized = [false] #assume unsized. may be sized but won't guess

    return TASOPT.aircraft(name, description,
        pari, parg, parm, para, pare, sized)

end

a = read_mdl()
print(a)