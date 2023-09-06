using TOML
include("../misc/index.inc")
include("../misc/constants.jl")
include("../misc/units.jl")

data = TOML.parsefile("src/IO/input.toml")
default = TOML.parsefile("src/IO/default_input.toml")

"""
    read_input(k::String, dict; default_dict = default)

Reads the input from a given dictonary (typically parsed from a TOML file).
If requested input does not exist in dictonary, looks for value in default input
and stores default value into the given dictonary (primarily for later output/
saving as an aircraft model file)
"""
function read_input(k::String, dict=data, default_dict = default)
    get!(dict, k) do 
        if k in keys(default_dict)
            @info """$k not found in user specified input file. 
            Reading $k from default TASOPT input.\n$k = $(default_dict[k])"""
            default_dict[k]
        else
            error("Requested key/parameter is not supported. Check the default 
            input file to see all available input options.")
        end
    end
end

data = TOML.parsefile("input.toml")

nmisx = data["N_missions"]
pari = zeros(Int64, iitotal)

# Setup option variables
options = data["Options"]
if options["fueltype"] == "LH2"
    pari[iifuel] = 1
else
    pari[iifuel] = 2
end
pari[iifwing]  = options["fuelinwing"]
pari[iifwcen]  = options["fuelinwingcen"]
pari[iiHTsize] = options["HTsize"]
pari[iiVTsize] = options["VTsize"]
pari[iifclose] = options["closefuse"]
pari[iiengloc] = options["engloc"]
pari[iiBLIc]   = options["coreBLI"]
pari[iixwmove] = options["movewingbox"]
pari[iiwplan]
pari[iiengwgt]