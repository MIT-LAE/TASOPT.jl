using TOML
include("index.inc")

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