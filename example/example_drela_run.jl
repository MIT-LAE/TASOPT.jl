
using Pkg
Pkg.activate("../")
Pkg.instantiate()

# 1) Load TASOPT
include("../src/TASOPT.jl")

# 2) Include input file for desired aircraft/
#  can also be a saved model
include("737input_2.jl")

# 3) Size aircraft
time_wsize = @elapsed size_aircraft(35, 0, false, true, saveOD)
println("Time to size aircraft = $time_wsize s")
# 4) Visualize outputs
# Output resulting geometry of aircraft
geometry(parg)

# Show weight buildup of the aircraft
weight_buildup(parg)

# 5) Plot figures
pygui(true)
stickfig(parg, pari, parm)
# plt.savefig("../../example/Example.png") # BUG