#!!! compat
#DELETE ME. I NO LONGER WORK

# 1) Load TASOPT
using TASOPT
# 2) Include input file for desired aircraft/ can also be a saved model
# Here this file defines a B737 model
include("737input_2.jl")

# 3) Size aircraft
time_wsize = @elapsed size_aircraft!(B737)
println("Time to size aircraft = $time_wsize s")

# 4) Visualize outputs
# Output resulting geometry of aircraft
# TASOPT.geometry(B737.parg) #TODO: needs update to use

# Show weight buildup of the aircraft
# TASOPT.weight_buildup(B737.parg) #TODO: needs update to use

# 5) Plot figures
# pygui(true) #TODO: needs update to use (this and below)
# stickfig(parg, pari, parm)
# plt.savefig("../../example/Example.png") # BUG