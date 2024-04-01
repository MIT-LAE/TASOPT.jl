# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
include("../src/misc/index.inc")
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
#  load default model
# example_ac = load_default_model() # simply a synonym to read_aircraft_model()
# Alternatively you can load your desired input file 
ac = read_aircraft_model("../src/IO/experiment_input.toml") # MODIFY <path> appropriately

# 3) Size aircraft
# time_wsize = @elapsed size_aircraft!(ac)
#println("Time to size aircraft = $time_wsize s")

# 3.5) Modify the aircraft model
display(size(ac.pari))
display(size(ac.parg))
display(size(ac.parm))
display(size(ac.para))
display(size(ac.pare))


# 4) Visualize outputs
# Output resulting geometry of aircraft
# summary(example_ac)
# Or individually look at certain aspects:
# Show weight buildup of the aircraft:
# TASOPT.weight_buildup(example_ac) 
# # Show aerodynamics:
# TASOPT.aero(example_ac)
# # Geometry:
# TASOPT.geometry(example_ac)

# 5) Plot figures
# using PyPlot
# TASOPT.stickfig(example_ac)
# plt.savefig("Example.png")