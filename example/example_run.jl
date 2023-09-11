# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
#  load default model
example_ac = load_default_model() # simply a synonym to read_aircraft_model()
# Alternatively you can load your desired input file 
# example_ac = read_aircraft_model("../src/IO/input.toml") # MODIFY <path> appropriately

# 3) Size aircraft
saveOD = false
time_wsize = @elapsed size_aircraft(example_ac, 35, 0, false, true, saveOD)
println("Time to size aircraft = $time_wsize s")

# 4) Visualize outputs
# Output resulting geometry of aircraft
summary(example_ac)
# Or individually look at certain aspects:
# Show weight buildup of the aircraft:
# TASOPT.weight_buildup(example_ac)
# # Show aerodynamics:
# TASOPT.aero(example_ac)
# # Geometry:
# TASOPT.geometry(example_ac)

# 5) Plot figures
using PyPlot
TASOPT.stickfig(example_ac)
plt.savefig("Example.png")