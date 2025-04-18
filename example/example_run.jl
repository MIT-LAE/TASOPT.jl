# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOpt. 

# 1) Load TASOPT
using TASOpt
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
#  load default model
example_ac = load_default_model() # simply a synonym to read_aircraft_model()
# Alternatively you can load your desired input file 
# example_ac = read_aircraft_model("../src/IO/input.toml") # MODIFY <path> appropriately

# 3) Size aircraft
time_size_aircraft = @elapsed size_aircraft!(example_ac)
println("Time to size aircraft = $time_size_aircraft s")

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
using Plots
p = TASOPT.stickfig(example_ac)
savefig(p, "Example.png")