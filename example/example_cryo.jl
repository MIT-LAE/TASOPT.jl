# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
include(TASOPT.__TASOPTindices__)
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
#  load default model

example_ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/example_fuel_cell.toml")) # MODIFY <path> appropriately

# 3) Size aircraft
time_size_aircraft = @elapsed size_aircraft!(example_ac, iter=50)
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