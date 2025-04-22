# Loading and saving models

## TOML input files

TASOPT models can be loaded from and saved to human-readable TOML files via [`read_aircraft_model()`](@ref) and [`save_aircraft_model()`](@ref). These scripts are written to process files based on the `/example/defaults/default_input.toml`, as is done when `load_default_model()` is run; they are usable for simple cases that alter the default inputs. For larger deviations from the default (e.g., mission-specific aero or engine performance), customized scripts are recommended when not required.

Once it is loaded in this fashion, the `aircraft` can be manipulated and re-sized. Note: TOML-based saves **do not** store the performance of the sized model, so a model must be sized before data can be referenced usefully:

```julia
using TASOPT

include(__TASOPTindices__)  #import array indices from ./src/data_structs/index.inc, including igTmetal

filepath = "/path/to/your/aircraft.toml"
ac = read_aircraft_model(filepath)   #creates new aircraft using default input .toml
size_aircraft!(ac)

ac.parg[igTmetal] = 2000    #set max metal temp to 2000 K
size_aircraft!(ac)          #resizing after parameter change
```

## Quicksave and quickload to `.jld2`

Complete `aircraft` structs can also be serialized to JLD2 files by calling [`quicksave_aircraft()`](@ref). Though inscrutable to the human eye, these saves capture **ALL** elements of the `aircraft` and can be loaded via [`quickload_aircraft()`](@ref) to a fresh REPL in an already-sized state.

```julia
#load and size the default aircraft model
using TASOPT
ac2 = load_default_model()
size_aircraft!(ac2)

#quicksave aircraft
filepath2 = "/path/to/your/new/quicksave.toml"
quicksave_aircraft(ac2, filepath2)

#quickload aircraft, ready to pull results
ac2 = quickload_aircraft(filepath2)
```


---

## Output data to CSV

`aircraft` data can be saved to CSVs via [`output_csv()`](@ref), where each column will be a `par` array quantity and each row can be a different `aircraft` or case. This is especially useful for post-processing of parameter sweeps.

The data that is output can be customized by specifying individual `par` array quantities, flight points, or missions. By default, only the first cruise point of the design mission is output for a representative set of values (i.e., those in [`default_output_indices`](@ref)). 

For ease of use, the default behavior is to create or append to the file at the user-given `filepath`, though if the requested indices are changed, a new file will be created to avoid inconsistency. This file will also appendable without further input. Overwrite behavior can also be specified.

The following example shows the basic functionalities for a parameter sweep. See the [`output_csv()`](@ref) docstring for details.


```julia 
# Sweeping a parameter space and outputting each design
using TASOPT
include(__TASOPTindices__) #import par array indices, including igTmetal
ac = load_default_model()
filepath = "path/to/your/newfile.csv"

#set up sweep of max metal temp in engine
Tmetals = [1800, 1850, 1900, 1950] # Kelvin
for Tmetal in Tmetals
    
    ac.parg[igTmetal] = Tmetal
    size_aircraft!(ac)

    #output default quantities, appending each time
    output_csv(ac, filepath) 
    #output with more engine quantities and with data at all flight points
    #due to inconsistent headings at filepath, a new file will be created and appended to
    output_csv(ac, filepath, indices = output_indices_wEngine, includeMissions = true)  
end
```

---

## IO Functions

```@docs
read_aircraft_model()

save_aircraft_model()

quicksave_aircraft()

quickload_aircraft()

output_csv()

default_output_indices
```
