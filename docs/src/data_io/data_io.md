
---
# Data structure basics

## `aircraft` struct

## `par` arrays

Refer to the [struct reference page](@ref datastructs) for add'l details.

---
# Inputting to models

## `.toml`s and the default input
tweak default input or import then tweak for example, to _____

```julia
igTmetal = 239 #index for max metal temp (see src/misc/index.inc)
ac = load_default_model() #creates new aircraft using default input .toml
ac.parg[igTmetal] = 2000 #set max metal temp to 2000 K
```

## Reading externals

temp


---

# Exporting from models

## Quicksave and quickload via .toml

## Saving model to a readable .toml

## Output data to .csv
