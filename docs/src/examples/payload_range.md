# Example for a Payload-Range diagram

![PayloadRangePlot](../assets/PayloadRangeExample.png)

## Choosing a design mission

To plot a payload-range diagram with a fleet of missions you must first load a model with >=1 non design mission

Start by choosing a design mission. Your design mission should be what you want the second corner point in your Payload Range plot to be. Once you have a chosen a specific design range and payload weight (For eg: 3500 nmi and 195 pax) you can add it to the input toml file: `default_input.toml`

```toml
[Mission]
    N_missions = 5 # Number of missions to be modeled (first mission is the design mission)
    pax = 195       # Number of passengers in each mission
    range = "3500.0 nmi" # Design Range + second mission range
                            #["3000.0 nmi", "500.0 nmi", "2500.0 nmi", "3550.0 nmi", "3734.0 nmi"] # Design Range + second mission range
    weight_per_pax = "215.0 lbf" # Specify weight per passenger - 
                            # includes luggage [lbm or lbf or kg or N] 
```

## Julia script for Payload Range Diagram
Start the script importing `TASOPT.jl`, `PyPlot` and `index.inc` and then loading the default `aircraft` model.
```julia
# Import modules
using PyPlot
using TASOPT
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand
include("../src/misc/index.inc")
# import indices for calling parameters

# Load aircraft using default module
ac = TASOPT.read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/PRD_input.toml"))
time_wsize = @elapsed size_aircraft!(ac)
```

Initialize some variables for mission range and payloads

```julia
nmis = size(ac.parm)[2] # Get number of missions
MAX_PAYLOAD = 230 * 215 * 4.448222 # Starting Payload, #Passengers * 215 lb * lb_to_N
# Ranges:
Range_arr_nmi = LinRange(100, 1000, nmis-1) # Starting range array (in nmi)
Range_arr = Range_arr_nmi .* 1852.0 # Convert to SI units
Max_range = Range_arr[nmis-1]
curr_range = Range_arr[1]
prev_range = 0.1
```

## Initialize variables for mission iterations

```julia
NPSS = Base.Process
NPSS_PT = true
itermax = 15
initeng = 0
Ldebug = false
saveOD = false
SEC_B = false
x_range = []
y_payload = []
y_Wemtpy = []
```

## Main iteration loop

```julia
while MAX_PAYLOAD >= 0
    println("-----------------------")
    println("Starting new set of missions, \n Ranges: ", Range_arr, "\n Payload: ", MAX_PAYLOAD./(215*4.448222))
    try
        # Set non design mission ranges and payloads
        ac.parm[imRange, 2:5] .= Range_arr
        ac.parm[imWpay, 2:5] .= MAX_PAYLOAD
        # Analyze each mission
        for mi in 2:size(ac.parm)[2]
            prev_range = curr_range
            curr_range = ac.parm[imRange, mi]
            # Call TASOPT woper function
            @views TASOPT.woper(ac.pari,ac.parg,ac.parm[:,mi:mi],ac.para[:,:,mi:mi],ac.pare[:,:,mi:mi], ac.para[:,:,1:1],ac.pare[:,:,1:1], itermax,initeng, NPSS_PT, NPSS)
            # Check if mission fuel and MTO weight is greater than design fuel and MTO weight
            if (ac.parm[imWfuel,mi] - ac.parg[igWfuel] > 1) || (ac.parm[imWTO,mi] - ac.parg[igWMTO] > 1 )
                printstyled("Mission designed beyond capacity!", "\n"; color=:red)
                println([ac.parg[igWfuel], ac.parm[imWfuel,mi], ac.parg[igWMTO], ac.parm[imWTO,mi]])
                throw(UndefVarError([ac.parg[igWfuel], ac.parm[imWfuel,mi], ac.parg[igWMTO], ac.parm[imWTO,mi]]))
            end

            # Add to dataframes
            append!(x_range, curr_range)
            append!(y_payload, MAX_PAYLOAD)
            append!(y_WTO,  ac.parm[imWTO,mi])

            # If in section B of Payload Range plot, decrease payload by 1 passenger
            if SEC_B
                println("decreasing payload since Section B")
                MAX_PAYLOAD = max(0, MAX_PAYLOAD-(1*215*4.448222))
            end
        end
        println("PRD converged for all: Increasing Ranges")
        prev_range = curr_range
        Max_range = prev_range+ (500* 1852.0) #Increase range array by 500 nmi
        catch
            println("PRD failed: Decreasing Payload and range")
            prev_range = curr_range
            Max_range = prev_range+ ((Max_range-prev_range)*0.75)
            MAX_PAYLOAD = MAX_PAYLOAD-pax2N(5)
            SEC_B = true
        end
    Range_arr = LinRange(prev_range, Max_range, nmis-1)
    println("-----------------------")
end
```

## Plot Payload Range diagram

```julia
using PyPlot
figure()
y_OEW = y_Wemtpy .+ y_payload
plot(x_range ./ (1852.0*1000), y_OEW./ (4.448222* 1000), linestyle="-",  color="b", label="OEW + Payload ")

xlabel("Range (1000 nmi)")
ylabel("Weight (1000 lbs)")
title("Payload Range Plot")

legend()
grid()
savefig("./PayloadRangeExample.png")
```
