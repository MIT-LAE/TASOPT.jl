# Example for a Payload-Range diagram

![PayloadRangePlot](../assets/PayloadRangeExample.png)

## Choosing a design mission

To plot a payload-range diagram with a fleet of missions you must first load any aircraft model. Start by choosing a design mission. Your design mission should be what you want the second corner point in your Payload-Range plot to be. Once you have a chosen a specific design range and payload weight (For eg: 3500 nmi and 195 pax) you can add it to the input toml file for eg: `default_input.toml`

```toml
[Mission]
    N_missions = 2 # Number of missions to be modeled (first mission is the design mission)
    pax = 195       # Number of passengers in each mission
    range = "3500.0 nmi" # Design Range + second mission range
                            #["3000.0 nmi", "500.0 nmi", "2500.0 nmi", "3550.0 nmi", "3734.0 nmi"] # Design Range + second mission range
    weight_per_pax = "215.0 lbf" # Specify weight per passenger - 
                            # includes luggage [lbm or lbf or kg or N] 
```

## Julia script for Payload Range Diagram

Start the script importing `TASOPT.jl` and then loading the default `aircraft` model.

```julia
# Import modules
using PythonPlot
using TASOPT
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand
include(__TASOPTindices__)
# import indices for calling parameters

# Load aircraft using default module
ac = load_default_model() #Use default model for payload-range diagram
size_aircraft!(ac)
```

One way is to call the `PayloadRange` function:

```julia
TASOPT.PayloadRange(ac)
```
In this approach, only one mission needs to be specified. TASOPT will copy the parameters from the sizing mission (e.g., takeoff altitude and temperature), and vary the payload and range to produce a payload-range diagram. 

If you want a more customizable diagram, you may specify a second mission and use an approach such as the one below. First initialize some variables for mission range and payloads

```julia
# Make an array of ranges to plot
RangeArray = ac.parm[imRange] * LinRange(0.1,1.2,Rpts)
# Store constant values to compare later
Wmax = ac.parg[igWMTO]
Fuelmax = ac.parg[igWfmax]
Wempty = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay]
# Arrays for plotting
RangesToPlot = []
PayloadToPlot = []
maxPay = ac.parm[imWpay ]
```

## Main iteration loop

```julia
for Range = RangeArray
    # Make an array of payloads to plot
    Payloads = (maxPay) * LinRange(1, 0.1, Ppts)
    ac.parm[imRange] = Range
    for mWpay = Payloads
        println("Checking for Range (nmi): ",Range/1852.0, " and Pax = ", mWpay/(215*4.44822))
        ac.parm[imWpay ] = mWpay
        # Try fly_off_design after setting new range and payload
        try
            TASOPT.fly_off_design(ac, 2, saveOffDesign = true)
            # fly_off_design success: store maxPay, break loop
            mWfuel = ac.parm[imWfuel,2]
            WTO = Wempty + mWpay + mWfuel

            if WTO > Wmax || mWfuel > Fuelmax || WTO < 0.0 || mWfuel < 0.0 
                WTO = 0.0
                mWfuel = 0.0
                println("Max out error!")
                if mWpay == 0
                    println("Payload 0 and no convergence found")
                    maxPay = 0
                end
            else
                maxPay = mWpay
                println("Converged - moving to next range...")
                break
            end     
        catch
            println("Not Converged - moving to lower payload...")      
        end
    end
    append!(RangesToPlot, Range)
    if OEW
        append!(PayloadToPlot, maxPay+Wempty)
    else
        append!(PayloadToPlot, maxPay)
    end
end
```

## Plot Payload Range diagram

```julia
using PythonPlot
fig, ax = pyplot.subplots(figsize=(8,5), dpi = 300)
ax.plot(RangesToPlot ./ (1000*1852.0), PayloadToPlot./ (9.8*1000), linestyle="-",  color="b", label="Payload ")
ax.set_xlabel("Range (1000 nmi)")
ax.set_ylabel("Weight (1000 kg)")
ax.legend()
ax.set_title("Payload Range Plot")
ax.grid()

fig.savefig("./PayloadRangeExample.png")
```
