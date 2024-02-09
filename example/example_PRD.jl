"""
An example script to plot a payload-range diagram with a 
fleet of missions
"""

# 1. Import modules
using PyPlot
using TASOPT
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand
include("../src/misc/index.inc")
# import indices for calling parameters

function PayloadRange(ac_og, Rpts = 20, Ppts = 20, filename = "./example/PayloadRangeExample.png", OEW = false, itermax = 20.0, initeng = true)
    # Copy aircraft structure as we will be changing range and payloads
    ac = deepcopy(ac_og)
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

    for Range = RangeArray
        # Make an array of payloads to plot
        Payloads = (maxPay) * LinRange(1, 0.1, Ppts)
        ac.parm[imRange] = Range
        for mWpay = Payloads
            println("Checking for Range (nmi): ",Range/1852.0, " and Pax = ", mWpay/(215*4.44822))
            ac.parm[imWpay ] = mWpay
            # Try woper after setting new range and payload
            try
                @views TASOPT.woper(ac, itermax, initeng)
                # woper success: store maxPay, break loop
                WTO = Wempty + mWpay + ac.parm[imWfuel]
                mWfuel = ac.parm[imWfuel]

                # Compare with previously stored constants
                if WTO > Wmax || mWfuel > Fuelmax || WTO < 0.0 || mWfuel < 0.0 
                    WTO = 0.0
                    mWfuel = 0.0
                    println("Max out error!")
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
    println(RangesToPlot)
    println(PayloadToPlot)
    fig, ax = plt.subplots(figsize=(8,5), dpi = 300)
    ax.plot(RangesToPlot ./ (1000*1852.0), PayloadToPlot./ (9.8*1000), linestyle="-",  color="b", label="Payload ")
    ax.set_xlabel("Range (1000 nmi)")
    ax.set_ylabel("Weight (1000 kg)")
    ax.legend()
    ax.set_title("Payload Range Plot")
    ax.grid()

    fig.savefig(filename)
end

# Load default model
ac = TASOPT.read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/PRD_input.toml"))
size_aircraft!(ac)