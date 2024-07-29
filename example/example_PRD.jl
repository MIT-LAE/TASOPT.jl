"""
An example script to plot a payload-range diagram with a 
fleet of missions
"""

# 1. Import modules
using PythonPlot
using TASOPT
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand
include(__TASOPTindices__)
# import indices for calling parameters

function PayloadRange(ac; Rpts = 20, Ppts = 20, filename = "PayloadRangeDiagram.png", OEW = false)
    RangeArray = ac.parm[imRange,1] * LinRange(0.1,2,Rpts)
    maxPay = 0

    Wmax = ac.parg[igWMTO]
    Fuelmax = ac.parg[igWfmax]
    Wempty = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay]

    RangesToPlot = []
    PayloadToPlot = []
    maxPay = ac.parm[imWpay]

    for Range = RangeArray
        if maxPay == 0
            break
        else
            Payloads = (maxPay) * LinRange(1, 0, Ppts)
        end
        ac.parm[imRange,2] = Range
        for mWpay = Payloads
            println("Checking for Range (nmi): ",Range/1852.0, " and Pax = ", mWpay/(215*4.44822))
            ac.parm[imWpay,2] = mWpay
            try
                TASOPT.woper(ac, 2, saveOffDesign = true)
                # woper success: store maxPay, break loop
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