"""
An example script to plot a payload-range diagram with a 
fleet of missions
"""

# 1. Import modules
# using PythonPlot
using TASOPT
using Plots
# import indices for calling parameters
include(__TASOPTindices__)

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
            catch err
                println("Not Converged - moving to lower payload...")      
                rethrow(err)
            end
        end
        append!(RangesToPlot, Range)
        if OEW
            append!(PayloadToPlot, maxPay+Wempty)
        else
            append!(PayloadToPlot, maxPay)
        end
    end

    # Convert values for plotting
    dpi = 300
    ranges_kft = RangesToPlot ./ (1000 * 1852.0)
    payload_tons = PayloadToPlot ./ (9.8 * 1000)


    # Plot with all attributes set in plot()
    plot1 = plot(ranges_kft, payload_tons, 
        lw=2,                   # Line width
        line=:solid,            # Line style
        color=:blue,            # Line color
        label="Payload",        # Legend label
        xlabel="Range (1000 nmi)", 
        ylabel="Weight (1000 kg)", 
        title="Payload-Range Diagram: "*string(ac.name), 
        grid=true,              # Enable grid
        dpi = 300)
        # size=(8*dpi, 5*dpi))  # Size in pixels (width, height)

    # Save with specified DPI
    savefig(filename)

    display(plot1)
    println(RangesToPlot)
    println(PayloadToPlot)

end

# Load default model
ac = TASOPT.read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/PRD_input.toml"))
        # defaultfile = joinpath(TASOPT.__TASOPTroot__, "IO/default_input.toml"))
size_aircraft!(ac)
PayloadRange(ac)