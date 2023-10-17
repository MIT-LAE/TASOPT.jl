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

function km2nmi(km)
    return km./1852.0
end

function nmi2km(nmi)
    return nmi.*1852.0
end

function pax2N(pax)
    pax.*215*4.448222
end

function N2pax(N)
    N./(215*4.448222)
end

# Alogrithm:

# Input: Range (arbitrary), WMpay
# Size the aircraft
# For a certain range of Ranges
    # Try woper for each mission other than design
    # If converges then append to data struct
        # Increase range of ranges
        # Reduce Payload if in section B
    # Else if error/ non converge
        # Reduce payload
        # Activate section B
    # End loop if difference between min range and max range below threshold
    # OR payload <= 0


lbf_to_N = 4.448222
MAX_PAYLOAD = pax2N(230)

R_thresh = Inf

Range_arr_nmi = LinRange(100, 1000, 4)
Range_arr = Range_arr_nmi .* 1852.0
Max_range = Range_arr[4]

curr_range = Range_arr[1]
prev_range = 0.1

# 5. Initialze some variables
NPSS = Base.Process
NPSS_PT = true
itermax = 15
initeng = 0
Ldebug = false
saveOD = false

x_range = []
y_payload = []
y_WTO = []
y_Wfuel = []
y_Wemtpy = []


range_thresh = nmi2km(5) # 50 nmi

SEC_B = false

# Load default model
ac = TASOPT.read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/PRD_input.toml"))

time_wsize = @elapsed size_aircraft!(ac)


while MAX_PAYLOAD >= 0
    println("-----------------------")
    println("Starting new set of missions, \n Ranges: ", km2nmi(Range_arr), "\n Payload: ", N2pax(MAX_PAYLOAD))
    
    try
        ac.parm[imRange, 2:5] .= Range_arr
        ac.parm[imWpay, 2:5] .= MAX_PAYLOAD
        
        for mi in 2:size(ac.parm)[2]
            prev_range = curr_range
            curr_range = ac.parm[imRange, mi]
            println("Analysing Mission Number: ",mi)
            println("Range: ",curr_range)
            @views TASOPT.woper(ac.pari,ac.parg,ac.parm[:,mi:mi],ac.para[:,:,mi:mi],ac.pare[:,:,mi:mi], ac.para[:,:,1:1],ac.pare[:,:,1:1], itermax,initeng, NPSS_PT, NPSS)
            if (ac.parm[imWfuel,mi] - ac.parg[igWfuel] > 1) || (ac.parm[imWTO,mi] - ac.parg[igWMTO] > 1 )
                printstyled("Mission designed beyond capacity!", "\n"; color=:red)
                println([ac.parg[igWfuel], ac.parm[imWfuel,mi], ac.parg[igWMTO], ac.parm[imWTO,mi]])
                throw(UndefVarError([ac.parg[igWfuel], ac.parm[imWfuel,mi], ac.parg[igWMTO], ac.parm[imWTO,mi]]))
            end
            #---- max TO weight
            WMTO = ac.parg[igWMTO]
            # ---- zero-fuel weight for this mission
            Wzero = WMTO-ac.parg[igWfuel]-ac.parg[igWpay]
            
            append!(y_Wemtpy, Wzero)
            append!(x_range, curr_range)
            append!(y_payload, MAX_PAYLOAD)
            append!(y_Wfuel,  ac.parm[imWfuel,mi])
            append!(y_WTO,  ac.parm[imWTO,mi])
            if SEC_B
                println("decreasing payload since Section B")
                MAX_PAYLOAD = max(0, MAX_PAYLOAD-pax2N(1))
            end
        end
        println("PRD converged for all: Increasing Ranges")
        prev_range = curr_range
        Max_range = prev_range+(500* 1852.0)
    catch
        println("PRD failed: Decreasing Payload and range")
        prev_range = curr_range
        Max_range = prev_range+ ((Max_range-prev_range)*0.75)
        MAX_PAYLOAD = MAX_PAYLOAD-pax2N(5)
        SEC_B = true
    end
    Range_arr = LinRange(prev_range, Max_range, 4)
    println("-----------------------", abs(Max_range-prev_range))
end


using PyPlot
figure()
y_OEW = y_Wemtpy .+ y_payload
plot(x_range ./ (1000*1852.0), y_OEW./ (4.448222* 1000), linestyle="-",  color="b", label="OEW + Payload ")
# plot(x_range./ (1000), (y_WTO.-y_Wfuel)./ (4.448222* 1000), marker="o", linestyle="-", color="g", label="MTO-wf")
# plot(x_range./ (1000), (y_WTO)./ (4.448222* 1000), marker="o", linestyle="-", color="g", label="MTO")
# plot(x_range./ (1000), y_Wfuel./ (4.448222* 1000), marker="o", linestyle="-", color="k", label="Wfuel")
# plot(x_range./ (1000), y_payload./ (4.448222* 1000), linestyle="-", color="b", label="Payload weight")


xlabel("Range (1000 nmi)")
ylabel("Weight (1000 lbs)")
title("Payload Range Plot")
legend()
grid()

savefig("./example/PayloadRangeExample_last.png")
