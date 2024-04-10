# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
using PyPlot
using CSV
# using Tables
include("../src/misc/index.inc")
ft_to_m = 0.3048
# you can optionally define
# const tas = TASOPT 
# to use as a shorthand

# 2) Include input file for desired aircraft/
#  load default model
# example_ac = load_default_model() # simply a synonym to read_aircraft_model()
# Alternatively you can load your desired input file 
nameAircraftModel = "../src/IO/experiment_input.toml"
ac = read_aircraft_model(nameAircraftModel) # MODIFY <path> appropriately
saveName = "EthaJetA31PerBlend1500nmi"
# 2.5) Change fuel type
ac.pari[iifuel] = 322431 #(JetA:24 Ethanol:32 JetAEtha31%Blend: 322431)
ac.parg[igrhofuel] = 808.1 #(JetA:817.0 Ethanol:789.0 JetAEtha31%Blend: 808.1)

# 3) Find Optimal Flight Altitude
AltList = LinRange(2e4,5e4,10) #ft
AltRec = [] #ft
RanRec = [] #nmi
WMTORec = [] #Ton (metric)
WFuelRec = [] #Ton (metric)
WPayRec = [] #Ton (metric)
PFEIRec = [] #(J/J)
WTO_WTOmaxRec = [] #N/N
Wf_WfmaxRec = [] #N/N
areaWingRec = [] #m2
ARWingRec = [] #Aspect ratio
spanWingRec = [] #m
diaFanRec = [] #m
FnTotCRRec = [] #N
for Alt = AltList
    ac.para[iaalt, ipcruise1, 1] =  Alt * ft_to_m # Cruise Altitude
    try
        size_aircraft!(ac,iter=500)
        append!(AltRec,Alt)
        append!(RanRec,ac.parg[igRange]./1852.0)
        append!(WMTORec,ac.parg[igWMTO]./(9.81*1000))
        append!(WFuelRec,ac.parg[igWfuel]./(9.81*1000))
        append!(WPayRec,ac.parg[igWpay]./(9.81*1000))
        append!(PFEIRec,ac.parm[imPFEI,1])
        append!(WTO_WTOmaxRec,ac.parm[imWTO,1]/ac.parg[igWMTO])
        append!(Wf_WfmaxRec,ac.parg[igWfuel]/ac.parg[igWfmax])
        append!(areaWingRec,ac.parg[igS])
        append!(ARWingRec,ac.parg[igAR])
        append!(spanWingRec,ac.parg[igb])
        append!(diaFanRec,ac.parg[igdfan])
        append!(FnTotCRRec,ac.pare[ieFe,ipcruise1,1])
    catch
        println("Failed at :",Alt)
    end
end
#Post Process
WEmpRec = WMTORec - WFuelRec - WPayRec # Ton Metric
EFuelRec  = PFEIRec.*WPayRec*1000*9.81.*RanRec*1852.0 #Joul
#Collect output Data
outputTup = (AltRec=AltRec,RanRec=RanRec,WMTORec=WMTORec,WFuelRec=WFuelRec
             ,WPayRec=WPayRec,PFEIRec=PFEIRec,WTO_WTOmaxRec=WTO_WTOmaxRec
             ,Wf_WfmaxRec=Wf_WfmaxRec,areaWingRec=areaWingRec,ARWingRec=ARWingRec
             ,spanWingRec=spanWingRec,diaFanRec=diaFanRec,FnTotCRRec=FnTotCRRec,WEmpRec=WEmpRec,EFuelRec=EFuelRec)
CSV.write(saveName*".csv",  outputTup, writeheader=true)
#Plot out the operating conditions for the optimal point
maskFeasi = (WTO_WTOmaxRec.<=1) .& (Wf_WfmaxRec.<=1) #These are the feasible solution
idxPFEIBest = findmin(PFEIRec[maskFeasi])[2] #Find the location of the best PFEI in the feasible PFEI Range
AltFeasi = AltRec[maskFeasi]
AltBest = AltFeasi[idxPFEIBest]
println("Best Altitude at [ft]: ",AltBest)
#Resize at the best altitude for detailed data extraction
ac2 = read_aircraft_model(nameAircraftModel)
ac2.pari[iifuel] = ac.pari[iifuel]
ac2.parg[igrhofuel] = ac.parg[igrhofuel]
ac2.para[iaalt, ipcruise1, 1] =  AltBest * ft_to_m # Cruise Altitude
size_aircraft!(ac2,iter=500)
PFEIBest = ac2.parm[imPFEI,1]
println("PFEI Numerical Error Check: ",100*(minimum(PFEIRec[maskFeasi])-PFEIBest)/PFEIBest," %")
#Read out the data at the optimal conditions
##Create a mask to mask out unreported phases
maskRep = ac2.pare[ieFe,:,1].>0 #Reported Phase has non zero thrust
phases = ["ST","RO","TO","CB","B1","B2","B3","B4","B5","C1","C2","D1","D2","D3","D4","D5","Test"]
phases = phases[maskRep]
print("Reported Phases Are:", phases,"\n")
##Read out other parameters
###Aero Parameters
timeOptMiss = ac2.para[iatime,maskRep,1] #second
ranOptMiss = ac2.para[iaRange,maskRep,1] #meter
altOptMiss = ac2.para[iaalt,maskRep,1] #meter
machOptMiss = ac2.para[iaMach,maskRep,1]
weiOptMiss = ac2.para[iafracW,maskRep,1]*(ac.parm[imWTO,1]/9.81) #kg
gamOptMiss = ac2.para[iagamV,maskRep,1] #rad
LDROptMiss = ac2.para[iaCL,maskRep,1]./ac2.para[iaCD,maskRep,1] #Lift to drag ratio
###Engine Parameters
hfOptMiss = ac2.pare[iehfuel,maskRep,1] #J/kg equivalent heating value
TfuelOptMiss = ac2.pare[ieTfuel,maskRep,1] #K fuel temperature
Tt3OptMiss = ac2.pare[ieTt3,maskRep,1] #K
Pt3OptMiss = ac2.pare[iept3,maskRep,1] #Pa 
Tt4OptMiss = ac2.pare[ieTt4,maskRep,1] #K
Pt4OptMiss = ac2.pare[iept4,maskRep,1] #Pa
FnOptMiss = ac2.pare[ieFe,maskRep,1] #N Total Thrust for all engines
mdotfOptMiss = FnOptMiss.*ac2.pare[ieTSFC,maskRep,1]/9.81 #kg/s for all engines
#Output Additional Data at the optimal mission
outputTup = (Phase=phases,Time=timeOptMiss,Range=ranOptMiss,Altitude=altOptMiss,MachNumber=machOptMiss,Weight=weiOptMiss
            ,ClimbAngle=gamOptMiss,LiftDragRatio=LDROptMiss,HeatingValue=hfOptMiss,FuelTemp=TfuelOptMiss
            ,Tt3=Tt3OptMiss,Pt3=Pt3OptMiss,Tt4=Tt4OptMiss,Pt4=Pt4OptMiss,Thrust=FnOptMiss,mdotFuel=mdotfOptMiss)
CSV.write(saveName*"OptimalMission.csv",  outputTup, writeheader=true)
#Backup Code Below
# time_wsize = @elapsed size_aircraft!(ac,iter=500)
#println("Time to size aircraft = $time_wsize s")

# 3.5) Read out the size of each variable
# display(size(ac.pari))
# display(size(ac.parg))
# display(size(ac.parm))
# display(size(ac.para))
# display(size(ac.pare))

# 3.75) Read out the total weight and flight range
# println("flight range (nmi): " , ac.parg[igRange]./1852.0)
# println("WMTO (1000 kg):" , ac.parg[igWMTO]./(9.81*1000))
# println("PFEI:",ac.parm[imPFEI])
# println("OPR:",ac.pare[iept3]/ac.pare[iept2])
# println("LPCPR:",maximum(ac.pare[iepilc, :, 1]))
# println("WTO/WTOmax:",ac.parm[imWTO,1]/ac.parg[igWMTO])
# println("Wf/Wfmax:",ac.parg[igWfuel]/ac.parg[igWfmax])

#Indexing
#A[(A .>2) .& .!(A.>3)]

# 4) Visualize outputs
# Output resulting geometry of aircraft
# summary(example_ac)
# Or individually look at certain aspects:
# Show weight buildup of the aircraft:
# TASOPT.weight_buildup(example_ac) 
# # Show aerodynamics:
# TASOPT.aero(example_ac)
# # Geometry:
# TASOPT.geometry(example_ac)

# 5) Plot figures
# using PyPlot
# TASOPT.stickfig(example_ac)
# plt.savefig("Example.png")
# println(isless.(Wf_WfmaxRec,1))

# fig, ax = plt.subplots(figsize=(8,5), dpi = 300)
# ax2 = ax.twiny()
# ax.plot(PFEIRec , AltRec, linestyle="-",  color="b", label="PFEI")
# ax2.plot(Wf_WfmaxRec , AltRec , linestyle="--",  color="r", label="Wf/Wfmax")
# ax.set_xlabel("Passenger Fuel Emission Index [J/J]")
# ax2.set_xlabel("Fuel Weight over Maximum Fuel Tank Capacity [N/N]")
# ax.set_ylabel("Cruise Altitude [ft]")
# ax.legend()
# ax2.legend()
# ax.set_title("PFEI")
# ax.grid()
# fig.savefig("PFEI.png")