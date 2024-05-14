# This is an example file to load an aircraft model/ input file and 
# size an aircraft using TASOPT. 

# 1) Load TASOPT
using TASOPT
using PyPlot
using CSV
# using Tables
include("../src/misc/index.inc")
ft_to_m = 0.3048
# From Optimization
"""
               AR                  Alt(ft)             Cl                  Λ                   λs                  λt                   hboxo                hboxs                rcls                rclt                Tt4CR               iepihc              iepif
Ethanol 1500: [11.999856667089247, 34443.14280029319 , 0.5754043145206887, 27.08567980701492,  0.9380520909882628, 0.10017099443433408, 0.14950245863108022, 0.12845986215006178, 1.0781134610041747, 0.9975636667912358, 1707.5445791618208, 15.0,               0.5097104925337879] 0.7669812558976667
JetA    1500: [11.967887176329073, 32479.866286551907, 0.5956532812392801, 27.382147149331143, 0.9447567896238773, 0.10644959191186304, 0.14694888676914802, 0.12003039648977087, 1.0596807683444869, 0.9905991041425743, 1749.8757464045927, 14.988409779779435, 2.6302018540124124] 0.7129078730647352
JetAEtha1500: [11.991011931628014, 32315.845997648077, 0.5748237916456843, 27.972858575055668, 0.9366717483405769, 0.12589668790310904, 0.14905975904703933, 0.1295658503413098,  1.0000313446876867, 0.977401864414243,  1749.621251634486,  14.987524381849813, 3.8295309599571854] 0.73997964807403
Ethanol 3000: [10.127599907033918, 32997.08887497184,  0.5147234038045866, 27.689288150041946, 0.7738191840909827, 0.1001515244423182,  0.14913014728970564, 0.14391182104054678, 1.3510738250216734, 1.0,                1642.7577794418612, 12.233750138515445, 1.6619068159205148] 1.140152289017824
JetA    3000: [11.380434406953572, 33198.13022096391,  0.5410047569882184, 26.66666248716794,  0.9378606924235586, 0.10281538084027678, 0.14977259425119216, 0.13296003568212714, 1.0627223742965342, 0.996592019369539,  1737.0012718430448, 14.997879755556585, 0.30657837795564297] 0.80283263924495
JetAEtha3000: [11.613463065609654, 34769.903770435114, 0.5545477127074685, 27.140602449287496, 0.8854377722367844, 0.10242836542719086, 0.14004048320503837, 0.1406054334210777, 1.1578677693088162, 0.9911248968752646, 1726.0197142159793, 14.944693491120033, 0.015218023095181485] 0.8942542470234152
"""

# 2) Include input file for desired aircraft/
nameAircraftModel = "../src/IO/experiment_input_1500.toml"
ac = read_aircraft_model(nameAircraftModel) # MODIFY <path> appropriately
saveName = "JetA1500nmi"
x = [11.967887176329073, 32479.866286551907, 0.5956532812392801, 27.382147149331143, 0.9447567896238773, 0.10644959191186304, 0.14694888676914802, 0.12003039648977087, 1.0596807683444869, 0.9905991041425743, 1749.8757464045927, 14.988409779779435, 2.6302018540124124]
# 2.5) Change fuel type
ac.pari[iifuel] = 32 #(JetA:24 Ethanol:32 JetAEtha31%Blend: 322431)
ac.parg[igrhofuel] = 789.0 #(JetA:817.0 Ethanol:789.0 JetAEtha31%Blend: 805)
# 3) Set the parameters based on optimization result
ac.parg[igAR] = x[1] # Aspect Ratio 
ac.para[iaalt, ipcruise1, :] .=  x[2] * ft_to_m # Cruise Altitude
ac.para[iaCL, ipclimb1+1:ipdescentn-1, :] .= x[3] # CL
ac.parg[igsweep] = x[4] # Wing sweep 
ac.parg[iglambdas] = x[5] #inner_panel_taper_ratio
ac.parg[iglambdat] = x[6] #outer_panel_taper_ratio
ac.parg[ighboxo] = x[7] #root_thickness_to_chord
ac.parg[ighboxs] = x[8] #spanbreak_thickness_to_chord
ac.para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= x[9]   #  rcls    break/root cl ratio = cls/clo
ac.para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= x[10]   #  rclt    tip  /root cl ratio = clt/clo
ac.pare[ieTt4, ipcruise1:ipcruise2, :] .= x[11] # Tt4
ac.pare[iepihc, ipclimb1+1 : ipdescentn-1, :] .= x[12] # High Pressure Compressor Pressure Ratio
ac.pare[iepif, ipclimbn, :] .= x[13] #Fan PR 
ac.pare[iepilc, :, :] .= 3 # Low Pressure Compressure Pressure Ratio set to 3
# 3) Size aircraft
TASOPT.size_aircraft!(ac, iter =500, printiter=true)
# 4) Collect data
AltRec        = [x[2]] #ft
RanRec        = [ac.parg[igRange]./1852.0] #nmi
WMTORec       = [ac.parg[igWMTO]./(9.81*1000)] #Ton (metric)
WFuelRec      = [ac.parg[igWfuel]./(9.81*1000)] #Ton (metric)
WPayRec       = [ac.parg[igWpay]./(9.81*1000)] #Ton (metric)
PFEIRec       = [ac.parm[imPFEI,1]] #(J/J)
WTO_WTOmaxRec = [ac.parm[imWTO,1]/ac.parg[igWMTO]] #N/N
Wf_WfmaxRec   = [ac.parg[igWfuel]/ac.parg[igWfmax]] #N/N
areaWingRec   = [ac.parg[igS]] #m2
ARWingRec     = [ac.parg[igAR]] #Aspect ratio
spanWingRec   = [ac.parg[igb]] #m
diaFanRec     = [ac.parg[igdfan]] #m
FnTotCRRec    = [ac.pare[ieFe,ipcruise1,1]] #N
WEmpRec       = [WMTORec - WFuelRec - WPayRec] # Ton Metric
EFuelRec      = [PFEIRec.*WPayRec*1000*9.81.*RanRec*1852.0] #Joul
SweepRec      = [ac.parg[igsweep]] #deg
outputTup = (AltRec=AltRec,RanRec=RanRec,WMTORec=WMTORec,WFuelRec=WFuelRec
             ,WPayRec=WPayRec,PFEIRec=PFEIRec,WTO_WTOmaxRec=WTO_WTOmaxRec
             ,Wf_WfmaxRec=Wf_WfmaxRec,areaWingRec=areaWingRec,ARWingRec=ARWingRec
             ,spanWingRec=spanWingRec,diaFanRec=diaFanRec,FnTotCRRec=FnTotCRRec,WEmpRec=WEmpRec,EFuelRec=EFuelRec,SweepRec=SweepRec)
CSV.write(saveName*"MissDetail.csv",  outputTup, writeheader=true)
##Create a mask to mask out unreported phases
maskRep = ac.pare[ieFe,:,1].>0 #Reported Phase has non zero thrust
phases = ["ST","RO","TO","CB","B1","B2","B3","B4","B5","C1","C2","D1","D2","D3","D4","D5","Test"]
phases = phases[maskRep]
print("Reported Phases Are:", phases,"\n")
##Read out other parameters
###Aero Parameters
timeOptMiss = ac.para[iatime,maskRep,1] #second
ranOptMiss = ac.para[iaRange,maskRep,1] #meter
altOptMiss = ac.para[iaalt,maskRep,1] #meter
machOptMiss = ac.para[iaMach,maskRep,1]
weiOptMiss = ac.para[iafracW,maskRep,1]*(ac.parm[imWTO,1]/9.81) #kg
gamOptMiss = ac.para[iagamV,maskRep,1] #rad
LDROptMiss = ac.para[iaCL,maskRep,1]./ac.para[iaCD,maskRep,1] #Lift to drag ratio
###Engine Parameters
hfOptMiss = ac.pare[iehfuel,maskRep,1] #J/kg equivalent heating value
TfuelOptMiss = ac.pare[ieTfuel,maskRep,1] #K fuel temperature
Tt3OptMiss = ac.pare[ieTt3,maskRep,1] #K
Pt3OptMiss = ac.pare[iept3,maskRep,1] #Pa 
Tt4OptMiss = ac.pare[ieTt4,maskRep,1] #K
Pt4OptMiss = ac.pare[iept4,maskRep,1] #Pa
FnOptMiss = ac.pare[ieFe,maskRep,1] #N Total Thrust for all engines
mdotfOptMiss = ac.pare[iemcore,maskRep,1].*ac.pare[ieff,maskRep,1] #kg/s for all engines
Cpa = 0.5.*(ac.pare[iecpt3,maskRep,1].+ac.pare[iecpt4,maskRep,1])
ffbMiss   = (Cpa.*(Tt4OptMiss.-Tt3OptMiss))./(hfOptMiss.*ac.pare[ieetab,maskRep,1].-Cpa.*(Tt4OptMiss.-TfuelOptMiss))
mdot3OptMiss = mdotfOptMiss./ffbMiss #kg/s for all engines air flow into the combustor (exclude bypass cooling flow)
#Output Additional Data at the optimal mission
outputTup = (Phase=phases,Time=timeOptMiss,Range=ranOptMiss,Altitude=altOptMiss,MachNumber=machOptMiss,Weight=weiOptMiss
            ,ClimbAngle=gamOptMiss,LiftDragRatio=LDROptMiss,HeatingValue=hfOptMiss,FuelTemp=TfuelOptMiss
            ,Tt3=Tt3OptMiss,Pt3=Pt3OptMiss,Tt4=Tt4OptMiss,Pt4=Pt4OptMiss,Thrust=FnOptMiss,mdotFuel=mdotfOptMiss,mdot3=mdot3OptMiss)
CSV.write(saveName*"MissDetail2.csv",  outputTup, writeheader=true)
#Plot Plane
TASOPT.stickfig(ac, label_fs = 8)
plt.savefig(saveName*"MissDetail3.png")