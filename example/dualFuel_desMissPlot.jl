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
Ethanol 1500: [11.995905780175995, 30485.17133529457, 0.5861830080074717, 27.83229852894221, 0.9401002404881813, 0.10310257798367431, 0.14984823368989658, 0.1261417727099628, 1.0614520934688336, 0.9993743124675745, 1749.587188323032, 14.998204289919045, 1.2724983690705405]    PFEI = 0.6533086505204532
Jet A   1500: [11.999617064124072, 31416.337576476286, 0.5992372107935361, 26.611195142395566, 0.8727598739228938, 0.10516478417308059, 0.1323530500901155, 0.1160935275052842, 1.1422855535237826, 0.9248590012462885, 1756.4958341529, 14.46487481499291, 1.4952189910099465]      PFEI = 0.6210760601248294
EthaJet 1500: [11.783893284058767, 31554.206888006105, 0.5768529573646356, 26.957606731479828, 0.9148646182695386, 0.12762781006928817, 0.1465610470464057, 0.12212882582998263, 1.1239760831273382, 0.8218721514277865, 1750.0088592006487, 14.930006848517355, 1.8338981449427103] PFEI = 0.6295258244546431
JetEtha 1500: [11.993075682683394, 30979.683483492467, 0.587091284820852, 27.777903528882714, 0.9504410725294687, 0.10504792915253236, 0.14821718935648662, 0.12613307372499902, 1.0546095476477841, 0.9595866434331847, 1750.6145512216945, 14.999005531097268, 1.8291310255710174] PFEI = 0.6452195127719252
Ethanol 3000: [11.7471582687464, 26729.063309342208, 0.5867387034272569, 27.769578896529747, 0.9144118980218298, 0.10176872301167492, 0.1499263303390325, 0.13028897338891732, 1.0896505831736039, 0.9750851475660576, 1766.083818184854, 14.827156913557616, 1.7994431664346604]    PFEI = 0.7999648883784056 
Jet A   3000: [11.99111053156572, 29444.882223445613, 0.5803619852473172, 27.72271286593455, 0.8385219382082921, 0.10144815339994472, 0.1492477328487266, 0.1257171304842019, 1.2394756830674365, 0.8114091756231697, 1788.490437444581, 14.831623197640358, 1.3554437947040672]     PFEI = 0.6849830753329049
EthaJet 3000: [11.979612806397814, 27691.251529245616, 0.5820175975162574, 27.61760629203382, 0.9203163081458398, 0.10153809941108247, 0.1488154484814802, 0.12678228028941435, 1.0920297288543899, 0.9188292667969202, 1795.7910521806994, 14.969907140889223, 1.2058189258251035]  PFEI = 0.7075994185989157
JetEtha 3000: [11.912918744799512, 26728.98273750468, 0.580659393330651, 27.66348262225188, 0.8707065663358251, 0.10637170507852939, 0.14782584258305143, 0.1300634226427048, 1.133565476774562, 0.9032372597136605, 1773.4626849145566, 14.654752137545362, 1.403347378421579]      PFEI = 0.7671913129937321
"""
savedir = "Movie/"
if !isdir(savedir)
    # If it doesn't exist, create the "optimization" directory
    mkdir(savedir)
end
# 2) Include input file for desired aircraft/
nameAircraftModel = "../src/IO/experiment_input_3000.toml"
ac = read_aircraft_model(nameAircraftModel) # MODIFY <path> appropriately
saveName = savedir*"JetA3000nmi"
x = [11.380434406953572, 33198.13022096391,  0.5410047569882184, 26.66666248716794,  0.9378606924235586, 0.10281538084027678, 0.14977259425119216, 0.13296003568212714, 1.0627223742965342, 0.996592019369539,  1737.0012718430448, 14.997879755556585, 0.30657837795564297]
# 2.5) Change fuel type
ac.pari[iifuel] = 24 #(JetA:24 Ethanol:32 JetAEtha31%Blend: 322431)
ac.parg[igrhofuel] = 817.0 #(JetA:817.0 Ethanol:789.0 JetAEtha31%Blend: 805)
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
FnOptMiss = ac.pare[ieFe,maskRep,1] #N Total Thrust for each engines
mdotfOptMiss = ac.pare[iemcore,maskRep,1].*ac.pare[ieff,maskRep,1] #kg/s for eacg engines
Cpa = 0.5.*(ac.pare[iecpt3,maskRep,1].+ac.pare[iecpt4,maskRep,1])
ffbMiss   = (Cpa.*(Tt4OptMiss.-Tt3OptMiss))./(hfOptMiss.*ac.pare[ieetab,maskRep,1].-Cpa.*(Tt4OptMiss.-TfuelOptMiss))
mdot3OptMiss = mdotfOptMiss./ffbMiss #kg/s for each engines air flow into the combustor (exclude bypass cooling flow)
#Output Additional Data at the optimal mission
outputTup = (Phase=phases,Time=timeOptMiss,Range=ranOptMiss,Altitude=altOptMiss,MachNumber=machOptMiss,Weight=weiOptMiss
            ,ClimbAngle=gamOptMiss,LiftDragRatio=LDROptMiss,HeatingValue=hfOptMiss,FuelTemp=TfuelOptMiss
            ,Tt3=Tt3OptMiss,Pt3=Pt3OptMiss,Tt4=Tt4OptMiss,Pt4=Pt4OptMiss,Thrust=FnOptMiss,mdotFuel=mdotfOptMiss,mdot3=mdot3OptMiss)
CSV.write(saveName*"MissDetail2.csv",  outputTup, writeheader=true)
#Plot Plane
TASOPT.stickfig(ac, label_fs = 8)
plt.savefig(saveName*"MissDetail3.png")