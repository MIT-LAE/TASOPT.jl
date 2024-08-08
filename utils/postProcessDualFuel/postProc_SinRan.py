#Baseline postprocessing code for TASOPT Dual Fuel Flight Output
##Print out:
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
keyLst = ["JetA3000","Etha3000","JetAEtha3000","EthaJetA3000"]
ExtNam = "_230MissDetail"
labels = ["Jet Fuel Both Zones","Ethanol Both Zones","Ethanol Main Zone","Ethanol Pilot Zone"]
formLine = ["x","x","x","x"]
colorLine = ["k","orange","g","r"]
AltLst = []
WMTOLst = []
fracFCarLst = [] #Percentage of fuel that is store inside the cargo 
spanWingLst = []
diaFanLst  = [] #Fan Diameter
for idxKey, keyCur in enumerate(keyLst):
    desReadNam = keyCur+ExtNam+".csv"
    desRead = pd.read_csv(desReadNam)
    AltLst.append(desRead["AltRec"].values[0]) #[ft]
    WMTOLst.append(desRead["WMTORec"].values[0]) #[ton]
    fracFCarLst.append(100*(desRead["Wf_WfmaxRec"].values[0]-1)/desRead["Wf_WfmaxRec"].values[0]) #[%]
    spanWingLst.append(desRead["spanWingRec"].values[0]) #[m]
    diaFanLst.append(desRead["diaFanRec"].values[0]) #[m]
#Process Para into Percentage
AltLst = np.array(AltLst)
WMTOLst = np.array(WMTOLst)
fracFCarLst = np.array(fracFCarLst)
spanWingLst = np.array(spanWingLst)
diaFanLst = np.array(diaFanLst)
fAltLst = (AltLst-AltLst[0])*100/AltLst[0]
fWMTOLst = (WMTOLst-WMTOLst[0])*100/WMTOLst[0]
ffracFCarLst = (fracFCarLst-fracFCarLst[0])*100/fracFCarLst[0]
fspanWingLst = (spanWingLst-spanWingLst[0])*100/spanWingLst[0]
fdiaFanLst = (diaFanLst-diaFanLst[0])*100/diaFanLst[0]
# Create Movie Directory
if os.path.isdir('Movie' + '/') == False:
    os.mkdir('Movie' + '/')
#Plot MultiMiss Comparison
bigPhases = ("$\%\Delta Alt_{cruise}$","$\%\Delta WTO_{max}$","$\%\Delta frac_{fuel,cargo}$","$\%\Delta span_{wing}$","$\%\Delta Dia_{fan}$")
PropMiss = {}
for indKey in range(1,len(labels)):
    PropMiss[labels[indKey]] = (fAltLst[indKey],fWMTOLst[indKey],ffracFCarLst[indKey],fspanWingLst[indKey],fdiaFanLst[indKey])
x = np.arange(len(bigPhases))  # the label locations
width = 0.25  # the width of the bars
multiplier = 1
fig, ax = plt.subplots(layout='constrained',dpi = 300)
for attribute, measurement in PropMiss.items():
    offset = width * (multiplier-1)
    rects = ax.bar(x + offset, measurement, width, label=attribute, color=colorLine[multiplier])
    ax.bar_label(rects, padding=0,fontsize=9,fmt='%.1f')
    multiplier += 1
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('%Change based to Jet Fuel Both Zones Injection Case')
ax.set_xticks(x + width, bigPhases)
ax.legend()
ax.yaxis.grid()
ax.set_ylim(bottom=-10, top=60)
plt.savefig('Movie' + '/' + "PropMiss.jpg")