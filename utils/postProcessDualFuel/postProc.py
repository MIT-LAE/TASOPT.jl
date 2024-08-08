#Baseline postprocessing code for TASOPT Dual Fuel Flight Output
##Print out:
import os
import pandas as pd
import matplotlib.pyplot as plt
keyLst = ["JetA1500","Etha1500","JetAEtha1500","EthaJetA1500","JetA2250","Etha2250","JetAEtha2250","EthaJetA2250","JetA3000","Etha3000","JetAEtha3000","EthaJetA3000"]
ExtNam = "_230MissDetail"
labels = ["Jet Fuel Both Zones","Ethanol Both Zones","Ethanol Main Zone","Ethanol Pilot Zone"]
formLine = ["x","x","x","x","x","x","x","x","x","x","x","x"]
colorLine = ["k","orange","g","r","k","orange","g","r","k","orange","g","r","k","orange","g","r"]
PFEILst = []
RanLst = []
for idxKey, keyCur in enumerate(keyLst):
    desReadNam = keyCur+ExtNam+".csv"
    desRead = pd.read_csv(desReadNam)
    PFEILst.append(desRead["PFEIRec"].values[0])
    RanLst.append(desRead["RanRec"].values[0])

# Create Movie Directory
if os.path.isdir('Movie' + '/') == False:
    os.mkdir('Movie' + '/')
#Plot out PFEI
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxRan in range(len(RanLst)):
    if idxRan>3:
        ax.plot(RanLst[idxRan], PFEILst[idxRan],formLine[idxRan],color=colorLine[idxRan])
    else:
        ax.plot(RanLst[idxRan], PFEILst[idxRan],formLine[idxRan],color=colorLine[idxRan],label=labels[idxRan])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('PFEI [J/J]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.6, top=0.8)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "PFEI.jpg")
plt.close('all')

#Plot out PFEI Change
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxRan in range(3):
    PFEIbase = PFEILst[idxRan*4]
    for idxFuel in range(1,4):
        PFEIPerc = 100*(PFEILst[idxRan*4+idxFuel]-PFEIbase)/PFEIbase
        if (idxRan*4+idxFuel)>3:
             ax.plot(RanLst[idxRan*4+idxFuel], PFEIPerc,formLine[idxRan*4+idxFuel],color=colorLine[idxRan*4+idxFuel])           
        else:
             ax.plot(RanLst[idxRan*4+idxFuel], PFEIPerc,formLine[idxRan*4+idxFuel],color=colorLine[idxRan*4+idxFuel],label=labels[idxRan*4+idxFuel])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('PFEI Percentage Increase [%]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.0, top=20.0)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "fracPFEI.jpg")
plt.close('all')

print("End")
