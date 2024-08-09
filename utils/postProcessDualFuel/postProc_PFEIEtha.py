#Baseline postprocessing code for TASOPT Dual Fuel Flight Output
##Print out:
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PostProcACS_RQLDeFunUnMixSZ import EmisInteg
hf_Etha = 26810134. #[J/kg] heating value of ethanol with evaporation
ResKey    = [["JetA1500_230MissDetail","Etha1500_230MissDetail","JetAEtha1500_230MissDetail","EthaJetA1500_230MissDetail"],["JetA2250_230MissDetail","Etha2250_230MissDetail","JetAEtha2250_230MissDetail","EthaJetA2250_230MissDetail"],["JetA3000_230MissDetail","Etha3000_230MissDetail","JetAEtha3000_230MissDetail","EthaJetA3000_230MissDetail"]]
InpKey    = [["ACS_CFmFpCom_Fli1500_0606_JetAJetA","ACS_CFmFpCom_Fli1500_0606_EthaEtha","ACS_CFmFpCom_Fli1500_0606_JetAEtha","ACS_CFmFpCom_Fli1500_0606_EthaJetA"],["ACS_CFmFpCom_Fli2250_0606_JetAJetA","ACS_CFmFpCom_Fli2250_0606_EthaEtha","ACS_CFmFpCom_Fli2250_0606_JetAEtha","ACS_CFmFpCom_Fli2250_0606_EthaJetA"],["ACS_CFmFpCom_Fli3000_0606_JetAJetA","ACS_CFmFpCom_Fli3000_0606_EthaEtha","ACS_CFmFpCom_Fli3000_0606_JetAEtha","ACS_CFmFpCom_Fli3000_0606_EthaJetA"]]
labels    = ["Jet Fuel Both Zones","Ethanol Both Zones","Ethanol Main Zone","Ethanol Pilot Zone"]
IdxEtha   = [[0,0],[1,1],[0,1],[1,0]]
formLine  = ["x","x","x","x"]
colorLine = ["k","orange","g","r"]
RanLstLst         = []
fECost_EthaLstLst = []
for idxRes in range(len(ResKey)):
    fECost_EthaLst = []
    RanLst         = []
    for idxCas in range(len(ResKey[idxRes])):
        #Get Total Fuel Energy Burnt
        Res = pd.read_csv(ResKey[idxRes][idxCas]+".csv")
        PFEI_Cur = Res["PFEIRec"].values[0] #[J/J]
        Ran_Cur = Res["RanRec"].values[0]*1852.0 #[m]
        WPay_Cur = Res['WPayRec'].values[0]*9.81*1000 #[N]
        EFuel_Cur = PFEI_Cur*Ran_Cur*WPay_Cur #[J] Total Fuel Energy Burnt
        RanLst.append(Ran_Cur/1852.0) #[nmi]
        #Get Ethanol Energy Burnt
        Inp = pd.read_csv(InpKey[idxRes][idxCas]+".csv")
        IdxEtha_Cur = IdxEtha[idxCas]
        Weight_Cur = Inp["Weight"].values[:] #[kg] Weight of the aircraft
        EIEtha_Cur = (Inp['Wf[kg/s]'].values[:]*IdxEtha_Cur[0] + Inp['WAS[kg/s]'].values[:]*IdxEtha_Cur[1])/Inp['WfTot'].values[:] #[kg/kg] Consumption Index of Ethanol Per Unit Fuel Burn
        InteEmis = EmisInteg(Weight_Cur,EIEtha_Cur)
        mEtha_Cur = InteEmis["mPolTot"] #[kg]
        EEtha_Cur = hf_Etha*mEtha_Cur #[J] Total Ethanol Energy Burnt
        ##Compute Ethanol Energy Cost
        if idxCas == 0: #Baseline Pure Jet Fuel Case
            EFuel_Base = EFuel_Cur #[J]
            fECost_EthaLst.append(0)
            if EEtha_Cur != 0.0:
                print("Warning: First Case has Ethanol: ",EEtha_Cur," J")
        else:
            Del_EFuel = EFuel_Cur-EFuel_Base #[J] Extra Energy Burnt
            fECost_Etha = Del_EFuel/EEtha_Cur #[J/J] Energy cost per unit ethanol burnt
            fECost_EthaLst.append(fECost_Etha)
    fECost_EthaLstLst.append(fECost_EthaLst)
    RanLstLst.append(RanLst)
fECost_EthaLstLst=np.transpose(fECost_EthaLstLst)
RanLstLst=np.transpose(RanLstLst)

# Create Movie Directory
if os.path.isdir('Movie' + '/') == False:
    os.mkdir('Movie' + '/')
#Plot out PFEI
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(1,len(labels)):
    ax.plot(RanLstLst[idxCas], fECost_EthaLstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Energy Penalty from Ethanol [J/J]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.04, top=0.18)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EthaPFEI.jpg")
plt.close('all')