#Baseline postprocessing code for TASOPT Dual Fuel Flight Output
##Print out:
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PostProcACS_RQLDeFunUnMixSZ import EmisInteg
hf_Etha = 26810134 #[J/kg] heating value of ethanol with evaporation
hf_JetA = 43215118 #[J/kg] heating value of jet fuel with evaporation
ResKey    = [["JetA1500_230MissDetail","Etha1500_230MissDetail","JetAEtha1500_230DMissDetail","EthaJetA1500_230DMissDetail"],["JetA2250_230MissDetail","Etha2250_230MissDetail","JetAEtha2250_230MissDetail","EthaJetA2250_230MissDetail"],["JetA3000_230MissDetail","Etha3000_230MissDetail","JetAEtha3000_230MissDetail","EthaJetA3000_230MissDetail"]]
InpKey    = [["ACS_CFmFpCom_Fli1500_0606_JetAJetA","ACS_CFmFpCom_Fli1500_0606_EthaEtha","ACS_CFmFpCom_Fli1500_0606_JetAEthaD","ACS_CFmFpCom_Fli1500_0606_EthaJetAD"],["ACS_CFmFpCom_Fli2250_0606_JetAJetA","ACS_CFmFpCom_Fli2250_0606_EthaEtha","ACS_CFmFpCom_Fli2250_0606_JetAEtha","ACS_CFmFpCom_Fli2250_0606_EthaJetA"],["ACS_CFmFpCom_Fli3000_0606_JetAJetA","ACS_CFmFpCom_Fli3000_0606_EthaEtha","ACS_CFmFpCom_Fli3000_0606_JetAEtha","ACS_CFmFpCom_Fli3000_0606_EthaJetA"]]
labels    = ["Jet Fuel Both Zones","Ethanol Both Zones","Ethanol Main Zone","Ethanol Pilot Zone"]
IdxEtha   = [[0,0],[1,1],[0,1],[1,0]]
IdxJetA   = 1-np.array(IdxEtha)
formLine  = ["x","x","x","x"]
colorLine = ["C0","orange","g","r"]
RanLstLst    = []
PFEILstLst   = []
fPFEI_LstLst = []
for idxRes in range(len(ResKey)):
    RanLst    = []
    PFEILst   = []
    fPFEI_Lst = []
    for idxCas in range(len(ResKey[idxRes])):
        #For PFEI Calculation
        Res = pd.read_csv(ResKey[idxRes][idxCas]+".csv")
        Ran_Cur = Res["RanRec"].values[0]*1852.0 #[m]
        RanLst.append(Ran_Cur/1852.0) #[nmi]
        WPay_Cur = Res['WPayRec'].values[0]*9.81*1000 #[N]
        #Get Ethanol Energy Burnt
        Inp = pd.read_csv(InpKey[idxRes][idxCas]+".csv")
        IdxEtha_Cur = IdxEtha[idxCas]
        IdxJetA_Cur = IdxJetA[idxCas]
        Weight_Cur = Inp["Weight"].values[:] #[kg] Weight of the aircraft
        EIEtha_Cur = (Inp['Wf[kg/s]'].values[:]*IdxEtha_Cur[0] + Inp['WAS[kg/s]'].values[:]*IdxEtha_Cur[1])/Inp['WfTot'].values[:] #[kg/kg] Consumption Index of Ethanol Per Unit Fuel Burn
        EIJetA_Cur = (Inp['Wf[kg/s]'].values[:]*IdxJetA_Cur[0] + Inp['WAS[kg/s]'].values[:]*IdxJetA_Cur[1])/Inp['WfTot'].values[:] #[kg/kg] Consumption Index of Jet Fuel Per Unit Fuel Burn
        InteEmis = EmisInteg(Weight_Cur,EIEtha_Cur)
        EEtha_TO = hf_Etha*InteEmis["mPolTO"] #[J]
        EEtha_CL = hf_Etha*InteEmis["mPolCL"] #[J]
        EEtha_CR = hf_Etha*InteEmis["mPolCR"] #[J]
        EEtha_DE = hf_Etha*InteEmis["mPolDE"] #[J]
        EEtha_TT = hf_Etha*InteEmis["mPolTot"] #[J]
        InteEmis = EmisInteg(Weight_Cur,EIJetA_Cur)
        EJetA_TT = hf_JetA*InteEmis["mPolTot"] #[J]
        #Calculate PFEI
        PFEI_TT = (EEtha_TT+EJetA_TT)/(Ran_Cur*WPay_Cur) #[J/J]
        PFEILst.append(PFEI_TT)#[J/J]
        if idxCas == 0: #Baseline Pure Jet Fuel Case
            PFEI_TT_Base = PFEI_TT
            fPFEI_Lst.append(0)
        else:
            Del_PFEI = PFEI_TT-PFEI_TT_Base # Change of PFEI
            fPFEI_Lst.append(100*Del_PFEI/PFEI_TT_Base) #[%] Relative Change of PFEI
    RanLstLst.append(RanLst)
    PFEILstLst.append(PFEILst)
    fPFEI_LstLst.append(fPFEI_Lst)
RanLstLst             = np.transpose(RanLstLst)
PFEILstLst            = np.transpose(PFEILstLst)
fPFEI_LstLst          = np.transpose(fPFEI_LstLst)

# Create Movie Directory
if os.path.isdir('Movie' + '/') == False:
    os.mkdir('Movie' + '/')
#Plot out PFEI
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(len(labels)):
    ax.plot(RanLstLst[idxCas], PFEILstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('PFEI [J/J]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.6, top=0.85)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "PFEI.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(1,len(labels)):
    ax.plot(RanLstLst[idxCas], fPFEI_LstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('PFEI Percentage Increase [%]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.0, top=20.0)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "fracPFEI.jpg")
plt.close('all')
