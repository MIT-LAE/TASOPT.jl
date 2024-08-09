#Baseline postprocessing code for TASOPT Dual Fuel Flight Output
##Print out:
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PostProcACS_RQLDeFunUnMixSZ import EmisInteg
hf_Etha = 25761470. #[J/kg] heating value of ethanol with evaporation
ResKey    = [["JetA1500_230MissDetail","Etha1500_230MissDetail","JetAEtha1500_230MissDetail","EthaJetA1500_230MissDetail"],["JetA2250_230MissDetail","Etha2250_230MissDetail","JetAEtha2250_230MissDetail","EthaJetA2250_230MissDetail"],["JetA3000_230MissDetail","Etha3000_230MissDetail","JetAEtha3000_230MissDetail","EthaJetA3000_230MissDetail"]]
InpKey    = [["ACS_CFmFpCom_Fli1500_0606_JetAJetA","ACS_CFmFpCom_Fli1500_0606_EthaEtha","ACS_CFmFpCom_Fli1500_0606_JetAEtha","ACS_CFmFpCom_Fli1500_0606_EthaJetA"],["ACS_CFmFpCom_Fli2250_0606_JetAJetA","ACS_CFmFpCom_Fli2250_0606_EthaEtha","ACS_CFmFpCom_Fli2250_0606_JetAEtha","ACS_CFmFpCom_Fli2250_0606_EthaJetA"],["ACS_CFmFpCom_Fli3000_0606_JetAJetA","ACS_CFmFpCom_Fli3000_0606_EthaEtha","ACS_CFmFpCom_Fli3000_0606_JetAEtha","ACS_CFmFpCom_Fli3000_0606_EthaJetA"]]
labels    = ["Jet Fuel Both Zones","Ethanol Both Zones","Ethanol Main Zone","Ethanol Pilot Zone"]
IdxEtha   = [[0,0],[1,1],[0,1],[1,0]]
formLine  = ["x","x","x","x"]
colorLine = ["k","orange","g","r"]
RanLstLst         = []
PFEILstLst        = []
fECost_Etha_TO_LstLst = []
fECost_Etha_CL_LstLst = []
fECost_Etha_CR_LstLst = []
fECost_Etha_DE_LstLst = []
fECost_Etha_TT_LstLst = []
fracEFuel_CRLstLst = []
fracEFuel_DELstLst = []
for idxRes in range(len(ResKey)):
    fECost_Etha_TO_Lst = []
    fECost_Etha_CL_Lst = []
    fECost_Etha_CR_Lst = []
    fECost_Etha_DE_Lst = []
    fECost_Etha_TT_Lst = []
    fracEFuel_CRLst = []
    fracEFuel_DELst = []
    RanLst         = []
    PFEILst        = []
    for idxCas in range(len(ResKey[idxRes])):
        #For PFEI Calculation
        Res = pd.read_csv(ResKey[idxRes][idxCas]+".csv")
        Ran_Cur = Res["RanRec"].values[0]*1852.0 #[m]
        RanLst.append(Ran_Cur/1852.0) #[nmi]
        WPay_Cur = Res['WPayRec'].values[0]*9.81*1000 #[N]
        #Get Ethanol Energy Burnt
        Inp = pd.read_csv(InpKey[idxRes][idxCas]+".csv")
        IdxEtha_Cur = IdxEtha[idxCas]
        Weight_Cur = Inp["Weight"].values[:] #[kg] Weight of the aircraft
        EIEtha_Cur = (Inp['Wf[kg/s]'].values[:]*IdxEtha_Cur[0] + Inp['WAS[kg/s]'].values[:]*IdxEtha_Cur[1])/Inp['WfTot'].values[:] #[kg/kg] Consumption Index of Ethanol Per Unit Fuel Burn
        InteEmis = EmisInteg(Weight_Cur,EIEtha_Cur)
        EEtha_TO = hf_Etha*InteEmis["mPolTO"] #[J]
        EEtha_CL = hf_Etha*InteEmis["mPolCL"] #[J]
        EEtha_CR = hf_Etha*InteEmis["mPolCR"] #[J]
        EEtha_DE = hf_Etha*InteEmis["mPolDE"] #[J]
        EEtha_TT = hf_Etha*InteEmis["mPolTot"] #[J]
        #Get Total Fuel Energy at different phases
        hf_Cur = Inp["HeatingValue"].values[:] #[J/kg] Heating Value at different phase of flight
        InteEmis = EmisInteg(Weight_Cur,hf_Cur) #[J] Energy at different phase
        EFuel_TO = InteEmis["mPolTO"] #[J]Total Energy Burnt at takeoff
        EFuel_CL = InteEmis["mPolCL"] #[J]
        EFuel_CR = InteEmis["mPolCR"] #[J]
        EFuel_DE = InteEmis["mPolDE"] #[J]
        EFuel_TT = InteEmis["mPolTot"] #[J]
        fracEFuel_CRLst.append(100*EFuel_CR/(EFuel_TO+EFuel_CL+EFuel_CR+EFuel_DE)) #Percentage Energy in Cruise
        fracEFuel_DELst.append(100*EFuel_DE/(EFuel_TO+EFuel_CL+EFuel_CR+EFuel_DE)) #Percentage Energy in Cruise
        if idxCas == 0: #Baseline Pure Jet Fuel Case
            EFuel_TO_Base = EFuel_TO #[J]
            EFuel_CL_Base = EFuel_CL #[J]
            EFuel_CR_Base = EFuel_CR #[J]
            EFuel_DE_Base = EFuel_DE #[J]
            EFuel_TT_Base = EFuel_TT #[J]
            fECost_Etha_TO_Lst.append(0)
            fECost_Etha_CL_Lst.append(0)
            fECost_Etha_CR_Lst.append(0)
            fECost_Etha_DE_Lst.append(0)
            fECost_Etha_TT_Lst.append(0)
        else:
            Del_EFuel = EFuel_TO-EFuel_TO_Base #[J] Extra Energy Burnt
            fECost_Etha_TO_Lst.append(Del_EFuel/EEtha_TO) #[J/J] Energy cost per unit ethanol burnt
            Del_EFuel = EFuel_CL-EFuel_CL_Base #[J] Extra Energy Burnt
            fECost_Etha_CL_Lst.append(Del_EFuel/EEtha_CL) #[J/J] Energy cost per unit ethanol burnt
            Del_EFuel = EFuel_CR-EFuel_CR_Base #[J] Extra Energy Burnt
            fECost_Etha_CR_Lst.append(Del_EFuel/EEtha_CR) #[J/J] Energy cost per unit ethanol burnt
            Del_EFuel = EFuel_DE-EFuel_DE_Base #[J] Extra Energy Burnt
            fECost_Etha_DE_Lst.append(Del_EFuel/EEtha_DE) #[J/J] Energy cost per unit ethanol burnt
            Del_EFuel = EFuel_TT-EFuel_TT_Base #[J] Extra Energy Burnt
            fECost_Etha_TT_Lst.append(Del_EFuel/EEtha_TT) #[J/J] Energy cost per unit ethanol burnt
        #Calculate PFEI
        PFEILst.append(EFuel_TT/(Ran_Cur*WPay_Cur))#[J/J]
    fECost_Etha_TO_LstLst.append(fECost_Etha_TO_Lst)
    fECost_Etha_CL_LstLst.append(fECost_Etha_CL_Lst)
    fECost_Etha_CR_LstLst.append(fECost_Etha_CR_Lst)
    fECost_Etha_DE_LstLst.append(fECost_Etha_DE_Lst)
    fECost_Etha_TT_LstLst.append(fECost_Etha_TT_Lst)
    RanLstLst.append(RanLst)
    PFEILstLst.append(PFEILst)
    fracEFuel_CRLstLst.append(fracEFuel_CRLst)
    fracEFuel_DELstLst.append(fracEFuel_DELst)
fECost_Etha_TO_LstLst = np.transpose(fECost_Etha_TO_LstLst)
fECost_Etha_CL_LstLst = np.transpose(fECost_Etha_CL_LstLst)
fECost_Etha_CR_LstLst = np.transpose(fECost_Etha_CR_LstLst)
fECost_Etha_DE_LstLst = np.transpose(fECost_Etha_DE_LstLst)
fECost_Etha_TT_LstLst = np.transpose(fECost_Etha_TT_LstLst)
fracEFuel_CRLstLst    = np.transpose(fracEFuel_CRLstLst)
fracEFuel_DELstLst    = np.transpose(fracEFuel_DELstLst)
RanLstLst             = np.transpose(RanLstLst)
PFEILstLst            = np.transpose(PFEILstLst)

# Create Movie Directory
if os.path.isdir('Movie' + '/') == False:
    os.mkdir('Movie' + '/')
#Plot out PFEI
fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(1,len(labels)):
    ax.plot(RanLstLst[idxCas], fECost_Etha_TT_LstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Ethanol induced Energy Penalty [J/J]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.04, top=0.18)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EthaPFEI_TT.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(1,len(labels)):
    ax.plot(RanLstLst[idxCas], fECost_Etha_TO_LstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Ethanol induced Energy Penalty in Take-off [J/J]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.0, top=0.4)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EthaPFEI_TO.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(1,len(labels)):
    ax.plot(RanLstLst[idxCas], fECost_Etha_CL_LstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Ethanol induced Energy Penalty in Climb-out [J/J]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.0, top=0.4)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EthaPFEI_CL.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(1,len(labels)):
    ax.plot(RanLstLst[idxCas], fECost_Etha_CR_LstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Ethanol induced Energy Penalty in Cruise [J/J]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0.04, top=0.18)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EthaPFEI_CR.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(1,len(labels)):
    ax.plot(RanLstLst[idxCas], fECost_Etha_DE_LstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Ethanol induced Energy Penalty in Descent [J/J]')
ax.set_xlim(left=1400, right=3200)
# ax.set_ylim(bottom=0.04, top=0.18)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EthaPFEI_DE.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(len(labels)):
    ax.plot(RanLstLst[idxCas], fracEFuel_CRLstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Percentage Energy Spent on Cruise [%]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0, top=100)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EnergySpentCR.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(len(labels)):
    ax.plot(RanLstLst[idxCas], fracEFuel_DELstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('Percentage Energy Spent on Descent [%]')
ax.set_xlim(left=1400, right=3200)
ax.set_ylim(bottom=0, top=4)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "EnergySpentDE.jpg")
plt.close('all')

fig = plt.figure(dpi = 300)
ax  = fig.add_subplot(1,1,1)
for idxCas in range(len(labels)):
    ax.plot(RanLstLst[idxCas], PFEILstLst[idxCas], formLine[idxCas],color=colorLine[idxCas],label=labels[idxCas])
ax.set_xlabel('Range [nmi]')
ax.set_ylabel('PFEI [J/J]')
ax.set_xlim(left=1400, right=3200)
# ax.set_ylim(bottom=0.04, top=0.18)
ax.grid(True)
plt.legend()
plt.savefig('Movie' + '/' + "PFEI.jpg")
plt.close('all')