import cantera as ct
import numpy as np
#User inputs
fracfEneMai = 0.69 #0.69 #Fraction of energy in main stage
hfPil = 26810134 #J/kg 26810134 for Ethanol & 43215118 for JetA
hfMai = 43215118 #J/kg
fuelStrYPil = "C2H5OH:1"#"NC12H26:0.4049576676308857,IC8H18:0.014289267305708273,C7H8:0.1531467916175248,IC16H34:0.29909723934977994,DECALIN:0.12850903409610107"
fuelStrYMai = "NC12H26:0.4049576676308857,IC8H18:0.014289267305708273,C7H8:0.1531467916175248,IC16H34:0.29909723934977994,DECALIN:0.12850903409610107"#"C2H5OH:1" #Fuel String
reacMech = "CRECKButanolHLNOx.yaml"#"CRECKButanolHLNOx.yaml"#"CRECKButanolHLNOx.yaml"#"CRECKButanolHLNOx.yaml"#"CRECKButanolHLNOx.yaml" #Reaction Mech
Tstd = 298.15 #298.1008134#298.1008134 #Temperature in K
Pstd = 101325.0 #101325.0 #Pressure in Pa
dT = 0.001 #Temperature Perturbation on both side in K for derivative calculation
Tstd_sutherland = 273 #Base Temperature for sutherland's law
rhoPil = 789.0 #For total desity calculation (JetA:817.0 Ethanol:789.0)
rhoMai = 817.0
########Generate Temperature Range
#### Uniform Spread Temperature

nDim = 73 #Number of data points
Tsta = 200.00 #Start T in K
Tend = 2000.00 #
tList = np.linspace(Tsta,Tend,nDim)
print("\nnDim: ",nDim)
print("tList: \n",list(tList),"\n")

#### Customized Spread Temperature
"""
tList = np.array([
        175.00e0, 200.00e0, 225.00e0, 250.00e0, 275.00e0, 300.00e0,
        325.00e0, 350.00e0, 375.00e0, 400.00e0, 450.00e0, 500.00e0,
        550.00e0, 600.00e0, 650.00e0, 700.00e0, 750.00e0, 800.00e0,
        850.00e0, 900.00e0, 950.00e0, 1000.00e0, 1050.00e0, 1100.00e0,
        1150.00e0, 1200.00e0, 1250.00e0, 1300.00e0, 1350.00e0, 1400.00e0,
        1500.00e0, 1600.00e0, 1700.00e0, 1800.00e0, 1900.00e0, 2000.00e0,
        2100.00e0, 2200.00e0, 2300.00e0, 2400.00e0, 2500.00e0, 2600.00e0,
        2700.00e0, 2800.00e0, 2900.00e0, 3000.00e0, 3500.00e0, 4000.00e0,
        4500.00e0, 5000.00e0, 5500.00e0, 6000.00e0])
nDim = len(tList)
print("\nnDim: ",nDim)
print("tList: \n",list(tList),"\n")
"""
########User inputs End

#Compute Main Zone Fuel Mass Fraction
fracfMdotMai = (1+((1-fracfEneMai)*hfMai)/(fracfEneMai*hfPil))**(-1) #Mass Fraction of fuel to the main zone
print("Mass fraction of fuel to main zone is: ",fracfMdotMai)
#Blend the fuel
gas = ct.Solution(reacMech)
gas.TPY = Tstd,Pstd,fuelStrYPil
YPil = gas.Y*1.0
gas.TPY = Tstd,Pstd,fuelStrYMai
YMai = gas.Y*1.0
fuelY = fracfMdotMai*YMai + (1-fracfMdotMai)*YPil
print("Combined Fuel Mass Fraction: ",fuelY[fuelY>0])
#Generate LogTemperature Range
tlList = np.log(tList)
print("tLogList: \n",list(tlList),"\n")
#Initialize gas solution
gasB = ct.Solution(reacMech) #Base Temperature
gasU = ct.Solution(reacMech) #Upper Temperature (Derivative Calculation)
gasD = ct.Solution(reacMech) #Lower Temperature (Derivative Calculation)
gasB.TPY = Tstd,Pstd,fuelY
gasU.TPY = Tstd,Pstd,fuelY
gasD.TPY = Tstd,Pstd,fuelY
MW = gasB.mean_molecular_weight/1000 #kg/mol
r = 8.3144598/MW #J/kgK
print("r: [J/(kg*K)]",r)
hform = gasB.enthalpy_mass #J/kg
print("hform: [J/kg]",hform)
sform = gasB.entropy_mass #J/kgK
print("sform: [J/kgK]",sform)
#Use for loop to compute thermal properties
cp = []
cpt = []
h = []
s = []
for idxT, TCur in enumerate(tList):
    gasB.TP = TCur,Pstd
    cp.append(gasB.cp_mass) #Cp
    gasU.TP = (TCur+dT),Pstd
    gasD.TP = (TCur-dT),Pstd
    cpt.append((gasU.cp_mass-gasD.cp_mass)/(2*dT)) #dCp/dT
    h.append(gasB.enthalpy_mass - hform) #Sensible Ethanlpy (J/kg)
    s.append(gasB.entropy_mass - sform) #Sensible Entropy (J/kgK)
print("cp: \n",cp,"\n")
print("cpt: \n",cpt,"\n")
print("h: \n",h,"\n")
print("s: \n",s,"\n")
#Compute the element concentration
elemMoleCHON = np.array([gasB.elemental_mole_fraction("C"),gasB.elemental_mole_fraction("H"),gasB.elemental_mole_fraction("O"),gasB.elemental_mole_fraction("N")])
elemMoleCHON = elemMoleCHON/np.min(elemMoleCHON[elemMoleCHON>0])
print("Element Number C,H,O,N:",elemMoleCHON)
#Compute Transport Properties Sutherland's Law
gasB.TP = Tstd_sutherland,Pstd
T0 = gasB.T
mu_0 = gasB.viscosity
k_0 = gasB.thermal_conductivity
print("Sutherland Based: ","T0[K]: ",T0,"& mu_0[Pa-s]: ",mu_0,"& k_0[W/m/K]: ",k_0)
tList_Sutherland = np.linspace(300,2700,100)
S_mu_Lst = []
S_k_Lst = []
for idxT, TCur in enumerate(tList_Sutherland):
    gasB.TP = TCur,Pstd
    C_mu = gasB.viscosity/(mu_0*((TCur/T0)**(3/2)))
    C_k = gasB.thermal_conductivity/(k_0*((TCur/T0)**(3/2)))
    S_mu = (C_mu*TCur-T0)/(1-C_mu)
    S_k = (C_k*TCur-T0)/(1-C_k)
    S_mu_Lst.append(S_mu)
    S_k_Lst.append(S_k)
SmuAve = np.mean(S_mu_Lst)
SkAve = np.mean(S_k_Lst)
SmuSTD = np.std(S_mu_Lst)*100/SmuAve
SkSTD = np.std(S_k_Lst)*100/SkAve
print("Sutherland Based: ","SmuAve: ",SmuAve,"& SkAve: ",SkAve)
print("Sutherland Based: ","SmuSTD[%]: ",SmuSTD,"& SkSTD[%]: ",SkSTD)
#Compute the average density of the fuel
rhoTot = ((1/rhoPil)*(1-fracfMdotMai)+(1/rhoMai)*fracfMdotMai)**(-1) 
print("Average density: [kg/m3]: ",rhoTot)
    



    

