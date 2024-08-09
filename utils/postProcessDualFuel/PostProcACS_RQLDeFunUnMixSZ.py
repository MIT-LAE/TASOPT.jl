import os
from scipy import interpolate,integrate
os.environ["OMP_NUM_THREADS"] = "1"


def EmisInteg(MassPlane,EI):#Emissions Intergrator
    #MassPlane: Mass of the airplane in [kg]
    #Emissions Index: [g/kg]
    #Get The Mass of Fuel Burnt
    MassFuel = MassPlane[0]-MassPlane #kg Mass of fuel consumed
    #Cut out the function into different pieces
    mFuelTO = MassFuel[0:2] #ST and RO Phases
    mFuelCL = MassFuel[2:7] #B1 to B5
    mFuelCR = MassFuel[7:9] #C1 to C2
    mFuelDE = MassFuel[9:] #D1 to D5
    EITO = EI[0:2] #ST and RO Phases
    EICL = EI[2:7] #B1 to B5
    EICR = EI[7:9] #C1 to C2
    EIDE = EI[9:] #D1 to D5
    #Function of EI
    EITOFun = interpolate.interp1d(mFuelTO,EITO,kind="linear",fill_value="extrapolate") #EI(mFuel)
    EICLFun = interpolate.interp1d(mFuelCL,EICL,kind="cubic",fill_value="extrapolate") #EI(mFuel)
    EICRFun = interpolate.interp1d(mFuelCR,EICR,kind="linear",fill_value="extrapolate") #EI(mFuel)
    EIDEFun = interpolate.interp1d(mFuelDE,EIDE,kind="cubic",fill_value="extrapolate") #EI(mFuel)
    #Integrate through the Emissions produced in each phase
    mPolTO = integrate.quad(EITOFun,mFuelTO[0],mFuelTO[-1])[0] #[kg] Total pollutant produced
    mPolCL = integrate.quad(EICLFun,mFuelCL[0],mFuelCL[-1])[0] #[kg]
    mPolCR = integrate.quad(EICRFun,mFuelCR[0],mFuelCR[-1])[0] #[kg]
    mPolDE = integrate.quad(EIDEFun,mFuelDE[0],mFuelDE[-1])[0] #[kg]
    mPolTot = mPolTO+mPolCL+mPolCR+mPolDE #[kg]
    return {"mPolTO":mPolTO,"mPolCL":mPolCL,"mPolCR":mPolCR,"mPolDE":mPolDE,"mPolTot":mPolTot}