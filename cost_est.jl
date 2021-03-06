"""
cost_est estimates the program and airframe cost for the aircraft
aircraft cost model developed from DAPCA IV, as found in Raymer's
Aircraft Design (2006)

Inputs:
parg[,] geometric parameters
pare[,] engine parameters

prod_Q  estimated production quantity over 5 years

Outputs:
cost    total program cost estimate

"""

function CostEst(parg, pare, prod_Q)
    Wempty  = (parg[igWMTO] - parg[igWfuel] - parg[igWpay]) / 9.81 #empty weight (kg)
    # TODO - add adjustments based on Markish for component-level weight differences (compared to conv. aircraft)
    Vmax = maximum(pare[ieu0,:]) * 3.6 #maximum velocity (km/h)
    Neng = parg[igneng] #number of engines
    Npax = 180 #number of pax (TODO - make a parameter)
    T4 = pare[ieTt4] * 5/9 #turbine inlet temperature (K)
    Fmax = maximum(pare[ieFe,:]) / Neng / 1000 #maximum thrust per eninge (kN)
    Mmax = maximum(pare[ieM0,:]) #maximum flight Mach number

    PPI = 242.8/144.8 #aerospace producer price index ratio from 1999 to 2020 (from US BLS)
    FTA = 6 #number of flight test aircraft (assume upper end due to novelty of propulsion system tech)

    RE = 86 #engineering wrap rate
    RT = 88 #tooling wrap rate
    RQ = 81 #quality control wrap rate
    RM = 73 #manufacturing wrap rate

    HE = 5.18*Wempty^0.777*Vmax^0.894*prod_Q^0.163 #engineering hours
    HT = 7.22*Wempty^0.777*Vmax^0.696*prod_Q^0.263 #tooling hours
    HM = 10.5*Wempty^0.82*Vmax^0.484*prod_Q^0.641 #manufacturing hours
    HQ = 0.133*HM #quality control hours

    CD = 48.7*Wempty^0.630*Vmax^1.3 #development support cost_est
    CF = 1408*Wempty^0.325*Vmax^0.822*FTA^1.21 #flight test cost
    CM = 22.6*Wempty^0.921*Vmax^0.621*prod_Q^0.799 #manufacturing materials cost (assumes Al construction)

    CE = 2251*(9.66*Fmax + 243.25*Mmax + 1.74*T4 - 2228)*1.2 #engine cost 
    #TODO - adjust to incorporate electric elements of propulsion system
    cA = 4/3 #adjustment factor, such that avionics costs are around 25% of flyaway cost
    CI = 2500*Npax #interior cost (e.g. seats, floors, etc.)
    # source: Raymer (2006), Collinson (2002)

    #TODO - other adjustments specific to H2 (e.g. fuel system)

    cost = PPI*(HE*RE + HT*RT + HM*RM + HQ*RQ + CD + CF + CM + CI + Neng*CE)*cA
    # returns total program cost - divide by prod_Q to yield cost per aircraft

    return cost
end