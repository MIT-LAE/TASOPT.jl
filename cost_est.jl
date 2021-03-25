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

function CostEst(parg, pare, parm, parpt, prod_Q)
    Wmod = (1.49*parg[igWfuse] + 2.42*(parg[igWhtail]+parg[igWvtail]) + 0.82*parg[igWwing] + 0.40*(parg[igWtesys])
    + 0.12*parg[igWMTO]*(parg[igflgnose]+parg[igflgmain]) + 1.59*(parg[igWMTO]*parg[igfhpesys]+parg[igWftank])) / 9.81
    # a modified OEW (in kg) that assigns more weight to systems with a higher cost per pound
    # based on methods using in Markish (2000), with componenent weights determined from data for a 777-200
    Vmax = maximum(pare[ieu0,:]) * 3.6 #maximum velocity (km/h)
    Npax = parm[imWpay]/215/lbf_to_N #number of pax
    P_hp = parpt[ipt_Ptshaft]/745.7 #turboshaft gas turbine engine power [hp]

    PPI = 242.8/144.8 #aerospace producer price index ratio from 1999 to 2020 (from US BLS)
    FTA = 6 #number of flight test aircraft (assume upper end due to novelty of propulsion system tech)

    RE = 86 #engineering wrap rate
    RT = 88 #tooling wrap rate
    RQ = 81 #quality control wrap rate
    RM = 73 #manufacturing wrap rate

    HE = 5.18*Wmod^0.777*Vmax^0.894*prod_Q^0.163 #engineering hours
    HT = 7.22*Wmod^0.777*Vmax^0.696*prod_Q^0.263 #tooling hours
    HM = 10.5*Wmod^0.82*Vmax^0.484*prod_Q^0.641 #manufacturing hours
    HQ = 0.133*HM #quality control hours

    CD = 48.7*Wmod^0.630*Vmax^1.3 #development support cost_est
    CF = 1408*Wmod^0.325*Vmax^0.822*FTA^1.21 #flight test cost
    CM = 22.6*Wmod^0.921*Vmax^0.621*prod_Q^0.799 #manufacturing materials cost (assumes Al construction)

    cA = 4/3 #adjustment factor, such that avionics costs are around 25% of flyaway cost
    CI = 2500*Npax #interior cost (e.g. seats, floors, etc.)
    # source: Raymer (2006), Collinson (2002)

    CTS = parpt[ipt_nTshaft]*3310*P_hp^0.7758 # source: Loh (2002)

    cost_airframe = PPI*(HE*RE + HT*RT + HM*RM + HQ*RQ + CD + CF + CM + CI)*cA
    cost_prop = CTS #plus other unmodeled components (motors, generators, fans, tanks)
    cost = cost_airframe + cost_prop
    # returns total program cost - divide by prod_Q to yield cost per aircraft

    return cost
end