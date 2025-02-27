"""
!!! warning "Unused and Unvetted"
    This legacy function is not used elsewhere in the code but has been
    retained for reference and in case we decide to update it in the future.
    Note that it has not been updated to work with the new `struct`s, that it
    has not been vetted, and is not endorsed by the current dev team. 

    Read the reference before using it!

Honestly, not a clue.
"""
function CostVal(prod_Q)
    # baseline cost of 737-MAX9 with CFM LEAP-1B engines, from TASOPT model provided by Prashanth
    wfuse = 20183.2
    whtail = 1202.0
    wvtail = 793.3
    wwing = 11813.0
    weng = 7995.9
    wlg = 1012.6
    wsys = 920.5
    Vmax = 855
    Npax = 220 
    Fmax = 130.4 
    Mmax = 1
    Tin = 1804 

    Wmod = 1.49*wfuse + 2.42*(whtail+wvtail) + 0.82*wwing + 0.40*weng + 0.12*wlg + 1.59*wsys

    PPI = 242.8/144.8 #aerospace producer price index ratio from 1999 to 2020 (from US BLS)
    PPI_07 = 242.8/186.8 #for 2007 (tank and pcec model)
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

    #CTS_h2 = 3560*P_hp^0.7758 # source: Loh (2002), includes inflation from 1998#
    CTS = 2251*(9.66*Fmax+243.25*Mmax+1.74*Tin-2228)*1.2*PPI #birkler's model for gas turbine engines as published in Raymer
    CPE = 0
    CEC = 0

    cost_airframe = PPI*(HE*RE + HT*RT + HM*RM + HQ*RQ + CD + CF + CM + CI)*cA / prod_Q
    cost_prop = 2*(CTS + CEC + CPE)
    #cost_prop = CTS*2 + CEC # x2 is rough estimate from Epstein & O'Flarity for turboelectric systems 
    cost_tank = 0
    cost = cost_prop + cost_airframe + cost_tank 
    # returns cost per aircraft

    dev_cost = (HE*RE + HT*RT + CD + CF)*cA*PPI/ prod_Q
    prod_cost = (HM*RM + HQ*RQ + CM + CI)*cA*PPI/ prod_Q
    prop_cost = 2*CTS

    return dev_cost, prod_cost, prop_cost
end