#TODO: make this work with new structs or delete
"""
cost_est estimates the program and airframe cost for the aircraft
aircraft cost model developed from DAPCA IV, as found in Raymer's
Aircraft Design (2006)

!!! warning "Unused and Unvetted"
    This legacy function is not used elsewhere in the code but has been
    retained for reference and in case we decide to update it in the future.
    Note that (i) it has not been updated to work with the new `struct`s, (ii) it refers to 
    turboelectric aircraft parameters that have been removed from the model, (iii) it
    has not been vetted recently, and (iv) it is not endorsed by the current dev team. 

    Read the references before adapting it!

Inputs:
parg[,] geometric parameters
pare[,] engine parameters
prod_Q  estimated production quantity over 5 years

Outputs:
cost    total program cost estimate [all \$ in USD 2020]

"""
function CostEst(pari, parg, pare, parm, parpt, prod_Q)

    @warn "This legacy function is not used elsewhere in the code but has been
    retained for reference and in case we decide to update it in the future.
    Note that (i) it has not been updated to work with the new `struct`s, (ii) it refers to 
    turboelectric aircraft parameters that have been removed from the model, (iii) it
    has not been vetted recently, and (iv) it is not endorsed by the current dev team. 

    Read the references before adapting it!"

    conv_dev, conv_prod, conv_prop = CostVal(prod_Q) #737 MAX9 baseline - no longer used

    # Wmod = (1.28*parg[igWfuse] + 2.09*(parg[igWhtail]+parg[igWvtail]) + 0.71*parg[igWwing]
    # + 0.10*parg[igWMTO]*(parg[igflgnose]+parg[igflgmain]) + 1.37*(parg[igWMTO]*parg[igfhpesys])) / 9.81

    Wmod = (1.19*parg[igWfuse] + 2.62*(parg[igWhtail]+parg[igWvtail]) + 0.89*parg[igWwing]
    + 0.13*parg[igWMTO]*(parg[igflgnose]+parg[igflgmain]) + 1.72*(parg[igWMTO]*parg[igfhpesys])
    + 0.44*parg[igWtesys]) / 9.81
    println("Wmod frac = ", (parg[igWfuse] + (parg[igWhtail]+parg[igWvtail]) + parg[igWwing]
    + parg[igWMTO]*(parg[igflgnose]+parg[igflgmain]) + (parg[igWMTO]*parg[igfhpesys])
    + parg[igWtesys]) / (parg[igWMTO] - parg[igWfuel] - parg[igWpay]))

    # Wmod = (parg[igWMTO] - parg[igWfuel] - parg[igWpay])/9.81
    # Wmod = parg[igWMTO]/9.81
    # a modified OEW (in kg) that assigns more weight to systems with a higher cost per pound
    # based on methods using in Markish (2000), with componenent weights determined from data for a 777-200
    #println(Wmod)
    Vmax = maximum(pare[ieu0,:]) * 3.6 #maximum velocity (km/h)
    #println(Vmax)
    Npax =  parm[imWpay]/215/lbf_to_N #number of pax
    P_hp = parpt[ipt_Ptshaft]/745.7 #turboshaft gas turbine engine power [hp]
    P_MMBTU = parpt[ipt_Ptshaft]/293e3 #turboshaft gas turbine engine power [hp]
    P_kW = parpt[ipt_Ptshaft]/1000 #in kW

    # T3 = parpt[ipt_Tt3] * 5/9 #K
    # P3 = parpt[ipt_pt3] * 6.895 #kPa
    # a = 6.25528852e-08
    # b = -1.17064467e-04
    # c = 7.36953400e-02
    # d = -1.50392850e+01
    # EINOx = P3^0.4*(a*T3^3 + b*T3^2 + c*T3 +d) #for RQL combustor on CFM56
    # mdotf_max = maximum(pare[iemdotf,:])
    # mdot_nox = EINOx*mdotf_max
    # life_nox = mdot_nox*8000*3.6 #in tons, using lifetime value in Sorrels 2014

    PPI = 242.8/144.8 #aerospace producer price index ratio from 1999 to 2020 (from US BLS)
    PPI_07 = 242.8/186.8 #for 2007 (tank and pcec model)
    FTA = 6 #number of flight test aircraft (assume upper end due to novelty of propulsion system tech)
    dev_adj = 1.15 #accounts for added technical complexity compared to standard aircraft
    if (pari[iifuel] == 1)
        h2_adj = 1.15 # accounts for added cost of H2 compatibile components (e.g. fuel pumps) + aerospace materials
    else
        h2_adj = 1.0
    end
    # 10% for aerospace compatibility, +5% for H2 compatibility
    aero_adj = 1.10 #accounts for added cost for aerospace-compatible (i.e. lightweight) materials

    RE = 86 # engineering wrap rate [$/hr]
    RT = 88 # tooling wrap rate
    RQ = 81 # quality control wrap rate
    RM = 73 # manufacturing wrap rate

    HE = 5.18*Wmod^0.777*Vmax^0.894*prod_Q^0.163 #engineering hours
    HT = 7.22*Wmod^0.777*Vmax^0.696*prod_Q^0.263 #tooling hours
    HM = 10.5*Wmod^0.82*Vmax^0.484*prod_Q^0.641 #manufacturing hours
    HQ = 0.133*HM #quality control hours

    CD = 48.7*Wmod^0.630*Vmax^1.3 #development support cost_est [$]
    CF = 1408*Wmod^0.325*Vmax^0.822*FTA^1.21 #flight test cost
    CM = 22.6*Wmod^0.921*Vmax^0.621*prod_Q^0.799 #manufacturing materials cost (assumes Al construction)

    cA = 4/3 #adjustment factor, such that avionics costs are around 25% of flyaway cost
    CI = 2500*Npax #interior cost (e.g. seats, floors, etc.)
    # source: Raymer (2006), Collinson (2002)

    CTS = 3560*P_hp^0.7758*h2_adj # source: Loh (2002), includes inflation from 1998
    CPE = 273*P_kW # assumes 150 EUR/kW for motors, 75 EUR/kW for power electronics from Hoelzen (2018)
    CEC = PPI_07*187647*P_MMBTU^0.54*aero_adj # fit from Sorrels (2014)

    cost_airframe = PPI*(HE*RE + HT*RT + HM*RM + HQ*RQ + CD + CF + CM + CI)*cA / prod_Q
    cost_prop = parpt[ipt_nTshaft]*(CTS + CEC + CPE) * (1 + FTA/prod_Q)
    #cost_prop = CTS*2 + CEC # x2 is rough estimate from Epstein & O'Flarity for turboelectric systems 
    cost_tank = PPI_07*800e3*parg[igWftank]/31870*17.094*prod_Q^(-0.387) * (1 + FTA/prod_Q) 
    # scaled from 4000kg tank (as calculated by tankWmech) for tanker trucks, cost taken from Yang(2008)
    # prod_Q factor calculated as aggregate of functions used in Raymer
    cost = cost_prop + cost_airframe + cost_tank 
    # returns cost per aircraft

    cost_dev  = (HE*RE + HT*RT + CD + CF)*cA*PPI/ prod_Q * dev_adj
    cost_prod = (HM*RM + HQ*RQ + CM + CI)*cA*PPI/ prod_Q
    cost_ts   = parpt[ipt_nTshaft]*CTS * (prod_Q + FTA)/prod_Q
    cost_pe   = parpt[ipt_nTshaft]*CPE * (prod_Q + FTA)/prod_Q
    cost_pcec = parpt[ipt_nTshaft]*CEC * (prod_Q + FTA)/prod_Q

    # delta_dev = cost_dev - conv_dev
    # delta_prod = cost_prod + cost_tank - conv_prod #should tank cost be included in airframe production or propulsion system?
    # delta_prop = cost_ts + cost_pe + cost_pcec - conv_prop
    # delta_total = delta_dev + delta_prod + delta_prop


    println(cost_dev)
    println(cost_prod + cost_tank)
    println(cost_ts + cost_pe + cost_pcec)

    println("Wmod = $Wmod")
    println("Dev cost = ", cost_dev) #development cost
    println("Airframe manufac. cost = ",cost_prod) #airframe manufacturing cost
    println("CTS  = ",parpt[ipt_nTshaft]*CTS)
    println("CPE  = ",parpt[ipt_nTshaft]*CPE)
    println("CEC  = ",parpt[ipt_nTshaft]*CEC)
    println("Tank = ",cost_tank)
    println("h2 adj = ", h2_adj)
    println("---------------")
    println("Total cost millions of USD 2020: ",round(cost/1e6; digits=2) )
    return cost
end

