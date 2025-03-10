function fuel_cell_with_ducted_fan_weight!(ac)
    parg = ac.parg
    wing = ac.wing
    neng = parg[igneng]
    
    Weng, Wnace, Webare, W_HXs, Snace1 = ductedfanweight(ac)

    parg[igWeng] = Weng
    parg[igWebare] = Webare
    parg[igWnace] = Wnace
    parg[igWHXs] = W_HXs

    # set new nacelle area / reference area  fraction fSnace
    S = wing.layout.S

    Snace = Snace1 * neng
    fSnace = Snace / S
    parg[igfSnace] = fSnace

    # set new nacelle area / reference area  fraction fSnace
    Snace = Snace1 * neng
    fSnace = Snace / S
    parg[igfSnace] = fSnace
    lnace = parg[igdfan] * parg[igrSnace] * 0.15
    parg[iglnace] = lnace

    #TODO add weight of fuel cell
    #TODO add weight of electric motor and machines
end