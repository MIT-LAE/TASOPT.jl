"""
    fuel_cell_with_ducted_fan_weight!(ac)

Computes and updates the weight and nacelle properties of the fuel-cell-powered ducted fan system. 
This function estimates the weight of the engine, nacelle, heat exchangers, and bare engine structure, 
while also updating the nacelle area fraction and nacelle length.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object containing propulsion and geometry parameters

    **Output:**
    No direct outputs. The `ac` object is modified with:
    - `parg[igWeng]`  → total engine weight
    - `parg[igWebare]` → bare engine weight
    - `parg[igWnace]` → nacelle weight
    - `parg[igWHXs]` → heat exchanger weight
    - `parg[igfSnace]` → nacelle area fraction relative to wing reference area
    - `parg[iglnace]` → nacelle length based on fan diameter and nacelle shape factor

!!! warning "⚠️ Work in Progress"
    - Fuel cell weight is not yet included.
    - Electric motor and associated machines are not yet included.
"""
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