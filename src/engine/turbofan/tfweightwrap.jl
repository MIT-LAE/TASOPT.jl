"""
    tfweightwrap!(ac, HXs)

General function to estimate and store the weight of a turbofan engine.
      
!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object
    - `HXs`: vector with heat exchanger performance data

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine weights.
"""
function tfweightwrap!(ac, HXs)
    pari = ac.pari
    parg = ac.parg
    wing = ac.wing
    neng = parg[igneng]
    
    if pari[iiengtype] == 1 #turbofan TODO: replace with better flag
        Weng, Wnace, Webare, W_HXs, Snace1 = tfweight(ac, HXs)
        
    end
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
end