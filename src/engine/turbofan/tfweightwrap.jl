"""
    tfweightwrap!(ac)

General function to estimate and store the weight of a turbofan engine. 
This function is basically a wrapper on tfweight, going from the
basic aircraft inputs to those required by the function and storing the outputs.
      
!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine weights.
"""
function tfweightwrap!(ac)
    parg = ac.parg
    wing = ac.wing
    neng = parg[igneng]
    
    Weng, Wnace, Webare, W_HXs, Snace1 = tfweight(ac)
        
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