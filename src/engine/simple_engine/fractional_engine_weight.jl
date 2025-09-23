"""
    fractional_engine_weight!(ac)

Simple function to estimate the total weight of an aircraft engine based on input
ratios of engine weight to MTOW.
      
!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine weights.
"""
function fractional_engine_weight!(ac)
    parg = ac.parg
    parg[igWeng] = parg[igWMTO] * parg[igfeng] #Engine weight is MTOW times engine weight fraction

    #Simple estimate of nacelle length, following orig. Drela assumptions
    lnace = parg[igdfan] * parg[igrSnace] * 0.15
    parg[iglnace] = lnace
end