function engineweight!(ac, engine_type, HXs)
    parg = ac.parg
    neng = parg[igneng]
    if engine_type == "turbofan"

        Weng, Wnace, Webare, Snace1 = tfweight(ac, HXs)
        
        parg[igWeng] = Weng
        parg[igWebare] = Webare
        parg[igWnace] = Wnace
    
        # set new nacelle area / reference area  fraction fSnace
        S = parg[igS]
    
        Snace = Snace1 * neng
        fSnace = Snace / S
        parg[igfSnace] = fSnace
    elseif engine_type == "ducted_fan"
    end
end