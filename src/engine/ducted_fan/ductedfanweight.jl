function ductedfanweight(ac, HXs)
    Dfan   = ac.parg[igdfan]
    Nmech  = max(ac.pared[ieNf, :])
    fpylon = ac.parg[igfpylon]
    neng = ac.parg[igneng]
    rSnace = ac.parg[igrSnace]

    ARfan  = 3   #Blade aspeect ratio
    bladesolidity = 0.4 # Blade solidity c/s
    ktech = 0.5 
    Utip  = Dfan/2* (2 * pi * Nmech/60);
    # Sagerser 1971, NASA TM X-2406
    # Note: The term "weight" in Sagerser1971 is actually mass
    mfan = ktech*(135.0 * Dfan^2.7/sqrt(ARfan) * (bladesolidity/1.25)^0.3 * (Utip/350.0)^0.3)

    Snace1 = rSnace * 0.25 * pi * Dfan^2
    Ainlet = 0.4*Snace1
    Acowl  = 0.2*Snace1
    Aexh   = 0.4*Snace1

    Wnace = 4.45*(Ainlet/0.3048^2.0) * (2.5+0.0238*Dfan/0.0254) +
            4.45*(Acowl /0.3048^2.0) *  1.9 +
            4.45*(Aexh  /0.3048^2.0) * (2.5+0.0363*Dfan/0.0254)

    Wfan = (mfan*9.81 + Wnace*0.8)*(1+fpylon)
    
    Weng = (Wfan + Wnace) * neng

    Webare = Wfan * neng
    Wnac = Wnace * neng
    return Weng, Wnac, Webare, Snace1
end