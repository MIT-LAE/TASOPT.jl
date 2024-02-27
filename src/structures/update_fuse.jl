function update_fuse!(pari, parg)

    nftanks = pari[iinftanks] #Number of fuel tanks in fuselage
    # parg[igRfuse   ] = 90 * in_to_m 
    lftank = parg[iglftank] # Get length of tank from previous iteration
    lftoffset = 2.0*ft_to_m #1 ft buffer for front and back of tanks

    #Useful relative distances to conserve
    lcyl = parg[igdxcabin]
    dxeng2wbox = parg[igdxeng2wbox]
    dxapu2end = parg[igxend] - parg[igxapu]
    dxshell2conend = parg[igxconend ] - parg[igxshell2 ]
    dxshell2apu = parg[igxapu ] - parg[igxshell2 ]
    dxhbox2conend = parg[igxconend] - parg[igxhbox ]
    dxvbox2conend = parg[igxconend] - parg[igxvbox ]

    if parg[igxftankaft] == 0.0 #if there is not a rear tank
        dxcyl2shellaft = parg[igdxcyl2shellaft]
    else #if there is a rear tank
        dxcyl2shellaft = 0.0 #no need for offset between shell2 and blend2 since rear space cannot be used
    end

    #Update positions and fuselage length
    parg[igxblend2] = parg[igxblend1] + nftanks * (lftank + lftoffset) + lcyl
    
    parg[igxshell2 ] = parg[igxblend2] + dxcyl2shellaft

    parg[igxconend ] = parg[igxshell2] + dxshell2conend
    parg[igxapu    ] = parg[igxshell2] + dxshell2apu
    parg[igxend    ] = parg[igxapu] + dxapu2end
    parg[igxhpesys] = parg[igxconend] * 0.52484 #TODO: address this
    
    parg[igxhbox   ] = parg[igxconend ] - dxhbox2conend
    parg[igxvbox   ] = parg[igxconend ] - dxvbox2conend
    
    parg[igxeng    ] =  parg[igxwbox] - dxeng2wbox

end