function update_fuse!(pari, parg)

    nftanks = pari[iinftanks] #Number of fuel tanks in fuselage
    # parg[igRfuse   ] = 90 * in_to_m 
    lftank = parg[iglftank] # Get length of tank from previous iteration
    lftoffset = 2.0*ft_to_m #1 ft buffer for front and back of tanks

    #Useful relative distances to conserve
    lcyl = parg[iglblend2blend]
    leng2wbox = parg[igxwbox] - parg[igxeng]
    lcyl2shell_aft = parg[igxshell2 ] - parg[igxblend2]
    lapu2end = parg[igxend] - parg[igxapu]
    lshell2conend = parg[igxconend ] - parg[igxshell2 ]
    lshell2apu = parg[igxapu ] - parg[igxshell2 ]
    lhbox2conend = parg[igxconend] - parg[igxhbox ]
    lvbox2conend = parg[igxconend] - parg[igxvbox ]


    #Update positions and fuselage length
    parg[igxblend2] = parg[igxblend1] + nftanks * (lftank + lftoffset) + lcyl
       
    parg[igxshell2 ] = parg[igxblend2] + lcyl2shell_aft

    # parg[igxftankfront ]  = parg[igxblend1] + 1.0*ft_to_m + ltank/2 #(buffer)
    # parg[igxftankaft ]  = parg[igxblend1] + 1.0*ft_to_m + ltank + 1.0*ft_to_m + lcabin + 1.0*ft_to_m + ltank/2.0

    parg[igxconend ] = parg[igxshell2] + lshell2conend
    parg[igxapu    ] = parg[igxshell2] + lshell2apu
    parg[igxend    ] = parg[igxapu] + lapu2end
    parg[igxhpesys] = parg[igxconend] * 0.52484 #TODO: address this
    
    parg[igxhbox   ] = parg[igxconend ] - lhbox2conend
    parg[igxvbox   ] = parg[igxconend ] - lvbox2conend
    parg[igxeng    ] =  parg[igxwbox] - leng2wbox


    
end