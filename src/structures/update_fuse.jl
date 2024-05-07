export update_fuse!, update_fuse_for_pax!
"""
    update_fuse!(pari, parg)

Function to update the fuselage layout when there is a change in fuselage fuel tank length.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pari::Vector{Int64}`: vector with aircraft integer parameters
    - `parg::Vector{Float64}`: vector with aircraft geometric and mass parameters

    **Outputs:**
    No direct outputs; parameters in `parg` are modified.
"""
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
        dxcyl2shellaft = parg[igxshell2] - parg[igxblend2]
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

"""
    update_fuse_for_pax!(pari, parg, parm, fuse_tank)

Function to update the fuselage layout when the cabin length is not known a priori, for example if the radius is changed. 
It sizes the cabin for the design number of passengers.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pari::Vector{Int64}`: vector with aircraft integer parameters
    - `parg::Vector{Float64}`: vector with aircraft geometric and mass parameters
    - `parm::Array{Float64}`: array with mission parameters
    - `fuse_tank::struct`: structure of type `fuselage_tank` with cryogenic fuel tank parameters

    **Outputs:**
    No direct outputs; parameters in `parg` are modified.
"""
function update_fuse_for_pax!(pari, parg, parm, fuse_tank)

    seat_pitch = parg[igseatpitch]
    seat_width = parg[igseatwidth]
    aisle_halfwidth = parg[igaislehalfwidth]

    if pari[iidoubledeck] == 1 #if the aircraft is a double decker
        #passenger count to size cabin is half of the maximum
        paxsize = ceil(parg[igWpaymax]/parm[imWperpax,1] / 2) 
    else
        paxsize = parg[igWpaymax]/parm[imWperpax,1] #maximum number of passengers
    end
    #TODO this double deck model assumes that both decks have a width equal to the fuselage diameter; 
    #in reality, at least one deck must be narrower

    #Useful relative distances to conserve
    dxeng2wbox = parg[igdxeng2wbox] #Distance from engine to wingbox
    dxcyl2shellaft = parg[igxshell2] - parg[igxblend2] #Distance from blend2 to shell2
    dxapu2end = parg[igxend] - parg[igxapu] #Distance from APU to end
    dxshell2conend = parg[igxconend ] - parg[igxshell2 ] #Distance from shell2 to conend
    dxshell2apu = parg[igxapu ] - parg[igxshell2 ] #Distance from shell2 to APU
    dxhbox2conend = parg[igxconend] - parg[igxhbox ] #Distance from conend to xhbox
    dxvbox2conend = parg[igxconend] - parg[igxvbox ] #Distance from conend to xvbox
    #Fraction of cabin length at which wing is located
    wbox_cabin_frac =  (parg[igxwbox]- parg[igxblend1] )/(parg[igxblend2] - parg[igxblend1]) 

    #Find new cabin length
    lcyl, _, _ = place_cabin_seats(paxsize, parg[igRfuse], seat_pitch, seat_width, aisle_halfwidth) #Size for max pax count

    #When there is a fuel tank at the back of the fuselage, there is no offset between the end of the seat rows
    #and the start of the tank. For this reason, leave a 5ft offset at back
    if (pari[iifwing]  == 0) && ((fuse_tank.placement == "rear") || (fuse_tank.placement == "both"))
        lcyl = lcyl + 5.0 * ft_to_m #Make cabin longer to leave room in the back
        #TODO the hardcoded 5 ft is not elegant
    end

    #Update positions and fuselage length
    parg[igxblend2] = parg[igxblend1] + lcyl

    #Update wingbox position
    parg[igxwbox] = parg[igxblend1] + wbox_cabin_frac * lcyl
       
    #Update other lengths
    parg[igxshell2 ] = parg[igxblend2] + dxcyl2shellaft

    parg[igxconend ] = parg[igxshell2] + dxshell2conend
    parg[igxapu    ] = parg[igxshell2] + dxshell2apu
    parg[igxend    ] = parg[igxapu] + dxapu2end
    parg[igxhpesys] = parg[igxconend] * 0.52484 #TODO: address this
    
    parg[igxhbox   ] = parg[igxconend ] - dxhbox2conend
    parg[igxvbox   ] = parg[igxconend ] - dxvbox2conend
    
    parg[igxeng    ] =  parg[igxwbox] - dxeng2wbox #Move engine

    parg[igdxcabin] = lcyl #Store new cabin length
end