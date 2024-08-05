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
function update_fuse!(fuselage, pari, parg)

    nftanks = pari[iinftanks] #Number of fuel tanks in fuselage
    # parg[igRfuse   ] = 90 * in_to_m 
    lftank = parg[iglftank] # Get length of tank from previous iteration
    lftoffset = 2.0*ft_to_m #1 ft buffer for front and back of tanks

    #Useful relative distances to conserve
    lcyl = parg[igdxcabin]
    dxeng2wbox = parg[igdxeng2wbox]
    dxapu2end = fuselage.layout.x_end - fuselage.APU.x
    dxshell2conend =fuselage.layout.x_cone_end - fuselage.layout.x_pressure_shell_aft
    dxshell2apu = fuselage.APU.x - fuselage.layout.x_pressure_shell_aft
    dxhbox2conend = fuselage.layout.x_cone_end - parg[igxhbox ]
    dxvbox2conend = fuselage.layout.x_cone_end - parg[igxvbox ]

    if parg[igxftankaft] == 0.0 #if there is not a rear tank
        dxcyl2shellaft = fuselage.layout.x_pressure_shell_aft - fuselage.layout.x_end_cylinder
    else #if there is a rear tank
        dxcyl2shellaft = 0.0 #no need for offset between shell2 and blend2 since rear space cannot be used
    end

    #Update positions and fuselage length
    fuselage.layout.x_end_cylinder = fuselage.layout.x_start_cylinder + nftanks * (lftank + lftoffset) + lcyl
    
    fuselage.layout.x_pressure_shell_aft = fuselage.layout.x_end_cylinder + dxcyl2shellaft

    fuselage.layout.x_cone_end = fuselage.layout.x_pressure_shell_aft + dxshell2conend
    fuselage.APU.r = [fuselage.layout.x_pressure_shell_aft + dxshell2apu, 0.0, 0.0]
    fuselage.layout.x_end = fuselage.APU.x + dxapu2end
    fuselage.HPE_sys.r = [fuselage.layout.x_cone_end * 0.52484, 0.0,0.0]#TODO: address this
    
    parg[igxhbox   ] = fuselage.layout.x_cone_end - dxhbox2conend
    parg[igxvbox   ] = fuselage.layout.x_cone_end - dxvbox2conend
    
    parg[igxeng    ] =  wing.layout.box_x - dxeng2wbox

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
function update_fuse_for_pax!(pari, parg, parm, fuse, fuse_tank)

    seat_pitch = parg[igseatpitch]
    seat_width = parg[igseatwidth]
    aisle_halfwidth = parg[igaislehalfwidth]

    if fuse.n_decks == 2 #if the aircraft is a double decker
        #passenger count to size cabin is half of the maximum
        paxsize = ceil(parg[igWpaymax]/parm[imWperpax,1] / 2) 
    else
        paxsize = parg[igWpaymax]/parm[imWperpax,1] #maximum number of passengers
    end
    #TODO this double deck model assumes that both decks have a width equal to the fuselage diameter; 
    #in reality, at least one deck must be narrower

    #Useful relative distances to conserve
    dxeng2wbox = parg[igdxeng2wbox] #Distance from engine to wingbox
    dxcyl2shellaft = fuse.layout.x_pressure_shell_aft - fuse.layout.x_end_cylinder #Distance from blend2 to shell2
    dxapu2end = fuse.layout.x_end - fuse.APU.x #Distance from APU to end
    dxshell2conend = fuse.layout.x_cone_end - fuse.layout.x_pressure_shell_aft #Distance from shell2 to conend
    dxshell2apu = fuse.APU.x - fuse.layout.x_pressure_shell_aft #Distance from shell2 to APU
    dxhbox2conend = fuse.layout.x_cone_end - parg[igxhbox ] #Distance from conend to xhbox
    dxvbox2conend = fuse.layout.x_cone_end - parg[igxvbox ] #Distance from conend to xvbox
    #Fraction of cabin length at which wing is located
    wbox_cabin_frac =  (wing.layout.box_x- fuse.layout.x_start_cylinder )/(fuse.layout.x_end_cylinder - fuse.layout.x_start_cylinder) 

    #Find new cabin length
    wcabin = find_cabin_width(fuse.layout.radius, fuse.layout.bubble_lower_downward_shift, fuse.layout.bubble_center_y_offset, fuse.layout.n_webs, parg[igfloordist]) #Find cabin width
    lcyl, _, _ = place_cabin_seats(paxsize, wcabin, seat_pitch, seat_width, aisle_halfwidth) #Size for max pax count

    #When there is a fuel tank at the back of the fuselage, there is no offset between the end of the seat rows
    #and the start of the tank. For this reason, leave a 5ft offset at back
    if (pari[iifwing]  == 0) && ((fuse_tank.placement == "rear") || (fuse_tank.placement == "both"))
        lcyl = lcyl + 5.0 * ft_to_m #Make cabin longer to leave room in the back
        #TODO the hardcoded 5 ft is not elegant
    end

    #Update positions and fuselage length
    fuselage.layout.x_end_cylinder = fuselage.layout.x_start_cylinder + lcyl

    #Update wingbox position
    wing.layout.box_x = fuselage.layout.x_start_cylinder + wbox_cabin_frac * lcyl
       
    #Update other lengths
    fuselage.layout.x_pressure_shell_aft = fuselage.layout.x_end_cylinder + dxcyl2shellaft

    fuselage.layout.x_cone_end = fuselage.layout.x_pressure_shell_aft + dxshell2conend
    fuselage.APU.r = [fuselage.layout.x_pressure_shell_aft + dxshell2apu, 0.0,0.0]
    fuselage.layout.x_end = fuselage.APU.x + dxapu2end
    fuselage.HPE_sys.r = [fuselage.layout.x_cone_end * 0.52484, 0.0, 0.0] #TODO: address this
    
    parg[igxhbox   ] = fuselage.layout.x_cone_end - dxhbox2conend
    parg[igxvbox   ] = fuselage.layout.x_cone_end - dxvbox2conend
    
    parg[igxeng    ] =  wing.layout.box_x - dxeng2wbox #Move engine

    parg[igdxcabin] = lcyl #Store new cabin length
end

# [fuselage.layout.x_end_cylinder, parg[igxwbox], fuselage.layout.x_pressure_shell_aft, fuselage.layout.x_cone_end,
# fuselage.APU.x, fuselage.layout.x_end, fuselage.HPE_sys.x, parg[igxhbox], parg[igxvbox],parg[igxeng],parg[igdxcabin]]