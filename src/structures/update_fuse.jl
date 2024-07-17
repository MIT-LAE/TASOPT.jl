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
    Parameters in `parg` are modified. It also outputs:
    - `seats_per_row::Float64`: number of seats per row in main cabin (lower deck if double decker)
"""
function update_fuse_for_pax!(pari, parg, parm, fuse_tank)

    seat_pitch = parg[igseatpitch]
    seat_width = parg[igseatwidth]
    aisle_halfwidth = parg[igaislehalfwidth]
    h_seat = parg[igseatheight]

    Rfuse = parg[igRfuse]
    dRfuse = parg[igdRfuse]
    wfb = parg[igwfb]
    nfweb = parg[ignfweb]

    #Find cabin length by placing seats
    if pari[iidoubledeck] == 1 #if the aircraft is a double decker
        xopt, seats_per_row = optimize_double_decker_cabin(parg, parm) #Optimize the floor layout and passenger distributions

        lcyl, _ = find_double_decker_cabin_length(xopt, parg, parm) #Total length is maximum of the two

    else
        θ = find_floor_angles(false, Rfuse, dRfuse, h_seat = h_seat) #Find the floor angle
        paxsize = parg[igWpaymax]/parm[imWperpax,1] #maximum number of passengers
        w = find_cabin_width(Rfuse, wfb, nfweb, θ, h_seat) #Cabin width
        lcyl, _, seats_per_row = place_cabin_seats(paxsize, w, seat_pitch, seat_width, aisle_halfwidth) #Cabin length
    end

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

    return seats_per_row
end

"""
    find_minimum_radius_for_seats_per_row(seats_per_row, ac_base)

This function calculates the minimum radius required to have a desired number of seats per row in the main cabin.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `seats_per_row::Float64`: number of seats per row in main cabin (lower deck if double decker)
    - `ac_base::aircraft`: aircraft object

    **Outputs:**
    - `R::Float64`: minimum radius for desired number of seats per row (m)
"""
function find_minimum_radius_for_seats_per_row(seats_per_row, ac_base)
    ac = deepcopy(ac_base) #Copy input ac to avoid modifying it
    obj(x, grad) = x[1] + 1e3 * abs(check_seats_per_row_diff(seats_per_row, x, ac))  #Objective function is the radius plus a big penalty if constraint is not met

    initial_x = [4.0]
    opt = Opt(:GN_DIRECT, length(initial_x)) #Use a global optimizer that can handle equality constraints
    opt.lower_bounds = [0.0]
    opt.upper_bounds = [5.0]

    # opt_local = Opt(:GN_DIRECT, length(initial_x))
    opt.maxeval = 5000  # Set the max number of evaluations
    # opt.local_optimizer = opt_local

    opt.min_objective = obj

    #Apply the equality constraint that ensures that the number of seats per row is the desired one
    #equality_constraint!(opt, (x, grad) -> check_seats_per_row_diff(seats_per_row, x, ac), 1e-5) 
    
    (minf,xopt,ret) = NLopt.optimize(opt, initial_x) #Solve optimization problem
    R = xopt[1]
    return R
end

"""
    check_seats_per_row_diff(seats_per_row, x, ac)

This function returns the difference between the desired number of seats per row and the one corresponding to a 
given radius

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `seats_per_row::Float64`: desired number of seats per row in main cabin (lower deck if double decker)
    - `x::Vector{Float64}`: vector with one entry containing the fuselage radius (m)
    - `ac::aircraft`: aircraft object

    **Outputs:**
    - `diff::Float64`: difference between desired number of seats per row and that for the input radius
"""
function check_seats_per_row_diff(seats_per_row, x, ac)
    Rfuse = x[1]
    ac.parg[igRfuse] = Rfuse
    try #Sometimes update_fuse_for_pax may fail
        seats_per_row_rad = update_fuse_for_pax!(ac.pari, ac.parg, ac.parm, ac.fuse_tank)
        diff = seats_per_row_rad - seats_per_row
        #println("R = $Rfuse, s = $seats_per_row_rad")
        return diff
    catch
        #println("failed")
        return 1.0
    end
end