"""
    place_cabin_seats(pax, cabin_width, seat_pitch = 30.0*in_to_m, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m)

Function to calculate the seat arrangement in the cabin and, therefore, the required cabin
length.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pax::Float64`: design number of passengers
    - `cabin_width::Float64`: width of cabin (m).
    - `seat_width::Float64`: width of one seat (m).
    - `aisle_halfwidth::Float64`: half the width of an aisle (m).
    - `fuse_offset::Float64`: distance from outside of fuselage to edge of closest window seat (m).

    **Outputs:**
    - `lcabin::Float64`: cabin length (m).
    - `xseats::Vector{Float64}`: longitudinal coordinate of each row of seats, measured from front of cabin (m).
    - `seats_per_row::Float64`: number of seats per row.
"""
function place_cabin_seats(pax, cabin_width, seat_pitch = 30.0*in_to_m, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m)

    cabin_offset = 10 * ft_to_m #Distance to the front and back of seats
    #TODO the hardcoded 10 ft is not elegant

    #Calculate the maximum number of seats per row
    seats_per_row = 1
    Dmin = seats_per_row*seat_width + 2*aisle_halfwidth + 2*fuse_offset #Minimum diameter for one passenger
    while (cabin_width > Dmin) || (cabin_width ≈ Dmin) #While the required diameter is smaller than the fuselage diameter
        seats_per_row = seats_per_row + 1 #Add one more seat
        layout = seat_layouts[seats_per_row] #Find corresponding seat layour from seat dict
        Dmin = seats_per_row*seat_width + (length(layout) - 1)*2*aisle_halfwidth + 2*fuse_offset #New minimum diameter
    end
    
    seats_per_row = seats_per_row - 1 #Subtract 1 seat to find maximum number of seats per row such that cabin_width > Dmin

    rows = Int(ceil(pax / seats_per_row))

    if seats_per_row <= 10
        emergency_rows = [12, 13]
    else
        emergency_rows = [19, 20]
    end

    xseats = zeros(rows)'
    xseats[1] = cabin_offset
    for r in 2:rows
        emergency_exit = 0.0
        if (r in emergency_rows)
            emergency_exit = seat_pitch/2
        end
        xseats[r] = xseats[r-1] + seat_pitch + emergency_exit
    end

    lcabin = xseats[end]
    return lcabin, xseats, seats_per_row
end # function place_cabin_seats

"""
    arrange_seats(seats_per_row, cabin_width,
     seat_width = 19.0 * in_to_m, 
     aisle_halfwidth = 10.0 * in_to_m,
     fuse_offset = 6.0*in_to_m)

Helper function to arrange seats given a number of `seats_per_row`
and cabin width. Assumes default `seat_width = 19"` and `aisle_halfwidth = 10"`,
but can be supplied by the user.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `seats_per_row::Float64`:  number of seats per row.
    - `cabin_width::Float64`: width of cabin (m).
    - `seat_width::Float64`: width of one seat (m).
    - `aisle_halfwidth::Float64`: half the width of an aisle (m).
    - `fuse_offset::Float64`: distance from outside of fuselage to edge of closest window seat (m).

    **Outputs:**
    - `yseats::Vector{Float64}`: transverse coordinate of each column of seats, measured from centerline (m).
"""
function arrange_seats(seats_per_row, cabin_width,
     seat_width = 19.0 * in_to_m, 
     aisle_halfwidth = 10.0 * in_to_m, 
     fuse_offset = 6.0*in_to_m)

    #Seats
    # Conditions:
    # - No more than 2 seats between any seat and the aisle

    layout = seat_layouts[seats_per_row] #find seat layout from dictionary
    n_aisles = length(layout) - 1 #Number of aisles
    Dmin = seats_per_row*seat_width + n_aisles*2*aisle_halfwidth + 2*fuse_offset #Minimum required diameter

    exp_aisle_halfwidth = aisle_halfwidth + (cabin_width - Dmin)/(2*n_aisles) #Expanded aisle to take up all available space

    yseats = zeros(seats_per_row)
    yseats[1] = fuse_offset + seat_width/2 #First seat is window seat
    for i = 2:seats_per_row
        flag_aisle = aisle_flag(i, layout) #1 if there is an aise to the left of seat, 0 if not
        yseats[i] = yseats[i - 1] + seat_width + flag_aisle*2*exp_aisle_halfwidth
    end

    #Shift seat coordinates to start in centerline
    yseats = yseats .- cabin_width/2

    return yseats
end  # function arrange_seats

"""
    aisle_flag(idx, layout)

Helper function to find if there is an aisle to the left of a given seat.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `idx::Int64`:  seat index.
    - `layout::Vector{Int64}`: seat layout map.

    **Outputs:**
    - `flag::Float64`: 1.0 if there is an aise to the left of seat, 0.0 if not.
"""
function aisle_flag(idx, layout)
    #Use cumulative sum to find total number of seats to the left of a given aisle.
    #If the difference between the cumsum and the index is exactly 1, the seat has an aisle to the left.
    if 1 in (idx .- cumsum(layout))
        flag = 1.0
    else
        flag = 0.0
    end
    return flag
end

"""
    find_cabin_width(Rfuse::Float64, wfb::Float64, nfweb::Float64, θ::Float64)

This function can be used to calculate the width of the passenger cabin from the double-bubble parameters
and the floor angular position.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Rfuse::Float64`: fuselage exterior radius (m)
    - `wfb::Float64`: lateral shift of double bubble (m)
    - `nfweb::Float64`: number of vertical webs in fuselage
    - `θ::Float64`: angle of floor wrt upper bubble center (rad)

    **Outputs:**
    - `w::Float64`: width of cabin (m).
"""
function find_cabin_width(Rfuse::Float64, wfb::Float64, nfweb::Float64, θ::Float64)
    #Use trigonometry to find cabin width
    w = nfweb*2*wfb + 2*Rfuse*cos(θ)
    return w
end

"""
    find_floor_angles(fdoubledecker::Bool, Rfuse::Float64, dRfuse::Float64; θ1::Float64 = 0.0, h_seat::Float64 = 0.0, d_floor::Float64 = 0.0)

This function can be used to place the passenger decks inside the fuselage. It works for single deck or double decker
cabins. It returns the angular position of each deck with respect to the center of the upper bubble.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Rfuse::Bool`: flag to indicate whether aircraft is a double decker
    - `Rfuse::Float64`: fuselage exterior radius (m)
    - `dRfuse::Float64`: vertical shift of downward bubble (m)
    - `θ1::Float64`: required in some cases; angle of main floor wrt upper bubble center (rad)
    - `h_seat::Float64`: required in some cases; seat height (m)
    - `d_floor::Float64`: required in some cases; distance between double decks (m)

    **Outputs:**
    - `θ1::Float64`: angle of main floor wrt upper bubble center (rad)
    - `θ2::Float64`: returned when double decker; angle of upper floor wrt upper bubble center (rad)
"""
function find_floor_angles(fdoubledecker::Bool, Rfuse::Float64, dRfuse::Float64; θ1::Float64 = 0.0, h_seat::Float64 = 0.0, d_floor::Float64 = 0.0)
    #If there is a lower bubble, the floor must be placed at the intersection of the two
    if dRfuse > 0.0
        θ1 = -asin(dRfuse / (2*Rfuse))
        if fdoubledecker
            θ2 = asin((h_seat + d_floor - dRfuse/2) / Rfuse)
            return θ1, θ2
        else
            return θ1
        end  
    else #If there is not a lower bubble
        if ~fdoubledecker #If it has a single deck
            θ1 = -asin(h_seat / (2*Rfuse)) #This angle maximizes the cabin width
            return θ1
        else #If it is a double decker with no lower bubble, the main cabin could be anywhere => Use provided angle
            θ2 = asin((h_seat + d_floor + Rfuse * sin(θ1)) / Rfuse)  
            return θ1, θ2
        end
    end

end

"""
    find_double_decker_cabin_length(x::Vector{Float64}, parg, parm)

This function can be used to calculate the length of a double decker cabin with different number of 
passengers on each deck.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x::Vector{Float64}`: vector with optimization variables
    - `parg::Vector{Float64}`: vector with aircraft geometric and mass parameters
    - `parm::Array{Float64}`: array with mission parameters

    **Outputs:**
    - `maxl::Float64`: required cabin length (m)
"""
function find_double_decker_cabin_length(x::Vector{Float64}, parg, parm)
    
    seat_pitch = parg[igseatpitch]
    seat_width = parg[igseatwidth]
    aisle_halfwidth = parg[igaislehalfwidth]
    h_seat = parg[igseatheight]

    Rfuse = parg[igRfuse]
    dRfuse = parg[igdRfuse]
    wfb = parg[igwfb]
    nfweb = parg[ignfweb]
    d_floor = parg[igfloordist]

    paxsize = parg[igWpaymax]/parm[imWperpax,1] #maximum number of passengers #TODO replace with better measure

    paxmain = x[1]
    paxtop = paxsize - paxmain

    try #The calculation could fail for some inputs if an asin returns an error
        if dRfuse ≈ 0.0 #If the cross-section is circular, the main deck could be anywhere

            θ1 = x[2] #Main deck angle
            _, θ2 = find_floor_angles(true, Rfuse, dRfuse, θ1 = θ1, h_seat= h_seat, d_floor = d_floor)
        else #Main deck angle depends on double bubble params
            θ1, θ2 = find_floor_angles(true, Rfuse, dRfuse, h_seat = h_seat, d_floor = d_floor)
        end

        #Find width of each cabin
        w1 = find_cabin_width(Rfuse, wfb, nfweb, θ1)
        w2 = find_cabin_width(Rfuse, wfb, nfweb, θ2)

        #Find length of each cabin
        l1, _, _ = place_cabin_seats(paxmain, w1, seat_pitch, seat_width, aisle_halfwidth)
        l2, _, _ = place_cabin_seats(paxtop, w2, seat_pitch, seat_width, aisle_halfwidth)

        maxl = max(l1, l2) #Required length
        
        return maxl 
    catch
        return 1e6
    end
end

"""
    find_double_decker_cabin_length(x::Vector{Float64}, parg, parm)

This function can be used to optimize the passenger distribution across two decks in a double decker aircraft. 
If the cross-section is circular, it also optimizes the deck layouts.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `parg::Vector{Float64}`: vector with aircraft geometric and mass parameters
    - `parm::Array{Float64}`: array with mission parameters

    **Outputs:**
    - `xopt::Vector{Float64}`: vector with optimization results
"""
function optimize_double_decker_cabin(parg, parm)
    dRfuse = parg[igdRfuse]

    paxsize = parg[igWpaymax]/parm[imWperpax,1] #maximum number of passengers #TODO replace with better measure

    if dRfuse ≈ 0.0 #If the cross-section is circular, the main deck could be anywhere so optimize it
        initial_x = [paxsize/2, -pi/4]
        lower = [1.0, -pi/2]
        upper = [paxsize - 1.0, pi/2]
        
    else #Only optimize PAX distribution
        initial_x = [paxsize/2]
        lower = [1.0]
        upper = [paxsize - 1.0]
    end

    obj(x, grad) = find_double_decker_cabin_length(x, parg, parm) #Objective function is cabin length

    opt = Opt(:GN_DIRECT_L, length(initial_x)) #Use a global optimizer as it is only 1 or 2 variables
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000  # Set the maximum number of function evaluations

    opt.min_objective = obj
    
    (minf,xopt,ret) = NLopt.optimize(opt, initial_x) #Solve optimization problem

    return xopt
end