"""
    place_cabin_seats(pax, cabin_width; seat_pitch = 30.0*in_to_m, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m,
    front_seat_offset = 10 * ft_to_m)	

Function to calculate the seat arrangement in the cabin and, therefore, the required cabin
length.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pax::Float64`: design number of passengers
    - `cabin_width::Float64`: width of cabin (m).
    - `seat_pitch::Float64`: longitudinal distance between seats (m).
    - `seat_width::Float64`: width of one seat (m).
    - `aisle_halfwidth::Float64`: half the width of an aisle (m).
    - `fuse_offset::Float64`: distance from outside of fuselage to edge of closest window seat (m).
    - `front_seat_offset::Float64`: distance from front of cabin to first row of seats (m).

    **Outputs:**
    - `lcabin::Float64`: cabin length (m).
    - `xseats::Vector{Float64}`: longitudinal coordinate of each row of seats, measured from front of cabin (m).
    - `seats_per_row::Float64`: number of seats per row.
"""
function place_cabin_seats(pax, cabin_width; seat_pitch = 30.0*in_to_m, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m,
    front_seat_offset = 10 * ft_to_m)	

    seats_per_row = findSeatsAbreast(cabin_width, seat_width, aisle_halfwidth, fuse_offset)

    rows = Int(ceil(pax / seats_per_row))

    if seats_per_row <= 10
        emergency_rows = [12, 13]
    else
        emergency_rows = [19, 20]
    end

    xseats = zeros(rows)'
    xseats[1] = front_seat_offset
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
    findSeatsAbreast(cabin_width, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m)

Function to find the number of seats abreast that can fit in a given cabin width.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `cabin_width::Float64`: width of cabin (m).
    - `seat_width::Float64`: width of one seat (m).
    - `aisle_halfwidth::Float64`: half the width of an aisle (m).
    - `fuse_offset::Float64`: distance from outside of fuselage to edge of closest window seat (m).

    **Outputs:**
    - `seats_per_row::Float64`: number of seats per row.
"""
function findSeatsAbreast(cabin_width, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m)

    #Calculate the maximum number of seats per row
    seats_per_row = 1
    Dmin = seats_per_row*seat_width + 2*aisle_halfwidth + 2*fuse_offset #Minimum diameter for one passenger
    while (cabin_width > Dmin) || (cabin_width ≈ Dmin) #While the required diameter is smaller than the fuselage diameter
        seats_per_row = seats_per_row + 1 #Add one more seat
        layout = seat_layouts[seats_per_row] #Find corresponding seat layour from seat dict
        Dmin = seats_per_row*seat_width + (length(layout) - 1)*2*aisle_halfwidth + 2*fuse_offset #New minimum diameter
    end
    
    seats_per_row = seats_per_row - 1 #Subtract 1 seat to find maximum number of seats per row such that cabin_width > Dmin
    return seats_per_row
end

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
        has_leftaisle = aisle_flag(i, layout) #1 if there is an aise to the left of seat, 0 if not
        yseats[i] = yseats[i - 1] + seat_width + has_leftaisle*2*exp_aisle_halfwidth
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
    - `has_leftaisle::Bool`: true if there is an aisle to the left of seat, false if not.
"""
function aisle_flag(idx, layout)
    #Use cumulative sum to find total number of seats to the left of a given aisle.
    #If the difference between the cumsum and the index is exactly 1, the seat has an aisle to the left.
    if 1 in (idx .- cumsum(layout))
        has_leftaisle = 1.0
    else
        has_leftaisle = 0.0
    end
    return has_leftaisle
end

"""
    find_cabin_width(Rfuse::Float64, wfb::Float64, nfweb::Int64, θ::Float64, h_seat::Float64)

This function can be used to calculate the width of the passenger cabin from the double-bubble parameters
and the floor angular position.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Rfuse::Float64`: fuselage exterior radius (m)
    - `wfb::Float64`: lateral shift of double bubble (m)
    - `nfweb::Float64`: number of vertical webs in fuselage
    - `θ::Float64`: angle of floor wrt upper bubble center (rad)
    - `h_seat::Float64`: seat height (m)

    **Outputs:**
    - `w::Float64`: width of cabin (m).
"""
function find_cabin_width(Rfuse::Float64, wfb::Float64, nfweb::Int64, θ::Float64, h_seat::Float64)
    #Use trigonometry to find cabin width
    θseat = asin((h_seat + Rfuse * sin(θ)) / Rfuse)  
    cosθ = min(cos(θ), cos(θseat)) #For the effective cabin width, take the minimum of the widths at the floor and at the seat height
    w = nfweb*2*wfb + 2*Rfuse*cosθ
    return w
end

"""
    find_floor_angles(is_doubledecker::Bool, Rfuse::Float64, dRfuse::Float64; θ1::Float64 = 0.0, h_seat::Float64 = 0.0, d_floor::Float64 = 0.0)

This function can be used to place the passenger decks inside the fuselage. It works for single deck or double decker
cabins. It returns the angular position of each deck with respect to the center of the upper bubble.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `is_doubledecker::Bool`: flag to indicate whether aircraft is a double decker
    - `Rfuse::Float64`: fuselage exterior radius (m)
    - `dRfuse::Float64`: vertical shift of downward bubble (m)
    - `θ1::Float64`: required in some cases; angle of main floor wrt upper bubble center (rad)
    - `h_seat::Float64`: required in some cases; seat height (m)
    - `d_floor::Float64`: required in some cases; distance between double decks (m)

    **Outputs:**
    - `θ1::Float64`: angle of main floor wrt upper bubble center (rad)
    - `θ2::Float64`: returned when double decker; angle of upper floor wrt upper bubble center (rad)
"""
function find_floor_angles(is_doubledecker::Bool, Rfuse::Float64, dRfuse::Float64; θ1::Float64 = 0.0, h_seat::Float64 = 0.0, d_floor::Float64 = 0.0)
    if ~is_doubledecker #If it has a single deck
        θ1 = -asin(h_seat / (2*Rfuse)) #This angle maximizes the cabin width
        return θ1
    else #If it is a double decker with no lower bubble, the main cabin could be anywhere => Use provided angle
        θ2 = asin((Rfuse * sin(θ1) + d_floor) / Rfuse)  
        return θ1, θ2
    end
end

"""
    find_double_decker_cabin_length(x::Vector{Float64}, parg, fuse)

This function can be used to calculate the length of a double decker cabin with different number of 
passengers on each deck.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x::Vector{Float64}`: vector with optimization variables
    - `parg::Vector{Float64}`: vector with aircraft geometric and mass parameters
    - `fuse::Fuselage`: structure with fuselage parameters

    **Outputs:**
    - `maxl::Float64`: required cabin length (m)
    - `pax_per_row_main::Float64`: number of seats abreast in lower cabin
"""
function find_double_decker_cabin_length(x::Vector{Float64}, fuse)
    seat_pitch = fuse.cabin.seat_pitch
    seat_width = fuse.cabin.seat_width 
    aisle_halfwidth = fuse.cabin.aisle_halfwidth
    h_seat = fuse.cabin.seat_height 

    Rfuse = fuse.layout.radius
    dRfuse = fuse.layout.bubble_lower_downward_shift
    wfb = fuse.layout.bubble_center_y_offset
    nfweb = fuse.layout.n_webs
    d_floor = fuse.cabin.floor_distance

    paxsize = fuse.cabin.exit_limit #maximum number of passengers

    paxmain = x[1]
    paxtop = paxsize - paxmain

    try #The calculation could fail for some inputs if an asin returns an error

        θ1 = x[2] #Main deck angle
        _, θ2 = find_floor_angles(true, Rfuse, dRfuse, θ1 = θ1, h_seat= h_seat, d_floor = d_floor)

        #Find width of each cabin
        w1 = find_cabin_width(Rfuse, wfb, nfweb, θ1, h_seat)
        w2 = find_cabin_width(Rfuse, wfb, nfweb, θ2, h_seat)

        #Find length of each cabin
        l1, _, pax_per_row_main = place_cabin_seats(paxmain, w1, seat_pitch = seat_pitch, seat_width = seat_width, aisle_halfwidth = aisle_halfwidth)
        l2, _, _ = place_cabin_seats(paxtop, w2, seat_pitch = seat_pitch, seat_width = seat_width, aisle_halfwidth = aisle_halfwidth)

        maxl = max(l1, l2) #Required length
        return maxl, pax_per_row_main 
    catch
        return 1e6, 1e6
    end
end

function EvaluateCabinProps!(fuse)
    seat_pitch = fuse.cabin.seat_pitch
    seat_width = fuse.cabin.seat_width 
    aisle_halfwidth = fuse.cabin.aisle_halfwidth
    h_seat = fuse.cabin.seat_height 

    Rfuse = fuse.layout.radius
    dRfuse = fuse.layout.bubble_lower_downward_shift
    wfb = fuse.layout.bubble_center_y_offset
    nfweb = fuse.layout.n_webs

    if fuse.n_decks == 1
        θ = find_floor_angles(false, Rfuse, dRfuse, h_seat = h_seat) #Find the floor angle
        wcabin = find_cabin_width(Rfuse, wfb, nfweb, θ, h_seat) #Cabin width
        seats_per_row = findSeatsAbreast(wcabin, seat_width, aisle_halfwidth)

        fuse.cabin.floor_angle_main = θ
        fuse.cabin.cabin_width_main = wcabin
        fuse.cabin.seats_abreast_main = seats_per_row
    else #Floor angles must be provided already, either as input or via optimization
        θ1 = fuse.cabin.floor_angle_main #Main deck angle
        θ2 = fuse.cabin.floor_angle_top

        #Find width of each cabin
        w1 = find_cabin_width(Rfuse, wfb, nfweb, θ1, h_seat)
        w2 = find_cabin_width(Rfuse, wfb, nfweb, θ2, h_seat)

        #Find seats abreast in each cabin
        seats_per_row = findSeatsAbreast(w1, seat_width, aisle_halfwidth)
        seats_per_row_top = findSeatsAbreast(w2, seat_width, aisle_halfwidth)

        fuse.cabin.cabin_width_main = w1
        fuse.cabin.seats_abreast_main = seats_per_row
        fuse.cabin.cabin_width_top = w2
        fuse.cabin.seats_abreast_top = seats_per_row_top
    end

end

"""
    optimize_double_decker_cabin(fuse)

This function can be used to optimize the deck layouts and passenger distribution in a double decker aircraft. 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `parg::Vector{Float64}`: vector with aircraft geometric and mass parameters
    - `fuse::Fuselage`: structure with fuselage parameters

    **Outputs:**
    - `xopt::Vector{Float64}`: vector with optimization results
    - `seats_per_row_main::Float64`: number of seats abreast in lower cabin
"""
function optimize_double_decker_cabin(fuse)
    dRfuse = fuse.layout.bubble_lower_downward_shift

    paxsize = fuse.cabin.exit_limit #maximum number of passengers

    initial_x = [paxsize/2, -pi/4]
    lower = [1.0, -pi/2]
    upper = [paxsize - 1.0, pi/2]
        
    obj(x, grad) = find_double_decker_cabin_length(x, fuse)[1] #Objective function is cabin length

    opt = Opt(:GN_AGS, length(initial_x)) #Use a global optimizer as it is only 1 or 2 variables
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000  # Set the maximum number of function evaluations

    opt.min_objective = obj

    inequality_constraint!(opt, (x, grad) -> MinCargoHeightConst(x, fuse), 1e-5) #Minimum height of cargo hold
    inequality_constraint!(opt, (x, grad) -> MinCabinHeightConst(x, fuse), 1e-5) #Minimum height of upper cabin
    
    (minf,xopt,ret) = NLopt.optimize(opt, initial_x) #Solve optimization problem

    #Evaluate number of passengers per row in main cabin
    _, seats_per_row_main = find_double_decker_cabin_length(xopt, fuse)

    return xopt, seats_per_row_main
end

"""
    MinCargoHeightConst(x, fuse, minheight = 1.626, minwidth = 3.13)

This function evaluates a minimum height constraint on the cargo hold. It returns a number less than or equal 0 if
the constraint is met and a value greater than 0 if it is not.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x::Vector{Float64}`: vector with optimization variables
    - `fuse::Fuselage`: structure with fuselage parameters

    **Outputs:**
    - `constraint::Float64`: this is ≤0 if constraint is met and >0 if not
"""
function MinCargoHeightConst(x, fuse)
    #Extract parameters
    θ1 = x[2]
    Rfuse = fuse.layout.radius
    dRfuse = fuse.layout.bubble_lower_downward_shift

    #Find size of unit load device that must fit in cargo bay
    ULD = fuse.cabin.unit_load_device
    ULDdims = UnitLoadDeviceDimensions[ULD]
    minheight = ULDdims[1]
    minwidth = ULDdims[2] #Base width

    θcargo = -acos(minwidth/(2*Rfuse)) #Angle of cargo hold floor
    hmax = dRfuse + Rfuse * (sin(θ1) - sin(θcargo)) #Maximum height of cargo hold

    constraint = minheight/hmax - 1.0 #Constraint has to be negative if hcargo > minheight
    return constraint
end

"""
    MinCabinHeightConst(x, parg, fuse, minheight = 2.0)

This function evaluates a minimum height constraint on the upper cabin. It returns a number less than or equal 0 if
the constraint is met and a value greater than 0 if it is not.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x::Vector{Float64}`: vector with optimization variables
    - `parg::Vector{Float64}`: vector with aircraft geometric and mass parameters
    - `fuse::Fuselage`: structure with fuselage parameters

    **Outputs:**
    - `constraint::Float64`: this is ≤0 if constraint is met and >0 if not
"""
function MinCabinHeightConst(x, fuse)
    #Extract parameters
    θ1 = x[2] #Main deck angle
    Rfuse = fuse.layout.radius
    dRfuse = fuse.layout.bubble_lower_downward_shift
    d_floor = fuse.cabin.floor_distance

    minheight = fuse.cabin.min_top_cabin_height

    try #The calculation could fail for some inputs if an asin returns an error

        _, θ2 = find_floor_angles(true, Rfuse, dRfuse, θ1 = θ1, h_seat= 0.0, d_floor = d_floor) #The seat height does not need to be included

        hcabin = Rfuse * (1 - sin(θ2)) #Upper cabin height

        constraint = minheight/hcabin - 1.0 #Constraint has to be negative if hcabin > minheight

        return constraint
    catch
        return 1.0
    end
end