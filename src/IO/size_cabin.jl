"""
    place_cabin_seats(pax, Rfuse, seat_pitch = 30.0*in_to_m, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m)

Function to calculate the seat arrangement in the cabin and, therefore, the required cabin
length.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pax::Float64`: design number of passengers
    - `Rfuse::Float64`: fuselage radius (m).
    - `seat_width::Float64`: width of one seat (m).
    - `aisle_halfwidth::Float64`: half the width of an aisle (m).
    - `fuse_offset::Float64`: distance from outside of fuselage to edge of closest window seat (m).

    **Outputs:**
    - `lcabin::Float64`: cabin length (m).
    - `xseats::Vector{Float64}`: longitudinal coordinate of each row of seats, measured from front of cabin (m).
    - `seats_per_row::Float64`: number of seats per row.
"""
function place_cabin_seats(pax, Rfuse, seat_pitch = 30.0*in_to_m, 
    seat_width = 19.0*in_to_m, aisle_halfwidth = 10.0*in_to_m, fuse_offset = 6.0*in_to_m)

    cabin_offset = 10 * ft_to_m #Distance to the front and back of seats
    #TODO the hardcoded 10 ft is not elegant

    #Calculate the maximum number of seats per row
    seats_per_row = 1
    Dmin = seats_per_row*seat_width + 2*aisle_halfwidth + 2*fuse_offset #Minimum diameter for one passenger
    while (2*Rfuse > Dmin) || (2*Rfuse ≈ Dmin) #While the required diameter is smaller than the fuselage diameter
        seats_per_row = seats_per_row + 1 #Add one more seat
        layout = seat_layouts[seats_per_row] #Find corresponding seat layour from seat dict
        Dmin = seats_per_row*seat_width + (length(layout) - 1)*2*aisle_halfwidth + 2*fuse_offset #New minimum diameter
    end
    
    seats_per_row = seats_per_row - 1 #Subtract 1 seat to find maximum number of seats per row such that 2*Rfuse > Dmin

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
    arrange_seats(seats_per_row, Rfuse,
     seat_width = 19.0 * in_to_m, 
     aisle_halfwidth = 10.0 * in_to_m,
     fuse_offset = 6.0*in_to_m)

Helper function to arrange seats given a number of `seats_per_row`
and fuselage radius. Assumes default `seat_width = 19"` and `aisle_halfwidth = 10"`,
but can be supplied by the user.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `seats_per_row::Float64`:  number of seats per row.
    - `Rfuse::Float64`: fuselage radius (m).
    - `seat_width::Float64`: width of one seat (m).
    - `aisle_halfwidth::Float64`: half the width of an aisle (m).
    - `fuse_offset::Float64`: distance from outside of fuselage to edge of closest window seat (m).

    **Outputs:**
    - `yseats::Vector{Float64}`: transverse coordinate of each column of seats, measured from centerline (m).
"""
function arrange_seats(seats_per_row, Rfuse,
     seat_width = 19.0 * in_to_m, 
     aisle_halfwidth = 10.0 * in_to_m, 
     fuse_offset = 6.0*in_to_m)

    #Seats
    # Conditions:
    # - No more than 2 seats between any seat and the aisle

    layout = seat_layouts[seats_per_row] #find seat layout from dictionary
    n_aisles = length(layout) - 1 #Number of aisles
    Dmin = seats_per_row*seat_width + n_aisles*2*aisle_halfwidth + 2*fuse_offset #Minimum required diameter

    exp_aisle_halfwidth = aisle_halfwidth + (2*Rfuse - Dmin)/(2*n_aisles) #Expanded aisle to take up all available space

    yseats = zeros(seats_per_row)
    yseats[1] = fuse_offset + seat_width/2 #First seat is window seat
    for i = 2:seats_per_row
        flag_aisle = aisle_flag(i, layout) #1 if there is an aise to the left of seat, 0 if not
        yseats[i] = yseats[i - 1] + seat_width + flag_aisle*2*exp_aisle_halfwidth
    end

    #Shift seat coordinates to start in centerline
    yseats = yseats .- Rfuse

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