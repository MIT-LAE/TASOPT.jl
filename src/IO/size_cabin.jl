"""
    place_cabin_seats(pax, Rfuse)

Function to calculate the seat arrangement in the cabin and, therefore, the required cabin
length.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `pax::Float64`: design number of passengers
    - `Rfuse::Float64`: fuselage radius (m).

    **Outputs:**
    - `lcabin::Float64`: cabin length (m).
    - `xseats::Vector{Float64}`: longitudinal coordinate of each row of seats, measured from front of cabin (m).
    - `seats_per_row::Float64`: number of seats per row.
"""
function place_cabin_seats(pax, Rfuse)
    seat_pitch = 30.0 * in_to_m 
    seat_width = 19.0 * in_to_m
    aisle_halfwidth = 10.0 * in_to_m # per CFR § 25.815 
    cabin_offset = 10 * ft_to_m #Distance to the front and back of seats
    #TODO the hardcoded 10 ft is not elegant

    seats_per_row = Int(2*Rfuse÷ (seat_width + aisle_halfwidth/3))
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
     aisle_halfwidth = 10.0 * in_to_m)

Helper function to arrange seats given a number of `seats_per_row`
and fuselage radius. Assumes default `seat_width = 19"` and `aisle_halfwidth = 10"`,
but can be supplied by the user.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `seats_per_row::Float64`:  number of seats per row.
    - `Rfuse::Float64`: fuselage radius (m).
    - `seat_width::Float64`: width of one seat (m).
    - `aisle_halfwidth::Float64`: half the width of an aisle (m).

    **Outputs:**
    - `yseats::Vector{Float64}`: transverse coordinate of each column of seats, measured from centerline (m).
    - `symmetric_seats::Boolean`: flag for symmetry of seat arrangement
"""
function arrange_seats(seats_per_row, Rfuse,
     seat_width = 19.0 * in_to_m, 
     aisle_halfwidth = 10.0 * in_to_m)

    #Seats
    # Conditions:
    # - No more than 2 seats between any seat and the aisle
    seats_per_row % 2 == 0 ? symmetric_seats = true : symmetric_seats = false

    if symmetric_seats # seating can be symmetric
        half_seats_per_row = seats_per_row ÷ 2
        yseats = zeros(half_seats_per_row)

        if half_seats_per_row <= 3 #Single aisle
            yseats[1] = aisle_halfwidth + seat_width/2 #Aisle in the center
            for col = 2:half_seats_per_row
                yseats[col] = yseats[col-1] + seat_width #offset every seat by width
            end 
        else # twin aisle 
            #If symmetric no more than 2 seats next to each other at the 
            # centerline (I'm not evil enough to create a x-6-x seating arrangement even if "technically" allowed)
            yseats[1] = seat_width/2.0
            yseats[2] = yseats[1] + seat_width
            #Aisle
            half_seats_remaining = half_seats_per_row - 2
            if half_seats_remaining > 4
                @warn "Potentially trying to design a 3 aisle aircraft?
                Seating arrangement not (yet) automatically handled, so check carefully."
            end
            yseats[3] = yseats[2] + aisle_halfwidth*2 + seat_width
            for col = 4:half_seats_per_row
                yseats[col] = yseats[col-1] + seat_width
            end 
        end

    else
        @info "Asymmetric seating only deals with 3 or 5 seats at the moment"
        seating_excess_space = 2*Rfuse - seats_per_row*seat_width - 2*aisle_halfwidth 
        yseats = zeros(seats_per_row)
        # Start from edge and give some space based on the excess space available.
        ind = 1
        yseats[ind] = -Rfuse + seating_excess_space/2 + seat_width/2 
        ind+=1
        if seats_per_row > 3
            yseats[ind] = yseats[ind-1] + seat_width
            ind+=1
        end
        yseats[ind] = yseats[ind-1] + aisle_halfwidth*2 + seat_width
        ind+=1
        for col = ind:seats_per_row
            yseats[col] = yseats[col-1] + seat_width
        end 
    end
    return yseats, symmetric_seats
end  # function arrange_seats