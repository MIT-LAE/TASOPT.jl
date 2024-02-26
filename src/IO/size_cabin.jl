function find_cabin_length(pax, Rfuse)
    seat_pitch = 30.0 * in_to_m 
    seat_width = 19.0 * in_to_m
    aisle_halfwidth = 10.0 * in_to_m # per CFR § 25.815 
    cabin_offset = 10 * ft_to_m #Distance to the front and back of seats

    seats_per_row = Int(2*Rfuse÷ (seat_width + aisle_halfwidth/3))
    rows = Int(ceil(pax / seats_per_row))

    println("Seats per row = $seats_per_row, Total rows = $rows")

    if seats_per_row <= 10
        emergency_rows = [12, 13]
    else
        emergency_rows = [19, 20]
    end

    xseats = zeros(rows)'
    xseats[1] = 0
    for r in 2:rows
        emergency_exit = 0.0
        if (r in emergency_rows)
            emergency_exit = seat_pitch/2
        end
        xseats[r] = xseats[r-1] + seat_pitch + emergency_exit
    end

    lcabin = xseats[end] + cabin_offset
    return lcabin
end