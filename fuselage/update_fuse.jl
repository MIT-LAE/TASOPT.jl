function update_fuse!(parg)

    # parg[igRfuse   ] = 90 * in_to_m 
    ltank = parg[iglftankin] # Get length of tank from previous iteration
    parg[igbo] = 2*parg[igRfuse]
    # parg[igbo      ] = 2 * (111 * in_to_m) # 2 × wing centerbox halfspan
    # parg[igbo      ] = 2 * (sqrt(parg[igRfuse]^2 - parg[igzwing]^2) + 1.0 *ft_to_m)# 2 × wing centerbox halfspan
    # Calculate number of seats per row, number of rows and length of cabin
    seats_per_row = 2*Int(parg[igRfuse] ÷ (seat_width + aisle_halfwidth/3))
    rows = Int(ceil(pax / seats_per_row))
    lcabin = rows * seat_pitch 
    # println("Seats per row = $seats_per_row, rows = $rows, lcabin = $(lcabin/ft_to_m) ft")

    parg[igxnose   ] =   0.0 * ft_to_m
    parg[igxblend1 ] =  20.0 * ft_to_m
    parg[igxshell1 ] =  17.0 * ft_to_m
    parg[igxshell2 ] = parg[igxshell1] + lcabin + 20.0*ft_to_m + 2*seat_pitch # 2 ends* 10 ft (for galley (6ft) + lavatory (4ft length) ) + space for emergency_exit

    parg[igxftank ]  = parg[igxshell2] + ltank/2 + 1.0*ft_to_m #(buffer)
    parg[igxblend2 ] = parg[igxftank]  + ltank/2 

    ltshaft = 9.0 * ft_to_m # length of T406 ~ 6.5 ft + 2.5 ft margin
    lgen    = 5.0 * ft_to_m 
    parg[igxtshaft]  = parg[igxblend2] + ltshaft/2
    parg[igxgen   ]  = parg[igxblend2] + ltshaft + lgen/2
    parg[igxcat   ]  = parg[igxgen   ]

    parg[igxconend ] = parg[igxgen] + lgen/2 + 2.0*ft_to_m
    parg[igxend    ] = parg[igxconend] + 5.0*ft_to_m # 5 ft margin/ other things not accounted for

    parg[igxwbox   ] =  70.0 * ft_to_m  # x location of wing box
    parg[igxhbox   ] = parg[igxconend ] - 2*ft_to_m
    parg[igxvbox   ] = parg[igxconend ] - 2*ft_to_m

    parg[igxinv   ]  =  60.0 * ft_to_m
    parg[igxmot   ]  =  parg[igxwbox] # 57.0 * ft_to_m
    parg[igxfan   ]  =  parg[igxwbox] # 55.0 * ft_to_m
    
    parg[igxapu    ] = parg[igxconend] * ft_to_m # xapu      APU location
    
end