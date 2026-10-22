

"""
    stickfig(ac::aircraft; plot_obj = nothing, label_fs = 16, 
            annotate_text = true, annotate_length = true, 
            annotate_group = true, show_grid = false)

`stickfig` plots a "stick figure" airplane on the plot or subplot object `plot_obj` if provided 
(else a new plot is created). This is used in conjuction with other functions
like [`plot_details`](@ref) to create summary plots to track progress of optimization
or present results.
"""
function stickfig(ac::aircraft; plot_obj = nothing, label_fs = 16, 
    annotate_text = true, annotate_length = true, annotate_group = true, show_grid = false)

    #if aircraft is not sized, cannot plot
    if !ac.is_sized[1] 
        @warn "The aircraft ($(ac.name)) must be sized before being plotted. Skipping `stick_fig`..."
        return nothing
    end

    # Unpack aircraft components
    parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, engine = unpack_ac(ac,1) #imission = 1 

    # Wing
        co = wing.layout.root_chord
        cs = wing.layout.root_chord*wing.inboard.λ
        ct = wing.layout.root_chord*wing.outboard.λ

        sweep = wing.layout.sweep
        λs = wing.inboard.λ
        λt = wing.outboard.λ

        bo = wing.layout.root_span
        bs = wing.layout.break_span
        b  = wing.layout.span

        xax = wing.layout.spar_box_x_c
        xcLE = -xax
        xcTE = 1.0 - xax

        dx = wing.layout.box_x    
        etas = bs/b
        etao = bo/b
        cosL = cos(sweep*pi/180.0)
        tanL = tan(sweep*pi/180.0)

        xs = tanL*(bs-bo)/2.0
        xt = tanL*(b -bo)/2.0

        xw = zeros(6)
        yw = zeros(6)

        #X locations of wing vertices
        xw[1] =      co*xcLE + dx
        xw[2] = xs + cs*xcLE + dx
        xw[3] = xt + ct*xcLE + dx
        xw[4] = xt + ct*xcTE + dx
        xw[5] = xs + cs*xcTE + dx
        xw[6] =      co*xcTE + dx
        
        #Y locations of wing vertices
        yw[1] = bo/2.0
        yw[2] = bs/2.0
        yw[3] = b /2.0
        yw[4] = b /2.0
        yw[5] = bs/2.0
        yw[6] = bo/2.0

    # Fuse
        fuselage = ac.fuselage
        Rfuse = fuselage.layout.radius
        wfb   = fuselage.layout.bubble_center_y_offset

        anose    = fuselage.layout.nose_radius
        btail    = fuselage.layout.tail_radius

        xnose    = fuselage.layout.x_nose
        xend     = fuselage.layout.x_end
        xblend1  = fuselage.layout.x_start_cylinder
        xblend2  = fuselage.layout.x_end_cylinder
        xhtail   = htail.layout.x
        xvtail   = vtail.layout.x
        xwing    = wing.layout.x

        xhbox    = htail.layout.box_x

        lcyl = xblend2 - xblend1
        xtail = xvtail 
        
        hwidth = Rfuse + wfb
        
        nnose = 10
        ntail = 10

        xf = zeros(nnose + ntail + 1)
        yf = zeros(nnose + ntail + 1)

        if compare_strings(fuselage.layout.opt_tapers_to, "point")
            dytail = -hwidth 
        else # to "edge"
            dytail = -0.2*hwidth
        end

        for i = 1: nnose
            fraci = float(i-1)/float(nnose-1)
            fracx = cos(0.5*pi*fraci)

            k = i
            xf[k] = xblend1 + (xnose-xblend1)*fracx
            yf[k] = hwidth*(1.0 - fracx^anose)^(1.0/anose)
        end

        for i = 1: ntail
            fraci = float(i-1)/float(ntail-1)
            fracx = fraci

            k = i + nnose
            xf[k] = xblend2 + (xend-xblend2)*fracx
            yf[k] = hwidth + dytail*fracx^btail
        end

        k = k+1
        xf[k] = xf[k-1]
        yf[k] = 0.

    
    # Tail
        xh = zeros(6)
        yh = zeros(6)
        
        boh = htail.layout.root_span
        Sh  = htail.layout.S
        ARh = htail.layout.AR
        lambdah = htail.outboard.λ
        sweeph  = htail.layout.sweep

        bh = sqrt(Sh*ARh)
        coh = Sh/(boh + (bh-boh)*0.5*(1.0+lambdah))


        dx = xhbox
        tanLh = tan(sweeph*π/180.0)
        cth = coh*lambdah

        xaxh = 0.40

        xoLEh = coh*(    - xaxh) + dx
        xoTEh = coh*(1.0 - xaxh) + dx
        xtLEh = cth*(    - xaxh) + dx + 0.5*(bh - boh)*tanLh
        xtTEh = cth*(1.0 - xaxh) + dx + 0.5*(bh - boh)*tanLh

        yoLEh = 0.5*boh
        yoTEh = 0.5*boh
        ytLEh = 0.5*bh
        ytTEh = 0.5*bh


        if compare_strings(fuselage.layout.opt_tapers_to,"point")
            xcLEh = xoLEh
            xcTEh = xoTEh
            ycLEh = yoLEh
            ycTEh = yoTEh
        else #to "edge"
            xcLEh = coh*(    - xaxh) + dx + 0.5*(0. - boh)*tanLh
            xcTEh = coh*(1.0 - xaxh) + dx + 0.5*(0. - boh)*tanLh
            ycLEh = 0.0
            ycTEh = 0.0
        end

        if  (false)
            yoLEh = hwidth
            xoLEh = yoLEh*tanLh + xcLEh
            if (xoLEh > xblend2)
                fracx = min( (xoLEh-xblend2)/(xend-xblend2) , 1.0 )
                yoLEh = hwidth + dytail*fracx^btail
                xoLEh = yoLEh*tanLh + xcLEh
            end
            if (xoLEh > xblend2)
                fracx = min( (xoLEh-xblend2)/(xend-xblend2) , 1.0 )
                yoLEh = hwidth + dytail*fracx^btail
                xoLEh = yoLEh*tanLh + xcLEh
            end

            yoTEh = 0.0
            xoTEh = yoTEh*tanLh + xcTEh
            if (xoTEh > xblend2)
                fracx = min( (xoTEh-xblend2)/(xend-xblend2) , 1.0 )
                yoTEh = hwidth + dytail*fracx^btail
                xoTEh = yoTEh*tanLh + xcTEh
            end
            if (xoTEh > xblend2)
                fracx = min( (xoTEh-xblend2)/(xend-xblend2) , 1.0 )
                yoTEh = hwidth + dytail*fracx^btail
                xoTEh = yoTEh*tanLh + xcTEh
            end
        
        end
        


        xh[ 1] = xcLEh
        xh[ 2] = xoLEh
        xh[ 3] = xtLEh
        xh[ 4] = xtTEh
        xh[ 5] = xoTEh
        xh[ 6] = xcTEh
  
        yh[ 1] = ycLEh
        yh[ 2] = yoLEh
        yh[ 3] = ytLEh
        yh[ 4] = ytTEh
        yh[ 5] = yoTEh
        yh[ 6] = ycTEh
  
        #Initialize seat start x-position
        xseats0 = fuselage.layout.x_pressure_shell_fwd
        # Fuel tank
        ntank = 8
        Rtank = Rfuse - 0.1 # Account for clearance_fuse
        l = max(parg[iglftankin], parg[iglftank])
        nftanks = fuse_tank.tank_count #Number of fuel tanks
        ARtank = 2.0

        if nftanks != 0
            tank_placement = ac.fuse_tank.placement
            if tank_placement == "front"
                xtanks = [parg[igxftank]]
                xseats0 = xtanks[1] + l/2 + 1.0 * ft_to_m #move seats backwards
            elseif tank_placement == "rear"
                xtanks = [parg[igxftankaft]]
                xseats0 = fuselage.layout.x_pressure_shell_fwd
            elseif tank_placement == "both"
                xtanks = [parg[igxftank], parg[igxftankaft]]
                xseats0 = xtanks[1] + l/2 + 1.0 * ft_to_m #move seats backwards
            end
            
            xt = zeros(nftanks, ntank*2 )
            yt = zeros(nftanks, ntank*2 )
            for m = 1:nftanks
                xcyl0 = xtanks[m] - l/2 + Rtank/ARtank
                xcyl1 = xtanks[m] + l/2 - Rtank/ARtank
                
                for i = 1: ntank
                    fraci = float(i-1)/float(ntank-1)
                    fracx = cos(0.5*pi*fraci)

                    k = i
                    xt[m, k] = xcyl0 - Rtank/ARtank*fracx
                    yt[m, k] = sqrt(Rtank^2 * max((1 - ((xt[m, k]-xcyl0)/(Rtank/ARtank))^2), 0.0) )
                end
                # k = k+1
                # xt[k] = xcyl0 + parg[iglftank]
                # yt[k] = Rtank
                for i = 1: ntank
                    fraci = float(i-1)/float(ntank-1)
                    fracx = sin(0.5*pi*fraci)

                    k = i + ntank
                    xt[m, k] = xcyl1 + (xcyl1 + Rtank/ARtank - xcyl1)*fracx
                    yt[m, k] = sqrt(Rtank^2 * max((1 - ((xt[m, k]-xcyl1)/(Rtank/ARtank))^2), 0.0) )
                end
            end
        end

        # xt = LinRange(xcyl0 - Rfuse/ARtank , xcyl0, 20 )
        # yt = zeros(length(xt))
        # @. yt = sqrt(Rfuse^2 * max((1 - ((xt-xcyl0)/(Rfuse/ARtank))^2), 0.0) )

        xshell = zeros(ntank)
        yshell = zeros(ntank)
        AR = 3.0
        xshellcenter = fuselage.layout.x_pressure_shell_aft - Rfuse/AR
        for i = 1: ntank
            fraci = float(i-1)/float(ntank-1)
            fracx = sin(0.5*pi*fraci)

            k = i
            xshell[k] = xshellcenter + Rfuse/AR *fracx
            yshell[k] = sqrt(Rfuse^2 * max((1 - ((xshell[k]-xshellcenter)/(Rfuse/AR))^2), 0.0) )
        end

    if !(options.is_doubledecker) #Only show seats in single deck arrangements
        h_seat = fuselage.cabin.seat_height
        pax = parg[igWpay]/parm[imWperpax]
        Rfuse = fuselage.layout.radius
        dRfuse = fuselage.layout.bubble_lower_downward_shift
        wfb = fuselage.layout.bubble_center_y_offset
        nfweb = fuselage.layout.n_webs

        θ = find_floor_angles(false, Rfuse, dRfuse, h_seat = h_seat) #Find the floor angle
        wcabin = find_cabin_width(Rfuse, wfb, nfweb, θ, h_seat) #Cabin width

        seat_pitch = fuselage.cabin.seat_pitch
        seat_width = fuselage.cabin.seat_width 
        aisle_halfwidth = fuselage.cabin.aisle_halfwidth
        front_seat_offset = fuselage.cabin.front_seat_offset

        _, xseats, seats_per_row = place_cabin_seats(pax, wcabin, seat_pitch = seat_pitch, seat_width = seat_width, 
            aisle_halfwidth = aisle_halfwidth, front_seat_offset = front_seat_offset)

        xseats = xseats[:] .+ xseats0
        rows = length(xseats)

        println("Seats per row = $seats_per_row, Total rows = $rows")
        yseats = arrange_seats(seats_per_row, wcabin)
    end

    # Plot
    #if starting from scratch, generate plot_obj
    if isnothing(plot_obj)
        #conditional statement
        # make fig size larger if adding annotate_text
        figsize = annotate_text ? (800, 500) : (650, 500)

        plot_obj = plot(
            legend = false,
            size = figsize,
            dpi = 300,
            grid = show_grid,
            show=true,
            aspect_ratio=:equal,
            margin=4mm
            )
    end
    
    # Plot wing (upper and lower surfaces)
    plot!(plot_obj, xw, yw, label = "", color = :black)
    plot!(plot_obj, xw, -yw, label = "", color = :black)

    # Panel break (alpha transparency and line width)
    plot!(plot_obj, xw[[2, 5]], yw[[2, 5]], label = "", color = :black, linewidth = 1.0, alpha = 0.5)
    plot!(plot_obj, xw[[2, 5]], -yw[[2, 5]], label = "", color = :black, linewidth = 1.0, alpha = 0.5)
    
    # Plot Tail (conditionally plot fill based on tailz value)
    tailz = htail.layout.z > 0 ? :back : 1
    plot!(plot_obj, xh, yh, label = "", color = :white, alpha = 0.8, linecolor = :black, linewidth = 2.0, z_order = tailz)
    # plot!(xh, -yh, label = "", color = :white, alpha = 0.8, linecolor = :black, linewidth = 2.0, z_order = tailz)
    plot!(plot_obj, xh, -yh, fillrange=yh, label = "", color = :white, alpha = 0.8, linecolor = :black, linewidth = 2.0, z_order = tailz)
    
    # Tail box (fill between xvt and yvt for the tail volume)
    xvt = [-0.4, -0.15, 0.2, 0.6].*vtail.layout.root_chord .+ vtail.layout.box_x
    yvt = hcat([0.0 ones(length(xvt) - 2)' .*(vtail.layout.root_chord*vtail.outboard.cross_section.thickness_to_chord/2) 0.0])[:]
    plot!(plot_obj, xvt, -yvt, fillrange=yvt, label = "", color = :black, alpha = 0.8, linecolor = :black, z_order = :front)
    
    # Plot fuselage (fill with edge color)
    plot!(plot_obj, xf, -yf, label = "", color = :black, linewidth = 2.0, z_order = 5)
    plot!(plot_obj, xf, yf, label = "", color = :black, linewidth = 2.0, z_order = 5)
    plot!(plot_obj, xf, -yf, fillrange=yf, label = "", color = :white, edgecolor = :black, linewidth = 2.0, z_order = 5)

    # Tank plotting
    if !(ac.options.has_wing_fuel)
        for m in 1:nftanks
            # Tank outline
            plot!(plot_obj, xt[m,:], yt[m,:], color=:black, lw=1.5, z_order=10)
            plot!(plot_obj, xt[m,:], -yt[m,:], color=:black, lw=1.5, z_order=10)
            
            # Filled area between tank outline
            plot!(xt[m, :], yt[m, :], fill_between=(-yt[m, :], yt[m, :]), color=:red, alpha=0.1)
            
            # Fuel name
            fuelname = ""
            if options.ifuel == 11
                fuelname ="\$CH_4\$"
            elseif options.ifuel == 40
                fuelname ="\$LH_2\$"
            end
            annotate!(plot_obj, (xtanks[m], 0.0, text(fuelname, label_fs - 2.0, :center, :center)))
        end
    end

    # Xshell2 plotting
    plot!(plot_obj, xshell, yshell, color=:black, lw=1.5, z_order=10)
    plot!(plot_obj, xshell, -yshell, color=:black, lw=1.5, z_order=10)

    # Plot Engines
    D = parg[igdaftfan]
    lnace = parg[iglnaceaft]
    x = parg[igxtshaft]

    # Aft fan outline
    plot!(plot_obj, [x, x, x + lnace, x + lnace, x], [D / 8, D / 8 + D, D / 8 + D * 3 / 4, D / 8 + D / 4, D / 8],
        lw=1.5, color=:red, z_order=:front)
    plot!(plot_obj, [x, x, x + lnace, x + lnace, x], [-D / 8, -D / 8 - D, -D / 8 - D * 3 / 4, -D / 8 - D / 4, -D / 8],
        lw=1.5, color=:red, z_order=:front)

    D = parg[igdfan]
    neng = parg[igneng]
    lnace = parg[iglnace]
    ηs = wing.layout.ηs
    dy = 2 * D  # Space to leave near wing root and tip [m]

    if neng == 2
        yi = [ηs * b / 2]
    else
        yi = range(bo / 2 + dy, b / 2 * 3 / 4, length=Int(neng / 2))
    end

    xi = zeros(length(yi))
    ηi = yi / (b / 2)
    ηo = bo / b
    ci = zeros(length(yi))

    for (i, η) in enumerate(ηi)
        if η <= ηs
            ci[i] = co * (1 + (λs - 1) * (η - ηo) / (ηs - ηo))
        else
            ci[i] = co * (λs + (λt - λs) * (η - ηs) / (1 - ηs))
        end
    end

    tanL = tan(wing.layout.sweep*π/180.0)
    @. xi = tanL * (yi - bo/2) - 0.4ci + wing.layout.box_x - 1.0
    xlocations = vec(hcat([xi, xi, xi.+lnace, xi.+lnace, xi]...)) #hcat is to avoid this being an array of arrays
    ylocations = vec(hcat([yi.-D/2, yi.+D/2, yi.+D/3, yi.-D/3, yi.-D/2]...))
    
    # Plot engine locations
    plot!(plot_obj, xlocations, ylocations, color=:red, linewidth=1.5)

    # Plot NP and annotate
    scatter!(plot_obj, [parg[igxNP]], [0.0], color=:blue, marker=:circle, z_order=15, label="NP")
    annotate!(plot_obj, (parg[igxNP]+2, 0.5, text("NP", label_fs - 2.0, :center, :center, color=:blue)))

    # Annotate CG range
    plot!(plot_obj, [parg[igxCGfwd], parg[igxCGaft]], [0.0, 0.0], color=:black, lw=2.0, label="CG movement")  # Line between points
    scatter!(plot_obj, [parg[igxCGfwd], parg[igxCGaft]], [0.0, 0.0], marker=:vline, color=:black)                # End points
    annotate!(plot_obj, parg[igxCGfwd]-2, 0, text("CG", label_fs - 2.0, :center, :center))

    # Show seats (single-deck case)
    if !(options.is_doubledecker)
        xgrid = repeat(xseats, inner=length(yseats))
        ygrid = repeat(yseats, outer=length(xseats))

        scatter!(plot_obj, xgrid, ygrid,
            color=:gray, alpha=0.1, marker=:rect, 
            ms=2, z_order=:front, label="")
    end

    # Annotations
    if annotate_text
        annotate!(plot_obj, 
            maximum(xh)*1.11, 21, #position based on aircraft length (outside plotting area)
            text("PFEI = $(round(parm[imPFEI], digits = 3))\n" *
            "Mₘₐₓ = $(round(maximum(para[iaMach,:,:]), digits=3))\n" *
            "WMTO = $(round(parg[igWMTO] / 9.81 / 1000, digits=1)) t\n" *
            "Span = $(round(wing.layout.span, digits=1)) m\n" *
            "c₀ = $(round(wing.layout.root_chord, digits=1)) m\n" *
            "Λ = $(round(wing.layout.sweep, digits=1))°\n" *
            "Rfuse = $(round(fuselage.layout.radius,digits=2)) m\n" *
            "L/D = $(round(para[iaCL,ipcruise1,1] / para[iaCD,ipcruise1,1], digits=2))",
            label_fs, #fontsize
            halign=:left, valign=:top, color=:black),
            z_order=:front)

        
        plot!(plot_obj, right_margin=200px) #"left-aligns" the plot within the figure
    end

    # Span annotations
    groups, bmax = find_aerodrome_code(wing.layout.max_span)
    xcode = -2.0

    if annotate_group
        box_color = :orange
        plot!(plot_obj, [xcode, xcode], [-bmax / 2, bmax / 2], lw=5, alpha=0.2, color=box_color, label="")
        plot!(plot_obj, [xcode, 40.0], [bmax / 2, bmax / 2], lw=5, alpha=0.2, color=box_color, label="")
        plot!(plot_obj, [xcode, 40.0], [-bmax / 2, -bmax / 2], lw=5, alpha=0.2, color=box_color, label="")
        annotate!(plot_obj, 20, bmax / 2 + 1, text("ICAO Code $(groups[1])/ FAA Group $(groups[2])",
            label_fs, color=box_color, #fontsize
            halign=:center, valign=:center))
    end
    
    if annotate_length
        yloc = (-bmax / 2) * 1.11
        plot!(plot_obj, [0.0, xf[end]], [yloc, yloc], lw=1.5, color=:black, label="")       #line
        scatter!(plot_obj, [0.0, xf[end]], [yloc, yloc], yerror=0.03*yloc, markersize=0)    #end ticks
        annotate!(plot_obj, xf[end] *0.3, yloc+1.5, text("\$\\ell\$ = $(round(xf[end], digits=1)) m", #text
            halign=:center, valign=:center),
            label_fs) #fontsize
    end

    # Set plot limits and labels
    ylims!(plot_obj, -1.2 * bmax / 2, 1.2 * bmax / 2)
    xlims!(plot_obj, -5, maximum(xh)*1.1)
    xlabel!(plot_obj, "x [m]", fontsize=label_fs)
    ylabel!(plot_obj, "y [m]", fontsize=label_fs)

    return plot_obj
end

"""
    find_aerodrome_code(b)

`find_aerodrome_code` finds the airport code corresponding to a given aircraft wingspan. It returns 
the codes in the ICAO and FAA formats, as well as the maximum wingspan for these codes.
"""
function find_aerodrome_code(b::Float64)
    max_bs = sort(collect(keys(aerodrome_codes)))
    idx = 1
    for cand_maxb in max_bs
        if b >= cand_maxb 
            idx += 1
        end
    end
    maxb = max_bs[idx]
    groups = aerodrome_codes[maxb]
    return groups, Float64(maxb)
end

#TODO: re-do/re-apply label_bars() to use here
"""
    plot_drag_breakdown(ac::aircraft; ip = [ipclimb1, ipcruise1, ipcruisen, ipdescent1, ipdescentn], imission = 1, show_fractions = true, show_values = false)

Generates a stacked bar plot showing the breakdown of drag components for specified flight points of an aircraft.

# Arguments
- `ac::aircraft`: Aircraft object containing aerodynamic data
- `ip::Union{Int, AbstractVector{Int}, Colon}`: Mission point(s) to analyze. Can be a single integer, vector of integers, or `:` for all points
- `imission::Int`: Mission number to analyze (default: 1)
- `show_fractions::Bool`: If true, shows components as fractions of total drag (default: true)
- `show_values::Bool`: If true, displays numerical values on each bar segment (default: false). Wouldn't advise true for len(ip) > 5

# Returns
- `Plots.Plot`: A stacked bar plot object showing drag components:
  - CDi: Induced drag
  - CDnace: Nacelle drag
  - CDvtail: Vertical tail drag
  - CDhtail: Horizontal tail drag
  - CDwing: Wing drag
  - CDfuse: Fuselage drag

# Notes
- Values less than 0.001 are not labeled when `show_values = true`
- When `show_fractions = true`, components are normalized by total drag coefficient
"""
function plot_drag_breakdown(ac::aircraft; 
    ip::Union{Int, AbstractVector{Int}, Colon} = [ipclimb1, ipcruise1, ipcruisen, ipdescent1, ipdescentn], imission::Int = 1,
    show_fractions::Bool = true, show_values::Bool = false)
    
    #get aircraft aero data
    para = view(ac.para, :, :, imission)

    #if all points requested, generate ip
    if ip == Colon()
        ip = 1:1:iptotal
    end

    # make ip a vector for consistency
    ips = isa(ip, Int) ? [ip] : ip

    # fetch/arrange the data
    LoD = [para[iaCL, i] for i in ips]./ [para[iaCD, i] for i in ips] 
    CL = [para[iaCL, i] for i in ips]
    CD = [para[iaCD, i] for i in ips]
    CDfuse  = [para[iaCDfuse, i]  for i in ips]
    CDi     = [para[iaCDi, i]     for i in ips]
    CDwing  = [para[iaCDwing, i]  for i in ips]
    CDhtail = [para[iaCDhtail, i] for i in ips]
    CDvtail = [para[iaCDvtail, i] for i in ips]
    CDnace  = [para[iaCDnace, i]  for i in ips]

    if show_fractions
        # Compute fractional contributions
        CDfuse  = CDfuse  ./ CD
        CDi     = CDi     ./ CD
        CDwing  = CDwing  ./ CD
        CDhtail = CDhtail ./ CD
        CDvtail = CDvtail ./ CD
        CDnace  = CDnace  ./ CD
    end

    # Combine drags into a matrix (points x components)
    data = hcat(CDi, CDnace, CDvtail, CDhtail, CDwing, CDfuse)

    # X positions
    x = 0:length(ips)-1

    # Plot
    p = groupedbar(
        x, data,
        bar_position = :stack,
        label = ["CDi" "CDnace" "CDvtail" "CDhtail" "CDwing" "CDfuse"],
        # xlabel = "ip", 
        ylabel = show_fractions ? "Fraction of total drag, \$C_{D,( )} / C_D\$" : "Drag component, \$C_{D, ( )}\$",
        title = "Drag Breakdown",
        xticks = (x, ["$(ip_labels[i])" for i in ips]),
        xrotation=-35,
        legend = :outerright
        )

    # Add data labels on each box
    if show_values
        for (i, series) in enumerate(eachcol(data))
            for (j, value) in enumerate(series)
                if value > 0.001  # Only label if value is significant
                    # Calculate cumulative sum up to but not including current value
                    y_base = sum(data[j, 1:i-1])
                    y_top = show_fractions ? 1 : CD[j]
                    y_pos = y_top - y_base - value/2 # Position label in middle of current segment
                    annotate!(p, j-1, y_pos, 
                        text(@sprintf("%.5f", value), 9, :white, :center))
                end
            end
        end    
    end
    
    return p #plot object
end

"""
    plot_details(ac::aircraft)

`plot_details` combines a [`stickfig`](@ref) plot along with a mission summary,
weight and drag buildup stacked bar charts to present results.
"""
function plot_details(ac::aircraft; imission::Int=1)

    parg, parm, para, pare, options, fuselage, fuse_tank, wing, htail, vtail, engine, landing_gear = unpack_ac(ac, imission)

  ## Gather data to be plotted
    # Drag build-up for subplot 2
    LoD     = para[iaCL, ipcruise1]/ para[iaCD, ipcruise1]
    CL      = para[iaCL, ipcruise1]
    CD      = para[iaCD, ipcruise1]
    CDfuse  = para[iaCDfuse, ipcruise1]
    CDi     = para[iaCDi, ipcruise1]
    CDwing  = para[iaCDwing, ipcruise1]
    CDhtail = para[iaCDhtail, ipcruise1]
    CDvtail = para[iaCDvtail, ipcruise1]
    CDnace  = para[iaCDnace, ipcruise1]

    CDfusefrac  = CDfuse /CD
    CDifrac     = CDi    /CD
    CDwingfrac  = CDwing /CD
    CDhtailfrac = CDhtail/CD
    CDvtailfrac = CDvtail/CD
    CDnacefrac  = CDnace /CD

    # Weight build-up for subplot 2
    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Whpesys = parg[igWMTO] * fuselage.HPE_sys.W
    Wlgnose = ac.landing_gear.nose_gear.weight.W
    Wlgmain = ac.landing_gear.main_gear.weight.W
    Wtotadd = Whpesys + Wlgnose + Wlgmain
    
    Wpay  = parg[igWpay]
    Wfuel = parg[igWfuel]
    WMTO  = parg[igWMTO]

    Wwing  = wing.weight
    Wfuse  = fuselage.weight
    Wvtail = vtail.weight
    Whtail = htail.weight
    Weng = parg[igWeng]
    Wtesys = parg[igWtesys]
    Wftank = parg[igWftank]

    Wwingfrac = Wwing /WMTO
    Wfusefrac = Wfuse /WMTO
    Wvtailfrac = Wvtail/WMTO
    Whtailfrac = Whtail/WMTO
    Wtesysfrac = Wtesys/WMTO
    Wftankfrac = Wftank/WMTO
    Wtotaddfrac = Wtotadd/WMTO
    Wengfrac = Weng/WMTO

    Wemptyfrac = Wempty/WMTO
    Wfuelfrac  = Wfuel /WMTO
    Wpayfrac   = Wpay  /WMTO

    
    # For subplot 3
    h     = [para[iaalt,ipclimb1:ipcruisen]; 0.0]./ft_to_m./1000 # show in 1000s of ft.
    R     = [para[iaRange,ipclimb1:ipcruisen]; para[iaRange, ipdescentn]]./nmi_to_m
    deNOx = pare[iedeNOx, :]
    fracW = [para[iafracW, ipclimb1:ipcruisen]; para[iafracW, ipdescentn]]
    mdotf = pare[iemdotf, :]
    mdotH2O = pare[iemdotf, :].*9.0
    gamV = [para[iagamV, ipclimb1:ipcruisen]; para[iagamV, ipdescentn]]

    bar_width = 0.2

  ## Do the plotting
    #subplots are created and synthesized into one figure at the end

  #=
    Subplot 1: draw stick figure
    =#

    # stickfigout = stickfig(ac, plot_obj = fig[1])
    p1 = stickfig(ac)
    # plot!(stickfigout, subplot=1, layout=layout)
    #undo the right margin adjustment

    #=
    Subplot 2: bar graphs for drag, weight buildups
    =#
    # Bar plot 1: Drag Components (CD)
    p2 = groupedbar(
            [0],
            [CDifrac, CDnacefrac, CDvtailfrac, CDhtailfrac, CDwingfrac, CDfusefrac]',
            bar_position=:stack,
            width=bar_width,
            label=["CDi" "CDnace" "CDvtail" "CDhtail" "CDwing" "CDfuse"],
            ylabel="Fraction",
            title="Drag, Weight Breakdowns",
            )

    # # Bar plot 2: Weight Fractions (WMTO)
    groupedbar!(p2,
            [1],
            [Wpayfrac, Wfuelfrac, Wemptyfrac]',
            bar_position=:stack,
            width=bar_width,
            label=["Wpay" "Wfuel" "Wempty"],
            # xlabel="Weight Fractions",
            # ylabel="Fraction",
            # title="Weight Fractions (WMTO)",
            )

    # Add lines for Wempty fractions
    plot!(p2,
        [0.5, 2.5], [Wemptyfrac, Wemptyfrac],
        lw=2, color=:black, linestyle=:dash,
        label=" "
    )

    # Bar plot 3: Other Weight Fractions
    groupedbar!(p2,
        [2],
        [Wengfrac, Wtotaddfrac, Wftankfrac, Wtesysfrac, Wvtailfrac, Whtailfrac, Wwingfrac, Wfusefrac]',
        bar_position=:stack,
        width=bar_width,
        label=["Weng" "Wadd" "Wftank" "Wtesys" "Wvtail" "Whtail" "Wwing" "Wfuse"],
        # xlabel="Weight Fractions",
        # ylabel="Fraction",
        # title="Detailed Weight Fractions",
        legend = :outerright
    )

    xticks!(p2, [0, 1, 2], ["Drag","Total\nWeight","Empty\nWeight"],
            )

    #TODO: re-do/re-apply label_bars() to use here

    #=
    Subplot 3: mission profile
    =#

    # Plot altitude vs range
    p3 = plot(R, h, label="Altitude", xlabel="Range [nmi]", 
            title="Mission Profile",ylabel="Altitude [kft]", lw=2.)

    yaxis2 = twinx(p3)
    # Overlay climb angle on the secondary y-axis
    plot!(yaxis2, R, gamV, color=:red, label="Traj. Climb Angle, γ", ylabel="Angle [rads]", lw=2,legend=:bottomright)
    # Add a horizontal line at climb angle = 0.015
    hline!(yaxis2, [0.015], lw=2.0, color=:red, linestyle=:dash, label="γ = 0.015")

    # Define layout
    layout = @layout [A; B{0.6w} C]

    # Generate the figure
    fig = plot(p1, p2, p3, layout=layout, 
                size=(800, 900), dpi=300,
                margin=4mm)

## Send it back
    return fig
end

"""
Simple utility function to label bars in a stacked bar chart
"""
#TODO: Bring this functionality back. Right now, we're working with legends in plots
# `label_bars()` is the original function. `label_bars!()` is the first cut at refactoring but hasn't been tested
# function label_bars(a, Bararray, labels; val_multiplier = 1, fontsize = 8)
#     # for (i,bar) in enumerate(Bararray)
#     #     w, h = bar[0].get_width(), bar[0].get_height()
#     #     x, y = bar[0].get_x(), bar[0].get_y()
#     #     a.text(x+w, y+h/2, @sprintf("%7.3f", pyconvert(Float64, h)*val_multiplier), ha = "left", va = "center", fontsize = fontsize)
#     #     a.text(x-w/2, y+h/2, @sprintf("%7s", labels[i]), ha = "right", va = "center", fontsize = fontsize)
#     # end
# end
# """
# Utility function to label bars in a stacked bar chart using Plots.jl
# """
# function label_bars!(xvals, yvals, labels; val_multiplier = 1, fontsize = 8)
#     for (i, y) in enumerate(yvals)
#         y_pos = sum(y[1:end-1]) + y[end] / 2
#         x_pos = xvals[i]
#         annotate!(x_pos, y_pos, text(@sprintf("%.2f", y[end] * val_multiplier), fontsize))
#         annotate!(x_pos, y_pos - y[end] / 2, text(labels[i], fontsize))
#     end
# end

"""
Function to plot the comparison of 737-800 weights from TASOPT w 220 pax
"""
function plot737compare(ac::aircraft; weightdetail = true, fracs = false)
    #note for any dev: the layout treatment in plot_details is preferred
    #plot layout
    layout = @layout([A B])
    fig = plot(layout=layout, size=(800, 500), dpi=300, margin=4mm)
    
    #=
    Subplot 1: Mass comparison
    =#
    # Constants
    lbf_to_N = 4.44822

    # 737-800 Weights breakdown
    Wempty737 = 105757.5 * lbf_to_N
    Wpay737 = 47318.9 * lbf_to_N
    Wfuel737 = 49862.7 * lbf_to_N
    WMTO737 = 202939.1 * lbf_to_N
    Wtotadd737 = 13191.0 * lbf_to_N
    Wfuse737 = 44496.3 * lbf_to_N
    Wwing737 = 26043.3 * lbf_to_N
    Whtail737 = 2650.0 * lbf_to_N
    Wvtail737 = 1749.0 * lbf_to_N
    Weng737 = 17628.0 * lbf_to_N

    # ac weights
    parg = ac.parg
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    fuselage = ac.fuselage

    Wpay  = parg[igWpay]
    Wfuel = parg[igWfuel]
    WMTO  = parg[igWMTO]
    Wempty  = parg[igWMTO] - parg[igWfuel] - parg[igWpay]
    Wwing = wing.weight
    Wfuse = ac.fuselage.weight
    Whtail = htail.weight
    Wvtail = vtail.weight
    Weng = parg[igWeng]
        Whpesys = parg[igWMTO] * ac.fuselage.HPE_sys.W
        Wlgnose = ac.landing_gear.nose_gear.weight.W
        Wlgmain = ac.landing_gear.main_gear.weight.W
    Wtotadd = Whpesys + Wlgnose + Wlgmain
    Wftank = parg[igWftank]

    # Convert to fractions or tonnes
    if fracs
        weight_factor = WMTO
        weight_factor737 = WMTO737
    else
        weight_factor = 9.81 * 1000
        weight_factor737 = weight_factor
    end

    # Stacked bar chart
    bar_width = 0.8
    grp_labels = ["737-800"; ac.name]

    if weightdetail
        yvals_737 = [
            Wfuse737,
            Wwing737,
            Whtail737,
            Wvtail737,
            Weng737,
            Wtotadd737,
            Wfuel737,
            Wpay737,
        ]/weight_factor737
        yvals_ac = [
            Wfuse,
            Wwing,
            Whtail,
            Wvtail,
            Weng,
            Wtotadd,
            Wfuel,
            Wpay,
        ]/weight_factor
    
        col_labels = ["Fuselage" "Wing" "Horz.\nTail" "Vert.\nTail" "Engine" "Add'l" "Fuel" "Payload"]
    
        groupedbar!(fig,
            subplot=1,
            hcat(yvals_737, yvals_ac)',
            label = col_labels,
            width = bar_width,
            color = [:black :blue :orange :purple :green :red :cyan :yellow],
            legend = :best,
            grid = true,
            bar_position = :stack
        )
    else
        yvals_737 = [Wempty737, Wfuel737, Wpay737,]/weight_factor737
        yvals_ac = [Wempty,Wfuel,Wpay]/weight_factor

        col_labels = ["Empty" "Fuel" "Payload"]

        groupedbar!(fig,
            subplot=1,
            hcat(yvals_737, yvals_ac)',
            label = col_labels,
            width = bar_width,
            color = [:red :blue :green],
            legend = :best,
            grid = true,
            bar_position=:stack
        )

    end

    xlabel!(fig[1],"Aircraft")
    ylabel!(fig[1], fracs ? "Fraction of MTOW" : "Mass [tonnes]")
    xticks!(fig[1], (1:length(grp_labels)), grp_labels)
    
    
    #=
    Subplot 2: Drag comparison
    =#

    # 737 perf data
    CD737     = 0.03185
    CDi737    = 0.01114
    CDfuse737 = 0.00622
    CDwing737 = 0.00830
    CDhtail737= 0.00261
    CDvtail737= 0.00176
    CDnace737 = 0.00182
    

    # ac perf data
    para = ac.para
    CD      = para[iaCD, ipcruise1,1]
    CDfuse  = para[iaCDfuse, ipcruise1,1]
    CDi     = para[iaCDi, ipcruise1,1]
    CDwing  = para[iaCDwing, ipcruise1,1]
    CDhtail = para[iaCDhtail, ipcruise1,1]
    CDvtail = para[iaCDvtail, ipcruise1,1]
    CDnace  = para[iaCDnace, ipcruise1,1]

    drag_yvals_ac = [CDfuse, CDi, CDwing, CDhtail, CDvtail, CDnace]
    drag_yvals_737 = [CDfuse737, CDi737, CDwing737, CDhtail737, CDvtail737, CDnace737]

    # Convert to fractions if needed
    if fracs
        drag_yvals_ac = drag_yvals_ac/CD
        drag_yvals_737 = drag_yvals_737/CD737
    end

    drag_col_labels = ["Fuselage" "Induced" "Wing" "Horz.\nTail" "Vert.\nTail" "Nacelle"]
    
    groupedbar!(fig,
        subplot=2,
        hcat(drag_yvals_737, drag_yvals_ac)',
        label = drag_col_labels,
        width = bar_width,
        color = [:black :blue :orange :purple :green :red :cyan :yellow],
        legend = :best,
        grid = true,
        bar_position = :stack
    )

    xlabel!(fig[2],"Aircraft")
    ylabel!(fig[2], fracs ? "Fraction of Total Drag Coeff." : "Component Drag Coeff.")
    xticks!(fig[2], (1:length(grp_labels)), grp_labels)

    return fig
end

"""
    MomentShear(ac)

Plot moment and shear diagrams
"""
function MomentShear(ac::aircraft)
    parg = ac.parg
    wing = ac.wing

    co = wing.layout.root_chord
    cs = co*wing.inboard.λ
    ct = co*wing.outboard.λ
  
    bo = wing.layout.root_span
    bs = wing.layout.break_span
    b  = wing.layout.span

    etas = bs/b
    etao = bo/b

    λs = wing.inboard.λ
    λt = wing.outboard.λ

    Ss = wing.outboard.max_shear_load
    Ms = wing.outboard.moment
    So = wing.inboard.max_shear_load
    Mo = wing.inboard.moment
    
    etaRange =[LinRange(etao,etas,20) ; LinRange(etas,1,20)]
    c = zero(etaRange)
    S = zero(etaRange)
    M = zero(etaRange)

    for (i,eta) in enumerate(etaRange)
        if eta<etao
            c[i] = co
        elseif eta<etas
            c[i] = co*( 1 + (λs - 1)*(eta - etao)/(etas - etao))
            S[i] = So
            M[i] = Mo
        else
            c[i] = co*(λs + (λt -λs)*(eta - etas)/(   1 - etas))
            S[i] = Ss*(c[i]/cs)^2
            M[i] = Ms*(c[i]/cs)^3
        end
    end

    # Define layout
    layout = @layout [A; B]
    fig = plot(layout=layout, size=(800, 500), dpi=300, link=:x, margin=4mm)

    # Shear and moment distribution along wings
    # Shear
    plot!(fig, 
        etaRange, S, 
        ylabel="Shear [N]", 
        label="",
        subplot=1,
        title="Shear and Moment Distribution along Wing"
    )
    # Moment
    plot!(fig, 
        etaRange, M, 
        ylabel="Moment [N-m]", 
        xlabel="Normalized Span Coord. [-]", # Set x-label only on the bottom plot
        label="",           # Disable legend for individual plots
        subplot=2           # Target second subplot
    )

    # Add dividing vertical lines for strut and outer section 
    for i in 1:2
        vline!(fig, [etas, etao], ls=:dash, color=:black, 
                subplot=i, label="") # Add dashed lines at etas and etao
    end

    annotate!(fig[2], etao, 0,
            text(" ← Wing start/\n     Wing box end",
            halign=:left, valign=:bottom, 11,))
    annotate!(fig[2], etas, Mo,
            text(" ← Planform Break",
            halign=:left, valign=:top, 11,))

    return fig

end

# Evaluate each altitude independently: failed missions must not seed later trials.
function _payload_range_trial(ac, range, payload, altitudes;
    itermax = 35, initializes_engine = true, printTO = false,
    mission_constraints = mission_constraint_report, constraint_tol = 1e-6)

    best = nothing
    context = MissionConstraintContext(2, ac.parg[igWMTO], true)
    for altitude in altitudes
        trial = deepcopy(ac)
        trial.parm[imRange, 2] = range
        trial.parm[imWpay, 2] = payload
        trial.para[iaalt, ipcruise1, 2] = altitude
        converged = try
            fly_mission!(trial, 2; itermax, initializes_engine, printTO,
                         opt_prescribed_cruise_parameter = "altitude")
        catch err
            # Failed engine gas-flow/temperature solves can return nothing.
            numerical_failure = err isa DomainError || err isa OverflowError ||
                (err isa MethodError && err.f in (iterate, engine.get_thermo) &&
                 any(isnothing, err.args)) ||
                (err isa ErrorException && (startswith(err.msg, "ht5 < h5") ||
                 startswith(err.msg, "Oops, you're searching for a NaN!")))
            numerical_failure || rethrow()
            false
        end
        converged || continue
        fuel = trial.parm[imWfuel, 2]
        # Reject trajectories whose climb/descent leaves no room for cruise.
        cruise_range = trial.para[iaRange, ipcruisen, 2] - trial.para[iaRange, ipcruise1, 2]
        isfinite(fuel) && fuel >= 0 || continue
        constraints = mission_constraints(trial, 2, context)
        feasible = is_feasible(constraints, constraint_tol) && cruise_range >= 0 &&
                   all(>(0), trial.para[iagamV, ipclimb1:ipclimbn, 2])
        # Reserves are a fixed fraction, so minimum fuel load also minimizes burn.
        if best === nothing || (feasible && !best.feasible) ||
           (feasible == best.feasible && fuel < best.fuel)
            best = (; range, payload, fuel, altitude, feasible, constraints,
                    pfei = iszero(payload) ? NaN : trial.parm[imPFEI, 2])
        end
    end
    best === nothing && error("No converged payload-range mission at range $range m, payload $payload N.")
    return best
end

# Bracket any aircraft/mission limit, retaining the feasible side of the bisection.
function _payload_range_limit(mission, payload, short_range, design_range, range_tol)
    low = mission(short_range, payload)
    low.feasible || error("Verification mission violates constraints at payload $payload N: $(low.constraints.raw)")
    upper = max(design_range, 2short_range)
    high = mission(upper, payload)
    for _ = 1:30
        !high.feasible && break
        low = high
        upper *= 1.25
        high = mission(upper, payload)
    end
    !high.feasible || error("Unable to bracket a payload-range constraint limit.")
    while high.range - low.range > range_tol
        mid = mission((low.range + high.range) / 2, payload)
        if mid.feasible
            low = mid
        else
            high = mid
        end
    end
    return low
end

# Sample the level segment in range and each descending segment in payload.
function _payload_range_samples(points, origin_check, mission, limit, Rpts, Ppts)
    sampled = [origin_check]
    for r in range(origin_check.range, points[1].range; length=Rpts)[2:end-1]
        point = mission(r, points[1].payload)
        point.feasible || error("Sampled payload-range mission violates constraints: $(point.constraints.raw)")
        push!(sampled, point)
    end
    push!(sampled, points[1])
    for (left, right) in zip(points[1:end-1], points[2:end])
        for p in range(left.payload, right.payload; length=Ppts)[2:end-1]
            push!(sampled, limit(p))
        end
        push!(sampled, right)
    end
    return unique(sampled)
end

"""
    payload_range_data(ac; altitude_ft = 28000:1000:50000, range_tol = 100.0, ...)

Calculate the maximum-payload, full-fuel and zero-payload range limits, with a
zero-range plotting intercept. Each mission searches the altitude grid for the
lowest fuel burn among feasible missions, preserving design Mach and reserves.
The default constraint policies are the model-defined aircraft/mission reports
(span, fan geometry, fuel, MTOW, BFL, climb and engine temperatures), with a scaled
acceptance tolerance of `1e-6`. Supply other report functions for other aircraft.
Operating limits can truncate the usual MTOW/full-fuel corners.
`range_tol` is in metres; `short_range` is the verification mission range in metres.
The intercept has no mission fuel, PFEI or altitude (NaN); `origin_check` contains
its actual short-mission verification. The input aircraft is not modified.
Set `all_points=true` to sample the boundary: `Rpts` counts missions from
`short_range` to corner 1 (including both ends); `Ppts` counts payload levels
on each descending segment (including both corners). The zero-range intercept
is additional, and shared corners appear only once.
"""
function payload_range_data(ac_og::TASOPT.aircraft;
    all_points::Bool = false, Rpts::Int = 10, Ppts::Int = 10,
    altitude_ft = 28000:1000:50000, range_tol::Real = 100.0,
    short_range::Real = 500 * 1852.0, itermax::Int = 35,
    initializes_engine::Bool = true, printTO::Bool = false, Ldebug::Bool = false,
    aircraft_constraints = aircraft_constraint_report,
    mission_constraints = mission_constraint_report, constraint_tol::Real = 1e-6)

    ac_og.is_sized[1] || error("Size the aircraft before calculating its payload-range diagram.")
    Rpts >= 2 && Ppts >= 2 || throw(ArgumentError("Rpts and Ppts must each be at least 2."))
    isfinite(range_tol) && range_tol > 0 || throw(ArgumentError("range_tol must be positive and finite."))
    isfinite(short_range) && short_range > 0 || throw(ArgumentError("short_range must be positive and finite."))
    altitudes = convertDist.(Float64.(collect(altitude_ft)), "ft", "m")
    !isempty(altitudes) && all(h -> isfinite(h) && h > 0, altitudes) ||
        throw(ArgumentError("altitude_ft must contain positive, finite altitudes."))

    aircraft_report = aircraft_constraints(ac_og)
    is_feasible(aircraft_report, constraint_tol) ||
        error("Aircraft violates shared constraints: $(aircraft_report.raw)")

    ac = deepcopy(ac_og)
    ac.parm = repeat(ac.parm[:, 1:1], 1, 2)
    ac.para = repeat(ac.para[:, :, 1:1], 1, 1, 2)
    ac.pare = repeat(ac.pare[:, :, 1:1], 1, 1, 2)
    for HX in ac.engine.heat_exchangers
        HX.HXgas_mission = repeat(HX.HXgas_mission[:, 1:1], 1, 2)
    end

    Wmax, Fuelmax = ac.parg[igWMTO], ac.parg[igWfmax]
    Wempty = Wmax - ac.parg[igWfuel] - ac.parg[igWpay]
    max_payload = ac.parg[igWpaymax]
    all(isfinite, (Wmax, Fuelmax, Wempty, max_payload)) &&
        Wmax > Wempty >= 0 && Fuelmax > 0 && max_payload > 0 ||
        throw(ArgumentError("Aircraft must have positive usable weight and fuel capacity."))
    # If the intersection lies outside the payload interval, two corners coincide.
    full_fuel_payload = clamp(Wmax - Wempty - Fuelmax, 0, max_payload)
    mission(r, p) = _payload_range_trial(ac, r, p, altitudes;
                                        itermax, initializes_engine, printTO, mission_constraints, constraint_tol)
    origin_check = mission(short_range, max_payload)
    limit(p) = _payload_range_limit(mission, p, short_range, ac.parm[imRange, 1], range_tol)
    points = map(unique([max_payload, full_fuel_payload, 0.0])) do payload
        point = limit(payload)
        Ldebug && println("Payload-range point: ", point)
        point
    end
    if all_points
        points = _payload_range_samples(points, origin_check, mission, limit, Rpts, Ppts)
    end
    origin = merge(origin_check, (; range = 0.0, fuel = NaN, altitude = NaN, pfei = NaN))
    design_report = mission_constraints(ac, 1, MissionConstraintContext(1, Wmax, false))
    pushfirst!(points, origin)
    return (; range = getproperty.(points, :range), payload = getproperty.(points, :payload),
            fuel = getproperty.(points, :fuel), pfei = getproperty.(points, :pfei),
            altitude = getproperty.(points, :altitude),
            WTO = [Wempty + p.payload + p.fuel for p in points],
            constraints = getproperty.(points, :constraints), aircraft_report, design_report,
            design_feasible = is_feasible(design_report, constraint_tol),
            Wempty, Wmax, Fuelmax, origin_check)
end

"""
    PayloadRange(ac; plots_OEW = false, all_points = false, Rpts = 10, Ppts = 10, ...)

Plot the characteristic payload-range points, marking the original sizing mission
separately and labeling it when outside the constraint tolerance. Straight lines
connect the solved points. Set
`plots_OEW=true` to include empty weight on the payload axis (false: payload only).
The default marks only the four characteristic points and the design point.
Set `all_points=true` to solve and mark additional boundary points, with `Rpts`
range samples on the level segment and `Ppts` payload levels per descending segment.
Other keywords are passed
to [`payload_range_data`](@ref).
"""
function PayloadRange(ac::TASOPT.aircraft; plots_OEW::Bool = false,
                      filename::String = "", kwargs...)
    pr = payload_range_data(ac; kwargs...)
    return _payload_range_plot(ac, pr; plots_OEW, filename)
end

function _payload_range_plot(ac, pr; plots_OEW = false, filename = "")
    ranges = convertDist.(pr.range, "m", "nmi")
    offset = plots_OEW ? pr.Wempty : 0.0
    design_range = convertDist(ac.parm[imRange, 1], "m", "nmi")
    design_payload = (ac.parm[imWpay, 1] + offset) / (gee * 1000)
    ylabel = plots_OEW ? "OEW + Payload Weight (tonnes)" : "Payload Weight (tonnes)"

    fig = plot(ranges, (pr.payload .+ offset) ./ (gee * 1000);
               xlabel = "Range (nmi)", ylabel, title = "Payload-Range Diagram: " * ac.name,
               lw = 2, marker = :circle, color = :blue, label = "Feasible boundary points",
               size = (800, 500), dpi = 300, margin = 4mm, legend = :bottomleft)
    scatter!(fig, [design_range], [design_payload]; marker = :star5, ms = 8,
             color = pr.design_feasible ? :blue : :orange,
             label = pr.design_feasible ? "Design point" : "Design point (outside tolerance)")
    isempty(filename) || savefig(fig, filename)
    return fig
end

"""
    DragPolar(ac; CL_range = 0.2:0.05:0.8, 
              show_drag_components=false, show_airfoil_data=false, 
              title=nothing, legend=true, print_results=false)

Generates drag polar plots for a given aircraft model by sweeping over a range of lift coefficients and computing aerodynamic performance metrics.

This function calls [`aeroperf_sweep`](@ref) to evaluate drag, lift-to-drag ratio, and component breakdowns across `CL_range`.  
It produces two side-by-side plots:  
- **CL vs CD** (with optional drag component breakdowns).  
- **CL vs L/D** (with optional airfoil section data).  

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac`: Aircraft model object.
    - `CL_range`: Range of lift coefficients to sweep (default: `0.2:0.05:0.8`).
    - `Mach`: Mach number at which points are evaluated. If none specified, defaults to value at specified ip, imission (optl.).
    - `show_drag_components::Bool`: If `true`, overlays induced drag and component drag contributions (`CDi`, `CDwing`, `CDfuse`, `CDhtail`, `CDvtail`, `CDother`) on the CL–CD plot (default: `false`).
    - `show_airfoil_data::Bool`: If `true`, overlays airfoil section performance (`clpss`, `cdss`) on both plots (default: `false`).
    - `title::Union{String,Nothing}`: Custom plot title. If `nothing`, uses a default title with aircraft name (default: `nothing`).
    - `legend::Bool`: Show legend in plots (default: `true`).
    - `print_results::Bool`: Print detailed output from [`aeroperf_sweep`](@ref) (default: `false`).

    **Outputs:**
    - Returns a combined `Plots.Plot` object with two subplots:
        - **Left subplot:** CL vs CD (with optional drag components and airfoil data).
        - **Right subplot:** CL vs L/D (with optional airfoil data).

    Sample usage:

        ```julia
        f = DragPolar(ac; CL_range=0.2:0.05:0.8, 
                         show_drag_components=true, 
                         show_airfoil_data=false)
        display(f)
        ```

See also: [`aeroperf_sweep`](@ref), [`TASOPT.balance_aircraft!`](@ref), [`TASOPT.aerodynamics.aircraft_drag!`](@ref).

"""
function DragPolar(ac; CL_range = 0.2:0.05:0.8, Mach=nothing,
    ip = ipcruise1, imission = 1,
    show_drag_components=false, show_airfoil_data=false, 
    title=nothing, legend=true, print_results = false)

    #get results from aeroperf_sweep
    results = aeroperf_sweep(ac, CL_range; 
                            imission=imission, ip=ip, 
                            #defaults: rfuel=1, rpay=1, ξpay=0.5,
                            Mach=Mach,
                            print_results = print_results)
    #get airfoil database limits
    #TODO: un-hardcode the limits?
    cl_lims = [0.4, 0.9]
    #TODO: warn if Mach outside of limits?


  # == Plot CL vs CD and CL vs CDi, CDhtail, CDwing 
    p1 = plot(results.CDs, results.CLs, 
    xlabel= show_airfoil_data ? "\$C_D,  c_d\$" : "\$C_D\$", 
    ylabel=show_airfoil_data ? "\$C_L,  c_l\$" : "\$C_L\$", 
    label="Aircraft", 
    lw=2, marker=:o,color=:skyblue)
    
    if show_drag_components
        plot!(p1, results.CDis, results.CLs, label="CL-CDi", lw=2, marker=:diamond)
        plot!(p1, results.CDwings, results.CLs, label="CL-CDp_wing", lw=2, marker=:square)
        plot!(p1, results.CDfuses, results.CLs, label="CL-CDp_fuse", lw=2, marker=:pentagon)
        plot!(p1, results.CDhtails, results.CLs, label="CL-CDp_htail", lw=2, marker=:rtriangle)
        plot!(p1, results.CDvtails, results.CLs, label="CL-CDp_vtail", lw=2, marker=:utriangle)
        plot!(p1, results.CDothers, results.CLs, label="CL-CDp_misc", lw=2, marker=:xcross)
    end

    # Set y-limits with optional airfoil data included
    if show_airfoil_data
        min_CL = min(minimum(results.clpss), minimum(results.CLs), 0)
        max_CL = max(maximum(results.clpss), maximum(results.CLs))
    else
        min_CL = min(0, minimum(results.CLs))
        max_CL = maximum(results.CLs)
    end

    # Apply offsets
    min_CL -= 0.05 * abs(min_CL)
    max_CL *= 1.05

    ylims!(p1, (min_CL, max_CL))

    #secondary title via hacky annotation
    annot_text = @sprintf "Mach = %.2f, ip = %d, im = %d" results.Mach ip imission
    annotate!(xlims(p1)[2]+.003, max_CL + .045, text(annot_text, 12, :gray22))

  # == Create a second plot: CL vs L/D
    p2 = plot(results.LDs, results.CLs, 
        xlabel="\$L/D\$", 
        ylabel= show_airfoil_data ? "\$C_L,  c_l\$" : "\$C_L\$", 
        label="Aircraft", 
        lw=2, marker=:o, color=:skyblue)
    # annotate!(p2, (2.5, 1.4, text("Aircraft L/D", :skyblue, :left, 10)))
    ylims!(p2, (min_CL, max_CL))

  # == On both plots:
    #plot airfoil spanbreak section data if requested
    if show_airfoil_data
        clcdss = results.clpss./results.cdss
        plot!(p1, results.cdss, results.clpss, 
            label="Airfoil at η_s", lw=2, marker=:hexagon, color=:black)
        plot!(p2, clcdss, results.clpss,
            label="Airfoil at η_s", lw=2, marker=:hexagon, color=:black)
        
        #limits of airfoil section
        #plot1
        plot!(p1, [0; maximum(results.cdss)], 
            [cl_lims[1]; cl_lims[1]], color=:black, label=false, linestyle=:dash)
        plot!(p1, [0; maximum(results.cdss)], 
            [cl_lims[2]; cl_lims[2]], color=:black, label=false, linestyle=:dash)
        annotate!(p1, maximum(results.CDs), cl_lims[1] + 0.002, text("Airfoil database limits (\$c_l\$)", :black, 7, :right, :bottom))
        #plot2
        
        plot!(p2, [0; maximum(clcdss)], 
            [cl_lims[1]; cl_lims[1]], color=:black, label=false, linestyle=:dash)
        plot!(p2, [0; maximum(clcdss)], 
            [cl_lims[2]; cl_lims[2]], color=:black, label=false, linestyle=:dash)
    end


    # Mark with dotted y=constant lines where airfoil data is valid 
    # (i.e., only the lowest and highest CLs_inf_swept_wing where clpss is in [0.4, 0.9])
    inds = findall(x -> cl_lims[1] <= x <= cl_lims[2], results.clpss)
    if !isempty(inds)
        for i in (first(inds), last(inds))
            CL = results.CLs[i]
            plot!(p1, [0; maximum(results.CDs)], 
                [CL; CL], color=:deepskyblue, label=false, linestyle=:dash)
            plot!(p2, [0; maximum(results.LDs)], 
                [CL; CL], color=:deepskyblue, label=false, linestyle=:dash)
        end
        annotate!(p1, maximum(results.CDs), results.CLs[first(inds)] + 0.002, text("Airfoil database limits (\$C_L\$)", :gray, 7, :right, :bottom))
    end

    # Combine the two plots side by side with custom spacing
    l = @layout [a{0.6w} b{0.4w}]
    f = plot(p1, p2, layout=l, size=(600,500))

    default_title = @sprintf "TASOPT.jl Aeroperf. Sweep: \n %s" ac.name

    plot!(f,
        suptitle=isnothing(title) ? default_title : title,
        bottom_margin=2Plots.mm,
        left_margin=3Plots.mm,
        top_margin=12Plots.mm,
        right_margin=3Plots.mm,
        legend= show_drag_components || show_airfoil_data, #only show legend if we're having many lines
        )

    return f
end


import .aerodynamics: plot_airf
"""
    plot_airf(ac::TASOPT.aircraft)

Convenience for `plot_airf(airf::airfoil)`
"""
function plot_airf(ac::aircraft)
    airfoil = ac.wing.airsection
    return plot_airf(airfoil)
end

