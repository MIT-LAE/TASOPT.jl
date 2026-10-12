

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
        
        (; xhbox, boh, bh, coh, cth) = _htail_geometry(ac)   # shared with side_view / front_view
        sweeph  = htail.layout.sweep

        dx = xhbox
        tanLh = tan(sweeph*π/180.0)

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
  
        # Seat start x-position (pushed behind a front tank), see `_cabin_x_extent`
        xseats0, _, _ = _cabin_x_extent(ac)

        # Fuel tank(s): outline of the upper half, shared with side_view
        ntank = 8
        Rtank = Rfuse - 0.1 # Account for clearance_fuse
        l = max(parg[iglftankin], parg[iglftank])
        xtanks = _fuse_tank_centers(ac)
        nftanks = length(xtanks) #Number of fuel tanks

        xt = zeros(nftanks, ntank*2 )
        yt = zeros(nftanks, ntank*2 )
        for m = 1:nftanks
            xt[m, :], yt[m, :] = _fuse_tank_half_outline(xtanks[m], l, Rtank; n = ntank)
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
        pax = parg[igWpay]/parm[imWperpax]   # seat the design passenger count
        _, wcabin, seats_per_row = _deck_cabin_layout(ac)[1]   # floor angle -> cabin width -> seats abreast

        xseats = _deck_seat_rows(ac, xseats0; pax = pax)[1]
        rows = length(xseats)

        println("Seats per row = $seats_per_row, Total rows = $rows")
        yseats = arrange_seats(seats_per_row, wcabin, fuselage.cabin.seat_width, fuselage.cabin.aisle_halfwidth)
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
    (; xvbox, cov, tc_v) = _vtail_geometry(ac)   # shared with side_view / front_view
    xvt = [-0.4, -0.15, 0.2, 0.6].*cov .+ xvbox
    yvt = hcat([0.0 ones(length(xvt) - 2)' .*(cov*tc_v/2) 0.0])[:]
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

    # Plot Engines (same engine layout as side_view / front_view, see `_engine_layout`)
    @assert iseven(Int(parg[igneng])) "Number of engines must be even"
    eng_layout = _engine_layout(ac)
    D     = eng_layout.D
    lnace = eng_layout.lnace

    # Mirror the one-wing engine list to the other wing
    yi_w = [y for (y, _, _) in eng_layout.engines]
    xi_w = [x for (_, x, _) in eng_layout.engines]
    yi = vcat(-reverse(yi_w), yi_w)
    xi = vcat( reverse(xi_w), xi_w)

    # Plot each engine separately
    for i in eachindex(yi)
        xengine = [xi[i], xi[i], xi[i] + lnace, xi[i] + lnace, xi[i]]
        yengine = [yi[i] - D/2, yi[i] + D/2, yi[i] + D/3, yi[i] - D/3, yi[i] - D/2]
        plot!(plot_obj, xengine, yengine, color = :red, linewidth = 1.5, label = "", z_order = :front)
    end

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

"""
    PayloadRange(ac_og; Rpts, Ppts, plots_OEW, filename, itermax, initializes_engine, Ldebug)

Function to plot a payload range diagram for an aircraft

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac_og::aircraft`: Aircraft structure for payload range diagram.
    - `Rpts::Int64`: Density of ranges to be plot (Optional).
    - `Ppts::Int64`: Density of payloads to be plot (Optional).
    - `filename::String`: filename string for the plot to be stored (Optional).
    - `plots_OEW::Bool`: Whether to plot OEW+Payload (true) on the y-axis or just Payload (false, default) (Optional).
    - `itermax::Int64`: Max Iterations for fly_mission! loop (Optional).
    - `initializes_engine::Bool`: Use design case as initial guess for engine state if true (Optional)
    - `specifying_cruise::String`: option for whether cruise altitude or lift coefficient is specified. Options are "altitude" or "lift_coefficient"
    - `Ldebug::Bool`: verbosity flag. false by default, hiding outputs as PR sweeps progress (Optional).
"""
function PayloadRange(ac_og::TASOPT.aircraft; 
    Rpts::Integer = 20, Ppts::Integer = 21, plots_OEW::Bool = false,
    filename::String = "", 
    itermax::Int64 = 35, initializes_engine::Bool = true, opt_prescribed_cruise_parameter = "CL",
    Ldebug::Bool = false)

    if !ac_og.is_sized[1]
        error("Aircraft $(ac_og.name) not sized. Please size aircraft before calling `PayloadRange()`.")
    end

    #Duplicate design mission as second aircraft, which will be modified
    parm = cat(ac_og.parm[:,1], ac_og.parm[:,1], dims=2)
    pare = cat(ac_og.pare[:,:,1], ac_og.pare[:,:,1], dims=3)
    para = cat(ac_og.para[:,:,1], ac_og.para[:,:,1], dims=3)
    ac = aircraft(ac_og.name, ac_og.description,
    ac_og.options, ac_og.parg, parm, para, pare, [true], 
    ac_og.fuselage, ac_og.fuse_tank, ac_og.wing, ac_og.htail, ac_og.vtail, ac_og.engine, ac_og.landing_gear)

    for HX in ac.engine.heat_exchangers
        HX.HXgas_mission = cat(HX.HXgas_mission[:,1], HX.HXgas_mission[:,1], dims=2)
    end
    #Extract aircraft parameters
    Wmax = ac.parg[igWMTO]
    Fuelmax = ac.parg[igWfmax]
    Wempty = ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay]
    PFEI = 0.0

    RangesToPlot = []
    PayloadToPlot = []
    PFEIsToPlot = []
    maxPay = ac.parg[igWpaymax]
    Wpax = ac.parm[imWperpax, 1]
    sizingWpay = ac.parm[imWpay, 1] #sizing payload

    RangeArray =  ac.parm[imRange,1] * sort([LinRange(0.1,2.0,Rpts-1); 1.0]) #This ensures design range is always shown
    tolweight = 1.0 #One newton tolerance for weight checks

    for Range = RangeArray
        if maxPay == 0
            break
        else
            #This ensures design payload is always included
            Payloads = reverse(sort([(maxPay) * LinRange(0, 1, Ppts - 1); sizingWpay]))
        end
        ac.parm[imRange,2] = Range
        for mWpay = Payloads
            if Ldebug println("Checking for Range (nmi): ", convertDist(Range, "m", "nmi"), " and Pax = ", mWpay/convertForce(Wpax, "N", "lbf")) end
            #reset the aircraft to the design altitude and cruise CL since one is changed during fly_mission!() per opt_prescribed_cruise_parameter
            ac.para[iaalt,ipcruise1,2] = ac.para[iaalt,ipcruise1,1]
            ac.para[iaCL,ipcruise1,2] = ac.para[iaCL,ipcruise1,1]
            
            ac.parm[imWpay,2] = mWpay
            try
                fly_mission!(ac, 2; itermax = itermax, initializes_engine = initializes_engine, opt_prescribed_cruise_parameter = opt_prescribed_cruise_parameter)
                # fly_mission! success: store maxPay, break loop
                mWfuel = ac.parm[imWfuel,2]
                WTO = Wempty + mWpay + mWfuel

                # if weights are negative or above their max, point is infeasible
                if ((WTO - Wmax) > tolweight) || ((mWfuel - Fuelmax) > tolweight) || (WTO < 0.0) || (mWfuel < 0.0)
                    WTO = 0.0
                    mWfuel = 0.0

                    if Ldebug println("Max out error!") end

                    if mWpay == 0
                        if Ldebug println("Payload 0 and no convergence found") end
                        
                        maxPay = 0
                    end
                else
                    maxPay = mWpay
                    PFEI = ac.parm[imPFEI, 2]
                    if Ldebug println("Converged - moving to next range...") end
                    
                    break
                end     
            catch
                if Ldebug println("Not Converged - moving to lower payload...") end
            end
        end
        append!(RangesToPlot, Range)
        append!(PFEIsToPlot, PFEI)
        if plots_OEW
            append!(PayloadToPlot, maxPay+Wempty)
        else
            append!(PayloadToPlot, maxPay)
        end
    end

    # Convert values for plotting
    ranges_nmi = convertDist.(RangesToPlot, "m", "nmi")
    payload_tons = PayloadToPlot ./ (9.81 * 1000)

    # Standard PR Diagram
    plot1 = plot(ranges_nmi, payload_tons, 
        lw=2,                   # Line width
        line=:solid,            # Line style
        color=:blue,            # Line color
        xlabel="Range (nmi)", 
        ylabel= plots_OEW ? "OEW + Payload Weight (tonnes)" : "Payload Weight (tonnes)",
        title="Payload-Range Diagram: "*string(ac.name), 
        grid=true,              # Enable grid
        dpi = 300,
        margin=4mm,
        legend=false,
        label="")

    # Add the design point to plot1
    design_range = convertDist.(ac_og.parm[imRange,1], "m", "nmi")
    design_payload = (plots_OEW ? ac_og.parm[imWpay, 1] + Wempty : ac_og.parm[imWpay, 1]) / (9.81 * 1000)
    scatter!(plot1, [design_range], [design_payload], 
             color=:blue, marker=:star5, ms=8, label="Design Point",
             legend=:bottomleft)

    # PFEI plot
    plot2 = plot(ranges_nmi, PFEIsToPlot, 
        lw=2,                   # Line width
        line=:solid,            # Line style
        color=:green,            # Line color
        xlabel="Range (nmi)", 
        ylabel="PFEI at max payload (kJ/kg-km)", 
        title="Payload-Range Diagram: "*string(ac.name), 
        grid=true,              # Enable grid
        dpi = 300,
        margin=4mm,
        legend=false)

    # Add the design point to plot2
    design_PFEI = ac_og.parm[imPFEI, 1]
    scatter!(plot2, [design_range], [design_PFEI], color=:green, marker=:star5, ms=8, label="Design Point")

    layout = @layout [A; B]
    fig = plot(plot1, plot2, layout=layout, 
    size=(600, 800), 
    dpi=300, margin=4mm)
        
    if filename != ""
        savefig(fig, filename)
    end

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


# ──────────────────────────────────────────────────────────────────────────────
# Additional aircraft view functions: side_view and front_view
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# Shared helpers (used by both side_view and front_view so the two views draw
# identical seats / people)
# ──────────────────────────────────────────────────────────────────────────────

"""
    _seated_person_geometry(h_seat)

Heights (above the cabin floor, in metres) of the seat / seated-person icon used
by both [`side_view`](@ref) and [`front_view`](@ref).  `h_seat` is the total seat
height (`fuselage.cabin.seat_height`).
"""
function _seated_person_geometry(h_seat)
    cushion    = 0.45                       # seat cushion height above floor [m]
    backh      = max(h_seat - cushion, 0.0) # seat-back height above cushion [m]
    z_shoulder = cushion + 0.75*backh
    head_r     = 0.12*h_seat
    z_head_c   = z_shoulder + 1.05*head_r
    z_arm      = 0.85*z_shoulder
    return (; cushion, backh, h_seat, z_shoulder, head_r, z_head_c,
              z_head_top = z_head_c + head_r, z_arm)
end

"""
    _deck_floor_z(ac)

Vertical position (z, relative to the upper-bubble centreline) of the passenger
floor(s), returned as `Vector{Float64}` (main deck first, then upper deck if
`n_decks == 2`).  The angles stored in `fuselage.cabin` are used when available
(`z = Rfuse*sin(θ)`); otherwise they are recomputed with [`find_floor_angles`](@ref)
so that the plot never silently puts the floor on the centreline.  For a double
decker the upper floor is `floor_distance` above the main floor (as in
`find_floor_angles`).
"""
function _deck_floor_z(ac)
    fus   = ac.fuselage
    R     = fus.layout.radius
    dR    = fus.layout.bubble_lower_downward_shift
    h     = fus.cabin.seat_height
    θ1    = fus.cabin.floor_angle_main
    if fus.n_decks == 2
        if abs(θ1) < 1e-9 && abs(fus.cabin.floor_angle_top) < 1e-9
            @warn "Cabin floor angles not set for double-decker; plotting decks symmetrically about z = 0"
            θ1 = asin(-0.5*fus.cabin.floor_distance/R)
        end
        θ2 = asin(clamp(R*sin(θ1) + fus.cabin.floor_distance, -R, R)/R)
        return [R*sin(θ1), R*sin(θ2)]
    else
        if abs(θ1) < 1e-9
            θ1 = find_floor_angles(false, R, dR, h_seat = h)
        end
        return [R*sin(θ1)]
    end
end

"""
    _deck_cabin_layout(ac)

Per-deck cabin cross-section data `[(z_floor, cabin_width, seats_abreast), …]`
(main deck first).  Uses exactly the same calculation as the sizing code and
[`stickfig`](@ref): floor angle -> [`find_cabin_width`](@ref) (effective width at
floor / seat height) -> `findSeatsAbreast`.  The values stored in
`fuselage.cabin` are *not* relied upon, since they are only populated for some
configurations (they are 0 for single-deck aircraft).
"""
function _deck_cabin_layout(ac)
    fus  = ac.fuselage
    lay  = fus.layout
    cab  = fus.cabin
    R    = lay.radius
    wfb  = lay.bubble_center_y_offset
    nweb = lay.n_webs
    h    = cab.seat_height
    zs   = _deck_floor_z(ac)
    out  = Tuple{Float64,Float64,Int}[]
    for z in zs
        θ = asin(clamp(z/R, -1.0, 1.0))
        w = find_cabin_width(R, wfb, nweb, θ, h)
        n = structures.findSeatsAbreast(w, cab.seat_width, cab.aisle_halfwidth)
        push!(out, (z, w, n))
    end
    return out
end

"""
    _cabin_x_extent(ac)

Longitudinal extent of the passenger cabin, following the convention of
[`stickfig`](@ref): seat rows start at `x_pressure_shell_fwd`, but are pushed
back behind a front (or "both") cryogenic tank (tank end + 1 ft).  A rear tank
truncates the cabin floor at the tank's front end (- 1 ft).

Returns `(x_seat0, x_floor_fwd, x_floor_aft)`.
"""
function _cabin_x_extent(ac)
    parg    = ac.parg
    layout  = ac.fuselage.layout
    options = ac.options
    ftank   = ac.fuse_tank

    x_seat0     = layout.x_pressure_shell_fwd
    x_floor_fwd = layout.x_start_cylinder
    x_floor_aft = layout.x_end_cylinder

    if ftank.tank_count > 0 && !options.has_wing_fuel
        l_ft = max(parg[iglftankin], parg[iglftank])
        if ftank.placement in ("front", "both")
            x_seat0     = parg[igxftank] + l_ft/2 + 1.0*ft_to_m
            x_floor_fwd = x_seat0
        end
        if ftank.placement in ("rear", "both")
            x_floor_aft = parg[igxftankaft] - l_ft/2 - 1.0*ft_to_m
        end
    end
    return x_seat0, x_floor_fwd, x_floor_aft
end


"""
    _deck_seat_rows(ac, x_seat0)

Longitudinal position of every seat row on each deck, `Vector{Vector{Float64}}`
(main deck first), measured from the front of the cabin `x_seat0` (see
[`_cabin_x_extent`](@ref)).  Uses [`place_cabin_seats`](@ref) (includes the extra
legroom at the emergency-exit rows) for `pax` passengers, by default the
aircraft's `fuselage.cabin.exit_limit` (the cabin the fuselage was sized for);
[`stickfig`](@ref) passes the design passenger count instead.

  • Single deck : all `exit_limit` passengers on the one deck.
  • Double deck : `exit_limit` is split between the decks.  The split chosen by
    the cabin optimizer (`optimize_double_decker_cabin`) is not stored in the
    aircraft, so it is recovered here for the stored floor angles: the integer
    split that minimises the longer of the two cabins (ties -> the most even
    split), without re-running the optimizer.
"""
function _deck_seat_rows(ac, x_seat0; pax = ac.fuselage.cabin.exit_limit)
    fus   = ac.fuselage
    cab   = fus.cabin
    decks = _deck_cabin_layout(ac)
    kw    = (seat_pitch = cab.seat_pitch, seat_width = cab.seat_width,
             aisle_halfwidth = cab.aisle_halfwidth, front_seat_offset = cab.front_seat_offset)
    rows(pax, w) = place_cabin_seats(pax, w; kw...)[2]   # row x-positions (1 x nrows)

    if !ac.options.is_doubledecker || length(decks) < 2
        return [collect(rows(pax, decks[1][2])[:]) .+ x_seat0]
    end

    w1, w2 = decks[1][2], decks[2][2]
    best = (Inf, Inf, 1)                                 # (max length, |l1 - l2|, pax on main deck)
    for p1 in 1:(pax - 1)
        l1 = rows(p1, w1)[end]
        l2 = rows(pax - p1, w2)[end]
        cand = (max(l1, l2), abs(l1 - l2), p1)
        cand[1:2] < best[1:2] && (best = cand)
    end
    p1 = best[3]
    return [collect(rows(p1, w1)[:]) .+ x_seat0, collect(rows(pax - p1, w2)[:]) .+ x_seat0]
end

"""
    _fuse_tank_centers(ac)

Longitudinal centres of the in-fuselage fuel tank(s) (`[]` if there are none), by `fuse_tank.placement`
(`"front"`, `"rear"` or `"both"`).  Shared by [`stickfig`](@ref) and [`side_view`](@ref).
"""
function _fuse_tank_centers(ac)
    parg = ac.parg
    ftank = ac.fuse_tank
    ftank.tank_count == 0 && return Float64[]
    return ftank.placement == "front" ? [parg[igxftank]] :
           ftank.placement == "rear"  ? [parg[igxftankaft]] :
                                        [parg[igxftank], parg[igxftankaft]]
end

"""
    _fuse_tank_half_outline(xc, l, Rtank; ARtank = 2.0, n = 20)

Upper half `(x, z)` (`2n` points) of the outline of a tank of length `l` centred at `xc`: a cylinder of
radius `Rtank` with elliptical end caps of aspect ratio `ARtank`.  Mirror about 0 for the lower half.
"""
function _fuse_tank_half_outline(xc, l, Rtank; ARtank = 2.0, n = 20)
    xcyl0 = xc - l/2 + Rtank/ARtank
    xcyl1 = xc + l/2 - Rtank/ARtank
    x = zeros(2n);  z = zeros(2n)
    for i in 1:n
        frac = (i-1)/(n-1)
        x[i]   = xcyl0 - Rtank/ARtank*cos(0.5π*frac)
        x[i+n] = xcyl1 + Rtank/ARtank*sin(0.5π*frac)
    end
    for i in 1:n
        z[i]   = sqrt(max(Rtank^2 * (1 - ((x[i]  -xcyl0)/(Rtank/ARtank))^2), 0.0))
        z[i+n] = sqrt(max(Rtank^2 * (1 - ((x[i+n]-xcyl1)/(Rtank/ARtank))^2), 0.0))
    end
    return x, z
end

"""
    _bubble_centers(ac)

Lateral (y) positions of the bubble centres of the fuselage cross-section.  The
fuselage has `n_webs + 1` bubbles of radius `Rfuse`, spaced `2·wfb`
(`bubble_center_y_offset`) apart and centred on `y = 0`, i.e. the outermost
centre is at `±n_webs·wfb` (same convention as the structural model in `fusew!`).
A single-bubble fuselage (`n_webs = 0`) returns `[0.0]`.
"""
function _bubble_centers(ac)
    lay = ac.fuselage.layout
    n   = lay.n_webs
    return [(k - 1 - n/2) * 2 * lay.bubble_center_y_offset for k in 1:(n + 1)]
end

"""
    _inside_bubble(y, z, yc, R, dR)

`true` if the point `(y, z)` lies strictly inside the bubble centred at `yc`.  The upper half of the
bubble is a semicircle centred at `z = 0`; the lower half is the same semicircle shifted down by `dR`
(`bubble_lower_downward_shift`), joined by vertical sides (`dR = 0` gives a circle).  A small relative
tolerance means points on the boundary are *not* inside.
"""
function _inside_bubble(y, z, yc, R, dR)
    zc = z > 0 ? z : (z < -dR ? z + dR : 0.0)   # distance to the centre line of the stadium
    return (y - yc)^2 + zc^2 < R^2 * (1 - 1e-6)
end

"""
    _bubble_outline(yc, R, dR; n = 240)

Closed polyline `(y, z)` of the outline of one bubble (see [`_inside_bubble`](@ref)), counter-clockwise from
the right-most point, with `n` points per half.
"""
function _bubble_outline(yc, R, dR; n = 240)
    θu = range(0, π, length = n)          # upper half, centred at z = 0
    θl = range(π, 2π, length = n)         # lower half, centred at z = -dR
    y = [yc .+ R .* cos.(θu); yc .+ R .* cos.(θl); yc + R]
    z = [R .* sin.(θu); R .* sin.(θl) .- dR; 0.0]
    return y, z
end

"""
    _nonoverlapping_mask(y, z, k, centers, R, dR)

Logical array (`BitVector`) over the outline points `(y, z)` of bubble number `k`: `true` where the point is
part of the **outer** skin, i.e. **not** inside any of the other bubbles in `centers`; `false` where it is
buried inside a neighbouring bubble (the overlapping part).
"""
function _nonoverlapping_mask(y, z, k, centers, R, dR)
    mask = trues(length(y))
    for (j, yc) in enumerate(centers)
        j == k && continue
        for i in eachindex(y)
            _inside_bubble(y[i], z[i], yc, R, dR) && (mask[i] = false)
        end
    end
    return mask
end

"""
    _bubble_sections(ac; n = 240)

Outline of every bubble of the fuselage cross-section together with the logical array that marks its
non-overlapping part (see [`_nonoverlapping_mask`](@ref)).  Returns a `NamedTuple` with `centers`, `R`, `dR`
and, for each bubble, `y[k]`, `z[k]` and `mask[k]`.  The outline points where the mask changes value are
refined (by bisection) so that the visible arcs end exactly where they meet the neighbouring bubble.
"""
function _bubble_sections(ac; n = 240)
    lay = ac.fuselage.layout
    R   = lay.radius
    dR  = lay.bubble_lower_downward_shift
    centers = _bubble_centers(ac)

    ys = Vector{Vector{Float64}}(); zs = Vector{Vector{Float64}}(); masks = Vector{BitVector}()
    for (k, yc) in enumerate(centers)
        y, z = _bubble_outline(yc, R, dR; n)
        m    = _nonoverlapping_mask(y, z, k, centers, R, dR)
        inside_other(a, b) = any(j -> j != k && _inside_bubble(a, b, centers[j], R, dR), eachindex(centers))

        # insert the exact boundary point wherever the mask changes between consecutive points
        yn = Float64[]; zn = Float64[]
        for i in 1:(length(y) - 1)
            push!(yn, y[i]); push!(zn, z[i])
            if m[i] != m[i+1]
                (ya, za), (yb, zb) = m[i] ? ((y[i], z[i]), (y[i+1], z[i+1])) : ((y[i+1], z[i+1]), (y[i], z[i]))
                for _ in 1:50        # (ya, za) is outside the other bubbles, (yb, zb) inside
                    ym, zm = 0.5 * (ya + yb), 0.5 * (za + zb)
                    inside_other(ym, zm) ? ((yb, zb) = (ym, zm)) : ((ya, za) = (ym, zm))
                end
                push!(yn, ya); push!(zn, za)
            end
        end
        push!(yn, y[end]); push!(zn, z[end])
        push!(ys, yn); push!(zs, zn)
        push!(masks, _nonoverlapping_mask(yn, zn, k, centers, R, dR))
    end
    return (; centers, R, dR, y = ys, z = zs, mask = masks)
end

"""
    _fuselage_tail_end_z(ac)

Upper and lower z of the fuselage at its aft end, as a `NamedTuple` `(; top, bottom)`.

In the side view the tail cone up-sweeps: the crown line stays at the top of the
fuselage (`top = Rfuse`) while the keel line rises to meet it.  A pointed tail
(`tapers_to = "point"`) closes onto the crown line (`bottom = top`); a tail that
tapers to an edge keeps a residual height of `0.4·Rfuse` (as in the top view's
`0.2·` half-width convention), hanging below the crown line.
Used by [`side_view`](@ref) and [`_fuselage_top_z`](@ref).
"""
function _fuselage_tail_end_z(ac)
    layout = ac.fuselage.layout
    R      = layout.radius
    h_end  = compare_strings(layout.opt_tapers_to, "point") ? 0.0 : 0.4*R
    return (; top = R, bottom = R - h_end)
end

"""
    _fuselage_tail_profile(ac, ξ)

Upper and lower z of the fuselage in the tail cone at the normalized station
`ξ = (x - x_end_cylinder)/(x_end - x_end_cylinder) ∈ [0, 1]`, returned as `(z_top, z_bot)`.

Each surface blends from its value on the cylinder to its value at the tail end
with the same exponent `b_tail` that governs the tail radius, `1 - ξ^b_tail`:

    z(ξ) = z_cyl + (z_end - z_cyl) ξ^b_tail

with `z_cyl = +Rfuse` (crown) / `-(Rfuse + ΔR)` (keel) and `z_end` from
[`_fuselage_tail_end_z`](@ref).  Because `z_end,top = Rfuse` the crown line stays
level, whereas the keel rises by `Rfuse + ΔR + z_end,bottom`, so the tail cone
closes at the top of the fuselage instead of at its mid-height.
"""
function _fuselage_tail_profile(ac, ξ)
    layout = ac.fuselage.layout
    R      = layout.radius
    z_bot0 = -(R + layout.bubble_lower_downward_shift)
    zend   = _fuselage_tail_end_z(ac)
    ξb     = ξ^layout.tail_radius
    return (R + (zend.top - R)*ξb, z_bot0 + (zend.bottom - z_bot0)*ξb)
end

"""
    _fuselage_top_z(ac, x)

Height (z) of the top of the fuselage at longitudinal station `x`.  Used to
place the vertical tail root in both views.
"""
function _fuselage_top_z(ac, x)
    layout  = ac.fuselage.layout
    xblend2 = layout.x_end_cylinder
    xend    = layout.x_end
    ξ = clamp((x - xblend2)/(xend - xblend2), 0.0, 1.0)
    return _fuselage_tail_profile(ac, ξ)[1]
end

"""
    _htail_geometry(ac)

Horizontal-tail dimensions shared by both views: reference box position,
root/total span, root & tip chords, z-position and thickness-to-chord.
"""
function _htail_geometry(ac)
    htail = ac.htail
    Sh    = htail.layout.S
    boh   = htail.layout.root_span
    bh    = sqrt(Sh * htail.layout.AR)
    coh   = Sh / (boh + (bh - boh)*0.5*(1 + htail.outboard.λ))
    return (; xhbox = htail.layout.box_x, zh = htail.layout.z, boh, bh, coh,
              cth = coh * htail.outboard.λ,
              tc_h = max(htail.outboard.cross_section.thickness_to_chord, 0.01))
end

"""
    _vtail_geometry(ac)

Vertical-tail dimensions shared by both views, including the z of the root
(on top of the fuselage at the tail box position) and of the tip.
"""
function _vtail_geometry(ac)
    vtail = ac.vtail
    bov   = vtail.layout.root_span
    bv    = vtail.layout.span
    cov   = vtail.layout.root_chord
    xvbox = vtail.layout.box_x
    z_root = _fuselage_top_z(ac, xvbox)
    return (; xvbox, bov, bv, cov, ctv = cov * vtail.outboard.λ,
              sweepv = vtail.layout.sweep,
              tc_v = max(vtail.outboard.cross_section.thickness_to_chord, 0.01),
              z_root, z_tip = z_root + (bv - bov))
end

"""
    _engine_layout(ac)

Wing-mounted engine positions used by **both** [`side_view`](@ref) and
[`front_view`](@ref), returned as a `NamedTuple` with `D`, `lnace` and a vector
`engines` of `(y, x_LE, z_ctr)` for **one** wing (the other wing mirrors `y`):

  • `y`     : spanwise position
  • `x_LE`  : nacelle leading edge, `parg[igxeng]` for the inner engine and swept
              aft by `tan(Λ)·(y − y_inner)` for the others
  • `z_ctr` : nacelle centre height, following the wing dihedral

Spanwise positions:
  • 2-engine aircraft: single engine per wing at the structural break (ηs·b/2)
  • 4-engine aircraft: inner at ηs·b/2, outer at twice that distance from the centreline
                       (TEMPORARY APPROXIMATION, same as in `size_aircraft!` for
                       `wing_weights!`; independent of `y_critical_engines`, which is
                       only used for the OEI yaw moment)
  • Other even count : evenly spaced from just outside the root to 0.75·b/2
"""
function _engine_layout(ac)
    parg = ac.parg
    wing = ac.wing
    b    = wing.layout.span
    bo   = wing.layout.root_span
    zw   = wing.layout.z
    ηs   = wing.layout.ηs
    D    = parg[igdfan]
    neng_wing = Int(round(parg[igneng])) ÷ 2

    tanΓ   = tan(ac.landing_gear.wing_dihedral_angle)
    tanΛ   = tan(wing.layout.sweep * π / 180)
    y_root = bo / 2

    if neng_wing == 1
        ys = [ηs * b / 2]
    elseif neng_wing == 2
        y_eng1 = ηs * b / 2
        ys = [y_eng1, 2*y_eng1]
        if abs(ys[2] - ys[1]) < D
            @warn "Engine layout: the two engines on each wing overlap " *
                  "(y1 = $(round(ys[1], digits=2)) m, y2 = $(round(ys[2], digits=2)) m, " *
                  "nacelle D = $(round(D, digits=2)) m)."
        end
    else # neng_wing => 3
        ys = collect(range(bo/2 + 2.0*D, b/2*0.75, length = neng_wing))
    end

    engines = [(y, parg[igxeng] + tanΛ*(y - ys[1]), zw + tanΓ*(y - y_root) - 0.5*D) for y in ys]
    return (; D, lnace = parg[iglnace], engines)
end


"""
    side_view(ac::aircraft; plot_obj = nothing, label_fs = 14,
              show_grid = false, annotate_dims = true, show_people = false)

Draw a **side (profile) view** of the aircraft: fuselage outline, wing cross
section at the root, horizontal tail, vertical tail, wing-mounted and aft-fan
engines, and – when a cryogenic tank is present – the in-fuselage fuel tank(s).

The x-axis is the longitudinal fuselage axis (same orientation as
[`stickfig`](@ref)).  The z-axis is positive upward.  Wing and htail z-offsets
(`wing.layout.z`, `htail.layout.z`) are respected so that low-wing and T-tail
configurations render correctly.

# Arguments
- `ac::aircraft`       : sized aircraft object.
- `plot_obj`           : existing `Plots.Plot` to draw on; `nothing` (default) creates a new one.
- `label_fs::Int`      : font size for dimension annotations (default 14).
- `show_grid::Bool`    : draw background grid (default `false`).
- `annotate_dims::Bool`: add length / radius callouts (default `true`).
- `show_people::Bool`  : draw the seated passengers (head, torso, arms, legs) on the seats (default `false`).
                         The seats themselves are always drawn.

# Returns
`Plots.Plot` object.
"""
function side_view(ac::aircraft;
        plot_obj      = nothing,
        label_fs      = 14,
        show_grid     = false,
        annotate_dims = true,
        show_people   = false)

    # ── guard ─────────────────────────────────────────────────────────────
    if !ac.is_sized[1]
        @warn "Aircraft ($(ac.name)) must be sized before plotting. Skipping `side_view`…"
        return nothing
    end

    # ── unpack ────────────────────────────────────────────────────────────
    parg, parm, para, pare, options, fuselage, fuse_tank, wing, htail, vtail, engine =
        unpack_ac(ac, 1)

    layout   = fuselage.layout
    Rfuse    = layout.radius
    dRfuse   = layout.bubble_lower_downward_shift
    anose    = layout.nose_radius
    n_decks  = fuselage.n_decks

    xnose    = layout.x_nose
    xend     = layout.x_end
    xblend1  = layout.x_start_cylinder
    xblend2  = layout.x_end_cylinder

    # ── fuselage profile ──────────────────────────────────────────────────
    nnose = 30
    ntail_pts = 30

    # upper profile
    xu_nose = zeros(nnose);  zu_nose = zeros(nnose)
    for i in 1:nnose
        frac  = (i-1)/(nnose-1)
        fracx = cos(0.5π * frac)
        xu_nose[i] = xblend1 + (xnose - xblend1)*fracx
        zu_nose[i] = Rfuse * (1 - fracx^anose)^(1/anose)
    end
    xu_tail = zeros(ntail_pts);  zu_tail = zeros(ntail_pts)
    for i in 1:ntail_pts
        frac = (i-1)/(ntail_pts-1)
        xu_tail[i] = xblend2 + (xend - xblend2)*frac
        zu_tail[i] = _fuselage_tail_profile(ac, frac)[1]
    end
    xu = [xu_nose; [xblend1, xblend2]; xu_tail]
    zu = [zu_nose; [Rfuse,   Rfuse  ]; zu_tail]

    # lower profile
    z_bot = -(Rfuse + dRfuse)
    xl_nose = zeros(nnose);  zl_nose = zeros(nnose)
    for i in 1:nnose
        frac  = (i-1)/(nnose-1)
        fracx = cos(0.5π * frac)
        xl_nose[i] = xblend1 + (xnose - xblend1)*fracx
        zl_nose[i] = z_bot * (1 - fracx^anose)^(1/anose)
    end


    xl_tail = zeros(ntail_pts);  zl_tail = zeros(ntail_pts)
    for i in 1:ntail_pts
        frac = (i-1)/(ntail_pts-1)
        xl_tail[i] = xblend2 + (xend - xblend2)*frac
        zl_tail[i] = _fuselage_tail_profile(ac, frac)[2]
    end
    xl = [xl_nose; [xblend1, xblend2]; xl_tail]
    zl = [zl_nose; [z_bot,   z_bot  ]; zl_tail]

    # closed polygon for fill
    xf_poly = [xu; reverse(xl); [xu[1]]]
    zf_poly = [zu; reverse(zl); [zu[1]]]

    # ── wing root chord cross-section ─────────────────────────────────────
    co      = wing.layout.root_chord
    xax     = wing.layout.spar_box_x_c
    dx_wing = wing.layout.box_x
    zw      = wing.layout.z
    tc_wing = wing.inboard.cross_section.thickness_to_chord
    hw      = tc_wing * co * 0.5

    xwLE = (-xax)*co + dx_wing
    xwTE = (1 - xax)*co + dx_wing
    xw_box = [xwLE, xwLE, xwTE, xwTE, xwLE]
    zw_box = [zw - hw, zw + hw, zw + hw, zw - hw, zw - hw]

    # ── horizontal tail cross-section ─────────────────────────────────────
    (; xhbox, zh, coh, tc_h) = _htail_geometry(ac)   # shared with front_view
    hh    = tc_h * coh * 0.5

    xhLE = -0.40*coh + xhbox
    xhTE = (1 - 0.40)*coh + xhbox
    xh_box = [xhLE, xhLE, xhTE, xhTE, xhLE]
    zh_box = [zh - hh, zh + hh, zh + hh, zh - hh, zh - hh]

    # ── vertical tail planform (side view shows the fin outline) ──────────
    (; xvbox, bov, bv, cov, ctv, sweepv, z_root, z_tip) = _vtail_geometry(ac)   # shared with front_view
    tanLv  = tan(sweepv*π/180)

    xvLE_root = -0.40*cov + xvbox
    xvTE_root = (1 - 0.40)*cov + xvbox

    # leading / trailing edges at root and tip
    xvLE_tip = xvLE_root + tanLv*(bv - bov)
    xvTE_tip = xvLE_tip  + ctv
    zv_root  = z_root     # z of fuselage top at the vtail x-location (on the tail cone)
    zv_tip   = z_tip

    xv_poly = [xvLE_root, xvLE_tip, xvTE_tip, xvTE_root, xvLE_root]
    zv_poly = [zv_root,   zv_tip,   zv_tip,   zv_root,   zv_root  ]

    # ── wing-mounted engine(s) ────────────────────────────────────────────────
    # Same spanwise / vertical positions as `front_view` (via `_engine_layout`);
    # in the side view engines at different y appear at different x due to sweep.
    eng_layout = _engine_layout(ac)
    D_eng  = eng_layout.D
    lnace  = eng_layout.lnace
    # (x_LE, nacelle_length, z_centre) for each engine on one wing
    sv_eng_nacelles = [(x_e, lnace, z_e) for (_, x_e, z_e) in eng_layout.engines]

    # For z_bot_lim calculation keep a combined ze reference
    ze_poly_ref = reduce(vcat, [[z_e - 0.5*D_eng, z_e + 0.5*D_eng] for (_, _, z_e) in sv_eng_nacelles])


    # ── aft-mounted / BLI fan ─────────────────────────────────────────────
    D_aft     = parg[igdaftfan]
    lnace_aft = parg[iglnaceaft]
    x_aft     = parg[igxtshaft]

    # ── cryo fuel tanks ───────────────────────────────────────────────────
    nftanks = fuse_tank.tank_count
    all_tank_polys = Vector{Tuple{Vector{Float64}, Vector{Float64}}}()

    if nftanks > 0 && !options.has_wing_fuel
        l_ftank = max(parg[iglftankin], parg[iglftank])
        Rtank   = Rfuse - 0.1
        xtanks  = _fuse_tank_centers(ac)

        for xc in xtanks
            xt_top, zt_top = _fuse_tank_half_outline(xc, l_ftank, Rtank; n = 20)

            # closed polygon (top then bottom reversed)
            xt_poly = [xt_top; reverse(xt_top); xt_top[1]]
            zt_poly = [zt_top; -reverse(zt_top); zt_top[1]]
            push!(all_tank_polys, (xt_poly, zt_poly))
        end
    end

    # ── cabin floor(s) and seat rows (side view) ────────────────────────
    # Floor z = Rfuse*sin(θ) (θ from the sized fuselage, negative below the
    # bubble centre).  Seats/floor are pushed behind a front cryogenic tank
    # exactly as in `stickfig`.
    x_cab_fwd = layout.x_pressure_shell_fwd
    x_cab_aft = layout.x_pressure_shell_aft
    x_seat0_sv, x_floor_fwd_sv, x_floor_aft_sv = _cabin_x_extent(ac)

    h_seat_sv = fuselage.cabin.seat_height
    pgeom_sv  = _seated_person_geometry(h_seat_sv)
    seat_height_body = pgeom_sv.cushion
    head_r           = pgeom_sv.head_r
    backh_sv         = pgeom_sv.backh

    z_floors_sv  = _deck_floor_z(ac)
    z_floor_main = z_floors_sv[1]            # z of main cabin floor

    seat_pitch_sv        = fuselage.cabin.seat_pitch
    front_seat_offset_sv = fuselage.cabin.front_seat_offset

    decks_sv  = _deck_cabin_layout(ac)
    wcabin_sv = decks_sv[1][2]
    shoulder_w_sv = min(0.42, wcabin_sv * 0.06)  # narrow shoulder in side view

    # Seat rows per deck (exit limit based, as in `stickfig`; see `_deck_seat_rows`)
    xrows_sv = _deck_seat_rows(ac, x_seat0_sv)
    for xr in xrows_sv
        filter!(x -> x <= x_floor_aft_sv, xr)
    end

    sv_floor_polys  = Tuple{Vector{Float64},Vector{Float64}}[]  # (x,z) per floor line
    sv_seat_data    = Tuple{Float64,Float64,Float64}[]          # (x_row, z_floor, shoulder_w)

    # Main deck
    push!(sv_floor_polys, ([x_floor_fwd_sv, x_floor_aft_sv], [z_floor_main, z_floor_main]))
    for xr in xrows_sv[1]
        push!(sv_seat_data, (xr, z_floor_main, shoulder_w_sv))
    end

    # Upper deck (double decker): floor_angle_top is the *sized* angle, i.e.
    # z_top = z_main + floor_distance
    if n_decks == 2
        z_floor_top = z_floors_sv[2]
        push!(sv_floor_polys, ([x_floor_fwd_sv, x_floor_aft_sv], [z_floor_top, z_floor_top]))
        wcabin_top = decks_sv[2][2]
        shoulder_w_top = min(0.42, wcabin_top * 0.06)
        for xr in xrows_sv[2]
            push!(sv_seat_data, (xr, z_floor_top, shoulder_w_top))
        end
    end

    # ── create plot ───────────────────────────────────────────────────────
    if isnothing(plot_obj)
        plot_obj = plot(
            legend       = false,
            size         = (900, 420),
            dpi          = 300,
            grid         = show_grid,
            show         = true,
            aspect_ratio = :equal,
            margin       = 6mm,
        )
    end

    # Fuselage filled silhouette
    plot!(plot_obj, xf_poly, zf_poly;
          fillrange = minimum(zf_poly), fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 2.0, label = "")

    # Cryo tanks (drawn *inside* fuselage, so on top of fuselage fill)
    for (xt_p, zt_p) in all_tank_polys
        plot!(plot_obj, xt_p, zt_p;
              fillrange = minimum(zt_p), fillalpha = 0.40, fillcolor = :cyan,
              linecolor = :steelblue4, linewidth = 1.5, label = "")
    end
    # Fuel label
    if !isempty(all_tank_polys) && !options.has_wing_fuel
        fuelname = options.ifuel == 11 ? "CH₄" :
                   options.ifuel == 40 ? "LH₂" : ""
        for xc in _fuse_tank_centers(ac)
            annotate!(plot_obj, xc, 0.0,
                      text(fuelname, label_fs - 2, :center, :center, :steelblue4))
        end
    end

    # Vertical tail
    plot!(plot_obj, xv_poly, zv_poly;
          fillrange = zv_poly[1], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.5, label = "")
          
    # Fuselage outline (drawn again on top of tank fill for clean boundary)
    plot!(plot_obj, xu, zu; linecolor = :black, linewidth = 2.0, label = "")
    plot!(plot_obj, xl, zl; linecolor = :black, linewidth = 2.0, label = "")

    # Horizontal tail
    plot!(plot_obj, xh_box, zh_box;
          fillrange = zh_box[1], fillalpha = 0.45, fillcolor = :navy,
          linecolor = :black, linewidth = 1.5, label = "")


    # Aft fan (only when non-trivial)
    if D_aft > 0.01 && lnace_aft > 0.01
        xa_poly = [x_aft, x_aft, x_aft + lnace_aft, x_aft + lnace_aft, x_aft]
        za_poly = [D_aft/8, D_aft/8 + D_aft,
                   D_aft/8 + D_aft*0.75, D_aft/8 + D_aft*0.25,
                   D_aft/8]
        plot!(plot_obj, xa_poly, za_poly;
              linecolor = :red, linewidth = 1.5,
              fillrange = za_poly[1], fillalpha = 0.30, fillcolor = :tomato,
              label = "")
    end

    # ── cabin floor lines ─────────────────────────────────────────────────
    for (xfl, zfl) in sv_floor_polys
        plot!(plot_obj, xfl, zfl;
              linecolor = :black, linewidth = 1.8,
              linestyle = :solid, label = "")
    end

    # ── seats (always) and stick figures (if `show_people`): one per row per deck ──
    # Each figure: seat box + torso + head + arms + legs  (side profile)
    for (xr, zf, sw) in sv_seat_data
        # seat box: legroom is in *front* (lower x) of the seat
        seat_d  = seat_pitch_sv * 0.50   # seat depth (fore-aft)
        #seat_h  = h_seat_sv              # cushion height
        seat_h = seat_height_body        # cushion height
        x_seat_fwd = xr - seat_d         # front edge of seat
        x_seat_aft = xr                  # back of seat = row reference
        z_seat_top = zf + seat_h         # top of cushion
        backh = backh_sv
        z_shoulder = zf + pgeom_sv.z_shoulder   # shared with front_view

        # Seat box (cushion)
        plot!(plot_obj,
              [x_seat_fwd, x_seat_aft, x_seat_aft, x_seat_fwd, x_seat_fwd],
              [zf,         zf,         z_seat_top, z_seat_top, zf        ];
              linecolor = :black, linewidth = 1.0, label = "")

        # Seat back (vertical at aft edge)
        plot!(plot_obj,
              [x_seat_aft, x_seat_aft],
              [z_seat_top, z_seat_top + backh];
              linecolor = :black, linewidth = 1.0, label = "")

        show_people || continue   # seated people are optional (`show_people` keyword)

        # Person: torso (from seat top up)
        plot!(plot_obj,
              [xr - sw*0.3, xr - sw*0.3],
              [z_seat_top,  z_shoulder];
              linecolor = :steelblue, linewidth = 1.4, label = "")

        # Head circle
        z_head_c = zf + pgeom_sv.z_head_c       # shared with front_view
        θh = LinRange(0, 2π, 24)
        plot!(plot_obj,
              xr .- sw*0.3 .+ head_r .* cos.(θh),
              z_head_c .+ head_r .* sin.(θh);
              linecolor = :steelblue, linewidth = 1.0, label = "")

        # Arms (horizontal from mid-torso)
        z_arm = zf + pgeom_sv.z_arm             # shared with front_view (relative to the deck floor)
        plot!(plot_obj,
              [xr - sw*0.3 - sw*0.9, xr - sw*0.3 + sw*0.2],
              [z_arm, z_arm];
              linecolor = :steelblue, linewidth = 1.0, label = "")

        # Legs (from front of seat down to floor)
        plot!(plot_obj,
              [x_seat_fwd + seat_d*0.2, x_seat_fwd - seat_pitch_sv*0.15],
              [zf + seat_h*0.1, zf];
              linecolor = :steelblue, linewidth = 1.0, label = "")
    end

    # Wing root box
    plot!(plot_obj, xw_box, zw_box;
          fillrange = zw_box[1], fillalpha = 1, fillcolor = :royalblue,
          linecolor = :black, linewidth = 1.5, label = "")

    # Wing-mounted engine(s)
    for (xi_e, ln_e, z_e) in sv_eng_nacelles
        xe_poly = [xi_e,            xi_e,
                   xi_e + ln_e,     xi_e + ln_e, xi_e]
        ze_poly = [z_e - 0.5*D_eng, z_e + 0.5*D_eng,
                   z_e + D_eng*0.25, z_e - D_eng*0.25,
                   z_e - 0.5*D_eng]
        plot!(plot_obj, xe_poly, ze_poly;
              fillrange = ze_poly[1], fillalpha = 0.35, fillcolor = :tomato,
              linecolor = :red, linewidth = 1.5, label = "")
    end

    # ── dimension callouts ────────────────────────────────────────────────
    if annotate_dims
        z_ann = 1.1*maximum(zv_poly)
        plot!(plot_obj, [xnose, xend], [z_ann, z_ann]; lw = 1.5, color = :black)
        scatter!(plot_obj, [xnose, xend], [z_ann, z_ann];
                 yerror = abs(z_ann)*0.07, markersize = 0)
        annotate!(plot_obj, (xnose + xend)*0.5, z_ann*1.25,
                  text("ℓ = $(round(xend - xnose, digits=1)) m",
                       label_fs, :center, :top))

        # radius callout
        #xmid = (xblend1 + xblend2)*0.5
        #plot!(plot_obj, [xmid, xmid], [0.0, Rfuse];
        #      lw = 1.2, color = :dimgray, linestyle = :dash)
        #annotate!(plot_obj, xmid + 1.0, Rfuse*0.5,
        #          text("R=$(round(Rfuse, digits=2)) m",
        #               label_fs - 2, :left, :center, color = :dimgray))
    end

    # ── labels / limits ───────────────────────────────────────────────────
    xlabel!(plot_obj, "x [m]"; fontsize = label_fs)
    ylabel!(plot_obj, "z [m]"; fontsize = label_fs)
    title!(plot_obj, "Side View — $(ac.name)"; fontsize = label_fs)

    z_top = max(zv_tip * 1.15, Rfuse * 1.5)
    z_bot_lim = min((Rfuse + dRfuse) * -1.8, minimum(ze_poly_ref)*1.1)
    ylims!(plot_obj, z_bot_lim, z_top)
    xlims!(plot_obj, xnose - 2.0, xend + 2.0)

    return plot_obj
end


"""
    front_view(ac::aircraft; plot_obj = nothing, label_fs = 14,
               show_grid = false, annotate_dims = true, show_people = false)

Draw a **front (head-on) view** of the aircraft: fuselage cross-section circle,
wing silhouette (semi-span projected onto the y–z plane), horizontal tail,
vertical tail, and engine nacelle circles.

The y-axis is lateral (positive starboard); z is positive upward.

# Arguments
- `ac::aircraft`       : sized aircraft object.
- `plot_obj`           : existing `Plots.Plot` to draw on; `nothing` (default) creates a new one.
- `label_fs::Int`      : font size for annotation text (default 14).
- `show_grid::Bool`    : draw background grid (default `false`).
- `annotate_dims::Bool`: add span / radius callouts (default `true`).
- `show_people::Bool`  : draw the seated passengers (head, torso, arms) on the seats (default `false`).
                         The seats themselves are always drawn.

# Returns
`Plots.Plot` object.
"""
function front_view(ac::aircraft;
        plot_obj      = nothing,
        label_fs      = 14,
        show_grid     = false,
        annotate_dims = true,
        show_people   = false)

    # ── guard ─────────────────────────────────────────────────────────────
    if !ac.is_sized[1]
        @warn "Aircraft ($(ac.name)) must be sized before plotting. Skipping `front_view`…"
        return nothing
    end

    # ── unpack ────────────────────────────────────────────────────────────
    parg, parm, para, pare, options, fuselage, fuse_tank, wing, htail, vtail, engine =
        unpack_ac(ac, 1)

    layout  = fuselage.layout
    Rfuse   = layout.radius
    dRfuse  = layout.bubble_lower_downward_shift
    wfb     = layout.bubble_center_y_offset
    nfweb   = layout.n_webs
    n_decks = fuselage.n_decks

    # ── dihedral angle (from landing gear data) ────────────────────────────
    Γ_wing = ac.landing_gear.wing_dihedral_angle   # [rad], positive = tips up
    tanΓ   = tan(Γ_wing)

    # ── wing semi-span silhouette ──────────────────────────────────────────
    b    = wing.layout.span
    bo   = wing.layout.root_span
    bs   = wing.layout.break_span
    zw   = wing.layout.z
    co   = wing.layout.root_chord
    λs   = wing.inboard.λ
    λt   = wing.outboard.λ
    tc_o = wing.inboard.cross_section.thickness_to_chord
    tc_s = wing.outboard.cross_section.thickness_to_chord

    hw_root  = tc_o * co * 0.5
    hw_break = tc_o * co * λs * 0.5
    hw_tip   = tc_s * co * λt * 0.5

    y_root   = bo / 2
    y_break  = bs / 2
    y_tip    = b  / 2

    # Dihedral z-rise at each span station (measured from root attachment)
    dz_root  = tanΓ * (y_root  - y_root)   # = 0
    dz_break = tanΓ * (y_break - y_root)
    dz_tip   = tanΓ * (y_tip   - y_root)

    # Trapezoid: upper then lower edge from root → break → tip, with dihedral
    yw_r = [y_root,  y_break,  y_tip,     y_tip,      y_break,     y_root  ]
    zw_r = [zw + dz_root  + hw_root,
            zw + dz_break + hw_break,
            zw + dz_tip   + hw_tip,
            zw + dz_tip   - hw_tip,
            zw + dz_break - hw_break,
            zw + dz_root  - hw_root]
    yw_l = -yw_r
    zw_l =  zw_r

    # ── horizontal tail ───────────────────────────────────────────────────
    (; zh, boh, bh, coh, cth, tc_h) = _htail_geometry(ac)   # shared with side_view

    hh_root = tc_h * coh * 0.5
    hh_tip  = tc_h * cth * 0.5

    yhr = [boh/2, bh/2,   bh/2,   boh/2   ]
    zhr = [zh + hh_root, zh + hh_tip, zh - hh_tip, zh - hh_root]
    yhl = -yhr
    zhl =  zhr

    # ── vertical tail (symmetric about y=0) ──────────────────────────────
    (; cov, ctv, tc_v, z_root, z_tip) = _vtail_geometry(ac)   # shared with side_view
    z_fuse_vtail = z_root
    z_vtip       = z_tip

    hw_v_root = tc_v * cov * 0.5
    hw_v_tip  = tc_v * ctv * 0.5

    # Trapezoid (root thickness → tip thickness), closed outline for the white fill
    yv_fill = [ hw_v_root,  hw_v_tip, -hw_v_tip, -hw_v_root]
    zv_fill = [z_fuse_vtail, z_vtip,   z_vtip,    z_fuse_vtail]
    # Outline without the root edge (the root sits on the fuselage crown)
    yv_poly = yv_fill
    zv_poly = zv_fill

    # ── engines ───────────────────────────────────────────────────────────
    # Same engine layout as `side_view` (see `_engine_layout`)
    eng_layout = _engine_layout(ac)
    D_eng      = eng_layout.D
    yi_engs    = [y_e for (y_e, _, _) in eng_layout.engines]
    ηs         = wing.layout.ηs
    z_eng_ctr  = minimum(z_e for (_, _, z_e) in eng_layout.engines)  # lowest nacelle centre (for plot limits)

    # ── cabin floor(s) in front view ─────────────────────────────────────
    # Same floor heights, seat height and person geometry as the side view.
    z_floors_fv = _deck_floor_z(ac)
    decks_fv    = _deck_cabin_layout(ac)
    z_floor_fv, wcabin_fv, seats_abreast_fv = decks_fv[1]

    h_seat_fv     = fuselage.cabin.seat_height          # total seat height [m]
    pgeom_fv      = _seated_person_geometry(h_seat_fv)  # shared with side_view
    seat_width_fv = fuselage.cabin.seat_width
    aisle_hw_fv   = fuselage.cabin.aisle_halfwidth

    yseats_fv = arrange_seats(seats_abreast_fv, wcabin_fv, seat_width_fv, aisle_hw_fv)

    fv_floor_lines = Tuple{Float64,Float64,Float64,Float64}[]  # (y1,z1,y2,z2)
    fv_seat_rows   = Tuple{Vector{Float64},Float64}[]          # (yseats, z_floor)

    push!(fv_floor_lines, (-wcabin_fv/2, z_floor_fv, wcabin_fv/2, z_floor_fv))
    push!(fv_seat_rows,   (yseats_fv, z_floor_fv))

    if n_decks == 2
        z_floor_top_fv, wcabin_top_fv, seats_abreast_top_fv = decks_fv[2]
        yseats_top_fv = arrange_seats(seats_abreast_top_fv, wcabin_top_fv,
                                      seat_width_fv, aisle_hw_fv)
        push!(fv_floor_lines, (-wcabin_top_fv/2, z_floor_top_fv,
                                wcabin_top_fv/2, z_floor_top_fv))
        push!(fv_seat_rows, (yseats_top_fv, z_floor_top_fv))
    end

    # ── build plot ────────────────────────────────────────────────────────
    if isnothing(plot_obj)
        plot_obj = plot(
            legend       = false,
            size         = (700, 600),
            dpi          = 300,
            grid         = show_grid,
            show         = true,
            aspect_ratio = :equal,
            margin       = 6mm,
        )
    end

    # Right wing
    plot!(plot_obj, yw_r, zw_r;
          fillrange = zw_r[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.5, label = "")

    # Left wing
    plot!(plot_obj, yw_l, zw_l;
          fillrange = zw_l[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.5, label = "")

    # Between wings
    plot!(plot_obj, [yw_r[1], yw_l[1]], [zw_r[1], zw_l[1]];
          fillrange = zw_l[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.5, label = "")
    plot!(plot_obj, [yw_r[end], yw_l[end]], [zw_r[end], zw_l[end]];
          fillrange = zw_l[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.5, label = "")

    # Right horizontal tail
    plot!(plot_obj, yhr, zhr;
          fillrange = zhr[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.2, label = "")

    # Left horizontal tail
    plot!(plot_obj, yhl, zhl;
          fillrange = zhl[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.2, label = "")

    # Between horizontal tails
    plot!(plot_obj, [yhr[1], yhl[1]], [zhr[1], zhl[1]];
          fillrange = zhl[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.5, label = "")
    plot!(plot_obj, [zhl[end], yhl[end]], [zhr[end], zhl[end]];
          fillrange = zw_l[end], fillalpha = 1, fillcolor = :white,
          linecolor = :black, linewidth = 1.5, label = "")

    # Vertical tail
    # The fin is a *filled* closed shape so that it hides everything behind it
    # (e.g. the lines joining the two horizontal-tail halves); the outline is
    # then drawn on top of the fill along the two sides and the tip only.
    plot!(plot_obj, Shape(yv_fill, zv_fill);
          fillcolor = :white, fillalpha = 1, linealpha = 0, label = "")
    plot!(plot_obj, yv_poly, zv_poly;
          linecolor = :black, linewidth = 1.5, label = "")

    # Engine circles
    θ_e = LinRange(0, 2π, 60)
    for y_eng in yi_engs
        # Engine z follows wing dihedral at this spanwise station
        z_eng_y = zw + tanΓ * (y_eng - y_root) - 0.5*D_eng
        # starboard
        yec =  y_eng .+ (D_eng/2) .* cos.(θ_e)
        zec = z_eng_y .+ (D_eng/2) .* sin.(θ_e)
        plot!(plot_obj, yec, zec;
              fillrange = z_eng_y, fillalpha = 0.95, fillcolor = :tomato,
              linecolor = :red, linewidth = 1.5, label = "")
        # port
        plot!(plot_obj, -yec, zec;
              fillrange = z_eng_y, fillalpha = 0.95, fillcolor = :tomato,
              linecolor = :red, linewidth = 1.5, label = "")
    end



    # ── fuselage cross section ─────────────────────────────────────────────
    ncirc  = 120
    θ_circ = LinRange(0, 2π, ncirc)

    if nfweb > 0
        # Multi-bubble fuselage: draw only the non-overlapping (outer) parts of the bubbles
        # (see `_bubble_sections` / `_nonoverlapping_mask`) instead of full overlapping circles.
        secs = _bubble_sections(ac)
        # white fill of the whole cross-section (union of the bubbles) to hide the wing root behind it
        for k in eachindex(secs.y)
            plot!(plot_obj, Shape(secs.y[k], secs.z[k]);
                  fillcolor = :white, fillalpha = 1, linealpha = 0, label = "")
        end
        # outer skin: points of the outline that are inside a neighbouring bubble become NaN (gap in the line)
        for k in eachindex(secs.y)
            plot!(plot_obj, ifelse.(secs.mask[k], secs.y[k], NaN), ifelse.(secs.mask[k], secs.z[k], NaN);
                  linecolor = :black, linewidth = 2.5, label = "")
        end
    else

        # Single-bubble fuselage: account for the lower-bubble downward shift.
        lay = ac.fuselage.layout
        R   = lay.radius
        dR  = lay.bubble_lower_downward_shift

        yf, zf = _bubble_outline(0.0, R, dR; n = length(θ_circ) ÷ 2)

        plot!(plot_obj, Shape(yf, zf);
            fillcolor = :white, fillalpha = 1, linealpha = 0, label = "")

        plot!(plot_obj, yf, zf; fillcolor = :white,
            linecolor = :black, linewidth = 2.5, label = "")
    end

    # ── cabin floor lines (front view) ───────────────────────────────────
    for (y1, z1, y2, z2) in fv_floor_lines
        plot!(plot_obj, [y1, y2], [z1, z2];
              linecolor = :black, linewidth = 1.8, label = "")
    end

    # ── one row of seats (always) and stick figures (if `show_people`), front view ──
    # Front-on: seat cushion + back, torso, head and arms.  All heights come
    # from `_seated_person_geometry`, identical to the side view.
    pg = pgeom_fv
    for (yseats_deck, zf_deck) in fv_seat_rows
        z_seat_top_fv = zf_deck + pg.cushion
        z_shoulder_fv = zf_deck + pg.z_shoulder
        z_head_c_fv   = zf_deck + pg.z_head_c
        z_arm_fv      = zf_deck + pg.z_arm
        sw2           = seat_width_fv * 0.5
        arm_half      = seat_width_fv * 0.45
        θh_fv         = LinRange(0, 2π, 24)

        for ys in yseats_deck
            # Seat cushion box
            plot!(plot_obj,
                  [ys - sw2, ys + sw2, ys + sw2, ys - sw2, ys - sw2],
                  [zf_deck,  zf_deck,  z_seat_top_fv, z_seat_top_fv, zf_deck];
                  linecolor = :black, linewidth = 0.9, label = "")
            # Seat back (full seat height)
            plot!(plot_obj,
                  [ys - sw2, ys - sw2, ys + sw2, ys + sw2],
                  [z_seat_top_fv, zf_deck + pg.h_seat, zf_deck + pg.h_seat, z_seat_top_fv];
                  linecolor = :black, linewidth = 0.9, label = "")
            show_people || continue   # seated people are optional (`show_people` keyword)
            # Torso
            plot!(plot_obj, [ys, ys], [z_seat_top_fv, z_shoulder_fv];
                  linecolor = :steelblue, linewidth = 1.5, label = "")
            # Head circle
            plot!(plot_obj,
                  ys .+ pg.head_r .* cos.(θh_fv),
                  z_head_c_fv .+ pg.head_r .* sin.(θh_fv);
                  linecolor = :steelblue, linewidth = 1.0, label = "")
            # Arms
            plot!(plot_obj, [ys - arm_half, ys + arm_half], [z_arm_fv, z_arm_fv];
                  linecolor = :steelblue, linewidth = 1.0, label = "")
        end
    end

    # ── dimension callouts ────────────────────────────────────────────────
    if annotate_dims
        # Full span arrow
        z_ann = 1.1*maximum(zv_poly)
        plot!(plot_obj, [-b/2, b/2], [z_ann, z_ann]; lw = 1.5, color = :black)
        scatter!(plot_obj, [-b/2, b/2], [z_ann, z_ann];
                 yerror = abs(z_ann)*0.07, markersize = 0)
        annotate!(plot_obj, 0.0, z_ann * 1.28,
                  text("Span = $(round(b, digits=1)) m",
                       label_fs, :center, :top))

        # Fuselage radius
        #plot!(plot_obj, [0.0, Rfuse*cosd(45)], [0.0, Rfuse*cosd(45)];
        #      lw = 1.2, color = :dimgray, linestyle = :dash)
        #annotate!(plot_obj, Rfuse * 2, Rfuse * 1.5,
        #          text("R=$(round(Rfuse, digits=2)) m",
        #               label_fs - 2, :left, :center))
    end

    # ── labels / limits ───────────────────────────────────────────────────
    xlabel!(plot_obj, "y [m]"; fontsize = label_fs)
    ylabel!(plot_obj, "z [m]"; fontsize = label_fs)
    title!(plot_obj, "Front View — $(ac.name)"; fontsize = label_fs)

    y_lim = b/2 * 1.12
    z_top = max(z_vtip * 1.15, Rfuse * 1.5)
    z_bot_lim = min((Rfuse + dRfuse) * -1.9, (z_eng_ctr .- (D_eng/2))*1.1)
    ylims!(plot_obj, z_bot_lim, z_top)
    xlims!(plot_obj, -y_lim, y_lim)

    return plot_obj
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

