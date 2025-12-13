
"""
    plot_wake_system(ac; show_points=true, show_control=true, show_elements=true,
                    show_geometry=true, annotate=false, mirror=true,
                    title="Trefftz-Plane Wake System", P=nothing)

Plots the wake points in the Trefftz plane.
"""
function plot_wake_system(ac; show_points::Bool=true, show_control::Bool=true,
    show_elements::Bool=true, show_geometry::Bool=true, annotate::Bool=false,
    mirror::Bool=true, title::String="Trefftz-Plane Wake System", P=nothing)

    ws = ac.wake_system

    ys = aerodynamics.field_ys(ws)
    zs = aerodynamics.field_zs(ws)
    ycp = aerodynamics.ctrl_ys(ws)
    zcp = aerodynamics.ctrl_zs(ws)

    if P === nothing
        P = plot(xlabel = "y [m]", ylabel = "z [m]", aspect_ratio = :equal,
        legend = :topright, title = title, grid=true,
        framestyle=:box, size=(800, 600), dpi=300)
    end
    # Get indices for wing and tail
    trefftz_config = ac.options.trefftz_config
    i_w1 = aerodynamics.i_first_wing(trefftz_config)
    i_w2 = aerodynamics.i_last_wing(trefftz_config) - 1
    i_t1 = aerodynamics.i_first_tail(trefftz_config)
    i_t2 = aerodynamics.i_last_tail(trefftz_config) - 1
    wing_range = i_w1:i_w2
    tail_range = i_t1:i_t2
    # Plot wake elements as connected lines
    if show_elements
        plot!(P, ys[wing_range], zs[wing_range], lw=1.5, color=:gray60, alpha=0.6,
              label=show_points ? "" : "Wake elements", linestyle=:solid)
        plot!(P, ys[tail_range], zs[tail_range], lw=1.5, color=:gray60, alpha=0.6,
              label=show_points ? "" : "Wake elements", linestyle=:solid)
        if mirror
            plot!(P, -ys[wing_range], zs[wing_range], lw=1.5, color=:gray60, alpha=0.6,
                label=show_points ? "" : "Wake elements", linestyle=:solid)
            plot!(P, -ys[tail_range], zs[tail_range], lw=1.5, color=:gray60, alpha=0.6,
                label=show_points ? "" : "Wake elements", linestyle=:solid)
        end
    end

    # Plot wake vortex points
    if show_points
        scatter!(P, ys, zs, marker=(:circle, 3), color=:dodgerblue,
                label="Wake points", alpha=0.9)
        if mirror
            scatter!(P, -ys, zs, marker=(:circle, 3), color=:dodgerblue,
                    label="", alpha=0.9)
        end
    end

    # Plot control points
    if show_control
        scatter!(P, ycp, zcp, marker=(:x, 4), color=:crimson, label="Control points")
        if mirror
            scatter!(P, -ycp, zcp, marker=(:x, 4), color=:crimson, label="")
        end
    end

    if annotate
        for (i,(y,z)) in enumerate(zip(ys,zs))
            annotate!(P, y, z + 0.3, text(string(i), :black, 7, :center))
            if mirror
                annotate!(P, -y, z + 0.3, text(string(i), :black, 7, :center))
            end
        end
    end

    # Mark wing geometry reference locations
    if show_geometry
        wing = ac.wing
        bo = wing.layout.root_span
        bs = wing.layout.break_span
        b = wing.layout.span

        z_bottom = minimum(zs) - 0.5
        z_top = maximum(zs) + 0.5

        # Root location
        plot!(P, [bo/2, bo/2], [z_bottom, z_top], lw=1.5,
              color=:orange, linestyle=:dash, alpha=0.7, label="")
        if mirror
            plot!(P, [-bo/2, -bo/2], [z_bottom, z_top], lw=1.5,
                  color=:orange, linestyle=:dash, alpha=0.7, label="")
        end

        # Break location
        plot!(P, [bs/2, bs/2], [z_bottom, z_top], lw=1.5,
              color=:green, linestyle=:dash, alpha=0.7, label="")
        if mirror
            plot!(P, [-bs/2, -bs/2], [z_bottom, z_top], lw=1.5,
                  color=:green, linestyle=:dash, alpha=0.7, label="")
        end

        # Tip location
        plot!(P, [b/2, b/2], [z_bottom, z_top], lw=1.5,
              color=:purple, linestyle=:dash, alpha=0.7, label="")
        if mirror
            plot!(P, [-b/2, -b/2], [z_bottom, z_top], lw=1.5,
                  color=:purple, linestyle=:dash, alpha=0.7, label="")
        end

        # Add geometry labels at bottom
        label_z = z_bottom - 0.3
        annotate!(P, bo/2, label_z, text("Root", :orange, 9, :center, :top))
        annotate!(P, bs/2, label_z, text("Break", :green, 9, :center, :top))
        annotate!(P, b/2, label_z, text("Tip", :purple, 9, :center, :top))
        if mirror
            annotate!(P, -bo/2, label_z, text("Root", :orange, 9, :center, :top))
            annotate!(P, -bs/2, label_z, text("Break", :green, 9, :center, :top))
            annotate!(P, -b/2, label_z, text("Tip", :purple, 9, :center, :top))
        end
    end

    return P
end


"""
    plot_lift_distribution(ac; ip=ipcruise1, mirror=true, show_elliptical=false,
                          normalize=false, show_loading=false, P=nothing, annotate=false)

Visualization of the Trefftz-plane lift distribution.
"""
function plot_lift_distribution(ac; ip=ipcruise1, mirror::Bool=true,
    show_elliptical::Bool=false, normalize::Bool=false, show_loading::Bool=false,
    P=nothing, annotate=false)

    ws = ac.wake_system

    # Recompute circulation for the specific mission point. This is kind of clunky and not sure if it makes sense to store the state at every point. 
    para_view = @view ac.para[:, ip, 1]
    aerodynamics.induced_drag!(para_view, ac, ac.options.trefftz_config)

    ycp = aerodynamics.ctrl_ys(ws)
    geom = aerodynamics.TREFFTZ_GEOM
    gc = geom.gc[1:length(ycp)]

    # Get indices for wing and tail
    trefftz_config = ac.options.trefftz_config
    i_w1 = aerodynamics.i_first_wing(trefftz_config)
    i_w2 = aerodynamics.i_last_wing(trefftz_config) - 1
    i_t1 = aerodynamics.i_first_tail(trefftz_config)
    i_t2 = aerodynamics.i_last_tail(trefftz_config) - 1

    # Apply normalization if requested
    if normalize
        max_gc = maximum(abs.(gc))
        gc = gc ./ max_gc
        ylabel_str = "\$\\Gamma / \\Gamma_{max}\$ [-]"
    elseif show_loading
        # Multiply by density and velocity for actual loading
        pare = ac.pared
        rho = pare[ierho0, ip]
        V = pare[ieu0, ip]
        gc = gc .* rho .* V
        ylabel_str = "Lift Loading \$\\rho V \\Gamma\$ [N/m]"
    else
        ylabel_str = "Circulation \$\\Gamma\$ [\$m^2/s\$]"
    end

    if P === nothing
        P = plot(xlabel="Spanwise Location y [m]", ylabel=ylabel_str,
                legend=:topright, title="Trefftz-Plane Lift Distribution",
                grid=false, framestyle=:box, size=(800, 600), dpi=300)
    end

    # Plot wing circulation
    plot!(P, ycp[i_w1:i_w2], gc[i_w1:i_w2],
          lw=2.5, marker=(:circle, 4), color=:royalblue,
          markerstrokewidth=1, markerstrokecolor=:navy,
          label="Wing \$\\Gamma\$", alpha=0.9)

    # Plot tail circulation
    plot!(P, ycp[i_t1:i_t2], gc[i_t1:i_t2],
          lw=2.5, marker=(:diamond, 4), color=:crimson,
          markerstrokewidth=1, markerstrokecolor=:darkred,
          label="Tail \$\\Gamma\$", alpha=0.9)

    # Mirror if requested
    if mirror
        plot!(P, -ycp[i_w1:i_w2], gc[i_w1:i_w2],
              lw=2.5, marker=(:circle, 4), color=:royalblue,
              markerstrokewidth=1, markerstrokecolor=:navy,
              label="", alpha=0.9)
        plot!(P, -ycp[i_t1:i_t2], gc[i_t1:i_t2],
              lw=2.5, marker=(:diamond, 4), color=:crimson,
              markerstrokewidth=1, markerstrokecolor=:darkred,
              label="", alpha=0.9)
    end

    # Overlay elliptical loading for comparison
    if show_elliptical && !normalize
        wing = ac.wing
        b = wing.layout.span
        Γ_max = maximum(gc[i_w1:i_w2])

        if mirror
            # Full span from -b/2 to b/2
            y_ell = range(-b/2, b/2, length=100)
            Γ_elliptical = @. Γ_max * sqrt(1 - (2*y_ell/b)^2)
            plot!(P, y_ell, Γ_elliptical, lw=2, color=:gray,
                  linestyle=:dash, alpha=0.7, label="Elliptical (ideal)")
        else
            # Half span from 0 to b/2
            y_ell = range(0, b/2, length=100)
            Γ_elliptical = @. Γ_max * sqrt(1 - (2*y_ell/b)^2)
            plot!(P, y_ell, Γ_elliptical, lw=2, color=:gray,
                  linestyle=:dash, alpha=0.7, label="Elliptical (ideal)")
        end
    end

    # Mark wing geometry stations
    wing = ac.wing
    bo = wing.layout.root_span
    bs = wing.layout.break_span
    b = wing.layout.span

    max_gc = maximum(gc)
    min_gc = minimum(gc)
    y_range = max_gc - min_gc
    annotation_y = max_gc - 0.5 * y_range  # Place at midpoint

    for (yloc, name, clr) in [(bo/2, "Root", :orange),
                               (bs/2, "Break", :green),
                               (b/2, "Tip", :purple)]
        vline!(P, [yloc], color=clr, lw=1.5, linestyle=:dot, alpha=0.6, label="")
        # Create text with colored background box
        txt = text(name, clr, 10, :center, :top)
        annotate!(P, yloc, annotation_y, txt)
        if mirror
            vline!(P, [-yloc], color=clr, lw=1.5, linestyle=:dot, alpha=0.6, label="")
            annotate!(P, -yloc, annotation_y, txt)
        end
    end

    # Add point annotations if requested
    if annotate
        for i in i_w1:i_w2
            annotate!(P, ycp[i], gc[i], text(string(i), :black, 7, :left))
        end
    end

    return P
end


"""
    plot_circulation_comparison(ac; mirror=true)

Create a multi-panel comparison of bound vs wake circulation distributions.

# Arguments
- `ac::aircraft`: Aircraft object with computed wake system
- `mirror::Bool=true`: Show both sides of the aircraft

# Returns
- Plot object with two subplots showing bound and wake circulation
"""
function plot_circulation_comparison(ac; mirror::Bool=true)
    ws = ac.wake_system
    if ws === nothing
        error("ac.wake_system is empty. Build Trefftz geometry before plotting.")
    end

    # Get data
    ycp = aerodynamics.ctrl_ys(ws)
    ys = aerodynamics.field_ys(ws)
    geom = aerodynamics.TREFFTZ_GEOM
    gc = geom.gc[1:length(ycp)]

    # Get indices for wing and tail
    trefftz_config = ac.options.trefftz_config
    i_w1 = aerodynamics.i_first_wing(trefftz_config)
    i_w2 = aerodynamics.i_last_wing(trefftz_config) - 1
    i_t1 = aerodynamics.i_first_tail(trefftz_config)
    i_t2 = aerodynamics.i_last_tail(trefftz_config) - 1

    # Calculate wake circulation
    gw = zeros(length(ys))
    ifrst = [i_w1, i_t1]
    ilast = [i_w2 + 1, i_t2 + 1]  # +1 because wake has one more point
    aerodynamics.calculate_wake_circulation!(gw, gc, ifrst, ilast, 2)

    # Create layout
    l = @layout [a; b]

    # Plot 1: Bound circulation
    p1 = plot(xlabel="", ylabel="Bound \$\\Gamma_c\$ [\$m^2/s\$]",
              title="Circulation Distributions", grid=true, legend=:topright,
              framestyle=:box)

    # Wing bound circulation with lines
    plot!(p1, ycp[i_w1:i_w2], gc[i_w1:i_w2],
          lw=2.5, marker=(:circle, 4), color=:royalblue,
          markerstrokewidth=1, markerstrokecolor=:navy,
          label="Wing", alpha=0.9)

    # Tail bound circulation with lines
    plot!(p1, ycp[i_t1:i_t2], gc[i_t1:i_t2],
          lw=2.5, marker=(:diamond, 4), color=:crimson,
          markerstrokewidth=1, markerstrokecolor=:darkred,
          label="Tail", alpha=0.9)

    if mirror
        plot!(p1, -ycp[i_w1:i_w2], gc[i_w1:i_w2],
              lw=2.5, marker=(:circle, 4), color=:royalblue,
              markerstrokewidth=1, markerstrokecolor=:navy,
              label="", alpha=0.9)
        plot!(p1, -ycp[i_t1:i_t2], gc[i_t1:i_t2],
              lw=2.5, marker=(:diamond, 4), color=:crimson,
              markerstrokewidth=1, markerstrokecolor=:darkred,
              label="", alpha=0.9)
    end

    # Plot 2: Wake circulation (shed vorticity)
    p2 = plot(xlabel="Spanwise Location y [m]", ylabel="Wake \$\\Gamma_w\$ [\$m^2/s\$]",
              grid=true, legend=:topright, framestyle=:box)

    # Wing wake circulation with lines
    plot!(p2, ys[i_w1:(i_w2+1)], gw[i_w1:(i_w2+1)],
          lw=2.5, marker=(:square, 4), color=:royalblue,
          markerstrokewidth=1, markerstrokecolor=:navy,
          label="Wing", alpha=0.9)

    # Tail wake circulation with lines
    plot!(p2, ys[i_t1:(i_t2+1)], gw[i_t1:(i_t2+1)],
          lw=2.5, marker=(:diamond, 4), color=:crimson,
          markerstrokewidth=1, markerstrokecolor=:darkred,
          label="Tail", alpha=0.9)

    if mirror
        plot!(p2, -ys[i_w1:(i_w2+1)], gw[i_w1:(i_w2+1)],
              lw=2.5, marker=(:square, 4), color=:royalblue,
              markerstrokewidth=1, markerstrokecolor=:navy,
              label="", alpha=0.9)
        plot!(p2, -ys[i_t1:(i_t2+1)], gw[i_t1:(i_t2+1)],
              lw=2.5, marker=(:diamond, 4), color=:crimson,
              markerstrokewidth=1, markerstrokecolor=:darkred,
              label="", alpha=0.9)
    end

    # Add zero reference line
    hline!(p2, [0], color=:black, linestyle=:dash, alpha=0.3, label="")

    fig = plot(p1, p2, layout=l, size=(800, 700), dpi=300,
               link=:x, margin=4Plots.mm)

    return fig
end

"""
    plot_induced_velocities(ac; ip=ipcruise1, show_components=true, mirror=true)

Visualize the induced velocity field at control points in the Trefftz plane.

# Arguments
- `ac::aircraft`: Aircraft object with computed wake system
- `ip::Int=ipcruise1`: Mission point index
- `show_components::Bool=true`: Show velocity and drag density separately
- `mirror::Bool=true`: Show both sides of aircraft

# Returns
- Plot object showing induced velocities
"""
function plot_induced_velocities(ac; ip=ipcruise1, show_components=true, mirror=true)
    ws = ac.wake_system
    if ws === nothing
        error("ac.wake_system is empty. Build Trefftz geometry before plotting.")
    end

    # Recompute circulation for the specific mission point
    para_view = @view ac.para[:, ip, 1]
    aerodynamics.induced_drag!(para_view, ac, ac.options.trefftz_config)

    # Get geometry and circulation data
    ycp = aerodynamics.ctrl_ys(ws)
    geom = aerodynamics.TREFFTZ_GEOM
    gc = geom.gc[1:length(ycp)]

    # Get wake circulation
    ys = aerodynamics.field_ys(ws)
    gw = zeros(length(ys))
    trefftz_config = ac.options.trefftz_config
    ifrst = [aerodynamics.i_first_wing(trefftz_config),
             aerodynamics.i_first_tail(trefftz_config)]
    ilast = [aerodynamics.i_last_wing(trefftz_config),
             aerodynamics.i_last_tail(trefftz_config)]
    aerodynamics.calculate_wake_circulation!(gw, gc, ifrst, ilast, 2)

    # Compute induced velocities using influence matrix
    bref = ac.wing.layout.span
    vnc = zeros(length(ycp))
    @views vnc = ws.influence_matrix * gw
    vnc .*= (bref / (2.0π))

    # Get wing indices
    i_w1 = aerodynamics.i_first_wing(trefftz_config)
    i_w2 = aerodynamics.i_last_wing(trefftz_config) - 1

    if show_components
        # Create two-panel plot for velocity and drag density
        l = @layout [a; b]

        p1 = plot(xlabel="", ylabel="Normal Velocity \$v_n\$ [m/s]",
                  title="Induced Velocities at Control Points",
                  grid=true, legend=:topright, framestyle=:box)
        plot!(p1, ycp[i_w1:i_w2], vnc[i_w1:i_w2],
              lw=2.5, marker=(:circle, 4), color=:royalblue,
              markerstrokewidth=1, markerstrokecolor=:navy,
              label="Wing", alpha=0.9)
        if mirror
            plot!(p1, -ycp[i_w1:i_w2], vnc[i_w1:i_w2],
                  lw=2.5, marker=(:circle, 4), color=:royalblue,
                  markerstrokewidth=1, markerstrokecolor=:navy,
                  label="", alpha=0.9)
        end
        hline!(p1, [0], color=:black, linestyle=:dash, alpha=0.3, label="")

        # Plot induced drag per unit span (related to v_n * Γ)
        p2 = plot(xlabel="Spanwise Location y [m]",
                  ylabel="Induced Drag Density [N/m]",
                  grid=true, legend=:topright, framestyle=:box)
        drag_density = -gc[i_w1:i_w2] .* vnc[i_w1:i_w2]
        plot!(p2, ycp[i_w1:i_w2], drag_density,
              lw=2.5, marker=(:square, 4), color=:crimson,
              markerstrokewidth=1, markerstrokecolor=:darkred,
              label="Wing \$-\\Gamma v_n\$", alpha=0.9)
        if mirror
            plot!(p2, -ycp[i_w1:i_w2], drag_density,
                  lw=2.5, marker=(:square, 4), color=:crimson,
                  markerstrokewidth=1, markerstrokecolor=:darkred,
                  label="", alpha=0.9)
        end

        fig = plot(p1, p2, layout=l, size=(800, 700), dpi=300,
                   link=:x, margin=4Plots.mm)
        return fig
    else
        # Single plot
        P = plot(xlabel="Spanwise Location y [m]",
                ylabel="Normal Velocity \$v_n\$ [m/s]",
                title="Induced Normal Velocity at Control Points",
                grid=true, legend=:topright, framestyle=:box,
                size=(800, 600), dpi=300)
        plot!(P, ycp[i_w1:i_w2], vnc[i_w1:i_w2],
              lw=2.5, marker=(:circle, 4), color=:royalblue,
              markerstrokewidth=1, markerstrokecolor=:navy,
              label="Wing", alpha=0.9)
        if mirror
            plot!(P, -ycp[i_w1:i_w2], vnc[i_w1:i_w2],
                  lw=2.5, marker=(:circle, 4), color=:royalblue,
                  markerstrokewidth=1, markerstrokecolor=:navy,
                  label="", alpha=0.9)
        end
        hline!(P, [0], color=:black, linestyle=:dash, alpha=0.3, label="")
        return P
    end
end


"""
    plot_mission_lift_distribution(ac; ip_range=nothing, mirror=true,
                                   show_elliptical=true, display_plots=true)

Create lift distribution plots for multiple mission points to visualize loading changes
throughout the flight mission.

# Arguments
- `ac::aircraft`: Aircraft object with computed wake system
- `ip_range::Union{Nothing,AbstractRange}=nothing`: Range of mission points to plot.
  If nothing, plots key mission points: rotation, climb, cruise start/end, descent
- `mirror::Bool=true`: Show both sides of the aircraft
- `show_elliptical::Bool=true`: Overlay ideal elliptical loading
- `display_plots::Bool=true`: Automatically display each plot (opens separate windows)

# Returns
- Vector of plot objects, one for each mission point

# Example
```julia
# Plot all cruise points
plots = plot_mission_lift_distribution(ac, ip_range=ipcruise1:ipcruisen)

# Plot key mission segments
plots = plot_mission_lift_distribution(ac)  # Uses default key points

# Plot and save without displaying
plots = plot_mission_lift_distribution(ac, display_plots=false)
for (i, p) in enumerate(plots)
    savefig(p, "lift_dist_ip\$(i).png")
end
```
"""
function plot_mission_lift_distribution(ac; ip_range=nothing, mirror=true,
                                       show_elliptical=true, display_plots=true)

    # Define default key mission points if no range specified
    if ip_range === nothing
        ip_range = [iprotate, ipclimb1, ipclimbn, ipcruise1,
                    Int(round((ipcruise1+ipcruisen)/2)), ipcruisen,
                    ipdescent1, ipdescentn]
        # Filter to only existing points
        ip_range = filter(ip -> ip <= size(ac.para, 2), ip_range)
    end

    # Define mission point names for titles
    mission_names = Dict(
        iprotate => "Takeoff Rotation",
        ipclimb1 => "Start of Climb",
        ipclimbn => "End of Climb (TOC)",
        ipcruise1 => "Start of Cruise",
        ipcruisen => "End of Cruise",
        ipdescent1 => "Start of Descent",
        ipdescentn => "End of Descent"
    )

    plots = []

    for ip in ip_range
        # Get mission point name
        mission_name = get(mission_names, ip, "Mission Point $ip")

        # Get aerodynamic parameters for title
        CL = ac.para[iaCL, ip, 1]
        alt_m = ac.para[iaalt, ip, 1]
        alt_ft = alt_m * 3.28084
        Mach = ac.para[iaMach, ip, 1]

        # Create title with mission info
        title_str = @sprintf("%s (ip=%d)\nCL=%.3f, M=%.3f, Alt=%.0f ft",
                            mission_name, ip, CL, Mach, alt_ft)

        # Create plot for this mission point
        p = plot_lift_distribution(ac, ip=ip, mirror=mirror,
                                  show_elliptical=show_elliptical)
        plot!(p, title=title_str)

        push!(plots, p)

        # Display if requested
        if display_plots
            display(p)
        end
    end

    return plots
end


"""
    plot_mission_animation(ac; ip_range=nothing, mirror=true,
                          show_elliptical=true, fps=2, filename="mission_lift.gif")

Create an animated GIF showing lift distribution changes throughout the mission.

# Arguments
- `ac::aircraft`: Aircraft object with computed wake system
- `ip_range::Union{Nothing,AbstractRange}=nothing`: Range of mission points.
  If nothing, animates all points from climb to descent
- `mirror::Bool=true`: Show both sides of the aircraft
- `show_elliptical::Bool=true`: Overlay ideal elliptical loading
- `fps::Int=2`: Frames per second in output animation
- `filename::String="mission_lift.gif"`: Output filename

# Returns
- Animation object (can be displayed or saved)

# Example
```julia
# Animate entire cruise segment
anim = plot_mission_animation(ac, ip_range=ipcruise1:ipcruisen, fps=5)

# Animate climb through descent
anim = plot_mission_animation(ac, ip_range=ipclimb1:ipdescentn)

# Save with custom filename
anim = plot_mission_animation(ac, filename="my_mission.gif")
```
"""
function plot_mission_animation(ac; ip_range=nothing, mirror=true,
                               show_elliptical=true, fps=2,
                               filename="mission_lift.gif")

    # Define default range if not specified (climb through descent)
    if ip_range === nothing
        ip_range = ipclimb1:ipdescentn
    end

    # Get circulation range for consistent y-axis
    ycp = aerodynamics.ctrl_ys(ac.wake_system)
    trefftz_config = ac.options.trefftz_config
    i_w1 = aerodynamics.i_first_wing(trefftz_config)
    i_w2 = aerodynamics.i_last_wing(trefftz_config) - 1

    # Compute circulation range across all mission points for consistent scaling
    gc_min, gc_max = Inf, -Inf
    for ip in ip_range
        para_view = @view ac.para[:, ip, 1]
        aerodynamics.induced_drag!(para_view, ac, ac.options.trefftz_config)
        geom = aerodynamics.TREFFTZ_GEOM
        gc = geom.gc[1:length(ycp)]
        gc_min = min(gc_min, minimum(gc[i_w1:i_w2]))
        gc_max = max(gc_max, maximum(gc[i_w1:i_w2]))
    end

    # Add 10% padding
    y_padding = 0.1 * (gc_max - gc_min)
    ylims_fixed = (gc_min - y_padding, gc_max + y_padding)

    # Create animation
    anim = @animate for ip in ip_range
        # Get aerodynamic parameters
        CL = ac.para[iaCL, ip, 1]
        alt_m = ac.para[iaalt, ip, 1]
        alt_ft = alt_m * 3.28084
        Mach = ac.para[iaMach, ip, 1]
        CDi = ac.para[iaCDi, ip, 1]

        # Create title
        title_str = @sprintf("Mission Point %d: CL=%.3f, CDi=%.5f\nM=%.3f, Alt=%.0f ft",
                            ip, CL, CDi, Mach, alt_ft)

        # Plot with fixed y-limits
        p = plot_lift_distribution(ac, ip=ip, mirror=mirror,
                                  show_elliptical=show_elliptical)
        plot!(p, title=title_str, ylims=ylims_fixed)
    end

    # Save animation
    gif(anim, filename, fps=fps)

    return anim
end
