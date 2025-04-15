"""
    calculate_shell_geometry(fuse::Fuselage, Δp::AbstractFloat)

Calculates the geometry (distances, angles, areas and volumes)
of the primary fuselage shell structure. The shell consists of the fuselage
skin and the "smeared" out stringers. For multi-bubble fuselages,
e.g., the D8 like designs, this also calculates the geometry of the webs between
the bubbles.

"""
function calculate_shell_geometry!(fuse::Fuselage, Δp::AbstractFloat)

    layout = fuse.layout
    R = layout.cross_section.radius
    ΔR = layout.cross_section.bubble_lower_downward_shift

    θ_web, h_web, sin2θ, web_length = web_geometry(layout.cross_section)
    perimeter = get_perimeter(layout.cross_section)

    fuse.skin.thickness = Δp * R / fuse.skin.σ
    web_thickness = 2.0 * Δp * layout.bubble_center_y_offset / fuse.skin.σ

    # Cross sectional areas
    A_skin = perimeter * fuse.skin.thickness
    A_web = web_length * web_thickness
    A_fuse = area(layout.cross_section)

    l_nose, l_shell, l_floor = calculate_shell_lengths(layout)

    # Surface areas
    S_bulk = (2π + 4.0 * layout.n_webs * θ_web) * R^2
    S_nose = S_bulk * (0.333 + 0.667 * (l_nose / R)^1.6)^0.625

    shell_centroid = 0.5 * (layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)
    nose_centroid = 0.5 * (layout.x_nose + layout.x_pressure_shell_fwd)
    rear_bulkhead_centroid = layout.x_pressure_shell_aft + 0.5 * R

    # Shell volumes
    V_cylinder = A_skin * l_shell
    V_nose = S_nose * fuse.skin.thickness
    V_bulk = S_bulk * fuse.skin.thickness
    V_web = A_web * l_shell

    # Weights
    skin_weight_density = fuse.skin.ρ * gee
    W_nose = Weight(W = skin_weight_density * V_nose, x = nose_centroid)
    W_bulk = Weight(W = skin_weight_density * V_bulk, x = rear_bulkhead_centroid)
    W_cylinder = Weight(W = skin_weight_density * V_cylinder, x = shell_centroid)

    W_skin = W_nose + W_cylinder + W_bulk
    W_web = Weight(W = fuse.skin.ρ * gee * (V_web), x = shell_centroid)

    fuse.shell.weight =
        W_skin * (
            1.0 +
            fuse.weight_frac_stringers +
            fuse.weight_frac_frame +
            fuse.weight_frac_skin_addl
        ) + W_web

    insulation = size_insulation(layout, S_nose, S_bulk, fuse.insulation_W_per_area)

    return A_fuse, insulation

end  # function calculate_shell_geometry

function calculate_shell_lengths(layout::FuselageLayout)
    l_nose = layout.x_pressure_shell_fwd - layout.x_nose
    l_shell = layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd
    l_floor =
        layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd + 2.0 * layout.radius
    return l_nose, l_shell, l_floor
end


"""
    skin_thickness(xSection::AbstractCrossSection, Δp::AbstractFloat, σ::AbstractFloat)

TBW
"""
function skin_thickness(xSection::AbstractCrossSection, Δp::AbstractFloat, σ::AbstractFloat)
    return Δp * xSection.radius / σ
end

"""
    web_thickness(xSection::AbstractCrossSection, Δp::AbstractFloat, σ::AbstractFloat)

TBW
"""
function web_thickness(xSection::AbstractCrossSection, Δp::AbstractFloat, σ::AbstractFloat)
    return 2 * Δp * xSection.bubble_center_y_offset / σ
end

function web_geometry(cs::MultiBubble)
    # fuselage bubble subtended half-angle
    θ_web = asin(cs.bubble_center_y_offset / cs.radius)
    h_web = sqrt(cs.radius^2 - cs.bubble_center_y_offset^2)
    cosθ = h_web / cs.radius
    sinθ = cs.bubble_center_y_offset / cs.radius
    sin2θ = 2 * sinθ * cosθ

    effective_web_length = cs.n_webs * (2 * h_web + cs.bubble_lower_downward_shift)

    return θ_web, h_web, sin2θ, cosθ, effective_web_length
end

function web_geometry(cs::SingleBubble)
    θ_web = 0.0
    h_web = cs.radius
    sin2θ = 0.0
    cosθ = 1.0
    effective_web_length = 0.0
    return θ_web, h_web, sin2θ, cosθ, effective_web_length
end

"""
    Iy(cs::AbstractCrossSection)
Calculates the second moment of area about the y-axis (i.e., the horizontal 
axis w.r.t the cross-section)

"""
function Iy(cs::AbstractCrossSection)
    θ_web, h_web, sin2θ, _ = web_geometry(cs)
    n_webs = cs.n_webs
    R = cs.radius
    skin =
        (
            (π + n_webs * (2θ_web + sin2θ)) * R^2 +
            8 * n_webs * cosθ * ΔR / 2 * R +
            (2π + 4n_webs * θ_web) * (ΔR / 2)^2
        ) *
        R *
        t_shell

    web = 2 // 3 * n_webs * (h_web + ΔR / 2)^3 * t_web

    return skin + web

end  # function Iy

# """
# """
# function Iz(cs::AbstractCrossSection)

# end  # function Iz

"""
    $(TYPEDSIGNATURES)
"""
function area(cs::SingleBubble)
    R = cs.radius
    ΔR = cs.bubble_lower_downward_shift
    enclosed_area = π * R^2 + 2 * R * ΔR
    return enclosed_area
end  # function area

"""
$(TYPEDSIGNATURES)
"""
function area(cs::MultiBubble)
    R = cs.radius
    ΔR = cs.bubble_lower_downward_shift
    θ_web, h_web, sin2θ, web_length = web_geometry(cs)
    enclosed_area = (π + cs.n_webs * (2θ_web + sin2θ)) * R^2 + 2R * ΔR
    return enclosed_area
end # function area

"""
    get_perimeter(x::SingleBubble)

$(TYPEDSIGNATURES)

Returns the perimeter of a given cross-section
"""
function get_perimeter(cs::SingleBubble)
    return (2π * cs.radius) + (2 * cs.bubble_lower_downward_shift)
end  # function perimeter

"""
    get_perimeter(x::MultiBubble)

$(TYPEDSIGNATURES)

Returns the perimeter of a given cross-section
"""
function get_perimeter(cs::MultiBubble)
    θ_web, _, _, _ = web_geometry(cs)
    perimeter =
        (2π + 4.0 * θ_web * cs.n_webs) * cs.radius + (2 * cs.bubble_lower_downward_shift)
    return perimeter
end  # function perimeter

"""
    size_insulation(cs::AbstractCrossSection,
     insulated_fraction::Float64=0.55)

Calculates the insulation weight. Assumes, by default, that 55% of the fuselage 
cross-sectional perimeter is the cabin which is insulated and 45% is the cargo 
hold that is not insulated.
"""
function size_insulation(
    layout::FuselageLayout,
    S_nose::Float64,
    S_bulk::Float64,
    weight_per_area::Float64,
    insulated_fraction::Float64 = 0.55,
)

    p = get_perimeter(layout.cross_section)
    A = p * layout.l_shell + S_nose + S_bulk
    W = weight_per_area * insulated_fraction * A

    shell_centroid = 0.5 * (layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)
    insulation = Weight(W = W, x = shell_centroid)

    return insulation

end  # function size_insulation


"""
    size_windows(x_cabin, l_cabin, weight_per_length, n_decks)
    weight_per_length, n_decks, n_fuel_tanks, l_tank)

"""
function size_windows(x_cabin, l_cabin, weight_per_length, n_decks)
    windows = Weight(W = n_decks * weight_per_length * l_cabin, x = x_cabin)
    return windows
end  # function size_windows

"""
"""
function get_cabin_dimensions(
    layout::FuselageLayout,
    tank_placement::String,
    n_fuel_tanks,
    l_tank,
)
    if n_fuel_tanks == 1
        if tank_placement == "front"
            x_cabin =
                0.5 * (
                    layout.x_pressure_shell_fwd +
                    l_tank +
                    2.0 * ft_to_m +
                    layout.x_pressure_shell_aft
                )
        else #tank is at rear
            x_cabin =
                0.5 * (
                    layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft -
                    (l_tank + 2.0 * ft_to_m)
                )
        end
    else
        #two symmetric or zero tanks
        x_cabin = 0.5 * (layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)
    end

    #cabin length is smaller if there are fuel tanks
    l_cabin =
        layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd -
        n_fuel_tanks * (l_tank + 2.0 * ft_to_m)

    return x_cabin, l_cabin
end  # function get_cabin_centroid

