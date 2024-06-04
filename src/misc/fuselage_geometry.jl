
"""
    calculate_shell_geometry(fuse::Fuselage, Δp::AbstractFloat)

Calculates the geometry (distances, angles, areas and volumes)
of the primary fuselage shell structure. The shell consists of the fuselage
skin and the "smeared" out stringers. For multi-bubble fuselages,
e.g., the D8 like designs, this also calculates the geometry of the webs between
the bubbles.

"""
function calculate_shell_geometry(fuse::Fuselage, Δp::AbstractFloat)
      
    layout = fuse.layout
    R = layout.cross_section.radius
    ΔR = layout.cross_section.bubble_lower_downward_shift

    θ_web, h_web, sin2θ, web_length = web_geometry(layout.cross_section)
    perimeter = get_perimeter(layout.cross_section)

    fuse.skin.thickness = Δp*R/fuse.skin.σ
    fuse.web.thickness = 2.0*Δp*layout.bubble_center_y_offset/fuse.skin.σ

    # Effective nose length and pressure-vessel length
    l_nose  = layout.x_pressure_shell_fwd - layout.x_nose
    l_shell = layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd

    # Cross sectional areas
    A_skin = perimeter * fuse.skin.thickness
    A_web = web_length * fuse.web.thickness
    A_fuse = (π + layout.n_webs*(2θ_web + sin2θ))*R^2 + 2R*ΔR

    # Surface areas
    S_bulk = (2π + 4.0*layout.n_webs*θ_web)*R^2
    S_nose = S_bulk * (0.333 + 0.667*(l_nose/R)^1.6 )^0.625 

    # Shell volumes
    V_cylinder = A_skin*l_shell
    V_nose = S_nose * fuse.skin.thickness
    V_bulk = S_bulk * fuse.skin.thickness
    V_web = A_web*l_shell
      
end  # function calculate_shell_geometry


function web_geometry(x::MultiBubble)
    # fuselage bubble subtended half-angle
    θ_web = asin(x.bubble_center_y_offset/x.radius)
    h_web = sqrt(x.radius^2 - x.bubble_center_y_offset^2)
    cosθ = h_web/x.radius
    sinθ = x.bubble_center_y_offset/x.radius
    sin2θ = 2 * sinθ * cosθ

    effective_web_length = x.n_webs*(2*h_web + x.bubble_center_y_offset)

    return θ_web, h_web, sin2θ, effective_web_length
end

function web_geometry(x::SingleBubble)
    θ_web = 0.0
    h_web = x.radius
    sin2θ = 0.0
    effective_web_length = 0.0
    return θ_web, h_web, sin2θ, effective_web_length
end


"""
"""
function get_perimeter(x::SingleBubble)
    return (2π*x.radius) + (2*x.bubble_lower_downward_shift)
end  # function perimeter

"""
"""
function get_perimeter(x::MultiBubble)
    θ_web, _, _, _ = web_geometry(x)
    perimeter = (2π + 4.0*θ_web*x.n_webs)*x.radius + 
                (2*x.bubble_lower_downward_shift)
    return perimeter
end  # function perimeter