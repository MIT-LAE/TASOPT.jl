using DocStringExtensions

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

    # Cross sectional areas
    A_skin = perimeter * fuse.skin.thickness
    A_web = web_length * fuse.web.thickness
    A_fuse = (π + layout.n_webs*(2θ_web + sin2θ))*R^2 + 2R*ΔR

    # Surface areas
    S_bulk = (2π + 4.0*layout.n_webs*θ_web)*R^2
    S_nose = S_bulk * (0.333 + 0.667*(l_nose/R)^1.6 )^0.625 

    l_nose, l_shell, l_floor = calculate_shell_lengths(layout)
    # Shell volumes
    V_cylinder = A_skin*l_shell
    V_nose = S_nose * fuse.skin.thickness
    V_bulk = S_bulk * fuse.skin.thickness
    V_web = A_web*l_shell
      
end  # function calculate_shell_geometry

function calculate_shell_lengths(layout::FuselageLayout)
    l_nose  = layout.x_pressure_shell_fwd - layout.x_nose
    l_shell = layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd
    l_floor = layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd + 2.0 * layout.radius
    return l_nose, l_shell, l_floor
end


function skin_thickness(xSection::AbstractCrossSection, Δp::AbstractFloat, σ::AbstractFloat)
    return Δp * xSection.radius / σ
end

function web_thickness(xSection::AbstractCrossSection, Δp::AbstractFloat, σ::AbstractFloat)
    return 2 * Δp * xSection.bubble_center_y_offset / σ
end

function web_geometry(cs::MultiBubble)
    # fuselage bubble subtended half-angle
    θ_web = asin(cs.bubble_center_y_offset/cs.radius)
    h_web = sqrt(cs.radius^2 - cs.bubble_center_y_offset^2)
    cosθ = h_web/cs.radius
    sinθ = cs.bubble_center_y_offset/cs.radius
    sin2θ = 2 * sinθ * cosθ

    effective_web_length = cs.n_webs*(2*h_web + cs.bubble_lower_downward_shift)

    return θ_web, h_web, sin2θ, effective_web_length
end

function web_geometry(cs::SingleBubble)
    θ_web = 0.0
    h_web = cs.radius
    sin2θ = 0.0
    effective_web_length = 0.0
    return θ_web, h_web, sin2θ, effective_web_length
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
    skin = ((π + n_webs*(2θ_web + sin2θ))*R^2 + 8*n_webs*cosθ*ΔR/2*R +
            (2π + 4n_webs*θ_web)*(ΔR/2)^2)*R*t_shell

    web = 2//3 * n_webs*(h_web + ΔR/2)^3 * t_web

    return skin + web

end  # function Iy
# """
# """
# function Iz(cs::AbstractCrossSection)
    
# end  # function Iz

"""
    get_perimeter(x::SingleBubble)

$(TYPEDSIGNATURES)

Returns the perimeter of a given cross-section
"""
function get_perimeter(cs::SingleBubble)
    return (2π*cs.radius) + (2*cs.bubble_lower_downward_shift)
end  # function perimeter

"""
    get_perimeter(x::MultiBubble)

$(TYPEDSIGNATURES)

Returns the perimeter of a given cross-section
"""
function get_perimeter(cs::MultiBubble)
    θ_web, _, _, _ = web_geometry(cs)
    perimeter = (2π + 4.0*θ_web*cs.n_webs)*cs.radius + 
                (2*cs.bubble_lower_downward_shift)
    return perimeter
end  # function perimeter