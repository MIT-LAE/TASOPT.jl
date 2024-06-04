
"""
    calculate_shell_geometry(fuse::Fuselage)

"""
function calculate_shell_geometry(fuse::Fuselage, Δp::AbstractFloat)
      
      layout = fuse.layout
      R = layout.radius
      ΔR = layout.bubble_lower_downward_shift

      θ_web, h_web, sin2θ = web_geometry(layout.cross_section)

      perimeter = (2π + 4.0*θ_web*layout.n_webs)*R + 2ΔR

      fuse.skin.thickness = Δp*R/fuse.skin.σ
      fuse.web.thickness = 2.0*Δp*layout.bubble_center_y_offset/fuse.skin.σ

      # Cross sectional areas
      A_skin = perimeter*fuse.skin.thickness

      
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
function calculate_bubble_web_geometry(layout::FuselageLayout)
    R = layout.radius
    ΔR = layout.bubble_lower_downward_shift

    #[TODO] the following few lines can be cleaner if we knew whether this was a double bubble design or standard fuse
    center_to_web_distance = max(min(layout.bubble_center_y_offset, R), 0.0)
    # fuselage bubble subtended half-angle
    θ_web = asin(center_to_web_distance/R)
    h_web = sqrt(R^2 - layout.bubble_center_y_offset^2)
    
    cosθ = h_web/R
    sinθ = layout.bubble_center_y_offset/R
    sin2θ = 2.0*sinθ*cosθ 
end  # function calculate_web_geometry