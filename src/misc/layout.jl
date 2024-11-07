using DocStringExtensions
abstract type AbstractLayout end
abstract type AbstractCrossSection end
abstract type AbstractCabin end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a conventional
single-bubble type fuselage.

$TYPEDFIELDS
"""
@kwdef mutable struct SingleBubble <: AbstractCrossSection
    """Fuselage Radius [m]"""
    radius::Float64 = 1.0
    """Downward shift of lower bubbles (dRfuse) [m] """
    bubble_lower_downward_shift::Float64 = 0.0
    """Skin thickness [m]"""
    skin_thickness::Float64 = 0.0

    # SingleBubble(radius, bubble_lower_downward_shift, Δp, σ) = 
    # calc_skin_thickness(new(radius, bubble_lower_downward_shift), Δp, σ)

end

"""
"""
function calc_skin_thickness(cs::AbstractCrossSection, Δp, σ)
    cs.skin_thickness = Δp*cs.radius/100.0
    return cs
end  # function calc_skin_thickness


# function SingleBubble(radius, dz, Δp)

# end

"""
$TYPEDEF

Contains dimensions of the fuselage cross section for a multi-bubble type 
fuselage

$TYPEDFIELDS
"""
@kwdef mutable struct MultiBubble <: AbstractCrossSection
    """Fuselage Radius [m]"""
    radius::Float64 = 1.0
    """Downward shift of lower bubbles (dRfuse) [m] """
    bubble_lower_downward_shift::Float64 = 0.2
    """Y offset of bubble center [m]"""
    bubble_center_y_offset::Float64 = 0.2
    """Number of webs [-]"""
    n_webs::Int64 = 1
    """Skin thickness [m]"""
    skin_thickness::Float64 = 0.0
    """Web thickness [m]"""
    web_thickness::Float64 = 0.0
end

# Trying to access these properties form a SingleBubble just gives
# constants
function Base.getproperty(obj::SingleBubble, sym::Symbol)
    if sym === :n_webs
        return 0
    elseif sym === :bubble_center_y_offset
        return 0.0
    elseif sym === :web_thickness
        return 0.0
    elseif sym === :perimeter
        return get_perimeter(obj)
    elseif sym === :area
        return area(obj)
    else
        return getfield(obj, sym)
    end
end  # function Base.getproperty

#Return area and perimeter of multi bubble
function Base.getproperty(obj::MultiBubble, sym::Symbol)
    if sym === :perimeter
        return get_perimeter(obj)
    elseif sym === :area
        return area(obj)
    else
        return getfield(obj, sym)
    end
end  # function Base.getproperty

"""
$TYPEDEF

Fuselage Layout Structure:
Contains dimensions, heights, etc. to design a fuselage

$TYPEDFIELDS
"""
@kwdef mutable struct FuselageLayout <: AbstractLayout
    """Cross section definition"""
    cross_section::AbstractCrossSection = SingleBubble()
    """Thickness of webs """
    thickness_webs::Float64 = 0
    """X position of nose [m]"""
    x_nose::Float64  = 0
    """X position of pressure shell forward [m]"""
    x_pressure_shell_fwd::Float64 = 0
    """X position of pressure shell aft [m]"""
    x_pressure_shell_aft::Float64 = 0
    """X position of cylinder start [m]"""
    x_start_cylinder::Float64 = 0
    """X position of cylinder end [m]"""
    x_end_cylinder::Float64 = 0
    """X position of fuselage cone end [m]"""
    x_cone_end::Float64 = 0
    """X position of fuselage end [m]"""
    x_end::Float64 = 0
    """Tailcone taper (lambdac) [m]"""
    taper_tailcone::Float64 = 0# lambdac
    """Floor depth (depth of floor beams) [m]"""
    floor_depth::Float64 = 0
    """Nose Radius [m]"""
    nose_radius::Float64 = 0
    """Tail Radius [m]"""
    tail_radius::Float64 = 0
    """Taper fuselage to Point (0) or Edge (1)"""
    taper_fuse::Int64 = 0 # 0 = point ; 1 = edge
    """Length of cylindrical portion of cabin that contains payload [m]"""
    l_cabin_cylinder::Float64 = 0.0
end

# Helper function to be able to simplify 
function Base.getproperty(layout::FuselageLayout, sym::Symbol)
    cross_section = getfield(layout, :cross_section)
    
    if sym === :l_nose
        return getfield(layout, :x_pressure_shell_fwd) - getfield(layout, :x_nose)
    elseif sym === :l_shell
        return getfield(layout, :x_pressure_shell_aft) - 
               getfield(layout, :x_pressure_shell_fwd)
    elseif sym === :l_floor
        return getproperty(layout, :l_shell) + 2.0*getproperty(layout, :radius)
    
    elseif sym ∈ (:radius, :n_webs, :bubble_lower_downward_shift, :bubble_center_y_offset)
        return getproperty(cross_section, sym)
    else
        return getfield(layout, sym)
    end
end

@kwdef mutable struct WingLayout
    """Aspect Ratio [m]"""
    AR::Float64 = 0 
    """Sweep [degrees]"""
    sweep::Float64 = 0
    """ Wing Span [m]"""
    b::Float64 = 0 
    """Max Wing Span [m]"""
    b_max::Float64 = 0 
    """Wing chord at root """
    chord::Float64 = 0
    """Wing planform area (including fuselage carryover) [m^2]"""
    S::Float64 = 0 
    """Wing Box wing x to chord"""
    box_width::Float64 = 0 
    """web-height/hbox ratio"""
    hweb_to_hbox::Float64 = 0 #igrh
    """Spar box axis x/c location """
    spar_box_x_c::Float64 = 0
    """X position of wing box"""
    box_x::Float64 = 0 
    """X location of wing"""
    x::Float64 = 0
    """Z location of wing"""
    z::Float64 = 0
end

"""
    scaled_cross_section(cross_section::SingleBubble, R::Float64)

Calculates the geometric properties of a scaled single-bubble cross section.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `cross_section::SingleBubble`: unscaled fuselage cross-section
    - `R::Float64`: radius of geometrically-similar cross-section (m)

    **Outputs:**
    - `p::Float64`: perimeter
    - `A::Float64`: cross-sectional area (m^2)
"""
function scaled_cross_section(cross_section::SingleBubble, R::Float64)
    scaled_cs = deepcopy(cross_section) #Deepcopy to avoid modifying
    #Scale geometric parameters 
    R_Rprev = R/cross_section.radius

    #Scale geometric parameters
    scaled_cs.radius = R_Rprev * cross_section.radius #Change radius 
    scaled_cs.bubble_lower_downward_shift = R_Rprev * cross_section.bubble_lower_downward_shift #Change downward shift

    return get_perimeter(scaled_cs), area(scaled_cs)
end

"""
    scaled_cross_section(cross_section::MultiBubble, R::Float64)

Calculates the geometric properties of a scaled multi-bubble cross section.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `cross_section::MultiBubble`: unscaled fuselage cross-section
    - `R::Float64`: radius of geometrically-similar cross-section (m)

    **Outputs:**
    - `p::Float64`: perimeter
    - `A::Float64`: cross-sectional area (m^2)
"""
function scaled_cross_section(cross_section::MultiBubble, R::Float64)
    scaled_cs = deepcopy(cross_section) #Deepcopy to avoid modifying
    #Scale geometric parameters 
    R_Rprev = R/cross_section.radius #Radii ratio

    #Scale geometric parameters
    scaled_cs.radius = R_Rprev * cross_section.radius #Change radius 
    scaled_cs.bubble_lower_downward_shift = R_Rprev * cross_section.bubble_lower_downward_shift #Change downward shift
    scaled_cs.bubble_center_y_offset = R_Rprev * cross_section.bubble_center_y_offset

    return get_perimeter(scaled_cs), area(scaled_cs)
end

"""
$TYPEDEF

Fuselage Layout Structure:
Contains dimensions, heights, etc. to design a fuselage

$TYPEDFIELDS
"""
@kwdef mutable struct Cabin <: AbstractCabin
    """Longitudinal seat pitch [m]"""
    seat_pitch = 0.0
    """Transverse seat width [m]"""
    seat_width = 0.0
    """Seat height [m]"""
    seat_height = 0.0
    """Aisle half-width [m]"""
    aisle_halfwidth = 0.0
    """Distance between double decker floors [m]"""
    floor_distance = 0.0
    """Main cabin width [m]"""
    cabin_width_main::Float64 = 0.0
    """Top cabin width [m]"""
    cabin_width_top::Float64 = 0.0
    """Number of seats abreast in main cabin"""
    seats_abreast_main::Int64 = 0
    """Number of seats abreast in top cabin"""
    seats_abreast_top::Int64 = 0
    """Floor angle of main cabin [rad]"""
    floor_angle_main::Float64 = 0.0
    """Floor angle of top cabin [rad]"""
    floor_angle_top::Float64 = 0.0
end

"""
$TYPEDEF

Wing Layout Structure:
Holds the span, taper, and chord thickness of a single wing section.
Ratios are with respect to the normal-plane chord c⟂=c×cos(Λ).

Cross-section of wing box:                                                                                                            
                  ┌──────────────────────────────────────┐        
               ┌──┘               ▲                      └──┐     
           ▲┌──┘                  │                         └──┐  
 web height ││               spar box height                  ││  
           ▼└──┐                  │                         ┌──┘  
               └──┐               ▼                      ┌──┘     
                  └──────────────────────────────────────┘        
             ◄───────────────── box width ──────────────────────►   
   ◄───────────────────────────── c⟂ ──────────────────────────────────►                                                    
                                                          
$TYPEDFIELDS
"""
@kwdef mutable struct WingSectionLayout
    """Wing Section Span [m]"""
    b::Float64 = 0.0
    """Wing Section taper"""
    λ::Float64 = 0.0
    """Wing section's spar box height to perpendicular chord (c⟂) [-]"""
    thickness_to_chord::Float64 = 0.0
    """Wing section's spar box width to c⟂[-]"""
    width_to_chord::Float64 = 0.50 #Default values from TASOPT docs #TODO needs to be connected to box_width of WingLayout
    """Wing section's web height to max box height [-]"""
    web_to_box_height::Float64 = 0.75 #Default values from TASOPT docs
    """Sparbox cap normalized thickness (i.e., h_cap/c⟂) [-]"""
    t_cap::Float64 = 0.0
    """Sparbox web normalized thickness"""
end

@kwdef mutable struct WingCrossSection
    """Section Chord [m]"""
    c_perp::Float64 = 0.0
    """Wing section's spar box height to perpendicular chord (c⟂) [-]"""
    thickness_to_chord::Float64 = 0.0
    """Wing section's spar box width to c⟂[-]"""
    width_to_chord::Float64 = 0.50 #Default values from TASOPT docs #TODO needs to be connected to box_width of WingLayout
    """Wing section's web height to max box height [-]"""
    web_to_box_height::Float64 = 0.75 #Default values from TASOPT docs
    """Sparbox cap normalized thickness (i.e., h_cap/c⟂) [-]"""
    t_cap::Float64 = 0.0
    """Sparbox web normalized thickness"""
    t_web::Float64 = 0.0
    """Internal Area normalized by chord2"""
    A_internal::Float64 = 0.0
end

@kwdef mutable struct WingSection_ #_ just to differentiate temporarily from the other struct WingSection
    start_section::WingCrossSection
    end_section::WingCrossSection
    """Section length [m]"""
    b_section::Float64
    """Sweep [°]"""
    Λ::Float64 = 0.0 # need to let this be the wing overall Λ by default

end

"""
"""
function normalized_chord(η; λs = 0.8, λt = 0.7, ηo=0.0, ηs = 0.5)
    if 0.0 ≤ η < ηo
        1.0
    elseif ηo ≤ η < ηs
        1 + (λs - 1) * (η - ηo)/(ηs - ηo)
    elseif ηs ≤ η ≤ 1
        λs + (λt - λs)*(η - ηs)/(1 - ηs)
    else
        error("η should be 0≤η≤1")
    end
end  # function normalized_chord

"""
    get_average_sparbox_heights(section::WingSectionLayout) -> (h̄_avg, h̄_rms)

Calculates the average and root mean square (RMS) heights of a spar box for a given wing section layout.
These are used in [`surfw`](@ref) for further calculations
"""
function get_average_sparbox_heights(section::WingSectionLayout)
    A = 1 - section.web_to_box_height
    h̄ = section.thickness_to_chord
    h̄_avg = h̄ * (1 - A / 3.0)
    h̄_rms = sqrt(h̄^2 * (1 - 2 * A / 3 + A^2 / 5))
    return h̄_avg, h̄_rms
end  # function get_average_sparbox_heights

