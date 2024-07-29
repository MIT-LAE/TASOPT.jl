using DocStringExtensions
abstract type AbstractLayout end
abstract type AbstractCrossSection end

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

"""
$TYPEDEF

Wing Layout Structure:
Contains dimensions, heights, etc. to design a Wing

$TYPEDFIELDS
"""
@kwdef mutable struct WingLayout
    """Aspect Ratio [m]"""
    AR::Float64 = 0 
    """Sweep [degrees]"""
    sweep::Float64 = 0
    """Wing Span [m]"""
    b::Float64 = 0
    """Span of inner wing (break/"snag") [m]"""
    b_inner::Float64 = 0 
    """Max Wing Span [m]"""
    b_max::Float64 = 0 
    """Wing chord at root """
    chord::Float64 = 0
    """Outer or "tip" taper ratio of chord"""
    λt::Float64 = 0 
    """Inner or break/"snag" taper ratio of chord"""
    λs::Float64 = 0 
    """Span fraction of inner wing break ("snag")"""
    ηs::Float64 = 0 
    """Wing center box width [m]"""
    box_halfspan::Float64 = 0 
    """Wing planform area (including fuselage carryover) [m^2]"""
    S::Float64 = 0 
    """Wing Box wing x to chord"""
    box_width_chord::Float64 = 0 
    """Wing Root thickness to chord"""
    root_chord_thickness::Float64 = 0 
    """Spanbreak thickness to chord"""
    spanbreak_chord_thickness::Float64 = 0 
    """web-height/hbox ratio"""
    hweb_to_hbox::Float64 = 0 #igrh
    """Spar box axis x/c location """
    spar_box_x_c::Float64 = 0
    """X position of wing box"""
    x_wing_box::Float64 = 0 
    """X location of wing"""
    x::Float64 = 0
    """Z location of wing"""
    z::Float64 = 0
end

"""
$TYPEDEF

Tail Layout Structure:
Contains dimensions, heights, etc. to design a Tail

$TYPEDFIELDS
"""
@kwdef mutable struct TailLayout
    """Tail planform area (including fuselage carryover) [m^2]"""
    S::Float64 = 0 
    """Aspect Ratio [m]"""
    AR::Float64 = 0 
    """Tail Span [m]"""
    b::Float64 = 0
    """Wing center box width [m]"""
    box_halfspan::Float64 = 0 
    """Taper [degrees]"""
    λ::Float64 = 0 
    """Wing chord at root """
    chord::Float64 = 0
    """Sweep [degrees]"""
    sweep::Float64 = 0
    """Tail box width"""
    box_width::Float64 = 0 
    """Tail box height"""
    box_height::Float64 = 0 
    """Tail box x position"""
    box_x::Float64 = 0 
    """web-height/hbox ratio"""
    hweb_to_hbox::Float64 = 0 
    """Tail Z position"""
    z::Float64 = 0
    """Tail X position"""
    x::Float64 = 0
end
