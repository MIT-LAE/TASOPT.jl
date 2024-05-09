using DocStringExtensions
"""
$TYPEDEF

Fuselage Layout Structure:
Contains dimensions, heights, etc. to design a fuselage

$TYPEDFIELDS
"""
@kwdef mutable struct FuselageLayout
    """Fuselage Radius [m]"""
    radius::Float64 = 1.9558 # = ac.parg[igRfuse] #Rfuse
    """Downward shift of lower bubbles (dRfuse) [m] """
    bubble_lower_downward_shift::Float64 = 0.381#dRfuse 
    """Y offset of bubble center [m]"""
    bubble_center_y_offset::Float64 = 0 #wfb
    """Number of webs (for double bubble designs)"""
    n_webs::Float64 = 1 #nfwebs
    """Thickness of webs """
    thickness_webs::Float64 = 0 #nfwebs
    """X position of nose [m]"""
    x_nose::Float64  = 0# = ac.parg[igxnose] #xnose
    """X position of pressure shell forward [m]"""
    x_pressure_shell_fwd::Float64 = 5.1816# = ac.parg[igxshell1] #xshell1
    """X position of pressure shell aft [m]"""
    x_pressure_shell_aft::Float64 = 31.0896# = ac.parg[igxshell2] #xshell2
    """X position of cylinder start [m]"""
    x_start_cylinder::Float64 = 6.096# = ac.parg[igxblend1] #xblend1
    """X position of cylinder end [m]"""
    x_end_cylinder::Float64 = 29.5656# = ac.parg[igxblend2] #xblend2
    """X position of fuselage cone end [m]"""
    x_cone_end::Float64 = 35.6616# = ac.parg[igxend] #xconeend
    """X position of fuselage end [m]"""
    x_end::Float64 = 35.6616# = ac.parg[igxend] #xend
    """Tailcone taper (lambdac) [m]"""
    tailcone_taper_ratio::Float64 = 0.3# lambdac
    """Floor depth (depth of floor beams) [m]"""
    floor_depth::Float64 = 0.127
    """Nose Radius [m]"""
    nose_radius::Float64 = 1.65
    """Tail Radius [m]"""
    tail_radius::Float64 = 2.0
    """Taper fuselage to Point (0) or Edge (1)"""
    taper_fuse::Int64 = 1 # 0 = point ; 1 = edge
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
    """Outer or "tip" taper ratio of chord"""
    λt::Float64 = 0
    """Inner or break/"snag" taper ratio of chord"""
    λs::Float64 = 0
    """Span fraction of inner wing break ("snag")"""
    ηs::Float64 = 0
    """Wing center box width [m]"""
    box_width::Float64 = 0
    """Wing planform area (including fuselage carryover) [m^2]"""
    S::Float64 = 0
end

# function FuselageLayout(;default = true)
#     #TODO add read input
#     # if default

#     # else

#     # end
#     return FuselageLayout(fuse_radius = 3, )
# end