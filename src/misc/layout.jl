using DocStringExtensions
"""
$TYPEDEF

Fuselage Layout Structure:
Contains dimensions, heights, etc. to design a fuselage

$TYPEDFIELDS
"""
@kwdef mutable struct FuselageLayout
    """Fuselage Radius [m]"""
    radius::Float64 = 0 # = ac.parg[igRfuse] #Rfuse
    """Downward shift of lower bubbles (dRfuse) [m] """
    bubble_lower_downward_shift::Float64 = 0#dRfuse 
    """Y offset of bubble center [m]"""
    bubble_center_y_offset::Float64 = 0 #wfb
    """Number of webs (for double bubble designs)"""
    n_webs::Float64 = 0 #nfwebs
    """Thickness of webs """
    thickness_webs::Float64 = 0 #tfwebs
    """X position of nose [m]"""
    x_nose::Float64  = 0# = ac.parg[igxnose] #xnose
    """X position of pressure shell forward [m]"""
    x_pressure_shell_fwd::Float64 = 0# = ac.parg[igxshell1] #xshell1
    """X position of pressure shell aft [m]"""
    x_pressure_shell_aft::Float64 = 0# = ac.parg[igxshell2] #xshell2
    """X position of cylinder start [m]"""
    x_start_cylinder::Float64 = 0# = ac.parg[igxblend1] #xblend1
    """X position of cylinder end [m]"""
    x_end_cylinder::Float64 = 0# = ac.parg[igxblend2] #xblend2
    """X position of fuselage cone end [m]"""
    x_cone_end::Float64 = 0# = ac.parg[igxend] #xconeend
    """X position of fuselage end [m]"""
    x_end::Float64 = 0# = ac.parg[igxend] #xend
    """Tailcone taper (lambdac) [m]"""
    tailcone_taper_ratio::Float64 = 0# lambdac
    """Floor depth (depth of floor beams) [m]"""
    floor_depth::Float64 = 0
    """Nose Radius [m]"""
    nose_radius::Float64 = 0
    """Tail Radius [m]"""
    tail_radius::Float64 = 0
    """Taper fuselage to Point (0) or Edge (1)"""
    taper_fuse::Int64 = 0 # 0 = point ; 1 = edge
end


"""
$TYPEDEF

Wing Layout Structure:
Contains dimensions, heights, etc. to design a Wing

$TYPEDFIELDS
"""
@kwdef mutable struct WingLayout
    """Aspect Ratio [m]"""
    AR::Float64 = 0 # igAr
    """Sweep [degrees]"""
    sweep::Float64 = 0 # igsweep
    """Wing Span [m]"""
    b::Float64 = 0 # igb
    """Span of inner wing (break/"snag") [m]"""
    b_inner::Float64 = 0 # igbs
    """Max Wing Span [m]"""
    b_max::Float64 = 0 # igbmax

    """Wing chord at root """
    chord::Float64 = 0
    """Outer or "tip" taper ratio of chord"""
    λt::Float64 = 0 # iglambdat
    """Inner or break/"snag" taper ratio of chord"""
    λs::Float64 = 0 # iglambdas
    """Span fraction of inner wing break ("snag")"""
    ηs::Float64 = 0 # igetas
    """Wing center box width [m]"""
    box_halfspan::Float64 = 0 # igbo
    """Wing planform area (including fuselage carryover) [m^2]"""
    S::Float64 = 0 # igS

    box_width_chord::Float64 = 0 #igwbox

    root_chord_thickness::Float64 = 0 # hboxo

    spanbreak_chord_thickness::Float64 = 0 # hboxs

    hweb_to_hbox::Float64 = 0 #igrh

    spar_box_x_c::Float64 = 0 #igXaxis

    x_wing_box::Float64 = 0 # igxwbox

    x::Float64 = 0

    z::Float64 = 0

end

"""
$TYPEDEF

Tail Layout Structure:
Contains dimensions, heights, etc. to design a Tail

$TYPEDFIELDS
"""
@kwdef mutable struct TailLayout
    tail_volume::Float64 = 0 #igV
    AR::Float64 = 0
    S::Float64 = 0
    b::Float64 = 0
    box_halfspan::Float64 = 0 # igboh or igbov
    λ::Float64 = 0 # iglambdah or iglambdav
    chord::Float64 = 0 # igcoh or igcov
    sweep::Float64 = 0 # igsweeph or igsweepv
    box_width::Float64 = 0 # igwboxh or igwboxv
    box_height::Float64 = 0 # ighboxh or ighboxv
    box_x::Float64 = 0 # igxhbox or igxvbox
    hweb_to_hbox::Float64 = 0 # igrhh or igrhv
    z::Float64 = 0
    x::Float64 = 0
end