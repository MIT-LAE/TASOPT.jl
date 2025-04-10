abstract type AbstractWingSection end

"""
$TYPEDEF

Stores all the non-dimensional parameters 
in the normal-plane cross section of the wing (mainly related to the spar-box)
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
@kwdef mutable struct WingCrossSection
    """Wing section's spar box height to perpendicular chord (c⟂) [-]"""
    thickness_to_chord::Float64 = 0.0
    """Wing section's spar box width to c⟂[-]"""
    width_to_chord::Float64 = 0.50 #Default values from TASOPT docs 
    """Wing section's web height to max box height [-]"""
    web_to_box_height::Float64 = 0.75 #Default values from TASOPT docs
    """Sparbox cap normalized thickness (i.e., h_cap/c⟂) [-]"""
    t_cap::Float64 = 0.0
    """Sparbox web normalized thickness"""
    t_web::Float64 = 0.0
    """Internal Area normalized by chord2"""
    A_internal::Float64 = 0.0
end

"""
$TYPEDEF

Wing Section

$TYPEDFIELDS
"""
@kwdef mutable struct WingSection <: AbstractWingSection
    """Wing Section layout"""
    cross_section::WingCrossSection = WingCrossSection()
    """Starting section chord [m]"""
    co::Float64 = 0.0
    """Taper ratio λ = c_end/c_start"""
    λ::Float64 = 0.0
    """Wing Section webs"""
    webs::StructuralMember = StructuralMember(material = StructuralAlloy("TASOPT-Al"))
    """Wing Section caps"""
    caps::StructuralMember = StructuralMember(material = StructuralAlloy("TASOPT-Al"))
    """Bending Stiffness EI matrix [N m^2]"""
    EI::Matrix{Float64} = zeros(Float64, 2, 2)
    """Torsional Stiffness GJ [N m^2]"""
    GJ::Float64 = 0
    """Max shear load [N]"""
    max_shear_load::Float64 = 0
    """Moment [N m]"""
    moment::Float64 = 0
    """Weight [N]"""
    weight::Float64 = 0 # make it a weight
    """Wing root moment contribution from wing weight section of engine [N m]"""
    dyW::Float64 = 0
end

@kwdef mutable struct TailSection <: AbstractWingSection
    """Tail Section layout"""
    cross_section::WingCrossSection = WingCrossSection()
    """Taper ratio λ = c_end/c_start"""
    λ::Float64 = 0.0
    """Bending Stiffness EI matrix [N m^2]"""
    EI::Matrix{Float64} = zeros(Float64, 2, 2)
    """Torsional Stiffness GJ [N m^2]"""
    GJ::Float64 = 0
    """Caps Thickness [m]"""
    thickness_cap::Float64 = 0
    """Webs Thickness [m]"""
    thickness_web::Float64 = 0
end

"""
Wing layout is a structure that defines the wing planform.
See `WingSection` and `WingCrossSection` as well.
"""
@kwdef mutable struct WingLayout
    """Aspect Ratio [m]"""
    AR::Float64 = 0
    """Sweep [degrees]"""
    sweep::Float64 = 0
    """ Wing start [m]"""
    root_span::Float64 = 0.0
    """ Wing Span [m]"""
    span::Float64 = 0
    """Max Wing Span used as a constraint for optimization [m]"""
    max_span::Float64 = 0
    """Span break location [-]"""
    ηs::Float64 = 0.0
    """Wing chord at root """
    root_chord::Float64 = 0
    """Wing planform area (including fuselage carryover) [m^2]"""
    S::Float64 = 0
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
"""
function Base.getproperty(obj::WingLayout, sym::Symbol)
    if sym === :ηo
        getfield(obj, :root_span)/getfield(obj, :span)
    elseif sym === :break_span
        getfield(obj, :ηs)*getfield(obj, :span)
    else
        getfield(obj, sym)
    end
    
end  # function Base.get_property


"""
    normalized_chord(η; λs = 0.8, λt = 0.7, ηo=0.0, ηs = 0.5)

Piecewise function to give the normalized chord C = c/co; co is the root chord.
"""
function normalized_chord(η; λs = 0.8, λt = 0.7, ηo = 0.0, ηs = 0.5)
    if 0.0 ≤ η < ηo
        1.0
    elseif ηo ≤ η < ηs
        1 + (λs - 1) * (η - ηo) / (ηs - ηo)
    elseif ηs ≤ η ≤ 1
        λs + (λt - λs) * (η - ηs) / (1 - ηs)
    else
        error("η should be 0≤η≤1")
    end
end  # function normalized_chord

"""
    get_average_sparbox_heights(section::WingCrossSection) -> (h̄_avg, h̄_rms)

Calculates the average and root mean square (RMS) heights of a spar box for a given wing section layout.
These are used in [`wing_weights!`](@ref) for further calculations
"""
function get_average_sparbox_heights(section::WingCrossSection)
    A = 1 - section.web_to_box_height
    h̄ = section.thickness_to_chord
    h̄_avg = h̄ * (1 - A / 3.0)
    h̄_rms = sqrt(h̄^2 * (1 - 2 * A / 3 + A^2 / 5))
    return h̄_avg, h̄_rms
end  # function get_average_sparbox_heights
