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

end

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

Contains dimensions of the fuselage cross section for an elliptical type 
fuselage

$TYPEDFIELDS
"""
@kwdef mutable struct EllipticalCrossSection <:AbstractCrossSection
    """Half Width [m]"""
    a::Float64 = 2.0
    """Half height [m]"""
    b::Float64 = 1.8
    """Skin thickness [m]"""
    skin_thickness::Float64 = 0.0
end

# =================
# Area calculations
# =================
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
$(TYPEDSIGNATURES)
"""
function area(cs::EllipticalCrossSection)
    return π * cs.a * cs.b
end

# ==========
# Perimeters
# ==========
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
    get_perimeter(cs::EllipticalCrossSection)

Use [Ramanujam's approximation](https://books.google.com/books?id=oSioAM4wORMC&pg=PA39#v=onepage&q&f=false)
for the circumference of an Ellipse.
"""
function get_perimeter(cs::EllipticalCrossSection)
    a = cs.a
    b = cs.b
    h = ((a - b)/(a + b))^2
    C = π * (a + b) * (1 + (3h / (10 + sqrt(4 - 3h))))
    return C
end

# =============
# Web geometry
# =============
"""
    web_geometry(cs::MultiBubble)

Calculates the geometric properties of the fuselage web
"""
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

function web_geometry(cs::EllipticalCrossSection)
    θ_web = 0.0
    h_web = cs.b
    sin2θ = 0.0
    cosθ = 1.0
    effective_web_length = 0.0
    return θ_web, h_web, sin2θ, cosθ, effective_web_length
end