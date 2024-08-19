using DocStringExtensions

"""
$TYPEDEF

Tail

$TYPEDFIELDS
"""
@kwdef mutable struct Tail
    """Tail Layout """
    layout::WingLayout = WingLayout()
    """Tail Section """
    outboard::TailSection = TailSection()
    """Tail Weight [N] """
    weight::Float64 = 0
    """Tail Added Weight Fraction"""
    weight_fraction_added::Float64 = 0
    """Tail Max CL """
    CL_max::Float64 = 0
    """Tail Volume [m^3] """
    volume::Float64 = 0
    """Tail Sizing factor: 1=set Sh via specified Vh, 2=et Sh via CLhCGfwd at max-forward CG during landing """
    size::Int64 = 0
    """Tail Downwash factor dε/dα """
    downwash_factor::Float64 = 0
    """Tail max fwd CG (only used if HTsize == "maxforwardCG") """
    CL_max_fwd_CG::Float64 = 0
    """Tail Minimum static margin"""
    SM_min::Float64 = 0
    """Max Tail down load. Tail download param at max load case"""
    CL_CLmax::Float64 = 0
    """Number of Tails"""
    ntails::Float64 = 0
    """Move wingbox factor. 0="fix" wing position ,1=move wing to get CLh="CLhspec" in cruise, 2= move wing to get min static margin = "SMmin"  """
    move_wingbox::Int64 = 0
    
end