@kwdef mutable struct Tail
    layout::TailLayout = TailLayout()

    weight::Float64 = 0
    dxW::Float64 = 0
    EI_bending::Float64 = 0 #igEIch or igEIcv
    EI_normal::Float64 = 0 #igEInh or igEInv
    GJ::Float64 = 0 #igGJh or igGJv
    thickness_cap::Float64 = 0 #igtbcaph or igtbcapv
    thickness_web::Float64 = 0 #igtbwebh or igtbwebv
    added_weight_fraction::Float64 = 0
    CL_max::Float64 = 0
    volume::Float64 = 0
    size::Float64 = 0

    downwash_factor::Float64 = 0
    CL_max_fwd_CG::Float64 = 0
    SM_min::Float64 = 0
    CL_CLmax::Float64 = 0
    ntails::Float64 = 0
    move_wingbox::Float64 = 0
    
end