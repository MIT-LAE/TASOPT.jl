@kwdef mutable struct Tail
    layout::TailLayout = TailLayout()

    weight::Float64 = 0
    dxW::Float64 = 0
    EI_bending::Float64 = 0 #igEIch or igEIcv
    EI_normal::Float64 = 0 #igEInh or igEInv
    GJ::Float64 = 0 #igGJh or igGJv
    thickness_cap::Float64 = 0 #igtbcaph or igtbcapv
    thickness_web::Float64 = 0 #igtbwebh or igtbwebv
end