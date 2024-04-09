@kwdef mutable struct StructuralMember
    weight::Float64 = 0 #OR Union{Float64, Nothing}
    σ::Float64 = 0
    EIh::Float64 = 0
    EIv::Float64 = 0
    GJ::Float64 = 0
    thickness::Float64 = 0
    ρ::Float64 = 0
    xh::Float64 = 0
    xv::Float64 = 0
end