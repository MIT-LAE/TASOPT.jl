@kwdef mutable struct StructuralMember
    material::StructuralAlloy = StructuralAlloy("TASOPT-Al")
    weight::Float64 = 0 #OR Union{Float64, Nothing}
    σ::Float64 = material.σmax
    EIh::Float64 = 0
    EIv::Float64 = 0
    GJ::Float64 = 0
    thickness::Float64 = 0
    ρ::Float64 = material.ρ
    x::Float64 = 0
    # xv::Float64 = 0
    #Material = StructuralAlloy
end