using Dierckx

function airfun(cl::Float64, τ::Float64, Ma::Float64, 
        AMa::Vector{Float64}, Acl::Vector{Float64}, Aτ::Vector{Float64}, ARe::Float64,
        A::Array{Float64},
        ∂cdf_∂M::Array{Spline1D, 2}, ∂cdp_∂M::Array{Spline1D, 2}, ∂cm_∂M::Array{Spline1D, 2})
    
    nAcl::Int  = length(Acl)
    nAtau::Int = length(Aτ)
    
    k = 1:nAtau
    j = 1:nAcl
    @inbounds cdf_clτ = evaluate.(∂cdf_∂M[j,k], Ma) #cdf as a function of cl and τ
    @inbounds cdp_clτ = evaluate.(∂cdp_∂M[j,k], Ma) #cdp as a function of cl and τ
    @inbounds  cm_clτ = evaluate.( ∂cm_∂M[j,k], Ma) #cm as a function of cl and τ

    ∂cdf_∂cl∂τ = Spline2D(Acl, Aτ, cdf_clτ; kx=3, ky=3, s=0.0)
    ∂cdp_∂cl∂τ = Spline2D(Acl, Aτ, cdp_clτ; kx=3, ky=3, s=0.0)
     ∂cm_∂cl∂τ = Spline2D(Acl, Aτ,  cm_clτ; kx=3, ky=3, s=0.0)

    cdf = evaluate(∂cdf_∂cl∂τ, cl, τ)
    cdp = evaluate(∂cdp_∂cl∂τ, cl, τ)
     cm = evaluate( ∂cm_∂cl∂τ, cl, τ)
    cdw = 0.0
    
    return cdf, cdp, cdw,  cm

end #function airfun
