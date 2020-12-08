using Dierckx

function airfun(Ma, cl, τ, AMa, Acl, Aτ, ARe, A)
    
    ∂cdf_∂M  = Array{Spline1D, 2}(undef, nAcl, nAtau)
    ∂cdp_∂M  = Array{Spline1D, 2}(undef, nAcl, nAtau)
     ∂cm_∂M  = Array{Spline1D, 2}(undef, nAcl, nAtau)
    
    for k = 1:nAtau
        for j = 1:nAcl
            ∂cdf_∂M[j, k] = Spline1D(AMa, A[:,j,k,1], k = 3)
            ∂cdp_∂M[j, k] = Spline1D(AMa, A[:,j,k,2], k = 3)
             ∂cm_∂M[j, k] = Spline1D(AMa, A[:,j,k,3], k = 3)
        end
    end
    k = 1:nAtau
    j = 1:nAcl
    cdf_clτ = evaluate.(∂cdf_∂M[j,k], Ma) #cdf as a function of cl and τ
    cdp_clτ = evaluate.(∂cdp_∂M[j,k], Ma) #cdp as a function of cl and τ
     cm_clτ = evaluate.( ∂cm_∂M[j,k], Ma) #cm as a function of cl and τ

    ∂cdf_∂cl = Array{Spline1D, 1}(undef, nAtau)
    ∂cdp_∂cl = Array{Spline1D, 1}(undef, nAtau)
     ∂cm_∂cl = Array{Spline1D, 1}(undef, nAtau)
    
    for k = 1:nAtau
         ∂cdf_∂cl[k] = Spline1D(Acl, cdf_clτ[:,k], k = 3)
         ∂cdp_∂cl[k] = Spline1D(Acl, cdp_clτ[:,k], k = 3)
          ∂cm_∂cl[k] = Spline1D(Acl,  cm_clτ[:,k], k = 3)
    end

    cdf_τ = evaluate.(∂cdf_∂cl[k], cl)
    cdp_τ = evaluate.(∂cdp_∂cl[k], cl)
     cm_τ = evaluate.( ∂cm_∂cl[k], cl)

    ∂cdf_∂τ = Spline1D(Atau, cdf_τ, k = 3)
    ∂cdp_∂τ = Spline1D(Atau, cdp_τ, k = 3)
     ∂cm_∂τ = Spline1D(Atau,  cm_τ, k = 3)

    cdf = evaluate(∂cdf_∂τ, τ)
    cdp = evaluate(∂cdp_∂τ, τ)
     cm = evaluate( ∂cm_∂τ, τ)

    return cdf, cdp, cm

end #function airfun
