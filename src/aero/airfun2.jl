using Base: NamedTuple_typename
@views function airfun(cl, τ, Mach, 
                Acl, Aτ, AMa,
                A,
                A_M,
                A_τ,
                A_cl,
                A_M_τ,
                A_M_cl,
                A_cl_τ,
                A_M_cl_τ)

    nAfun::Int = 3

    Ai = zeros(2,2)
    Ai_cl = zeros(2,2)
    Ai_tau = zeros(2,2)
    Ai_cl_tau = zeros(2,2)
    Aij = zeros(2)
    Aij_tau = zeros(2)
    Aijk = zeros(nAfun) # for the three vars cdf, cdp and cm
    
    # io, im, dMa , tMa  = findsegment(Mach, AMa)
    jo, jm, dcl , tcl  = findsegment(cl, Acl)
    ko, km, dtau, ttau = findsegment(τ, Aτ)
    
    for l = 1:nAfun

        for jd = 1:2
            for kd = 1:2
                j = jo + jd-2
                k = ko + kd-2

                       Ai[jd,kd] = SEVAL(Mach,      A[:, j, k, l],      A_M[:, j, k, l], AMa)
                    Ai_cl[jd,kd] = SEVAL(Mach,   A_cl[:, j, k, l],   A_M_cl[:, j, k, l], AMa)
                   Ai_tau[jd,kd] = SEVAL(Mach,    A_τ[:, j, k, l],    A_M_τ[:, j, k, l], AMa)
                Ai_cl_tau[jd,kd] = SEVAL(Mach, A_cl_τ[:, j, k, l], A_M_cl_τ[:, j, k, l], AMa)

                # fxm = dMa*A_M[im,j,k,l] - A[io,j,k,l] + A[im,j,k,l]
                # fxo = dMa*A_M[io,j,k,l] - A[io,j,k,l] + A[im,j,k,l]
                # Ai[jd, kd] = tMa*A[io,j,k,l] + (1 - tMa)*A[im,j,k,l] + tMa*(1 - tMa)*((1-tMa)*fxm - tMa*fxo)
            end
        end

        for kd = 1:2
            Δ = Ai[2, kd] - Ai[1, kd]
            fxm = dcl*Ai_cl[1, kd] - Δ
            fxo = dcl*Ai_cl[2, kd] - Δ
            Aij[kd] = tcl*Ai[2,kd] + (1-tcl)*Ai[1,kd] + tcl*(1-tcl)*((1-tcl)*fxm - tcl*fxo)

            Δ = Ai_tau[2, kd] - Ai_tau[1, kd]
            fxm = dcl*Ai_cl_tau[1, kd] - Δ
            fxo = dcl*Ai_cl_tau[2, kd] - Δ
            Aij_tau[kd] = tcl*Ai_tau[2,kd] + (1-tcl)*Ai_tau[1,kd] + tcl*(1-tcl)*((1-tcl)*fxm - tcl*fxo)
        end

        Δ = Aij[2] - Aij[1]
        fxm = dtau*Aij_tau[1] - Δ
        fxo = dtau*Aij_tau[2] - Δ
        Aijk[l] = ttau*Aij[2] + (1-ttau)*Aij[1] + ttau*(1-ttau)*((1-ttau)*fxm - ttau*fxo)

    end
    cdf = Aijk[1]
    cdp = Aijk[2]
    cm  = Aijk[3]
    cdw = 0.0

    #Add quadratic penalties for exceeding database cl or h/c limits
    # clmin, clmax = Acl[[1, end]]
    if cl < Acl[1]
        cdp = cdp + 1.0*(cl-Acl[1])^2
    elseif cl > Acl[end]
        cdp = cdp + 1.0*(cl-Acl[end])^2
    end

    # τmin, τmax = Aτ[[1, end]]
    if τ < Aτ[1]
        cdp = cdp + 25.0*(τ - Aτ[1])^2
    elseif τ > Aτ[end]
        cdp = cdp + 25.0*(τ - Aτ[end])^2
    end

    return cdf, cdp, cdw, cm

end


"""
    findsegment(x::T, xarr::Vector{T})

Uses bisection to find the right interval of the array `xarr` where x lies.
Returns im and io s.t. `xarr[im] < x < xarr[io]`.

Additionally returns the interval `dx = xarr[io] - xarr[im]`
"""
function findsegment(x::Float64, xarr::AbstractArray{Float64})

    io::Int = length(xarr)

    if x ≤ xarr[1]
        im = 1
        io = 2
    elseif x ≥ xarr[end]
        im = io - 1
    else
        im = searchsortedlast(xarr, x)
        io = im + 1
    end

    dx = xarr[io] - xarr[im]
    t = (x - xarr[im])/dx

    return io, im, dx, t
end
