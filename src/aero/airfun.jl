export airfun

"""
    airfun(cl, τ, Mach, air::airfoil)

Looks up airfoil performance data at specified conditions, as precomputed and found in `./src/airfoil_data/`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**      
      - `cl::Float64`: Airfoil section lift coefficient.
      - `τ::Float64`: Airfoil section thickness-to-chord ratio.
      - `Mach::Float64`: Mach number.
      - `air::TASOPT.aerodynamics.airfoil`: `airfoil` structure with performance data.

      **Outputs:**
      - `cdf::Float64`: Airfoil section skin friction drag.
      - `cdp::Float64`: Airfoil section pressure drag.
      - `cdw::Float64`: Airfoil section wave drag (unused and assumed 0 here; placeholder left for future implementation).
      - `cm::Float64`: Airfoil section pitching moment.
"""
@views function airfun(cl, τ, Mach, air::airfoil)

    return airfun(cl, τ, Mach, 
                    air.cl, air.τ, air.Ma,
                    air.A,
                    air.A_M,
                    air.A_τ,
                    air.A_cl,
                    air.A_M_τ,
                    air.A_M_cl,
                    air.A_cl_τ,
                    air.A_M_cl_τ)

end

@views function airfun(cl, τ, Mach, 
                Acl,
                Aτ,
                AMa,
                A,
                A_M,
                A_τ,
                A_cl,
                A_M_τ,
                A_M_cl,
                A_cl_τ,
                A_M_cl_τ)

    nAfun::Int = 3
    Ai = @MMatrix zeros(2,2)
    Ai_cl = @MMatrix zeros(2,2)
    Ai_tau = @MMatrix zeros(2,2)
    Ai_cl_tau = @MMatrix zeros(2,2)
    Aij = @MVector zeros(2)
    Aij_tau = @MVector zeros(2)
    Aijk = @MVector zeros(nAfun) # for the three vars cdf, cdp and cm
    
    io, im, dMa , tMa  = findsegment(Mach, AMa)
    jo, jm, dcl , tcl  = findsegment(cl, Acl)
    ko, km, dtau, ttau = findsegment(τ, Aτ)
    
    @inbounds for l = 1:nAfun
        @inbounds for jd = 1:2
            @inbounds for kd = 1:2
                j = jo + jd-2
                k = ko + kd-2

                Ai[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A, :, j, k, l),
                    dYdX=view(A_M, :, j, k, l))
                
                Ai_cl[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_cl, :, j, k, l),
                    dYdX=view(A_M_cl, :, j, k, l))
                
                Ai_tau[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_τ, :, j, k, l),
                    dYdX=view(A_M_τ, :, j, k, l))
                
                Ai_cl_tau[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_cl_τ, :, j, k, l),
                    dYdX=view(A_M_cl_τ, :, j, k, l))

            end
        end

        @inbounds for kd = 1:2
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
function findsegment(x::Float64, xarr)

    if isnan(x)
        error("Oops, you're searching for a NaN! Go fix your bug!")
    end
    
    im::Int = find_bisection(x, xarr)
    io = im+1

    dx = xarr[io] - xarr[im]
    t = (x - xarr[im])/dx

    return io, im, dx, t
end

function find_bisection(x::T, X::AbstractVector{T}) where T
    ilow::Int = 1
    i::Int = length(X)
    while (i-ilow > 1)
        imid = (i+ilow)÷2
        if x<X[imid]
            i = imid
        else
            ilow = imid
        end
    end
    return i-1
end  # function find_bisection

function interpolate(;ix1, ix2, dx, t,
    Y, dYdX)

    yim = Y[ix1] 
    yi = Y[ix2]
    Δ = yi - yim
    fxm = dx * dYdX[ix1] - Δ
    fxo = dx * dYdX[ix2] - Δ

    return eval_spline(t = t, yim = yim, yi = yi, fxm = fxm, fxo = fxo)
end

function eval_spline(;t, yim, yi, fxm, fxo)
    tmi = 1.0 - t
    return t*yi + tmi*yim + t*tmi*(tmi*fxm - t*fxo)
end