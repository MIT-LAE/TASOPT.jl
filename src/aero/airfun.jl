export airfun
const nAfun::Int = 3
const Ai = @MMatrix zeros(2,2)
const Ai_cl = @MMatrix zeros(2,2)
const Ai_tau = @MMatrix zeros(2,2)
const Ai_cl_tau = @MMatrix zeros(2,2)
const Aij = @MVector zeros(2)
const Aij_tau = @MVector zeros(2)
const Aijk = @MVector zeros(nAfun) # for the three vars cdf, cdp and cm

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
            Aij[kd] = interpolate(ix1=1, ix2=2, 
                dx=dcl, t=tcl, Y=view(Ai, :, kd), dYdX=view(Ai_cl, :, kd))
            
            Aij_tau[kd] = interpolate(ix1=1, ix2=2,
                dx=dcl, t = tcl, Y=view(Ai_tau, :, kd), dYdX=view(Ai_cl_tau,:,kd))
        end

        Aijk[l] = interpolate(ix1=1, ix2=2,
            dx=dtau, t=ttau, Y=view(Aij, :), dYdX=view(Aij_tau, :))
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
"""
    find_bisection(x::T, X::AbstractVector{T}) where T

Simple, fast bisection search in a **sorted** array. This does not sort the array for you.
"""
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

"""
    interpolate(;ix1, ix2, dx, t, Y, dYdX)

Convenience function to evaluate the symmetric form of the cubic spline equation, given
the x interval defined by the indices (ix1, ix2),  the interval (dx), normalized x location (t),
the Y array to be interpolated and the slopes (dYdX) at the interval end points.

Note: This is intentionally defined as a kw only argument list to encourage anyone working 
with this deep part of the code knows what they are doing. See [`eval_spline`](@ref).
"""
function interpolate(;ix1, ix2, dx, t, Y, dYdX)

    yim = Y[ix1] 
    yi = Y[ix2]
    Δ = yi - yim
    fxm = dx * dYdX[ix1] - Δ
    fxo = dx * dYdX[ix2] - Δ

    return eval_spline(t = t, yim = yim, yi = yi, fxm = fxm, fxo = fxo)
end

"""
    eval_spline(;t, yim, yi, fxm, fxo)

Given normalized interpolation location t, y₋ and y, calculates the
spline interpolation of the form    
    (1 - t)*y₁ + t*y₂ + t*(1 - t)*((1 - t)fxm + t*fx),
    where t = (x - x₁)/Δx.
"""
function eval_spline(;t, yim, yi, fxm, fxo)
    tmi = 1.0 - t
    return t*yi + tmi*yim + t*tmi*(tmi*fxm - t*fxo)
end