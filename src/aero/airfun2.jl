export airfun

"""
    airfun(cl, τ, Mach, air::airfoil)

Looks up airfoil performance data at specified conditions, as precomputed and found in `./src/air/C.air`.

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
                Acl::AbstractVector{Float64},
                Aτ::AbstractVector{Float64},
                AMa::AbstractVector{Float64},
                A::AbstractArray{Float64},
                A_M::AbstractArray{Float64},
                A_τ::AbstractArray{Float64},
                A_cl::AbstractArray{Float64},
                A_M_τ::AbstractArray{Float64},
                A_M_cl::AbstractArray{Float64},
                A_cl_τ::AbstractArray{Float64},
                A_M_cl_τ::AbstractArray{Float64})

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
    
    @inbounds for l = 1:nAfun
        @inbounds for jd = 1:2
            @inbounds for kd = 1:2
                j = jo + jd-2
                k = ko + kd-2

                @views        Ai[jd,kd] = SEVAL(Mach,      A[:, j, k, l],      A_M[:, j, k, l], AMa)
                @views     Ai_cl[jd,kd] = SEVAL(Mach,   A_cl[:, j, k, l],   A_M_cl[:, j, k, l], AMa)
                @views    Ai_tau[jd,kd] = SEVAL(Mach,    A_τ[:, j, k, l],    A_M_τ[:, j, k, l], AMa)
                @views Ai_cl_tau[jd,kd] = SEVAL(Mach, A_cl_τ[:, j, k, l], A_M_cl_τ[:, j, k, l], AMa)

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
function findsegment(x::Float64, xarr::AbstractArray{Float64})

    if isnan(x)
        error("Oops, you're searching for a NaN! Go fix your bug!")
    end
    
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
