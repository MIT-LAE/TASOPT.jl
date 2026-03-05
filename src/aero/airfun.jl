export airfun
const nAfun::Int = 6
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
      - `cdw::Float64`: Airfoil section wave drag (unused and assumed 0 here; placeholder left for future implementation).*
      - `cm::Float64`: Airfoil section pitching moment.
      - `alpha::Float64`: Airfoil section angle of attack (perpendicular to leading edge).*
      *assumed 0.0 for original `C.air` files
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
   
#=
 Explanation of variables, thanks Claude
    Variable	    Purpose
    io, im	        Upper/lower indices in Mach dimension
    jo, jm	        Upper/lower indices in CL dimension
    ko, km	        Upper/lower indices in τ dimension
    dMa, dcl, dtau	Interval sizes between indexed points [i.e., x2 - x1 in an interp]
    tMa, tcl, ttau	Normalized positions (0-1) within intervals [i.e., interp at (x_sample - x1) / (x2-x1)]
    A	            Core 4D database of aerodynamic coefficients
    A_M, A_τ, A_cl	First derivative matrices (by Mach, τ, CL respectively)
    A_M_τ, A_M_cl, 	Mixed partial derivatives
    A_cl_τ, A_M_cl_τ
    Ai, Ai_cl,      Intermediate interpolated values at Mach level
    Ai_tau, Ai_cl_tau	
    Aij, Aij_tau	Interpolated values at CL level
    Aijk	        Final interpolated value at target point
=#

    #TODO: rename the indices. 
    #find indices of known data points adjacent to query point
    io, im, dMa , tMa  = findsegment(Mach, AMa)
    jo, jm, dcl , tcl  = findsegment(cl, Acl)
    ko, km, dtau, ttau = findsegment(τ, Aτ)

    #Extrapolating in CL 
    #if any of the 4 elements in A[i,j,k,1] for combinations of {jo, jm} x {ko, km} 
    # are NaN (i.e., infeasible per MSES), use the next valid segments (i.e., lower cl points at given tau)
    # to perform the interpolation (thus, an extrapolation). 
    # drag performance is penalized later in this fxn to discourage such extrapolation (e.g., for optimization)

    # Check if any points are NaN, making extrapolation necessary
    any_nan = any(isnan.(A[im:io, jo-1:jo, ko-1:ko,1]))
    jupdate = nothing #initialize index for updated cl segment; subs for jo in case of NaNs
    if any_nan   
        # Find last valid cl segment at current tau to use instead of current segment
        for j_valid in jo-1:-1:1 #step backwards to lower cls
            #if all points relevant to tricubic interp are valid (i.e., none are invalid)...
            if !any(isnan.(A[im:io, j_valid, ko-1:ko, 1])) 
                jupdate = j_valid

                #update reference cl values for tricubic eval
                dcl = Acl[jupdate] - Acl[jupdate-1]
                tcl = (cl - Acl[jupdate-1])/dcl
                jo = jupdate
                break
            end
        end
    end
    
    #tricubic interpolation
    @inbounds for l = 1:nAfun #for every fxn to interpolate
        @inbounds for jd = 1:2  #loop over adjacent cl values
            @inbounds for kd = 1:2 #loop over adjacent tau values
                j = jo + jd-2
                k = ko + kd-2
                
                #interpolate vals at target Mach (along Mach)
                # values
                Ai[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A, :, j, k, l),
                    dYdX=view(A_M, :, j, k, l))
                
                #cl derivatives
                Ai_cl[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_cl, :, j, k, l),
                    dYdX=view(A_M_cl, :, j, k, l))
                
                #tau derivatives
                Ai_tau[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_τ, :, j, k, l),
                    dYdX=view(A_M_τ, :, j, k, l))
                
                #cross derivatives
                Ai_cl_tau[jd, kd] = interpolate(ix1=im, ix2=io,
                    dx=dMa, t=tMa, Y=view(A_cl_τ, :, j, k, l),
                    dYdX=view(A_M_cl_τ, :, j, k, l))
            end
        end

        #interpolate vals at target cl (along cl)
        @inbounds for kd = 1:2
            #values
            Aij[kd] = interpolate(ix1=1, ix2=2, 
                dx=dcl, t=tcl, Y=view(Ai, :, kd), dYdX=view(Ai_cl, :, kd))
            #tau derivatives
            Aij_tau[kd] = interpolate(ix1=1, ix2=2,
                dx=dcl, t = tcl, Y=view(Ai_tau, :, kd), dYdX=view(Ai_cl_tau,:,kd))
        end

        #interpolate vals at target tau (along tau)
        Aijk[l] = interpolate(ix1=1, ix2=2,
            dx=dtau, t=ttau, Y=view(Aij, :), dYdX=view(Aij_tau, :))
    end

    # extract data to return
    # updated airfoil database contains:
    # ["CD", "CDp", "CDv", "CDw", "CM", "alpha"]
    # (note that cd = cdv + cdw [viscous + wave] = cdp + cdf [pressure + friction])
    cdp = Aijk[2]           #pressure drag
    cdf = Aijk[1] - cdp     #friction drag
    cdw  = Aijk[4]*0          #wave drag #TODO: make nonzero after verification of csv import
    cm = Aijk[5]            #pitching moment
    alpha = Aijk[6]         #angle of attack, aoa

    #Add quadratic penalties for exceeding database cl or h/c limits
    clmax = isnothing(jupdate) ? Acl[end] : Acl[jupdate] #adjusts edge of database cl limits to account for NaNs 
    # clmax = Acl[end] #to turn off NaN penalty
    if cl < Acl[1]
        cdp = cdp + 1.0*(cl-Acl[1])^2
    elseif cl > clmax
        cdp = cdp + 1.0*(cl-clmax)^2
    end

    # τmin, τmax = Aτ[[1, end]]
    if τ < Aτ[1]
        cdp = cdp + 25.0*(τ - Aτ[1])^2
    elseif τ > Aτ[end]
        cdp = cdp + 25.0*(τ - Aτ[end])^2
    end

    return cdf, cdp, cdw, cm, alpha
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
        if x<=X[imid]
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

    yim = Y[ix1]           # y-value at lower endpoint
    yi = Y[ix2]            # y-value at upper endpoint
    Δ = yi - yim           # difference: Δ = y₂ - y₁
    fxm = dx * dYdX[ix1] - Δ   # "curvature" at lower point
    fxo = dx * dYdX[ix2] - Δ   # "curvature" at upper point

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