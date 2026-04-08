export airfun, airfoil_cl_limits
const nAfun::Int = 6
const Ai = @MMatrix zeros(2,2)
const Ai_cl = @MMatrix zeros(2,2)
const Ai_toc = @MMatrix zeros(2,2)
const Ai_cl_toc = @MMatrix zeros(2,2)
const Aij = @MVector zeros(2)
const Aij_toc = @MVector zeros(2)
const Aijk = @MVector zeros(nAfun) # for the three vars cdf, cdp and cm

"""
    airfun(cl, toc, Mach, air::airfoil)

Looks up airfoil performance data at specified conditions, as precomputed and found in `./src/airfoil_data/`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**      
      - `cl::Float64`: Airfoil section lift coefficient.
      - `toc::Float64`: Airfoil section thickness-to-chord ratio.
      - `Mach::Float64`: Mach number.
      - `air::TASOPT.aerodynamics.airfoil`: `airfoil` structure with performance data.

      **Outputs:**
      - `cdf::Float64`: Airfoil section skin friction drag.
      - `cdp::Float64`: Airfoil section pressure drag.
      - `cdw::Float64`: Airfoil section wave drag (unused and assumed 0 here; placeholder left for future implementation).*
      - `cm::Float64`: Airfoil section pitching moment.
      - `alpha::Float64`: Airfoil section angle of attack (perpendicular to leading edge).*
"""
@views function airfun(cl, toc, Mach, air::airfoil)

    return airfun(cl, toc, Mach, 
                    air.cl, air.toc, air.Ma,
                    air.A,
                    air.A_M,
                    air.A_toc,
                    air.A_cl,
                    air.A_M_toc,
                    air.A_M_cl,
                    air.A_cl_toc,
                    air.A_M_cl_toc)

end

@views function airfun(cl, toc, Mach, 
                Acl,
                Atoc,
                AMa,
                A,
                A_M,
                A_toc,
                A_cl,
                A_M_toc,
                A_M_cl,
                A_cl_toc,
                A_M_cl_toc)
   
#=
 Explanation of variables, thanks Claude
    Variable	    Purpose
    i2, i1	        Upper/lower indices in Mach dimension
    j2, j1	        Upper/lower indices in CL dimension
    k2, k1	        Upper/lower indices in toc dimension
    dMa, dcl, dtoc	Interval sizes between indexed points [i.e., x2 - x1 in an interp]
    tMa, tcl, ttoc	Normalized positions (0-1) within intervals [i.e., interp at (x_sample - x1) / (x2-x1)]
    A	            Core 4D database of aerodynamic coefficients
    A_M, A_toc, A_cl	First derivative matrices (by Mach, toc, CL respectively)
    A_M_toc, A_M_cl, 	Mixed partial derivatives
    A_cl_toc, A_M_cl_toc
    Ai, Ai_cl,      Intermediate interpolated values at Mach level
    Ai_toc, Ai_cl_toc	
    Aij, Aij_toc	Interpolated values at CL level
    Aijk	        Final interpolated value at target point
=#

    #TODO: rename the indices. 
    #find indices of known data points adjacent to query point
    i2, i1, dMa , tMa  = findsegment(Mach, AMa)
    j2, j1, dcl , tcl  = findsegment(cl, Acl)
    k2, k1, dtoc, ttoc = findsegment(toc, Atoc)

    #Extrapolating in CL 
    #if any of the 4 elements in A[i,j,k,1] for combinations of {j2, j1} x {k2, k1} 
    # are NaN (i.e., infeasible per MSES), use the next valid segments (i.e., lower cl points at given toc)
    # to perform the interpolation (thus, an extrapolation). 
    # drag performance is penalized later in this fxn to discourage such extrapolation (e.g., for optimization)

    # Check if any points are NaN, making extrapolation necessary
    any_nan = any(isnan.(A[i1:i2, j1:j2, k1:k2,1]))
    jupdate = nothing #initialize index for updated cl segment; subs for j2 in case of NaNs
    if any_nan   
        # Find last valid cl segment at current toc to use instead of current segment
        for j_valid in j2-1:-1:1 #step backwards to lower cls
            #if all points relevant to tricubic interp are valid (i.e., none are invalid)...
            if !any(isnan.(A[i1:i2, j_valid, k1:k2, 1])) 
                jupdate = j_valid

                #update reference cl values for tricubic eval
                dcl = Acl[jupdate] - Acl[jupdate-1]
                tcl = (cl - Acl[jupdate-1])/dcl
                j2 = jupdate
                break
            end
        end
    end
    
    #tricubic interpolation
    @inbounds for l = 1:nAfun #for every fxn to interpolate
        @inbounds for jd = 1:2  #loop over adjacent cl values
            @inbounds for kd = 1:2 #loop over adjacent toc values
                j = j2 + jd-2
                k = k2 + kd-2
                
                #interpolate vals at target Mach (along Mach)
                # values
                Ai[jd, kd] = interpolate(ix1=i1, ix2=i2,
                    dx=dMa, t=tMa, Y=view(A, :, j, k, l),
                    dYdX=view(A_M, :, j, k, l))
                
                #cl derivatives
                Ai_cl[jd, kd] = interpolate(ix1=i1, ix2=i2,
                    dx=dMa, t=tMa, Y=view(A_cl, :, j, k, l),
                    dYdX=view(A_M_cl, :, j, k, l))
                
                #toc derivatives
                Ai_toc[jd, kd] = interpolate(ix1=i1, ix2=i2,
                    dx=dMa, t=tMa, Y=view(A_toc, :, j, k, l),
                    dYdX=view(A_M_toc, :, j, k, l))
                
                #cross derivatives
                Ai_cl_toc[jd, kd] = interpolate(ix1=i1, ix2=i2,
                    dx=dMa, t=tMa, Y=view(A_cl_toc, :, j, k, l),
                    dYdX=view(A_M_cl_toc, :, j, k, l))
            end
        end

        #interpolate vals at target cl (along cl)
        @inbounds for kd = 1:2
            #values
            Aij[kd] = interpolate(ix1=1, ix2=2, 
                dx=dcl, t=tcl, Y=view(Ai, :, kd), dYdX=view(Ai_cl, :, kd))
            #toc derivatives
            Aij_toc[kd] = interpolate(ix1=1, ix2=2,
                dx=dcl, t = tcl, Y=view(Ai_toc, :, kd), dYdX=view(Ai_cl_toc,:,kd))
        end

        #interpolate vals at target toc (along toc)
        Aijk[l] = interpolate(ix1=1, ix2=2,
            dx=dtoc, t=ttoc, Y=view(Aij, :), dYdX=view(Aij_toc, :))
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

    # tocmin, tocmax = Atoc[[1, end]]
    if toc < Atoc[1]
        cdp = cdp + 25.0*(toc - Atoc[1])^2
    elseif toc > Atoc[end]
        cdp = cdp + 25.0*(toc - Atoc[end])^2
    end

    return cdf, cdp, cdw, cm, alpha

"""
    airfoil_cl_limits(airf::airfoil, Mach_perp::Float64, toc::Float64)

Returns the valid `(cl_min, cl_max)` range of the airfoil database at the given
perpendicular Mach number and thickness-to-chord ratio.

Validity is determined by checking for `NaN` entries in the CD column of the data
array `A`, which are used to mark infeasible (e.g. stalled) operating points in the
MSES database. The bracketing Mach and toc grid cells are found with
`searchsortedlast`; inputs outside the grid are clamped to the nearest cell so the
function does not error on extrapolation queries.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `airf::airfoil`: Airfoil database struct (e.g. `ac.wing.airsection`).
      - `Mach_perp::Float64`: Perpendicular Mach number (= aircraft Mach × cos(sweep)).
      - `toc::Float64`: Thickness-to-chord ratio at the section of interest.

      **Outputs:**
      - `cl_min::Float64`: Lowest cl in the database with valid (non-NaN) data at `(Mach_perp, toc)`.
      - `cl_max::Float64`: Highest cl in the database with valid (non-NaN) data at `(Mach_perp, toc)`.
"""
function airfoil_cl_limits(airf::airfoil, Mach_perp::Float64, toc::Float64)
    i1 = clamp(searchsortedlast(airf.Ma,  Mach_perp), 1, length(airf.Ma)  - 1)
    k1 = clamp(searchsortedlast(airf.toc, toc),       1, length(airf.toc) - 1)
    valid = j -> !any(isnan, @view airf.A[i1:i1+1, j, k1:k1+1, 1])
    jmin  = findfirst(valid, eachindex(airf.cl))
    jmax  = findlast( valid, eachindex(airf.cl))
    return (isnothing(jmin) ? airf.cl[1]   : airf.cl[jmin],
            isnothing(jmax) ? airf.cl[end] : airf.cl[jmax])
end


"""
    findsegment(x::T, xarr::Vector{T})

Uses bisection to find the right interval of the array `xarr` where x lies.
Returns i1 and i2 s.t. `xarr[i1] < x < xarr[i2]`.

Additionally returns the interval `dx = xarr[i2] - xarr[i1]`
"""
function findsegment(x::Float64, xarr)

    if isnan(x)
        error("Oops, you're searching for a NaN! Go fix your bug!")
    end
    
    i1::Int = find_bisection(x, xarr)
    i2 = i1+1

    dx = xarr[i2] - xarr[i1]
    t = (x - xarr[i1])/dx

    return i2, i1, dx, t
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