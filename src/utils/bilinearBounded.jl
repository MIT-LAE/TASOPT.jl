"""
        boundingIndices(x, array)

Gets the indices of the two array elements 1 and 2 such that x1 < x < x2 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x`:     `x` value to query
    - `array`: Array to search

    **Outputs:**
    - `iLow`: Index of lower bound (-1 if x less than all array elements)
    - `iUpp`: Index of upper bound (-1 if x greater than all array elements)
"""
function boundingIndices(x, array)
    indUpper = searchsortedfirst(array, x)

    if indUpper == length(array)
        return (indUpper, -1)
    elseif indUpper == 1
        return (-1, indUpper)
    else
        return indUpper-1, indUpper
    end
end # boundingIndices


"""
        bilinearBoundedLookup(x, y, xGrid, yGrid, lookupGrid)

Implementation of bilinear interpolation on a bounded grid. Used by `NcTblMap` and `ecTblMap`
to calculate off-design engine performance information. Custom implementation is required due
need for warning user operating at out-of-bounds conditions. Data is assumed to be defined on
a convex hull.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x`:          `x` position of query
    - `y`:          `y` position of query
    - `xGrid`:      1D-array of on-grid `x` values
    - `yGrid`:      1D-array of on-grid `y` values
    - `lookupGrid`: 2D-matrix of function values

    **Outputs:**
    - `val`: Bilinearly interpolated result (-1.0 if out-of-bounds or for impossible input)
"""
function bilinearBoundedLookup(x, y, xGrid, yGrid, lookupGrid)
    xIndLow, xIndUpp = boundingIndices(x, xGrid)
    yIndLow, yIndUpp = boundingIndices(y, yGrid)

    # println(xGrid[xIndLow], " ", xGrid[xIndUpp])

    if (xIndLow == -1 || xIndUpp == -1) || (yIndLow == -1 || yIndUpp == -1)
        # error()
        # @warn 
        # TODO: raise error for outside of interpolation grid
        return -1.0
    end

    # ---- Get function values at bounding points
    p1 = lookupGrid[yIndLow, xIndLow]
    p2 = lookupGrid[yIndUpp, xIndLow]
    p3 = lookupGrid[yIndUpp, xIndUpp]
    p4 = lookupGrid[yIndLow, xIndUpp]

    xLow, xUpp = xGrid[xIndLow], xGrid[xIndUpp]
    yLow, yUpp = yGrid[yIndLow], yGrid[yIndUpp]

    dx = xUpp - xLow
    dy = yUpp - yLow

    wxLow = (xUpp-x) / dx
    wxUpp = (x-xLow) / dx
    wyLow = (yUpp-y) / dy
    wyUpp = (y-yLow) / dy

    # ---- Is each point in bounds? (true if in; false if out)
    p1in = !(abs(p1+1.0) < 1.0e-6)
    p2in = !(abs(p2+1.0) < 1.0e-6)
    p3in = !(abs(p3+1.0) < 1.0e-6)
    p4in = !(abs(p4+1.0) < 1.0e-6)
    boundSum = p1in + p2in + p3in + p4in

    if boundSum == 0 # ---- Point lies outside map
        return -1.0
    end

    if boundSum == 1 # ---- Single point in bounds; return funciton value at that point
        if p1in
            return p1
        elseif p2in
            return p2
        elseif p3in
            return p3
        elseif p4in
            return p4
        else
            # Something's gone wrong
            return -1.0
        end
    end 

    if boundSum == 2 # ---- Two points in bounds
        
        # ---- Prepare useful quantities
        if p1in && p2in # Left vertical
            return p1 * wyLow + p2 * wyUpp
        elseif p4in && p3in
            return p4 * wyLow + p3 * wyUpp
        elseif p2in && p3in
            return p2 * wxLow + p3 * wxUpp
        elseif p1in && p4in
            return p1 * wxLow + p2 * wxUpp
        else # ---- Diagonal points in bounds; this should never be able to happen
            return -1.0
        end

    elseif boundSum == 3 # ---- Three points in bounds
        # ---- Map each three point configuration to (corner, vertical, horizontal) indices on
        # ---- the filtered arrays
        if (p1in && p2in) && p3in # Upper Left Cornered
            fx = p2 * wxLow + p3 * wxUpp
            fy = p2 * wyUpp + p1 * wyLow

            fxy = fx * wyUpp + p1 * wyLow
            fyx = fy * wxLow + p3 * wxUpp
        elseif (p2in && p3in) && p4in # Upper Right Cornered
            fx = p3 * wxUpp + p2 * wxLow
            fy = p3 * wyUpp + p4 * wyLow

            fxy = fx * wyUpp + p4 * wyLow
            fyx = fy * wxUpp + p2 * wxLow
        elseif (p1in && p2in) && p4in # Lower Left Cornered
            fx = p1 * wxLow + p4 * wxUpp
            fy = p1 * wyLow + p2 * wyUpp

            fxy = fx * wyLow + p2 * wyUpp
            fyx = fy * wxLow + p4 * wxUpp
        elseif (p1in && p3in) && p4in # Lower Right Cornered
            fx = p4 * wxUpp + p1 * wxLow
            fy = p4 * wyLow + p3 * wyUpp

            fxy = fx * wyLow + p3 * wyUpp
            fyx = fy * wxUpp + p1 * wxLow
        else
            return -1.0
            # TODO raise error - this shouldn't be possible
        end

        return (fxy + fyx) / 2

    elseif boundSum == 4 # ---- Four points in bounds
        # ---- Calculate the linear interpolation of function value along the `y`-axes
        fy12 = p1 * wyLow + p2 * wyUpp
        fy34 = p4 * wyLow + p3 * wyUpp

        # ---- Perform the bilinear interpolation
        fapprox = fy12 * wxLow + fy34 * wxUpp
        return fapprox
    
    else # ---- Something has gone very wrong
        # TODO: Raise error/warning
        return -1.0
    end

end # bilinearBoundedLookup


"""
        bilinearDerivatives(x, y, xGrid, yGrid, lookupGrid; eps=1.0e-6)

Calculates the nominal value and x/y derivatives of a function defined on a 2D
grid using bilinear interpolation.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x`:          `x` position of query
    - `y`:          `y` position of query
    - `xGrid`:      1D-array of on-grid `x` values
    - `yGrid`:      1D-array of on-grid `y` values
    - `lookupGrid`: 2D-matrix of function values

    **Outputs:**
    - `f`: Bilinearly interpolated result (-1.0 if out-of-bounds or for impossible input)
    - `f_?`: Central derivatives of f at ?*(1+eps) and ?*(1-eps)

"""
function bilinearDerivatives(xNom, yNom, xGrid, yGrid, lookupGrid; eps=1.0e-3)
    xUpp, xLow = xNom * (1+eps), xNom * (1-eps)
    yUpp, yLow = yNom * (1+eps), yNom * (1-eps)

    f      = bilinearBoundedLookup(xNom, yNom, xGrid, yGrid, lookupGrid)

    if (f + 1.0 < 1.0e-6)
        return -1.0, -1.0, -1.0
    else
        f_xUpp = bilinearBoundedLookup(xUpp, yNom, xGrid, yGrid, lookupGrid)
        f_xLow = bilinearBoundedLookup(xLow, yNom, xGrid, yGrid, lookupGrid)
        # println(xNom, " ", xUpp, " ", f, " ", f_xUpp)

        f_yUpp = bilinearBoundedLookup(xNom, yUpp, xGrid, yGrid, lookupGrid)
        f_yLow = bilinearBoundedLookup(xNom, yLow, xGrid, yGrid, lookupGrid)

        if f_xUpp == -1.0 || f_xLow == -1.0
            f_x = -1.0
        else
            f_x = (f_xUpp - f_xLow) / (2 * xNom * eps)
        end

        if f_yUpp == -1.0 || f_yLow == -1.0
            f_y = -1.0
        else
            f_y = (f_yUpp - f_yLow) / (2 * yNom * eps)
        end

        return f, f_x, f_y
    end
end # bilinearDerivatives