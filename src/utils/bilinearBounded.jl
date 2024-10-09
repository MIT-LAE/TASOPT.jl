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
        # TODO: raise error for outside of interpolation grid
        return -1.0
    end

    lookupIndices = [[xIndLow, yIndLow], [xIndLow, yIndUpp], [xIndUpp, yIndUpp], [xIndUpp, yIndLow]]

    # ---- Get function values at bounding points
    p1 = lookupGrid[yIndLow, xIndLow]
    p2 = lookupGrid[yIndUpp, xIndLow]
    p3 = lookupGrid[yIndUpp, xIndUpp]
    p4 = lookupGrid[yIndLow, xIndUpp]
    pts = [p1, p2, p3, p4]

    # ---- Is each point in bounds? (true if in; false if out)
    p1in = !(abs(p1+1.0) < 1.0e-6)
    p2in = !(abs(p2+1.0) < 1.0e-6)
    p3in = !(abs(p3+1.0) < 1.0e-6)
    p4in = !(abs(p4+1.0) < 1.0e-6)
    ptsIn = [p1in, p2in, p3in, p4in]
    boundSum = p1in + p2in + p3in + p4in
    # println(boundSum)

    filterIndIn = findall(x->x, ptsIn)
    filterPts = pts[filterIndIn]
    filterLookup = lookupIndices[filterIndIn]
    # println(filterLookup)

    if boundSum == 0 # ---- Point lies outside map
        return -1.0
    end

    if boundSum == 1 # ---- Single point in bounds; return funciton value at that point
        return filterPts[1]
    end 

    if boundSum == 2 # ---- Two points in bounds
        # ---- Prepare useful quantities
        xInd1 = filterLookup[1][1]
        xInd2 = filterLookup[2][1]
        yInd1 = filterLookup[1][2]
        yInd2 = filterLookup[2][2]

        x1 = xGrid[xInd1]
        x2 = xGrid[xInd2]
        y1 = yGrid[yInd1]
        y2 = yGrid[yInd2]

        dx = x2 - x1
        dy = y2 - y1

        # ---- Calculate weights based on point orientation
        if filterInd == [1,2] || filterInd == [3,4] # ---- Vertical points in bounds
            weight2 = abs((y - y1)/dy)
            weight1 = abs((y - y2)/dy)
        elseif filterInd == [2,3] || filterInd == [1,4] # ---- Horizontal points in bounds
            weight2 = abs((x - x1)/dx)
            weight1 = abs((x - x2)/dx)
        else # ---- Diagonal points in bounds; this can never happen on a convex hull of points
            return -1.0
        end

        # ---- Return the linear interpolation of the two points
        return filterPts[1] * weight1 + filterPts[2] * weight2

    elseif boundSum == 3 # ---- Three points in bounds
        # ---- Map each three point configuration to (corner, vertical, horizontal) indices on
        # ---- the filtered arrays
        if filterIndIn == [1,2,3]
            cornerInd = 2
            vertInd = 1
            horInd = 3
        elseif filterIndIn == [2,3,4]
            cornerInd = 2
            vertInd = 3
            horInd = 1
        elseif filterIndIn == [1,2,4]
            cornerInd = 1
            vertInd = 2
            horInd = 3
        elseif filterIndIn == [1,3,4]
            cornerInd = 3
            vertInd = 2
            horInd = 1
        else
            return -1.0
            # TODO raise error - this shouldn't be possible
        end

        # ---- Get x and y values for each point
        xCorner, yCorner = xGrid[ filterLookup[cornerInd][1] ], yGrid[ filterLookup[cornerInd][2] ]
        yVert            =                                      yGrid[ filterLookup[vertInd][2]   ]
        xHor             = xGrid[ filterLookup[horInd][1]    ]

        # ---- Calculate the length and area infinitesimals
        dx = abs(xHor  - xCorner)
        dy = abs(yVert - yCorner)

        # ---- Get the function value at each point
        fCorner, fVert, fHor = filterPts[cornerInd], filterPts[vertInd], filterPts[horInd]

        # Calculate the linear interpolation of function value along each axis
        fx = abs(x - xHor)  / dx * fCorner + abs(x - xCorner) / dx * fHor
        fy = abs(y - yVert) / dy * fCorner + abs(y - yCorner) / dy * fVert

        # Calculate the linear interpolation of function and axis interpolation along the other axis
        fxy = abs(y - yVert) / dy * fx + abs(y - yCorner) / dy * fVert
        fyx = abs(x - xHor ) / dx * fy + abs(x - xCorner) / dx * xHor

        # Average the two semi-bilinear interpolations to approximate the true function value
        return (fxy + fyx) / 2

    elseif boundSum == 4 # ---- Four points in bounds
        x1, y1 = xGrid[filterLookup[1][1]], yGrid[filterLookup[1][2]]
        x2, y2 = xGrid[filterLookup[2][1]], yGrid[filterLookup[2][2]]
        x3, y3 = xGrid[filterLookup[3][1]], yGrid[filterLookup[3][2]]
        x4, y4 = xGrid[filterLookup[4][1]], yGrid[filterLookup[4][2]]
        
        # ---- Get differentials assuming a regular grid
        dx = abs(x3 - x2)
        dy = abs(y2 - y1)

        # ---- Calculate the linear interpolation of function value along the `y`-axes
        fy12 = abs(y - y2) / dy * p1 + abs(y - y1) / dy * p2
        fy34 = abs(y - y4) / dy * p3 + abs(y - y3) / dy * p4

        # ---- Perform the bilinear interpolation
        fapprox = (abs(x-x3) / dx) * fy12 + (abs(x-x1) / dx) * fy34
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