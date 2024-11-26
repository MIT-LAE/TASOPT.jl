function regGridBoundingIndices(x::Float64, dx::Float64, N::Int, xLow::Float64)
    fracInd = (x - xLow) / dx + 1
    remainder = fracInd % 1
    clamp = 0
    if fracInd > N
        # TODO: Add upper OOB warning
        clamp = 1
        return N-1, N, clamp
    end
    if fracInd < 1
        # TODO: Add lower OOB warning
        clamp = -1
        return 1, 2, clamp
    end
    
    if isapprox(remainder, 0.)
        iLow, iHigh = Int(round(fracInd)), Int(round(fracInd))
    else
        iLow, iHigh = Int(floor(fracInd)), Int(ceil(fracInd))
    end

    return iLow, iHigh, clamp
end # regGridBoundingIndices


"""
        bilinearBoundedLookup(x, y, xGrid, yGrid, nominal)

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
    - `nominal`: 2D-matrix of function values

    **Outputs:**
    - `val`: Bilinearly interpolated result (-1.0 if out-of-bounds or for impossible input)
"""
function bilinearBoundedLookup(x::Float64, y::Float64, dx::Float64, dy::Float64, Nx::Int, Ny::Int, 
                               xGrid::AbstractArray{Float64}, yGrid::AbstractArray{Float64}, 
                               nominal::AbstractArray{Float64}, dqdx::AbstractArray{Float64}, 
                               dqdy::AbstractArray{Float64})
    xLow, yLow = xGrid[1], yGrid[1]
    xUpp, yUpp = xGrid[Nx], yGrid[Ny]

    xIndLow, xIndUpp, xClamp = regGridBoundingIndices(x, dx, Nx, xLow)
    yIndLow, yIndUpp, yClamp = regGridBoundingIndices(y, dy, Ny, yLow)

    # ---- Clamp the x and y values to the grid bounding points to correspond to the function value clamping
    if xClamp == 1
        x = xUpp
    elseif xClamp == -1
        x = xLow
    end

    if yClamp == 1
        y = yUpp
    elseif yClamp == -1
        y = yLow
    end

    # ---- Get function values at bounding points
    p1 = nominal[xIndLow, yIndLow]
    p2 = nominal[xIndUpp, yIndLow]
    p3 = nominal[xIndUpp, yIndUpp]
    p4 = nominal[xIndLow, yIndUpp]

    p1x = dqdx[xIndLow, yIndLow]
    p2x = dqdx[xIndUpp, yIndLow]
    p3x = dqdx[xIndUpp, yIndUpp]
    p4x = dqdx[xIndLow, yIndUpp]

    p1y = dqdy[xIndLow, yIndLow]
    p2y = dqdy[xIndUpp, yIndLow]
    p3y = dqdy[xIndUpp, yIndUpp]
    p4y = dqdy[xIndLow, yIndUpp]

    xLow, xUpp = xGrid[xIndLow], xGrid[xIndUpp]
    yLow, yUpp = yGrid[yIndLow], yGrid[yIndUpp]

    wxLow = (xUpp-x) / dx
    wxUpp = (x-xLow) / dx
    wyLow = (yUpp-y) / dy
    wyUpp = (y-yLow) / dy

    # ---- Handle cases where a grid point is matched exactly
    # ---- Both x and y grids matched
    if xIndLow == xIndUpp && yIndLow == yIndUpp
        fapprox  = p1
        dxApprox = p1x
        dyApprox = p1y
    # ---- Matching only the x grid point; interpolate along y axis
    elseif xIndLow == xIndUpp && yIndLow != yIndUpp
        fapprox  = p1  * wyLow + p2  * wyUpp
        dxApprox = p1x * wyLow + p2x * wyUpp
        dyApprox = p1y * wyLow + p2y * wyUpp
    # ---- Matching only the y grid point; interpolate along x axis
    elseif xIndLow != xIndUpp && yIndLow == yIndUpp
        fapprox  = p1  * wxLow + p4  * wxUpp
        dxApprox = p1x * wxLow + p4x * wxUpp
        dyApprox = p1y * wxLow + p4y * wxUpp
    # ---- Not matching a grid point
    else
        # ---- Calculate the linear interpolation of function value along the `y`-axes
        fy12  = p1  * wyLow + p2  * wyUpp
        fy34  = p4  * wyLow + p3  * wyUpp

        dxy12 = p1x * wyLow + p2x * wyUpp
        dxy34 = p4x * wyLow + p3x * wyUpp

        dyy12 = p1y * wyLow + p2y * wyUpp
        dyy34 = p4y * wyLow + p3y * wyUpp

        # ---- Perform the bilinear interpolation
        fapprox  = fy12  * wxLow + fy34  * wxUpp
        dxApprox = dxy12 * wxLow + dxy34 * wxUpp
        dyApprox = dyy12 * wxLow + dyy34 * wxUpp
    end
    
    return fapprox::Float64, dxApprox::Float64, dyApprox::Float64

end # bilinearBoundedLookup
