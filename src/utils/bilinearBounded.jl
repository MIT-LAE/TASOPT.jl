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
        return (indUpper-1, indUpper)
    end
end # boundingIndices


"""
        bilinearBoundedLookup(x, y, xGrid, yGrid, lookupGrid)

Implementation of bilinear interpolation on a bounded grid. Used by `NcTblMap` and `ecTblMap`
to calculate off-design engine performance information. Custom implementation is required due
need for warning user operating at out-of-bounds conditions

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `x`:          `x` position of query
    - `y`:          `y` position of query
    - `xGrid`:      1D-array of on-grid `x` values
    - `yGrid`:      1D-array of on-grid `y` values
    - `lookupGrid`: 2D-array of function values

    **Outputs:**
    - `val`: Bilinearly interpolated result (-1.0 if out-of-bounds)
"""
function bilinearBoundedLookup(x, y, xGrid, yGrid, lookupGrid)
    (xIndLow, xIndUpp) = boundingIndices(x, xGrid)
    (yIndLow, yIndUpp) = boundingIndices(y, yGrid)

    if (xIndLow == -1 || xIndUpp == -1) || (yIndLow == -1 || yIndUpp == -1)
        # raise error for out of map operation
        return -1.0
    end

end # bilinearBoundedLookup