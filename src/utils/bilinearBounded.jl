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
    minX, minXind = findmin(xi->abs(xi-x), a)

end # bilinearBoundedLookup