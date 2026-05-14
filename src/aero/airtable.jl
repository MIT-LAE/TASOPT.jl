using DataFrames, CSV
"""
    airtable(fname)

Reads an airfoil database CSV and returns a populated [`airfoil`](@ref) struct, including
cubic-spline derivative tables at the grid knots for tricubic interpolation.

# Required CSV schema
The file must contain a rectilinear `(Mach, CL, toc)` grid (i.e., like `C_airfoil.csv`) with columns:

| column | description |
|--------|-------------|
| `toc`            | thickness/chord ratio |
| `Re`             | Reynolds number (must be the same value on every row) |
| `Mach`           | freestream Mach |
| `alpha`          | section angle of attack, **degrees** |
| `CL`             | section lift coefficient |
| `CD`             | total section drag |
| `CM`             | section pitching moment |
| `CDw`, `CDv`, `CDp` | wave / viscous / pressure drag components |
| `file`           | provenance tag (source MSES `polar`/`blade` file) |
| `Missing`        | `true`/`false` flag for infeasible rows (the row's drag/moment columns should already be `NaN`) |
| `wasInterpolated`, `InterpAxis` | interpolation provenance (informational) |

The 6-tuple of interpolated quantities `[CD, CDp, CDv, CDw, CM, alpha]` is the **fixed**
4th-axis order of `airfoil.A` and is consumed positionally by [`airfun`](@ref).

# Assumptions
- **Single Reynolds number** — errors out if multiple `Re` values are present. Re-scaling is
  applied separately (see [`wing_profiledrag_direct`](@ref) and [`wing_profiledrag_scaled`](@ref)).
- **Rectilinear grid** in `(Mach, CL, toc)`; missing combinations emit a warning.
- **Degree → radian** conversion of `alpha` happens here so the rest of the codebase
  works in radians.
- **`NaN` propagation** — `Missing == true` rows should carry `NaN` in the drag/moment
  columns; these are passed through to `airfoil.A` and handled downstream by [`spline`](@ref)
  and [`airfun`](@ref).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `fname::String`: Path to a CSV file matching the schema above.

    **Outputs:**
    - `airf::TASOPT.aerodynamics.airfoil`: `struct` with airfoil performance characteristics
      and derivative arrays for interpolation in `airfun`.

"""
function airtable(fname)

# Read airfoil data
df = DataFrame(CSV.File(fname))

    cols_to_interp = ["CD", "CDp", "CDv", "CDw", "CM", "alpha"]

    AMa = sort(unique(df.Mach))     # array of Mach numbers in database
    Acl = sort(unique(df.CL))       # array of cl           in database
    Atoc = sort(unique(df.toc))       # array of t/c          in database
    
    nAMa::Int = length(AMa)
    nAcl::Int = length(Acl)
    nAtoc::Int = length(Atoc)
    nAfun::Int = length(cols_to_interp)

    A   = zeros(nAMa, nAcl, nAtoc, nAfun)

    # Read in Reynolds number
    ARe = unique(df.Re) 
    if length(ARe) != 1
        @error("Only single-Reynolds Number airfoil CSV databases are supported at this time. \nTry to pare it down: $(fname)")
    else
        ARe = Float64(ARe[1])
    end

    #Read in airfoil data from df. not the most efficient, but doesn't need to be
    @inbounds for k = 1:nAtoc
        @inbounds for j = 1:nAcl
            @inbounds for i = 1:nAMa
                # Find the row matching current Mach, CL, and t/c
                mask = (df.Mach .== AMa[i]) .& (df.CL .== Acl[j]) .& (df.toc .== Atoc[k])
                row_idx = findfirst(mask)
                
                if !isnothing(row_idx)
                    # Extract values for each column in cols_to_interp
                    for (idx, col) in enumerate(cols_to_interp)
                        A[i, j, k, idx] = df[row_idx, col]
                    end
                else
                    @warn "No data found for Mach=$(AMa[i]), CL=$(Acl[j]), t/c=$(Atoc[k])"
                end
            end
        end
    end

    #convert alpha / aoa from degrees to radians for internal consistency
    A[:,:,:,6] = deg2rad.(A[:,:,:,6])
    
#assemble matrices of derivatives
# ---------------------
     A_M = zeros(nAMa, nAcl, nAtoc, nAfun)
     A_toc = zeros(nAMa, nAcl, nAtoc, nAfun)
    A_cl = zeros(nAMa, nAcl, nAtoc, nAfun)
   A_M_toc = zeros(nAMa, nAcl, nAtoc, nAfun)
  A_M_cl = zeros(nAMa, nAcl, nAtoc, nAfun)
  A_cl_toc = zeros(nAMa, nAcl, nAtoc, nAfun)
A_M_cl_toc = zeros(nAMa, nAcl, nAtoc, nAfun)

for l = 1:nAfun
    for k = 1:nAtoc
        for j = 1:nAcl
            A_M[:, j, k, l] = spline(AMa, A[:, j, k, l]) # This gets the derivatives ∂cdf/∂M, ∂cdp/∂M, ∂cm/∂M
        end
    end

    for k = 1:nAtoc
        for i = 1:nAMa
            # for j = 1:nAcl
                f   =   A[i, :, k, l]
                f_M = A_M[i, :, k, l]
                  A_cl[i, :, k, l] = spline(Acl, f)
                A_M_cl[i, :, k, l] = spline(Acl, f_M)
        end
    end

    for j = 1:nAcl
        for i = 1:nAMa
            f   = A[i, j, :, l]
            f_M = A_M[i, j, :, l]
              A_toc[i, j, :, l] = spline(Atoc, f)
            A_M_toc[i, j, :, l] = spline(Atoc, f_M)
        end
    end

    for i = 1:nAMa
        for j = 1:nAcl
            f   = A_cl[i, j, :, l]
            f_M = A_M_cl[i, j, :, l]
              A_cl_toc[i, j, :, l] = spline(Atoc, f)
            A_M_cl_toc[i, j, :, l] = spline(Atoc, f_M)
        end
        for k = 1:nAtoc
            f   = A_toc[i, :, k, l]
            f_M = A_M_toc[i, :, k, l]
              A_cl_toc[i, :, k, l] = 0.5*(  A_cl_toc[i, :, k, l] + spline(Acl, f))
            A_M_cl_toc[i, :, k, l] = 0.5*(A_M_cl_toc[i, :, k, l] + spline(Acl, f_M))
        end
    end


end

return  airfoil(AMa, Acl, Atoc, ARe,
        A,
        A_M,
        A_toc,
        A_cl,
        A_M_toc,
        A_M_cl,
        A_cl_toc,
        A_M_cl_toc)
end

