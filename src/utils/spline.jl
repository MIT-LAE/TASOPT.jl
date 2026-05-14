using LinearAlgebra

"""
    spline(S, X)

Calculate cubic-spline derivatives `dX/dS` at the knots of `S` using natural end
conditions (zero third derivative over the first and last intervals).

To evaluate the spline at some value of `S`, use `SEVAL` and/or `DEVAL`.

# NaN handling — assumption
The airfoil database (see [`airtable`](@ref)) marks infeasible operating points with
`NaN`. This function tolerates `NaN` entries in `X` **only if the valid (non-NaN) entries
form a single contiguous block** at one or both ends of the array — i.e. no internal
holes.

- If `X` is all `NaN`, returns a `NaN`-filled `XS` of the same length.
- If `X` has a single contiguous valid block, the spline is computed on that block and
  `XS` is `NaN` at the invalid positions.
- If a `NaN` is *surrounded* by valid data (an internal hole), this function **errors**.
  Callers are expected to ensure the database / input was constructed with this in mind.

- Inputs

  S        independent variable array (input)
  X        dependent variable array   (input)

- Outputs

  XS       dX/dS array                (calculated)
"""
function spline(S, X)
# Pre-processing for NaNs (i.e., missing or infeasible regions)
N = length(S)
nan_mask = isnan.(X) # Find NaN locations

# If any NaNs in X
if any(nan_mask)
    # Find indices for valid data
    valid_indices = findall(.!nan_mask)

    # if all NaN, return all NaN derivatives
    if isempty(valid_indices)
        return fill(NaN, N)
    end
    
    first_valid = first(valid_indices)
    last_valid = last(valid_indices)
    
    # Check for internal holes (NaN between valid data)
    if last_valid - first_valid + 1 != length(valid_indices)
        error("NaN holes detected in airfoil: cannot compute spline with internal missing values. " *
                "Valid indices: $(first_valid):$(last_valid), but only $(length(valid_indices)) valid points found.")
    end
    
    # Extract valid segment
    S_valid = S[first_valid:last_valid]
    X_valid = X[first_valid:last_valid]
    
    # Recursively call spline on valid segment
    XS_valid = spline(S_valid, X_valid)
    
    # Initialize output with NaNs
    XS = fill(NaN, N)
    
    # Insert computed derivatives into valid positions
    XS[first_valid:last_valid] = XS_valid
    
    return XS

else #no NaNs
# Create the N x N system to solve for spline parameters --> dX/dS
    N = length(S)
    B = zeros(N-1) # sub-diagonal
    A = zeros(N)   # diagonal
    C = zeros(N-1) # sup-diagonal
    D = zeros(N)   # RHS

    for i = 2:N - 1
        ΔS₋ = S[i    ] - S[i - 1]
        ΔS₊ = S[i + 1] - S[i    ]
        ΔX₊ = X[i + 1] - X[i    ]
        ΔX₋ = X[i    ] - X[i - 1]

        B[i-1] = ΔS₊
        A[i] = 2.0 * (ΔS₋ + ΔS₊)
        C[i] = ΔS₋
        D[i] = 3.0 * (ΔX₊ * ΔS₋/ΔS₊ + ΔX₋ * ΔS₊/ΔS₋)
    end

# ---- set zero 3rd derivative end conditions
    A[1] = 1.0
    C[1] = 1.0
    D[1] = 2.0 * (X[2] - X[1]) / (S[2] - S[1])

    B[N-1] = 1.0
    A[N] = 1.0
    D[N] = 2.0 * (X[N] - X[N - 1]) / (S[N] - S[N - 1])

    if (N == 2)
# ----- if only two points are present, specify zero 2nd derivative instead
# -     (straight line interpolation will result)
        B[N-1] = 1.0
        A[N] = 2.0
        D[N] = 3.0 * (X[N] - X[N - 1]) / (S[N] - S[N - 1])
    end

    # @show A
    # @show B
    # @show C
# ---- solve for derivative array XS
    SYS = Tridiagonal(B, A, C)
    XS = SYS \ D    
    # @show XS
    return XS
end #if-else
end # SPLINE


function SEVAL(SS::H, X::AbstractVector{H}, 
               XS::AbstractVector{H}, S::AbstractVector{H}) where H
# --------------------------------------------------
#     Calculates X(SS)                             |
#     XS array must have been calculated by SPLINE |
#                                                  |
#     X is dependent variable, S is independent    |
#     - Inputs                                     |
#        XS = dX/dS                                |
#        SS independent point to interpolate at    |
#        X is dependent array                      |
#        S is independent array                    |
#                                                  |
#     - Outputs                                    |
#        XX is interpolated value                  |
# --------------------------------------------------
    if SS ≤ S[1]
        im = 1
        i  = 2
    elseif SS ≥ S[end]
        i = length(S)
        im = i - 1
    else
        i::Int = searchsortedlast(S, SS) + 1
        im::Int = i - 1
    end

    ΔS = S[i] - S[im]
    ΔX = X[i] - X[im]
    T  = (SS - S[im]) / ΔS

    CX1 = ΔS*XS[im] - ΔX
    CX2 = ΔS*XS[i] - ΔX

    XX = T * X[i] + (1.0 - T) * X[im] + (T - T * T) * ((1.0 - T) * CX1 - T * CX2)

    return XX 
end # SEVAL
