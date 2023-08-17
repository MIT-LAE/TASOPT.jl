using LinearAlgebra

"""
Calculates spline coefficients for X(S).          
Natural end conditions are used (zero 3rd        
derivative over first, last intervals).         
                                                  
To evaluate the spline at some value of S,       
use SEVAL and/or DEVAL.                          
- Inputs

  S        independent variable array (input)      
  X        dependent variable array   (input) 
  
- Outputs

  XS       dX/dS array                (calculated) 
    
                                                      
"""
function spline(S, X)
 
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

        B[i] = ΔS₊
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

# ---- solve for derivative array XS
    SYS = Tridiagonal(B, A, C)
    XS = SYS \ D    

    return XS
end # SPLINE


function SEVAL(SS, X, XS, S)
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
    i_low = 1
    i = length(S)

    # Use bisection to find right segment where SS is.
    # i stores the closest index S[i] > SS
    while (i - i_low > 1)
        i_mid = (i + i_low) ÷ 2
        if (SS < S[i_mid])
            i = i_mid
        else
            i_low = i_mid
        end
    end

    ΔS = S[i] - S[i - 1]
    ΔX = X[i] - X[i - 1]
    T  = (SS - S[i - 1]) / ΔS

    CX1 = ΔS*XS[i - 1] - ΔX
    CX2 = ΔS*XS[i    ] - ΔX

    XX = T * X[i] + (1.0 - T) * X[i - 1] + (T - T * T) * ((1.0 - T) * CX1 - T * CX2)

    return XX 
end # SEVAL
