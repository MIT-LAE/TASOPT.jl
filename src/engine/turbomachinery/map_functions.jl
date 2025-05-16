using Interpolations, NLsolve

"""
    MapDefaults

Holds the default reference conditions for a compressor map, used for scaling and normalization.

!!! details "💾 Fields"
    - `Nc::Float64`: normalized speed at design point.
    - `Rline::Float64`: R-line value at design point.
    - `Wc::Float64`: corrected mass flow at design point.
    - `PR::Float64`: pressure ratio at design point.
    - `polyeff::Float64`: maximum polytropic efficiency in the map.
"""
struct MapDefaults
    Nc::Float64
    Rline::Float64
    Wc::Float64
    PR::Float64 
    polyeff::Float64
end

"""
    CompressorMap

Contains all the map data and interpolations needed to compute compressor performance, including design defaults,
interpolated surfaces, and extrapolated lookup tables.

!!! details "💾 Fields"
    - `defaults::MapDefaults`: reference design conditions for scaling.
    - `RlineMap::Vector{Float64}`: grid of R-line values.
    - `NcMap::Vector{Float64}`: grid of normalized speeds.
    - `WcMap::Matrix{Float64}`: corrected mass flow values.
    - `PRMap::Matrix{Float64}`: pressure ratio values.
    - `effMap::Matrix{Float64}`: isentropic efficiency values.
    - `polyeffMap::Matrix{Float64}`: polytropic efficiency values.
    - `itp_Wc::GriddedInterpolation`: interpolant for `Wc(N, R)`.
    - `itp_PR::GriddedInterpolation`: interpolant for `PR(N, R)`.
    - `itp_polyeff::GriddedInterpolation`: interpolant for polytropic efficiency.
"""
struct CompressorMap
    defaults::MapDefaults
    RlineMap::Vector{Float64}
    NcMap::Vector{Float64}
    WcMap::Matrix{Float64}
    PRMap::Matrix{Float64}
    effMap::Matrix{Float64}
    polyeffMap::Matrix{Float64}
    itp_Wc::Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}
    itp_PR::Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}
    itp_polyeff::Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}
end

"""
    poly_efficiency_map_from_isentropic(effMap, PRMap)

Computes the polytropic efficiency map from an isentropic efficiency map and pressure ratio map.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `effMap::Matrix{Float64}`: isentropic efficiency at each point in the map.
    - `PRMap::Matrix{Float64}`: pressure ratio at each point in the map.

    **Outputs:**
    - `polyeff_Map::Matrix{Float64}`: computed polytropic efficiency at each point.

!!! note
    Uses γ = 1.4. NaNs are replaced with 0.0 in the output map.
"""
function poly_efficiency_map_from_isentropic(effMap, PRMap)
    γ = 1.4
    polyeff_Map = (γ - 1)/γ * log.(PRMap) ./ log.((PRMap.^((γ-1)/γ) .- 1.0) ./ effMap .+ 1)
    polyeff_Map[isnan.(polyeff_Map)] .= 0.0 #Replace NaN values with 0.0
    return polyeff_Map
end

"""
    find_NR_inverse_with_derivatives(itp_Wc, itp_PR, Wc_target, PR_target; Ng=0.5, Rg=2.0)

Finds the normalized speed `N` and R-line `R` corresponding to a target corrected mass flow `Wc_target` and pressure ratio `PR_target`,
using a nonlinear solver with Jacobian information from interpolation gradients.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `itp_Wc::GriddedInterpolation`: interpolation of corrected mass flow over (N, R).
    - `itp_PR::GriddedInterpolation`: interpolation of pressure ratio over (N, R).
    - `Wc_target::Float64`: target corrected mass flow.
    - `PR_target::Float64`: target pressure ratio.
    - `Ng::Float64`: initial guess for normalized speed (optional).
    - `Rg::Float64`: initial guess for R-line (optional).

    **Outputs:**
    - `N::Float64`: normalized speed at the matched point.
    - `R::Float64`: R-line value at the matched point.
    - `dN_dw::Float64`: ∂N/∂Wc
    - `dN_dpr::Float64`: ∂N/∂PR
    - `dR_dw::Float64`: ∂R/∂Wc
    - `dR_dpr::Float64`: ∂R/∂PR
"""
function find_NR_inverse_with_derivatives(itp_Wc::Interpolations.GriddedInterpolation, itp_PR::Interpolations.GriddedInterpolation, 
                                    Wc_target::Float64, PR_target::Float64; Ng::Float64 = 0.5, Rg::Float64 = 2.0)
    # Define the system of equations: 
    function residuals!(F::Vector{Float64}, x::Vector{Float64})
        x[1] = clamp(x[1], 1e-4, 1.9999) #Extrapolated speed map goes from 0 to 2.0
        x[2] = clamp(x[2], 1e-4, 3.9999) #Extrapolated Rline map goes from 0 to 4.0
        # Return the residuals for both equations
        F[1] = itp_Wc(x...) - Wc_target 
        F[2] = itp_PR(x...) - PR_target
    end

    # Define the Jacobian of the system (partial derivatives)
    function jacobian!(J::Matrix{Float64}, x::Vector{Float64})
        x[1] = clamp(x[1], 1e-4, 1.9999) #Extrapolated speed map goes from 0 to 2.0
        x[2] = clamp(x[2], 1e-4, 3.9999) #Extrapolated Rline map goes from 0 to 4.0
        # Compute the partial derivatives of W and PR with respect to N and R
        dw_dN, dw_dR = Interpolations.gradient(itp_Wc, x[1], x[2])
        dpr_dN, dpr_dR = Interpolations.gradient(itp_PR, x[1], x[2])
        
        # Return the Jacobian matrix
        J[1,1] = dw_dN
        J[1,2] = dw_dR
        J[2,1] = dpr_dN
        J[2,2] = dpr_dR
    end

    # Solve the system of equations using root finding (non-linear solver)
    sol = nlsolve(residuals!, jacobian!, [Ng, Rg], factor = 1.0, iterations = 100)

    if ~converged(sol) #Try again from a different initial guess if the first one fails
        sol = nlsolve(residuals!, jacobian!, [1.0, 2.0], method = :newton, iterations = 150)
    end

    # Extract the solution: the x and y corresponding to the given Wc_target and PR_target
    N_found, R_found = sol.zero

    # Compute the Jacobian at the found solution
    jac = zeros(2, 2)
    jacobian!(jac, [N_found, R_found])

    # Calculate the derivatives (inverse of the Jacobian matrix)
    jac_inv = inv(jac)

    # The derivatives of N, R with respect to Wc and PR are the components of the inverse Jacobian
    dN_dw, dN_dpr = jac_inv[1, :]
    dR_dw, dR_dpr = jac_inv[2, :]

    return N_found, R_found, dN_dw, dN_dpr, dR_dw, dR_dpr
end

"""
    calculate_compressor_speed_and_efficiency(map, pratio, mb, piD, mbD, NbD, epol0; Ng=0.5, Rg=2.0)

Calculates corrected speed and polytropic efficiency for a compressor, along with derivatives to pressure ratio and mass flow.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `map::CompressorMap`: structure containing the compressor map and interpolations.
    - `pratio::Float64`: compressor pressure ratio.
    - `mb::Float64`: current mass flow rate.
    - `piD::Float64`: design pressure ratio.
    - `mbD::Float64`: design mass flow rate.
    - `NbD::Float64`: design rotational speed.
    - `epol0::Float64`: maximum polytropic efficiency.
    - `Ng::Float64`: initial guess for normalized speed (optional).
    - `Rg::Float64`: initial guess for R-line (optional).

    **Outputs:**
    - `Nb::Float64`: corrected compressor speed.
    - `epol::Float64`: polytropic efficiency.
    - `dNb_dpi::Float64`: derivative of corrected speed w.r.t. `piD`.
    - `dNb_dmb::Float64`: derivative of corrected speed w.r.t. `mb`.
    - `depol_dpi::Float64`: derivative of efficiency w.r.t. `piD`.
    - `depol_dmb::Float64`: derivative of efficiency w.r.t. `mb`.
    - `N::Float64`: matched normalized speed.
    - `R::Float64`: matched R-line.
"""
function calculate_compressor_speed_and_efficiency(map::CompressorMap, pratio::Float64, mb::Float64, piD::Float64, 
                                                    mbD::Float64, NbD::Float64, epol0::Float64; Ng::Float64 = 0.5, Rg::Float64 = 2.0)
    #Cache conversions
    dWc_dmb = map.defaults.Wc/mbD
    dNb_dN = NbD/map.defaults.Nc
    depol_dep = epol0 / map.defaults.polyeff
    
    #Calculate objective Wc and PR
    Wc = mb * dWc_dmb
    PR = 1.0 + (pratio-1.0)/(piD-1.0) * (map.defaults.PR-1.0)
    dPR_dpi = (map.defaults.PR-1.0)/(piD-1.0)

    #Calculate speed and Rline for the map
    N, R, dN_dw, dN_dpr, dR_dw, dR_dpr = find_NR_inverse_with_derivatives(map.itp_Wc, map.itp_PR, Wc, PR, Ng = Ng, Rg = Rg)

    #Compute corrected speed and derivatives
    Nb = N * dNb_dN
    dNb_dmb = dN_dw * dWc_dmb * dNb_dN
    dNb_dpi = dN_dpr * dPR_dpi * dNb_dN

    #Compute efficiency and derivatives
    ep = map.itp_polyeff(N, R)
    dep_dN, dep_dR = Interpolations.gradient(map.itp_polyeff, N, R)

    #Rescale efficiency to design point
    epol = ep * depol_dep
    depol_dN = dep_dN * depol_dep
    depol_dR = dep_dR * depol_dep

    #Convert to adjusted values
    depol_dw = depol_dN * dN_dw + depol_dR * dR_dw
    depol_dpr = depol_dN * dN_dpr + depol_dR * dR_dpr

    depol_dmb = depol_dw * dWc_dmb
    depol_dpi = depol_dpr * dPR_dpi
    return Nb, epol, dNb_dpi, dNb_dmb, depol_dpi, depol_dmb, N, R
end

"""
    create_extrapolated_maps(NcMap, RlineMap, WcMap, PRMap, polyeff_Map)

Creates extended 2D interpolants for compressor map data, with extrapolation buffers to avoid interpolation edge errors.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `NcMap::Vector{Float64}`: normalized speed grid.
    - `RlineMap::Vector{Float64}`: R-line grid.
    - `WcMap::Matrix{Float64}`: corrected mass flow data.
    - `PRMap::Matrix{Float64}`: pressure ratio data.
    - `polyeff_Map::Matrix{Float64}`: polytropic efficiency data.

    **Outputs:**
    - `itp_Wc`: interpolated `Wc(N, R)` map.
    - `itp_PR`: interpolated `PR(N, R)` map.
    - `itp_polyeff`: interpolated polytropic efficiency map.
"""
function create_extrapolated_maps(NcMap, RlineMap, WcMap, PRMap, polyeff_Map)
    # Extend the Rline and Nc axes with padding values
    mRlineMap = [0; RlineMap; 4]
    mNcMap = [0; NcMap; 2]

    large_fac = 1.5
    small_fac = 1.01

    # Determine a max value for extrapolated Wc
    Wcmax = large_fac * maximum(WcMap)

    # Pad Wc map along Rline (left = zeros, right = slightly above max Rline)
    lowRWc = zeros(size(WcMap)[1])
    highRWc = small_fac * WcMap[:, end]
    mWcMap = [lowRWc WcMap highRWc]

    # Pad Wc map along Nc (top = zeros, bottom = Wcmax)
    lowNWc = zeros(1, size(mWcMap)[2])
    highNWc = Wcmax * ones(1, size(mWcMap)[2])
    mWcMap = vcat(lowNWc, mWcMap, highNWc)

    # Ensure bottom-left corner is 0 to avoid NaNs
    mWcMap[end, 1] = 0.0

    # Extend PR map similarly
    PRmax = large_fac * maximum(PRMap)
    lowRPR = small_fac * PRMap[:, 1]                   # Slightly extrapolate on the left
    highRPR = ones(size(PRMap)[1])                # Constant 1.0 padding on the right
    mPRMap = [lowRPR PRMap highRPR]

    lowNPR = 0.99*ones(1, size(mPRMap)[2])             # 1.0 padding on top
    highNPR = PRmax * ones(1, size(mPRMap)[2])    # max value padding on bottom
    mPRMap = vcat(lowNPR, mPRMap, highNPR)

    # Ensure bottom-right corner is 0.99 for robustness
    mPRMap[end, end] = 0.99

    # Pad efficiency map with zeros
    lowReff = zeros(size(polyeff_Map)[1])
    highReff = zeros(size(polyeff_Map)[1])
    mpolyeff_Map = [lowReff polyeff_Map highReff]

    lowNeff = zeros(1, size(mpolyeff_Map)[2])
    highNeff = zeros(1, size(mpolyeff_Map)[2])
    mpolyeff_Map = vcat(lowNeff, mpolyeff_Map, highNeff)

    # Create linear interpolants over the extended grids
    itp_Wc = interpolate((mNcMap, mRlineMap), mWcMap, Gridded(Linear()))
    itp_PR = interpolate((mNcMap, mRlineMap), mPRMap, Gridded(Linear()))
    itp_polyeff = interpolate((mNcMap, mRlineMap), mpolyeff_Map, Gridded(Linear()))

    return itp_Wc, itp_PR, itp_polyeff
end

# Function to compute interpolator gradients by finite differences; doesn't seem to 
# make much of a difference
# function interpolator_gradient(itp::Interpolations.GriddedInterpolation, x::Float64, y::Float64)
#     eps = 1e-6
#     # Compute the gradient using finite differences
#     df_dx = (itp(x + eps, y) - itp(x - eps, y)) / (2*eps)
#     df_dy = (itp(x, y + eps) - itp(x, y - eps)) / (2*eps)
#     return df_dx, df_dy
# end