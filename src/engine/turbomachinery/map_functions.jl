using Interpolations, NLsolve
#using BenchmarkTools

struct MapDefaults
    Nc::Float64
    Rline::Float64
    Wc::Float64
    PR::Float64 
end

struct CompressorMap
    defaults::MapDefaults
    RlineMap::Vector{Float64}
    NcMap::Vector{Float64}
    WcMap::Matrix{Float64}
    PRMap::Matrix{Float64}
    effMap::Matrix{Float64}
    polyeffMap::Matrix{Float64}
    itp_Wc::Interpolations.GriddedInterpolation
    itp_PR::Interpolations.GriddedInterpolation
    itp_polyeff::Interpolations.GriddedInterpolation
end

function poly_efficiency_map_from_isentropic(effMap, PRMap)
    γ = 1.4
    polyeff_Map = (γ - 1)/γ * log.(PRMap) ./ log.((PRMap.^((γ-1)/γ) .- 1.0) ./ effMap .+ 1)
    polyeff_Map[isnan.(polyeff_Map)] .=0.0 #Replace NaN values with 0.0
    return polyeff_Map
end

function find_NR_inverse_with_derivatives(itp_Wc::Interpolations.GriddedInterpolation, itp_PR::Interpolations.GriddedInterpolation, 
                                    Wc_target::Float64, PR_target::Float64; Ng::Float64 = 0.5, Rg::Float64 = 2.0)

    # Define the system of equations: 
    function residuals!(F::Vector{Float64}, x::Vector{Float64})
        x[1] = clamp(x[1], 0.1, 3.9)
        x[2] = clamp(x[2], 0.1, 3.9)
        # Return the residuals for both equations
        F[1] = itp_Wc(x...) - Wc_target 
        F[2] = itp_PR(x...) - PR_target
    end

    # Define the Jacobian of the system (partial derivatives)
    function jacobian!(J::Matrix{Float64}, x::Vector{Float64})
        x[1] = clamp(x[1], 0.1, 3.9)
        x[2] = clamp(x[2], 0.1, 3.9)
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
    sol = nlsolve(residuals!, jacobian!, [Ng, Rg], iterations = 100)

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

function calculate_compressor_speed_and_efficiency(map::CompressorMap, pratio::Float64, mb::Float64, piD::Float64, 
                                                    mbD::Float64, NbD::Float64; Ng::Float64 = 0.5, Rg::Float64 = 2.0)
    Wc = mb/mbD * map.defaults.Wc
    PR = pratio/piD * map.defaults.PR

    #Calculate speed and Rline for the map
    N, R, dN_dw, dN_dpr, dR_dw, dR_dpr = find_NR_inverse_with_derivatives(map.itp_Wc, map.itp_PR, Wc, PR, Ng = Ng, Rg = Rg)

    #Compute corrected speed and derivatives
    Nb = N * NbD/map.defaults.Nc
    dNb_dmb = dN_dw * map.defaults.Wc/mbD * NbD/map.defaults.Nc
    dNb_dpi = dN_dpr * map.defaults.PR/piD * NbD/map.defaults.Nc

    #Compute efficiency and derivatives
    epol = map.itp_polyeff(N, R)
    depol_dN, depol_dR = Interpolations.gradient(map.itp_polyeff, N, R)

    #Convert to adjusted values
    depol_dw = depol_dN * dN_dw + depol_dR * dR_dw
    depol_dpr = depol_dN * dN_dpr + depol_dR * dR_dpr

    depol_dmb = depol_dw * map.defaults.Wc/mbD
    depol_dpi = depol_dpr * map.defaults.PR/piD
    return Nb, epol, dNb_dpi, dNb_dmb, depol_dpi, depol_dmb, N, R
end

function create_extrapolated_maps(NcMap, RlineMap, WcMap, PRMap, polyeff_Map)
    mRlineMap = [0; RlineMap; 4]
    mNcMap = [0; NcMap; 4]
    Wcmax = 2*maximum(WcMap)
    lowRWc = zeros(size(WcMap)[1])
    highRWc = 1.01*WcMap[:,end]
    mWcMap = [lowRWc WcMap highRWc]

    lowNWc = zeros(1, size(mWcMap)[2])
    highNWc = Wcmax * ones(1, size(mWcMap)[2])
    mWcMap = vcat(lowNWc, mWcMap, highNWc)
    mWcMap[end,1] = 0.0

    PRmax = 2*maximum(PRMap)
    lowRPR = 1.01*PRMap[:,1]
    highRPR = ones(size(PRMap)[1])
    mPRMap = [lowRPR PRMap highRPR]

    lowNPR = ones(1, size(mPRMap)[2])
    highNPR = PRmax * ones(1, size(mPRMap)[2])
    mPRMap = vcat(lowNPR, mPRMap, highNPR)
    mPRMap[end,end] = 1.0

    lowReff = zeros(size(polyeff_Map)[1])
    highReff = zeros(size(polyeff_Map)[1])
    mpolyeff_Map = [lowReff polyeff_Map highReff]

    lowNeff = zeros(1, size(mpolyeff_Map)[2])
    highNeff = zeros(1, size(mpolyeff_Map)[2])
    mpolyeff_Map = vcat(lowNeff, mpolyeff_Map, highNeff)

    itp_Wc = interpolate((mNcMap, mRlineMap), mWcMap, Gridded(Linear()))
    itp_PR = interpolate((mNcMap, mRlineMap), mPRMap, Gridded(Linear()))
    itp_polyeff = interpolate((mNcMap, mRlineMap), mpolyeff_Map, Gridded(Linear()))

    return itp_Wc, itp_PR, itp_polyeff
end