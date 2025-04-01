using Interpolations, NLsolve
using BenchmarkTools
using TASOPT

struct Defaults
    Rline::Float64
    PR::Float64 
    Nc::Float64
end

struct CompressorMap
    defaults::Defaults
    RlineMap::Vector{Float64}
    NcMap::Vector{Float64}
    WcMap::Matrix{Float64}
    PRMap::Matrix{Float64}
    effMap::Matrix{Float64}
end

RlineMap = [1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000]
NcMap = [0.300, 0.400, 0.500, 0.600, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 1.000, 1.050, 1.100, 1.150]

WcMap = [
    121.797 150.895 179.422 207.275 234.359 260.582 285.861 310.120 333.291 355.315 369.552;
    194.867 227.417 258.872 289.101 317.981 345.407 371.284 395.535 418.097 438.921 457.973;
    265.640 302.320 337.109 369.834 400.351 428.543 454.320 477.621 498.412 516.684 532.453;
    330.650 373.304 412.750 448.759 481.163 509.851 534.771 555.926 573.366 587.188 597.526;
    380.373 433.991 481.896 523.736 559.314 588.580 611.621 628.642 639.948 645.923 647.172;
    408.374 467.955 520.262 564.886 601.649 630.583 651.910 666.008 673.384 674.865 674.865;
    432.494 499.118 556.491 604.137 641.921 670.010 688.838 699.045 701.536 701.536 701.536;
    470.522 540.248 599.144 646.786 683.165 708.637 723.856 729.705 729.841 729.841 729.841;
    527.146 594.695 650.583 694.581 726.881 748.023 758.831 760.796 760.796 760.796 760.796;
    593.025 653.568 702.700 740.438 767.140 783.445 790.213 790.533 790.533 790.533 790.533;
    643.809 697.196 739.568 771.126 792.350 803.950 806.892 806.892 806.892 806.892 806.892;
    684.120 729.283 764.846 791.078 808.441 817.554 819.416 819.416 819.416 819.416 819.416;
    723.881 760.550 789.192 810.094 823.667 830.421 831.443 831.443 831.443 831.443 831.443;
    762.811 790.037 811.271 826.753 836.786 841.724 842.410 842.410 842.410 842.410 842.410
]

effMap = [
    0.6931 0.7672 0.8306 0.8802 0.9112 0.9167 0.8887 0.8138 0.6634 0.3791 0.0000;
    0.7418 0.8033 0.8552 0.8953 0.9199 0.9241 0.9030 0.8489 0.7461 0.5672 0.2591;
    0.7517 0.8112 0.8613 0.8997 0.9233 0.9277 0.9089 0.8607 0.7712 0.6201 0.3720;
    0.7322 0.7985 0.8543 0.8973 0.9240 0.9297 0.9103 0.8596 0.7650 0.6063 0.3482;
    0.6641 0.7519 0.8264 0.8844 0.9213 0.9308 0.9071 0.8412 0.7151 0.4988 0.1381;
    0.6357 0.7325 0.8148 0.8789 0.9199 0.9310 0.9059 0.8350 0.6992 0.4665 0.0800;
    0.6109 0.7155 0.8044 0.8738 0.9183 0.9307 0.9044 0.8294 0.6858 0.4403 0.0347;
    0.6261 0.7262 0.8108 0.8763 0.9179 0.9294 0.9052 0.8376 0.7105 0.4991 0.1625;
    0.6875 0.7679 0.8348 0.8857 0.9173 0.9253 0.9067 0.8573 0.7685 0.6279 0.4175;
    0.7557 0.8118 0.8576 0.8916 0.9119 0.9161 0.9030 0.8710 0.8162 0.7337 0.6169;
    0.7947 0.8341 0.8656 0.8885 0.9014 0.9030 0.8926 0.8695 0.8315 0.7761 0.6998;
    0.8153 0.8449 0.8684 0.8851 0.8940 0.8942 0.8853 0.8669 0.8376 0.7958 0.7392;
    0.8289 0.8503 0.8669 0.8783 0.8838 0.8829 0.8753 0.8608 0.8386 0.8076 0.7666;
    0.8385 0.8532 0.8643 0.8715 0.8745 0.8728 0.8664 0.8552 0.8387 0.8163 0.7873
]

PRmap = [
    1.0546 1.0558 1.0553 1.0532 1.0494 1.0440 1.0372 1.0290 1.0196 1.0089 1.0000;
    1.1002 1.1010 1.0994 1.0955 1.0892 1.0807 1.0702 1.0578 1.0436 1.0278 1.0102;
    1.1593 1.1606 1.1584 1.1527 1.1434 1.1307 1.1149 1.0964 1.0753 1.0517 1.0258;
    1.2313 1.2363 1.2353 1.2283 1.2154 1.1968 1.1730 1.1449 1.1125 1.0764 1.0366;
    1.3042 1.3232 1.3306 1.3261 1.3100 1.2825 1.2450 1.1992 1.1457 1.0854 1.0193;
    1.3479 1.3775 1.3917 1.3901 1.3728 1.3402 1.2944 1.2376 1.1708 1.0956 1.0133;
    1.3935 1.4356 1.4581 1.4602 1.4420 1.4040 1.3488 1.2797 1.1982 1.1065 1.0068;
    1.4654 1.5127 1.5381 1.5409 1.5209 1.4789 1.4178 1.3413 1.2513 1.1501 1.0406;
    1.5830 1.6221 1.6408 1.6386 1.6155 1.5723 1.5116 1.4367 1.3492 1.2512 1.1448;
    1.7258 1.7494 1.7565 1.7469 1.7208 1.6787 1.6229 1.5555 1.4778 1.3912 1.2969;
    1.8381 1.8472 1.8432 1.8260 1.7960 1.7537 1.7006 1.6386 1.5684 1.4910 1.4073;
    1.9316 1.9312 1.9197 1.8973 1.8642 1.8209 1.7687 1.7091 1.6427 1.5702 1.4923;
    2.0237 2.0149 1.9970 1.9704 1.9352 1.8918 1.8414 1.7850 1.7231 1.6562 1.5848;
    2.1043 2.0885 2.0659 2.0366 2.0008 1.9588 1.9114 1.8595 1.8033 1.7433 1.6797
]

mRlineMap = [0.0; RlineMap; 4.0]
mNcMap = [0.0; NcMap; 10.0]

Wcmax = 1e4
lowRWc = zeros(size(WcMap)[1])
highRWc = Wcmax * ones(size(WcMap)[1])
mWcMap = [lowRWc WcMap highRWc]

lowNWc = zeros(1, size(mWcMap)[2])
highNWc = Wcmax * ones(1, size(mWcMap)[2])
mWcMap = vcat(lowNWc, mWcMap, highNWc)

PRmax = 3.0
lowRPR = ones(size(PRmap)[1])
highRPR = PRmax * ones(size(PRmap)[1])
mPRMap = [lowRPR PRmap highRPR]

lowNPR = ones(1, size(mPRMap)[2])
highNPR = PRmax * ones(1, size(mPRMap)[2])
mPRMap = vcat(lowNPR, mPRMap, highNPR)

function find_xy_inverse_with_derivatives(itp_W, itp_Z, w_target, z_target, xg = 0.5, yg = 2.0)

    # Define the system of equations: Z(x, y) = z_target, W(x, y) = w_target
    function residuals(p)
        # Return the residuals for both equations
        return [itp_W(p...) - w_target, itp_Z(p...) - z_target]
    end

    # Define the Jacobian of the system (partial derivatives)
    function jacobian(p)
        # Compute the partial derivatives of W and Z with respect to x and y
        dw_dx, dw_dy = gradient(itp_W, p[1], p[2])
        dz_dx, dz_dy = gradient(itp_Z, p[1], p[2])
        
        # Return the Jacobian matrix
        return [dw_dx dw_dy; dz_dx dz_dy]
    end

    # Solve the system of equations using root finding (non-linear solver)
    sol = nlsolve(residuals, jacobian, [xg, yg])

    # Extract the solution: the x and y corresponding to the given w_target and z_target
    x_found, y_found = sol.zero

    # Compute the Jacobian at the found solution
    jac = jacobian([x_found, y_found])

    # Calculate the derivatives (inverse of the Jacobian matrix)
    jac_inv = inv(jac)

    # The derivatives of x, y with respect to w and z are the components of the inverse Jacobian
    dx_dw, dx_dz = jac_inv[1, :]
    dy_dw, dy_dz = jac_inv[2, :]

    return x_found, y_found, dx_dw, dx_dz, dy_dw, dy_dz
end

Wcobj = 1000
PRobj = 1.2
#eps = 1e-4

# Create interpolation objects for W and Z
itp_W = interpolate((mNcMap, mRlineMap), mWcMap, Gridded(Linear()))
itp_PR = interpolate((mNcMap, mRlineMap), mPRmap, Gridded(Linear()))

(x, y, dx_dw, dx_dz, dy_dw, dy_dz) = find_xy_inverse_with_derivatives(itp_W, itp_PR, Wcobj, PRobj)
#(xw, yw, _, _, _, _) = find_xy_inverse_with_derivatives(itp_W, itp_PR, Wcobj + eps, PRobj)
#(xp, yp, _, _, _, _) = find_xy_inverse_with_derivatives(itp_W, itp_PR, Wcobj, PRobj + eps)

#@benchmark find_xy_inverse_with_derivatives(itp_W, itp_PR, Wcobj, PRobj)
#@benchmark TASOPT.engine.Ncmap(PRobj, Wcobj, 1.8, 900, 1, TASOPT.engine.Cmapf)