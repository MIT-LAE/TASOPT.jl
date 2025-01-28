"""

"""
function check_simple_bound(val, bounds::SVector{2, Float64})
    if val < bounds[1] || val > bounds[2]
        return false
    else
        return true
    end
end


"""
ind: 1   2   3   4   5   6   7   8   9
x = [pf, pl, ph, mf, ml, mh, Tb, Pc, Mi]
"""
function is_in_TR(x, fan, lpc, hpc, Tbb, Pcb, Mib)
    if (!check_simple_bound(x[7], Tbb) || !check_simple_bound(x[8], Pcb)) || !check_simple_bound(x[9], Mib)
        # print("outside simple bound   ")
        return false
    end

    if !check_in_bounds_reggrid(x[4], x[1], fan)
        # print("outside fan map   ")
        return false
    end
    if !check_in_bounds_reggrid(x[5], x[2], lpc)
        # print("outside lpc map   ")
        return false
    end
    if !check_in_bounds_reggrid(x[6], x[3], hpc)
        # print("outside hpc map   ")
        return false
    end

    return true
end


"""
"""
function relax_trust_region(xNew, rlx, xOld, dxOld, fan, lpc, hpc; 
                            rlx_fac=1.414, trialMax=50)
    trials = 0

    startsOutside = !is_in_TR(xOld, fan, lpc, hpc, Tb_bounds, Pc_bounds, Mi_bounds)
    if startsOutside
        return rlx, xNew
    end

    while !is_in_TR(xNew, fan, lpc, hpc, Tb_bounds, Pc_bounds, Mi_bounds)
        trials = trials + 1
        if trials > trialMax
            break
        end

        rlx = rlx / rlx_fac
        xNew = xOld - rlx * dxOld
    end

    return rlx, xNew
end # relax_trust_region


"""
"""
function sol_norm(x; ndims=9)
    res = 0.0
    for i in 1:ndims
        res = res + x[i]*x[i]
    end 

    return res
end


"""
"""
function conv_criteria(x, dx; ndims=9)
    wn = rel_tol .* x .+ abs_tol

    conv = 0.0
    for i in 1:ndims
        conv = conv + (dx[i] * dx[i]) / (wn[i] * wn[i])
    end

    return conv
end