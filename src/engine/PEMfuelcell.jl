function dy_dx(y, params, x)
    u, p = params
    return [dλ_dz_membrane(y[1], p[1], p[2])]
end