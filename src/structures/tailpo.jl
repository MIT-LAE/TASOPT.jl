"""
    tailpo(S, AR, λa, qne, CLmax)

Calculate stabilizer span, root chord and root loading based on the 
never-exceed dynamic pressure, maximum CL, sweep, and aspect ratio.

"""
function tailpo(S, AR, λa, qne, CLmax)

    b  = sqrt(S*AR)
    co = S/(0.5*b*(1.0+λa))
    po = qne*S*CLmax/b * 2.0/(1.0 + λa)
    return b, co, po
end