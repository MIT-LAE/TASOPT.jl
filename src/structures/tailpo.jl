"""
    tailpo(S, AR, λa, qne, CLmax)

Calculates stabilizer span, root chord, and root loading based on the 
never-exceed dynamic pressure, maximum CL, sweep, and aspect ratio.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `S::Float64`: Stabilizer area.
    - `AR::Float64`: Stabilizer aspect ratio.
    - `λa::Float64`: Stabilizer taper ratio (tip chord / root chord).
    - `qne::Float64`: Never-exceed dynamic pressure.
    - `CLmax::Float64`: Maximum coefficient of lift.
    
    **Outputs:**
    - `b::Float64`: Stabilizer wingspan.
    - `co::Float64`: Stabilizer root chord length.
    - `po::Float64`: Stabilizer root loading.

See [here](@ref wingtail) or Section 2.3.2 and 2.9.6 of TASOPT docs.
"""
function tailpo(S, AR, λa, qne, CLmax)

    b  = sqrt(S*AR)
    co = S/(0.5*b*(1.0+λa))
    po = qne*S*CLmax/b * 2.0/(1.0 + λa)
    return b, co, po
end