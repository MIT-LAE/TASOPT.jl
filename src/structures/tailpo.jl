"""
    tailpo!(tail,S,qne; t_fac = 1.0)

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

See [Geometry](@ref geometry) or Section 2.3.2 and 2.9.6 of the [TASOPT Technical Description](@ref dreladocs).
"""
function tailpo!(tail, S, qne; t_fac = 1.0)

    b  = sqrt(S*tail.layout.AR*t_fac)
    tail.layout.root_chord = S/(0.5*b*(1.0+tail.outboard.λ))
    po = qne*S*tail.CL_max/b * 2.0/(1.0 + tail.outboard.λ)
    tail.outboard.co = tail.layout.root_chord*tail.inboard.λ
    tail.inboard.co = tail.layout.root_chord
    return po,b
end