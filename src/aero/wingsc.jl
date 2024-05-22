"""
    wingsc(W, CL, qinf, AR, ηsi, bo, λt, λs)

Sizes wing area, span, root chord from `q`, `CL`, `W`, `AR` at given point (taken as start-of-cruise in `wsize`).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `W::Float64`: Aircraft weight.
    - `CL::Float64`: Lift coefficient.
    - `qinf::Float64`: Freestream dynamic head.
    - `AR::Float64`: Wing aspect ratio.
    - `ηsi::Float64`: Span fraction of inner wing break ("snag").
    - `bo::Float64`: Wing center box width.
    - `λt::Float64`: Outer or "tip" taper ratio of chord.
    - `λs::Float64`: Inner or break/"snag" taper ratio of chord.

    **Outputs:**
    - `S::Float64`: Wing planform area (including fuselage carryover).
    - `b::Float64`: Wing span.
    - `bs::Float64`: Span of inner wing (break/"snag").
    - `co::Float64`: Chord at wing root, "center."

See Sections 2.5 and 3.4.1 of the [TASOPT Technical Desc](@ref dreladocs).
"""
function wingsc!(W,CL,qinf,wing)
    AR = wing.layout.AR
    ηsi = wing.layout.ηs
    bo = wing.layout.box_halfspan 
    λt = wing.layout.λt
    λs = wing.layout.λs
    wing.layout.S = W/(qinf*CL)
    wing.layout.b = sqrt(wing.layout.S*AR)

    wing.layout.b_inner = max( wing.layout.b*ηsi , bo )

    ηo = bo/wing.layout.b
    ηs = wing.layout.b_inner/wing.layout.b

    Kc = ηo +
    0.5*(1.0    +λs)*(ηs-ηo) +
    0.5*(λs+λt)*(1.0 -ηs)

    wing.layout.chord = wing.layout.S/(Kc*wing.layout.b)
end # wingsc


#"""
#  Sets wing area, AR, root chord 
#  to be consistent with q,CL,weight,span
#"""
#      function wingAc(W,CL,qinf,b,ηsi,bo,bs,λt,λs)
#
#      S = W/(qinf*CL)
#
#      ηo = bo/b
#      ηs = bs/b
#
#      Kc = ηo +
#	 0.5*(1.0    +λs)*(ηs-ηo) +
#	 0.5*(λs+λt)*(1.0 -ηs)
#
#      co = S/(Kc*b)
#
#      AR = b^2 / S
#
#      return  S,AR,co
#      end # wingAc
#


