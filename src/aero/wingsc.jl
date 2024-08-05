"""
    wingsc(W, CL, qinf, AR, ηsi, bo, λt, λs)

Sizes wing area, span, root chord from `q`, `CL`, `W`, `AR` at given point (taken as start-of-cruise in `wsize`).

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `W::Float64`: Aircraft weight.
    - `CL::Float64`: Lift coefficient.
    - `qinf::Float64`: Freestream dynamic head.
    - `wing::TASOPT.structures.Wing`: Wing structure 

See Sections 2.5 and 3.4.1 of the [TASOPT Technical Desc](@ref dreladocs).
"""
function wingsc!(W,CL,qinf,wing)
    wing.layout.S = W/(qinf*CL)
    wing.outboard.layout.b= sqrt(wing.layout.S*wing.layout.AR)

    wing.inboard.layout.b = max( wing.outboard.layout.b*wing.ηs , wing.layout.box_halfspan  )

    ηo = wing.layout.box_halfspan /wing.outboard.layout.b
    ηs = wing.inboard.layout.b/wing.outboard.layout.b

    Kc = ηo +
    0.5*(1.0    +wing.inboard.layout.λ)*(ηs-ηo) +
    0.5*(wing.inboard.layout.λ+ wing.outboard.layout.λ)*(1.0 -ηs)

    wing.layout.chord = wing.layout.S/(Kc*wing.outboard.layout.b)
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


