"""
      surfdx(b, bs, bo, λt, λs, sweep)

Calculates area centroid x-offset due to sweep
and the mean aerodynamic chord (normalized by root chord, `co`)

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `b::Float64`: Wingspan
    - `bs::Float64`: Spanwise location of the start of the taper
    - `bo::Float64`: Spanwise location of the root chord
    - `λt::Float64`: Tip chord ratio (tip chord / root chord)
    - `λs::Float64`: Start chord ratio (start chord / root chord).
    - `sweep::Float64`: Sweep angle in degrees.

    **Outputs:**
    - `dx::Float64`: Area centroid x-offset due to sweep
    - `macco::Float64`: Mean aerodynamic chord normalized by root chord.

See [Geometry](@ref geometry) or Section 2.5.1  of the [TASOPT Technical Description](@ref dreladocs).
"""

function surfdx(b,bs,bo,λt,λs,sweep)

      tanL = tan(deg2rad(sweep))

      ηo = bo/b
      ηs = bs/b

#---- 2 Int c dy /(co b)  =  S/(co b)  =  Kc
      Kc = ηo +
	 0.5*(1.0    +λs)*(ηs-ηo) +
	 0.5*(λs+λt)*(1.0 -ηs)

#---- 2 Int c (x-xo) dy /(co b^2) 
      Kcx = (1.0     + 2.0*λs)*(ηs-ηo)^2 / 12.0 +
	 (λs + 2.0*λt)*(1.0 -ηs)^2 / 12.0 +
	 (λs +     λt)*(1.0 -ηs)*(ηs-ηo) / 4.0

#---- 2 Int c^2 dy / (co^2 b)
      Kcc = ηo +
	 (1.0        + λs         + λs^2)*(ηs-ηo)/3.0 +
	 (λs^2 + λs*λt + λt^2)*(1.0 -ηs)/3.0

      dx    = Kcx/Kc * b*tanL
      macco = Kcc/Kc

      return dx,macco

end # surfdx

"""
      surfdx!(wing, b, bs, parg)

surfdx wrapper for Wing

"""
function surfdx!(wing::Wing; b::Float64 = 0.0, bs::Float64 = 0.0, parg::Vector{Float64} = Float64[])
      if isempty(parg)
          dx, _ = surfdx(b, bs,
                      wing.outboard.layout.b,
                      wing.outboard.layout.λ,
                      wing.inboard.layout.λ,
                      wing.layout.sweep)
      else
          dx, macco = surfdx(wing.layout.b,
                          wing.inboard.layout.b,
                          wing.outboard.layout.b,
                          wing.outboard.layout.λ,
                          wing.inboard.layout.λ,
                          wing.layout.sweep)
          parg[igcma] = macco * wing.layout.chord
      end
      wing.layout.x = wing.layout.box_x + dx
end

"""
      surfdx!(wing, b, λs)

surfdx wrapper for Tail

"""
function surfdx!(tail::Tail, b::Float64, λs::Float64)
      dx, _ = surfdx(b,
                  tail.outboard.layout.b,
                  tail.outboard.layout.b,
                  tail.outboard.layout.λ,
                  λs,
                  tail.layout.sweep)
      tail.layout.x = tail.layout.box_x + dx         
end