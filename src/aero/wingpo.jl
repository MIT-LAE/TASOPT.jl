"""
    wingpo(b, bs, bo, 
        λt, λs, γt, γs, 
        AR, N, W, Lhtail, fLo, fLt)

Computes wing root ("center") loading ``p_o`` to balance the net load.

```math
N*W - L_{h tail} \\times 2*∫p(η) dy + 2ΔL₀ + 2ΔLₜ = N*W - (L_{htail}).
```

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `b::Float64`: Wing span.
    - `bs::Float64`: Span of inner wing section.
    - `bo::Float64`:  Span of wing box (span at wing root).
    - `λt::Float64`, `λs::Float64` : Wing chord taper ratios at tip and break ("snag"), respectively.
    - `γt::Float64`,`γs::Float64` : Wing lift distribution "taper" ratios for outer and inner wing sections, respectively.
    - `AR::Float64`, `N::Float64`, `W::Float64`, `Lhtail` : Aspect ratio, Load factor, weight and H-tail lift.
    - `fLo::Float64`, `fLt::Float64` : Wing root and tip load adjustment factors.

    **Outputs:**
    - `po::Float64`: Wing root loading magnitude.

See Section 2.6.2 of the [TASOPT Technical Desc](@ref dreladocs).
"""
# function wingpo(b, bs, bo,
#                λt, λs, γt, γs,
#                AR, N, W, Lhtail, fLo, fLt)

# (wing.outboard.layout.b, wing.inboard.layout.b, wing.layout.box_halfspan,
#                     wing.outboard.layout.λ, wing.inboard.layout.λ, γt, γs,
#                     wing.layout.AR, Nlift, BW, Lhtail, fLo, fLt)
function wingpo(wing, rclt, rcls, N, W, Lhtail, fLo, fLt)
    
    γt, γs = wing.outboard.layout.λ * rclt, wing.inboard.layout.λ * rcls
        
    ηo = wing.layout.box_halfspan/wing.outboard.layout.b #calculate non-dim. span locations eta
    ηs = wing.inboard.layout.b/wing.outboard.layout.b
    
    Kc = ηo +
    0.5*(1.0    +wing.inboard.layout.λ)*(ηs-ηo) +
    0.5*(wing.inboard.layout.λ+wing.outboard.layout.λ)*(1.0 -ηs)

    Ko = 1.0/(wing.layout.AR*Kc)

    Kp = ηo +
    0.5*(1.0    +γs )*(ηs-ηo) +
    0.5*(γs +γt )*(1.0 -ηs) +
    fLo*ηo + 2.0*fLt*Ko*γt*wing.outboard.layout.λ

    po = (N*W - Lhtail)/(Kp*wing.outboard.layout.b)

    return po
end # wingpo


"""
    wingcl(b,bs,bo,
        λt,λs,γt,γs,
        sweep,AR,CL,CLhtail,fLo,fLt,
        duo,dus,dut)

Calculates section cl at  eta = ηo,ηs,1

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `b::Float64`, `bs::Float64`, `bo::Float64`: Span of wing, inner wing section (up to "snag"), and wing root, respectively.
    - `λt::Float64`, `λs::Float64`: Wing chord taper ratios at tip and break ("snag"), respectively.
    - `γt::Float64`, `γs::Float64`: Wing lift distribution "taper" ratios for outer and inner wing sections, respectively.
    - `sweep::Float64`: Wing sweep for reference axis.
    - `AR::Float64`: Wing aspect ratio, ``AR = b^2/S``.
    - `CL::Float64`, `CLhtail::Float64`: Overall lift coefficient of wing and horizontal tail, respectively.
    - `fLo::Float64`, `fLt::Float64`: Correction factors for lift of wingbox and tip.
    - `duo::Float64`, `dus::Float64`, `dut::Float64`: Velocity-change fractions at wing root, break ("snag"), and tip due to fuselage flow.

    **Outputs:**
    - `clo::Float64`, `cls::Float64`, `clt::Float64`: Section lift coefficient at root, wing break ("snag"), and tip.

See Sections 2.5 and 2.6 of the [TASOPT Technical Desc](@ref dreladocs). Called by `cdsum!`.
"""
function wingcl(b,bs,bo,
                λt,λs,γt,γs,
                sweep,AR,CL,CLhtail,fLo,fLt,
                duo,dus,dut)

      cosL = cos(sweep*pi/180.0)

      ηo = bo/b
      ηs = bs/b
     
      Kc = ηo +
	 0.5*(1.0    +λs)*(ηs-ηo) +
	 0.5*(λs+λt)*(1.0 -ηs) #surface area factor S = co*bo*K

      Ko = 1.0/(AR*Kc) #root chord def'n factor

      Kp = ηo + #planform loading factor
	 0.5*(1.0    +γs )*(ηs-ηo) +
	 0.5*(γs +γt )*(1.0 -ηs) +
	 fLo*ηo + 2.0*fLt*Ko*γt*λt

      cl1 = (CL-CLhtail)/cosL^2 * (Kc/Kp)

      clo = cl1                / (1.0+duo)^2
      cls = cl1*γs/λs / (1.0+dus)^2
      clt = cl1*γt/λt / (1.0+dut)^2

      return clo, cls, clt

end # wingcl


