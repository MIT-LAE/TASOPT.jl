"""
    wingpo(b, bs, bo, 
        λt, λs, γt, γs, 
        AR, N, W, Lhtail)

Computes wing root ("center") loading ``p_o`` to balance the net load.

```math
N*W - L_{h tail} \\times 2*∫p(η) dy + 2ΔL₀ + 2ΔLₜ = N*W - (L_{htail}).
```

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.wing`: Wing structure.
    - `rclt::Float64`: .
    - `rcls::Float64`: .
    - `N::Float64`: 
    - `W::Float64`: 
    - `Lhtail::Float64`: 

    **Outputs:**
    - `po::Float64`: Wing root loading magnitude.

See Section 2.6.2 of the [TASOPT Technical Desc](@ref dreladocs).
"""
function wingpo(wing, rclt, rcls, N, W, Lhtail)
    
    γt, γs = wing.outboard.layout.λ * rclt, wing.inboard.layout.λ * rcls
        
    ηo = wing.layout.root_span/wing.layout.span #calculate non-dim. span locations eta
    ηs = wing.inboard.layout.b/wing.layout.span
    
    Kc = ηo +
    0.5*(1.0    +wing.inboard.layout.λ)*(ηs-ηo) +
    0.5*(wing.inboard.layout.λ+wing.outboard.layout.λ)*(1.0 -ηs)

    Ko = 1.0/(wing.layout.AR*Kc)

    Kp = ηo +
    0.5*(1.0    +γs )*(ηs-ηo) +
    0.5*(γs +γt )*(1.0 -ηs) +
    wing.fuse_lift_carryover*ηo + 2.0*wing.tip_lift_loss*Ko*γt*wing.outboard.layout.λ

    po = (N*W - Lhtail)/(Kp*wing.layout.span)

    return po
end # wingpo


"""
    wingcl(wing,gammat,gammas,
            CL,CLhtail,
	        fduo,fdus,fdut)

Calculates section cl at  eta = ηo,ηs,1

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `Wing::TASOPT.Wing`: Wing Structure
    - `γt::Float64`, `γs::Float64`: Wing lift distribution "taper" ratios for outer and inner wing sections, respectively.
    - `CL::Float64`, `CLhtail::Float64`: Overall lift coefficient of wing and horizontal tail, respectively.
    - `duo::Float64`, `dus::Float64`, `dut::Float64`: Velocity-change fractions at wing root, break ("snag"), and tip due to fuselage flow.

    **Outputs:**
    - `clo::Float64`, `cls::Float64`, `clt::Float64`: Section lift coefficient at root, wing break ("snag"), and tip.

See Sections 2.5 and 2.6 of the [TASOPT Technical Desc](@ref dreladocs). Called by `cdsum!`.
"""
function wingcl(wing,γt,γs,
                CL,CLhtail,
                duo,dus,dut)

      cosL = cosd(wing.layout.sweep)

      ηo = wing.layout.root_span/wing.layout.span
      ηs = wing.inboard.layout.b/wing.layout.span
     
      Kc = ηo +
	 0.5*(1.0    +wing.inboard.layout.λ)*(ηs-ηo) +
	 0.5*(wing.inboard.layout.λ+wing.outboard.layout.λ)*(1.0 -ηs) #surface area factor S = co*bo*K

      Ko = 1.0/(wing.layout.AR*Kc) #root chord def'n factor

      Kp = ηo + #planform loading factor
	 0.5*(1.0    +γs )*(ηs-ηo) +
	 0.5*(γs +γt )*(1.0 -ηs) +
	 wing.fuse_lift_carryover*ηo + 2.0*wing.tip_lift_loss*Ko*γt*wing.outboard.layout.λ

      cl1 = (CL-CLhtail)/cosL^2 * (Kc/Kp)

      clo = cl1                / (1.0+duo)^2
      cls = cl1*γs/wing.inboard.layout.λ / (1.0+dus)^2
      clt = cl1*γt/wing.outboard.layout.λ / (1.0+dut)^2

      return clo, cls, clt

end # wingcl


