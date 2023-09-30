"""
    wingpo(b, bs, bo, 
        λt, λs, γt, γs, 
        AR, N, W, Lhtail, fLo, fLt)

Calculates wing root loading po to balance net load 

```math
N*W - L_{h tail} \\times 2*∫p(η) dy + 2ΔL₀ + 2ΔLₜ = N*W - (L_{htail}).
```


# Inputs
- `b::Float64`: span in [m]
- `bs::Float64`: panel break location [m]
- `bo::Float64`:  wing root location [m]
- `λt::Float64`, `λs::Float64` : inner and outer taper ratios.
- `γt::Float64`,`γs::Float64` : inner and outer local ``c_l`` factors - ``\\gamma_t = r_{c_l,t}\\times \\lambda_t`` etc.
- `AR::Float64`, `N::Float64`, `W::Float64`, `Lhtail` : Aspect ratio, Load factor, weight and H-tail lift.
- `fLo::Float64`, `flt::Float64` : wing root and tip load adjustment factors.

# Outputs
- `po::Float64`: wing's root loading magnitude.

See Eqn. 154 to 160 of TASOPT docs.
"""
function wingpo(b, bs, bo,
               λt, λs, γt, γs,
               AR, N, W, Lhtail, fLo, fLt)
      
      ηo = bo/b #calculate non-dim. span locations eta
      ηs = bs/b
     
      Kc = ηo +
	 0.5*(1.0    +λs)*(ηs-ηo) +
	 0.5*(λs+λt)*(1.0 -ηs)

      Ko = 1.0/(AR*Kc)

      Kp = ηo +
	 0.5*(1.0    +γs )*(ηs-ηo) +
	 0.5*(γs +γt )*(1.0 -ηs) +
	 fLo*ηo + 2.0*fLt*Ko*γt*λt

      po = (N*W - Lhtail)/(Kp*b)

      return po
end # wingpo


"""
Calculates section cl at  eta = ηo,ηs,1
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
	 0.5*(λs+λt)*(1.0 -ηs)

      Ko = 1.0/(AR*Kc)

      Kp = ηo +
	 0.5*(1.0    +γs )*(ηs-ηo) +
	 0.5*(γs +γt )*(1.0 -ηs) +
	 fLo*ηo + 2.0*fLt*Ko*γt*λt

      cl1 = (CL-CLhtail)/cosL^2 * (Kc/Kp)

      clo = cl1                / (1.0+duo)^2
      cls = cl1*γs/λs / (1.0+dus)^2
      clt = cl1*γt/λt / (1.0+dut)^2

      return clo, cls, clt

end # wingcl


