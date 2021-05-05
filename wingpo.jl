"""
Calculates wing root loading po to balance 
net load  N*W - Lhtail
2*∫̃p(η) dy + 2*ΔL₀ + 2ΔLₜ = N*W - (Lhtail) Eqn. 154 to 160 of TASOPT docs

## Inputs: 
- `b`, `bs`, `bo`: span, panel break location, wing root location
- `λt`, `λs` : inner and outer taper ratios
- `γt`,`γs` : inner and outer local cl factors - γt = rclt*λt etc.
- `AR`, `N`, `W`, `Lhtail` : Aspect ratio, Load factor, weight and H-tail lift
- `fLo`, `flt` : wing root and tip load adjustment factors

## Outputs:

- `po` : wing's root loading magnitude
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


