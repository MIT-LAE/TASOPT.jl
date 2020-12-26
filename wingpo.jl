"""
Calculates wing root loading po to balance 
net load  N*W - Lhtail

## Inputs: 
- `b`, `bs`, `bo`: span, panel break location, wing root location
- `λt`, `λs` : inner and outer taper ratios
- `γt`,`γs` : inner and outer local cl factors - γt = rclt*λt etc.
- `AR`, `N`, `W`, `Lhtail` : Aspect ratio, Load factor, weight and H-tail lift
- `fLo`, `flt` : wing root and tip load adjustment factors

## Outputs:

- `po` : wing's root loading magnitude
"""
function wingpo(b,bs,bo,
               λt,λs,γt,γs,
               AR,N,W,Lhtail,fLo,fLt)
      
      etao = bo/b #calculate non-dim. span locations eta
      etas = bs/b
     
      Kc = etao +
	 0.5*(1.0    +λs)*(etas-etao) +
	 0.5*(λs+λt)*(1.0 -etas)

      Ko = 1.0/(AR*Kc)

      Kp = etao +
	 0.5*(1.0    +γs )*(etas-etao) +
	 0.5*(γs +γt )*(1.0 -etas) +
	 fLo*etao + 2.0*fLt*Ko*γt*λt

      po = (N*W - Lhtail)/(Kp*b)

      return po
end # wingpo


"""
Calculates section cl at  eta = etao,etas,1
"""
function wingcl(b,bs,bo,
                λt,λs,γt,γs,
                sweep,AR,CL,CLhtail,fLo,fLt,
                duo,dus,dut)

      cosL = cos(sweep*pi/180.0)

      etao = bo/b
      etas = bs/b
     
      Kc = etao +
	 0.5*(1.0    +λs)*(etas-etao) +
	 0.5*(λs+λt)*(1.0 -etas)

      Ko = 1.0/(AR*Kc)

      Kp = etao +
	 0.5*(1.0    +γs )*(etas-etao) +
	 0.5*(γs +γt )*(1.0 -etas) +
	 fLo*etao + 2.0*fLt*Ko*γt*λt

      cl1 = (CL-CLhtail)/cosL^2 * (Kc/Kp)

      clo = cl1                / (1.0+duo)^2
      cls = cl1*γs/λs / (1.0+dus)^2
      clt = cl1*γt/λt / (1.0+dut)^2

      return clo, cls, clt

end # wingcl


