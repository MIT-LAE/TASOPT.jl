"""
Calculates wing root loading po to balance 
net load  N*W - Lhtail

## Inputs: 
- `b`, `bs`, `bo`: span, panel break location, wing root location
- `lambdat`, `lambdas` : inner and outer taper ratios
- `gammat`,`gammas` : inner and outer local cl factors - gammat = rclt*lambdat etc.
- `AR`, `N`, `W`, `Lhtail` : Aspect ratio, Load factor, weight and H-tail lift
- `fLo`, `flt` : wing root and tip load adjustment factors

## Outputs:

- `po` : wing's root loading magnitude
"""
function wingpo(b,bs,bo,
               lambdat,lambdas,gammat,gammas,
               AR,N,W,Lhtail,fLo,fLt)
      
      etao = bo/b #calculate non-dim. span locations eta
      etas = bs/b
     
      Kc = etao +
	 0.5*(1.0    +lambdas)*(etas-etao) +
	 0.5*(lambdas+lambdat)*(1.0 -etas)

      Ko = 1.0/(AR*Kc)

      Kp = etao +
	 0.5*(1.0    +gammas )*(etas-etao) +
	 0.5*(gammas +gammat )*(1.0 -etas) +
	 fLo*etao + 2.0*fLt*Ko*gammat*lambdat

      po = (N*W - Lhtail)/(Kp*b)

      return po
end # wingpo


"""
Calculates section cl at  eta = etao,etas,1
"""
function wingcl(b,bs,bo,
                lambdat,lambdas,gammat,gammas,
                sweep,AR,CL,CLhtail,fLo,fLt,
                duo,dus,dut)

      cosL = cos(sweep*pi/180.0)

      etao = bo/b
      etas = bs/b
     
      Kc = etao +
	 0.5*(1.0    +lambdas)*(etas-etao) +
	 0.5*(lambdas+lambdat)*(1.0 -etas)

      Ko = 1.0/(AR*Kc)

      Kp = etao +
	 0.5*(1.0    +gammas )*(etas-etao) +
	 0.5*(gammas +gammat )*(1.0 -etas) +
	 fLo*etao + 2.0*fLt*Ko*gammat*lambdat

      cl1 = (CL-CLhtail)/cosL^2 * (Kc/Kp)

      clo = cl1                / (1.0+duo)^2
      cls = cl1*gammas/lambdas / (1.0+dus)^2
      clt = cl1*gammat/lambdat / (1.0+dut)^2

      return clo, cls, clt

end # wingcl


