"""
    surfcm(b, bs, bo, sweep, Xaxis, 
		λt, λs, γt, γs, 
		AR, fLo, fLt, cmpo, cmps, cmpt)

Calculates components of wing CM about wing root axis:

``C_M = C_{M,0} + C_{M,1} (C_L - C_{L,htail})``

``ΔC_{m, wing} = ΔC_{m, 0} + dCₘ/dCL × (C_L - C_{L,h})``


See eqn 172 of TASOPT docs.
See also [`surfcd`](@ref) and [`surfcd2`](@ref).
"""
function surfcm(b,bs,bo, sweep, Xaxis,
                       λt,λs,γt,γs,
                       AR,fLo,fLt,cmpo,cmps,cmpt)

      cosL = cos(sweep*pi/180.0)
      tanL = tan(sweep*pi/180.0)

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

      C1 = (1.0 + 0.5*(λs+γs ) + λs*γs)*(ηs-ηo) +
	     (λs*γs +
	 0.5*(λs*γt + γs *λt) + λt*γt)*(1.0 -ηs)

      C2 =     (1.0    + 2.0*γs)*(ηs-ηo)^2 +
	 (γs + 2.0*γt)*(1.0 -ηs)^2 +
	 3.0*(γs +     γt)*(ηs-ηo)*(1.0-ηs)

      C3 = ( cmpo*(3.0            + 2.0*λs         + λs^2) +
	         cmps*(3.0*λs^2 + 2.0*λs         + 1.0       ) ) * (ηs-ηo) +
	       ( cmps*(3.0*λs^2 + 2.0*λs*λt + λt^2) +
	         cmpt*(3.0*λt^2 + 2.0*λs*λt + λs^2) ) * (1.0-ηs)

      CM1 = (1.0/Kp) *
	 (  ηo*(1.0+fLo)*(Xaxis-0.25) +
	 (Xaxis-0.25)*cosL^2 * C1/3.0 -
	 (tanL/Ko) * C2/12.0 +
	 2.0*fLt*λt*γt *
	( Ko*λt*(Xaxis-0.25)*cosL^2 -
	 0.5*(1.0-ηo)*tanL) )

      CM0 = (cosL^4/Kc) * C3/12.0


      return CM0,CM1
end # surfcm
