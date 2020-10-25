"""
Calculates components of wing CM about wing root axis

    CM = CM0 + CM1*(CL-CLhtail)
"""
function surfcm(b,bs,bo, sweep, Xaxis,
                       lambdat,lambdas,gammat,gammas,
                       AR,fLo,fLt,cmpo,cmps,cmpt)

      cosL = cos(sweep*pi/180.0)
      tanL = tan(sweep*pi/180.0)

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

      C1 = (1.0 + 0.5*(lambdas+gammas ) + lambdas*gammas)*(etas-etao) +
	 (lambdas*gammas +
	 0.5*(lambdas*gammat +
	 gammas *lambdat) + lambdat*gammat)*(1.0 -etas)

      C2 =     (1.0    + 2.0*gammas)*(etas-etao)^2 +
	 (gammas + 2.0*gammat)*(1.0 -etas)^2 +
	 3.0*(gammas +     gammat)*(etas-etao)*(1.0-etas)

      C3 = ( cmpo*(3.0            + 2.0*lambdas         + lambdas^2) +
	 cmps*(3.0*lambdas^2 + 2.0*lambdas         + 1.0       ) ) *
	 (etas-etao) +
	 ( cmps*(3.0*lambdas^2 + 2.0*lambdas*lambdat + lambdat^2) +
	 cmpt*(3.0*lambdat^2 + 2.0*lambdas*lambdat + lambdas^2) ) *
	 (1.0-etas)

      CM1 = (1.0/Kp) *
	 (  etao*(1.0+fLo)*(Xaxis-0.25) +
	 (Xaxis-0.25)*cosL^2 * C1/3.0 -
	 (tanL/Ko) * C2/12.0 +
	 2.0*fLt*lambdat*gammat *
	( Ko*lambdat*(Xaxis-0.25)*cosL^2 -
	 0.5*(1.0-etao)*tanL) )

      CM0 = (cosL^4/Kc) * C3/12.0


      return CM0,CM1
end # surfcm
