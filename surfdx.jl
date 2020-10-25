"""
     Calculates area centroid x-offset due to sweep
     Calculates mean aerodynamic chord/co
"""
function surfdx(b,bs,bo,lambdat,lambdas,sweep)

      tanL = tan(sweep*pi/180.0)

      etao = bo/b
      etas = bs/b

#---- 2 Int c dy /(co b)  =  S/(co b)  =  Kc
      Kc = etao +
	 0.5*(1.0    +lambdas)*(etas-etao) +
	 0.5*(lambdas+lambdat)*(1.0 -etas)

#---- 2 Int c (x-xo) dy /(co b^2) 
      Kcx = (1.0     + 2.0*lambdas)*(etas-etao)^2 / 12.0 +
	 (lambdas + 2.0*lambdat)*(1.0 -etas)^2 / 12.0 +
	 (lambdas +     lambdat)*(1.0 -etas)*(etas-etao) / 4.0

#---- 2 Int c^2 dy / (co^2 b)
      Kcc = etao +
	 (1.0        + lambdas         + lambdas^2)*(etas-etao)/3.0 +
	 (lambdas^2 + lambdas*lambdat + lambdat^2)*(1.0 -etas)/3.0

      dx    = Kcx/Kc * b*tanL
      macco = Kcc/Kc
   
      return dx, macco
end # surfdx

