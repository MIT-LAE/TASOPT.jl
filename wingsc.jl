"""
     Sets wing area, span, root chord 
     to be consistent with q,CL,weight,AR
"""
      function wingsc(W,CL,qinf,AR,etasi,bo,lambdat,lambdas)

      S = W/(qinf*CL)
      b = sqrt(S*AR)

      bs = max( b*etasi , bo )

      etao = bo/b
      etas = bs/b

      Kc = etao +
	 0.5*(1.0    +lambdas)*(etas-etao) +
	 0.5*(lambdas+lambdat)*(1.0 -etas)

      co = S/(Kc*b)

      return  S,b,bs,co
      end # wingsc


#"""
#  Sets wing area, AR, root chord 
#  to be consistent with q,CL,weight,span
#"""
#      function wingAc(W,CL,qinf,b,etasi,bo,bs,lambdat,lambdas)
#
#      S = W/(qinf*CL)
#
#      etao = bo/b
#      etas = bs/b
#
#      Kc = etao +
#	 0.5*(1.0    +lambdas)*(etas-etao) +
#	 0.5*(lambdas+lambdat)*(1.0 -etas)
#
#      co = S/(Kc*b)
#
#      AR = b^2 / S
#
#      return  S,AR,co
#      end # wingAc
#


