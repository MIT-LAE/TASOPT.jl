"""
     Sets wing area, span, root chord 
     to be consistent with q,CL,weight,AR
"""
      function wingsc(W,CL,qinf,AR,ηsi,bo,λt,λs)

      S = W/(qinf*CL)
      b = sqrt(S*AR)

      bs = max( b*ηsi , bo )

      ηo = bo/b
      ηs = bs/b

      Kc = ηo +
	 0.5*(1.0    +λs)*(ηs-ηo) +
	 0.5*(λs+λt)*(1.0 -ηs)

      co = S/(Kc*b)

      return  S,b,bs,co
      end # wingsc


#"""
#  Sets wing area, AR, root chord 
#  to be consistent with q,CL,weight,span
#"""
#      function wingAc(W,CL,qinf,b,ηsi,bo,bs,λt,λs)
#
#      S = W/(qinf*CL)
#
#      ηo = bo/b
#      ηs = bs/b
#
#      Kc = ηo +
#	 0.5*(1.0    +λs)*(ηs-ηo) +
#	 0.5*(λs+λt)*(1.0 -ηs)
#
#      co = S/(Kc*b)
#
#      AR = b^2 / S
#
#      return  S,AR,co
#      end # wingAc
#


