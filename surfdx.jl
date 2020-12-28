"""
Calculates area centroid x-offset due to sweep
Calculates mean aerodynamic chord/co
"""
function surfdx(b,bs,bo,λt,λs,sweep)

      tanL = tan(sweep * π/180.0)

      ηo = bo/b
      ηs = bs/b

#---- 2 Int c dy /(co b)  =  S/(co b)  =  Kc
      Kc = ηo +
	 0.5*(1.0    +λs)*(ηs-ηo) +
	 0.5*(λs+λt)*(1.0 -ηs)

#---- 2 Int c (x-xo) dy /(co b^2) 
      Kcx = (1.0     + 2.0*λs)*(ηs-ηo)^2 / 12.0 +
	 (λs + 2.0*λt)*(1.0 -ηs)^2 / 12.0 +
	 (λs +     λt)*(1.0 -ηs)*(ηs-ηo) / 4.0

#---- 2 Int c^2 dy / (co^2 b)
      Kcc = ηo +
	 (1.0        + λs         + λs^2)*(ηs-ηo)/3.0 +
	 (λs^2 + λs*λt + λt^2)*(1.0 -ηs)/3.0

      dx    = Kcx/Kc * b*tanL
      macco = Kcc/Kc
   
      return dx, macco
end # surfdx

