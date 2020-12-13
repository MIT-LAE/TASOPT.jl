"""
tankWmech calculates weight of the tank that holds the LH2

Inputs:

-LD: L/D of airplane
-

Outputs:

- `EI`
"""
function tankWmech(gee, rhoFuel,
                      fstring,fframe,ffadd,deltap,
                      Wppinsul,
                      rMh,rMv,Lhmax,Lvmax,
                      bv,lambdav,
                      Rfuse,dRfuse,wfb,nfweb,lambdac,
                      xshell1,xshell2,
                      sigskin,sigbend, rhoskin,rhobend,
                      Eskin,Ebend,Gskin, m_airplane, range, lcv, eta, LD)

#--- effective pressure-vessel length
      lshell = xshell2 - xshell1
      xshell = 0.5*(xshell1+xshell2)


#--- Calculate Wfuel
     m_st = m_airplane * exp(range * gee / (lcv * eta * LD))
     m_fuel = m_st - m_airplane

#--- fuselage cross-section geometric parameters
      wfblim = max( min( wfb , Rfuse ) , 0.0 )
      thetafb = asin(wfblim/Rfuse)
      hfb = sqrt(Rfuse^2 - wfb^2)
      sin2t = 2.0*hfb*wfb/Rfuse^2
      cost  = hfb/Rfuse
      perim = (2.0*pi + 4.0*thetafb)*Rfuse + 2.0*dRfuse

#--- fuselage skin and center web thicknesses to withstand pressure load
      tskin =     deltap*Rfuse/sigskin
      tfweb = 2.0*deltap*wfb  /sigskin

#--- cross-sectional areas
      Askin = (2.0*pi+4.0*nfweb*thetafb)*Rfuse*tskin + 2.0*dRfuse*tskin
      Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      Atank = (pi + nfweb*(2.0*thetafb + sin2t))*Rfuse^2 + 2.0*Rfuse*dRfuse
#          + 2.0*(Rfuse+nfweb*wfb)*dRfuse
#--- component volumes and volume moments
      Vcyl  = Askin*lshell
      xVcyl  = Vcyl  * 0.5*(xshell1 + xshell2)
      #xVfweb = Vfweb * 0.5*(xshell1 + xshell2)

#--- weights and weight moments
      Wskin = rhoskin*gee*(Vcyl)
      #Wfweb = rhoskin*gee* Vfweb
      xWskin = rhoskin*gee*(xVcyl)
      #xWfweb = rhoskin*gee* xVfweb

      Wshell = Wskin*(1.0+fstring+fframe+ffadd)
      xWshell = xWskin*(1.0+fstring+fframe+ffadd)

#--- insulation weight
      Winsul = Wppinsul*(1.1*pi+2.0*thetafb)*Rfuse*lshell# + 0.55*(Snose+Sbulk))
      xWinsul = Winsul * 0.5*(xshell1 + xshell2)

#--- various weight moments
      #xWfix  = Wfix  * xfix
      xWfuel = Wfuel * 0.5*(xshell1 + xshell2)

#--- shell bending inertias
      tshell = tskin*(1.0 + rE*fstring*rhoskin/rhobend)

#--- bending stress to be matched
      sigMh = sigbend - rE*0.5*deltap*Rfuse/tshell
      sigMv = sigbend - rE*0.5*deltap*Rfuse/tshell

#--- rear bulkhead where tail structure attaches
      xbulk = xshell2

#--- overall tank weight and moment
      Wtank = Wshell + Winsul + Wfuel
             #Whbend + Wvbend + Wfuel

#--- pressurized tank volume
      #tankVol = Afuse*(lshell + 0.67*Rfuse)
      tankVol = Wfuel/rhoFuel

return  tskin, Wtank, tankVol
end # fusew
