"""
tankWmech calculates weight of the tank that holds the LH2

Inputs:

-LD is L/D of airplane
-xshell1, xshell2 are start and end coordinates of tank
-gee is gravitational acceleration
-Rfuse is fuselage radius, dRfuse is the subtraction factor that accounts for flatness of fuselage at bottom
-wfb, nfb are parameters for multiple-bubble configuration
-m_airplane is airplane mass
-range is specified for a given mission
-lcv is lower calorific value of fuel
-eta is overall efficiency of gas turbine/turboelectric powertrain etc.
-sigskin, sigbend, rhoskin, rhobend are material properties

Outputs:

- Wtank: tank weight
- tankVol: tank volume
"""
function tankWmech(gee, rhoFuel,
                      fstring,fframe,ffadd,deltap,
                      Rfuse,dRfuse,wfb,nfweb,
                      xshell1,xshell2,
                      sigskin,sigbend, Wppinsul, rhoskin,rhobend,
                      m_airplane, range, lcv, eta, LD)

#--- effective pressure-vessel length
      lshell = xshell2 - xshell1

#--- Calculate Wfuel
      m_st = m_airplane * exp(range * gee / (lcv * eta * LD))
      Wfuel = (m_st - m_airplane)*gee

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

#--- areas
      Askin = (2.0*pi+4.0*nfweb*thetafb)*Rfuse*tskin + 2.0*dRfuse*tskin
      Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      Atank = (pi + nfweb*(2.0*thetafb + sin2t))*Rfuse^2 + 2.0*Rfuse*dRfuse + 2.0*(Rfuse+nfweb*wfb)*dRfuse
#--- component volumes
      Vcyl  = Askin*lshell

#--- weights and weight moments
      Wtank = rhoskin*gee*(Vcyl)
      Wtank = Wtank*(1.0+fstring+fframe+ffadd)

#--- insulation weight
      Winsul = Wppinsul*(1.1*pi+2.0*thetafb)*Rfuse*lshell

#--- overall tank weight and moment
      Wtank = Wshell + Winsul + Wfuel

#--- pressurized tank volume
      #tankVol = Afuse*(lshell + 0.67*Rfuse)
      tankVol = Wfuel/rhoFuel

return  tskin, Wtank, tankVol
end # fusew
