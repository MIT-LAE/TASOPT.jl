"""
fuseW sizes the fuselage and calculates the component weights

Inputs:

- Geometry `xnose`, `Rfuse`, etc..
- Fixed weights `Wfix`, `Wpay`, `Wseat` etc...
- Material props `sigskin`, `rhoskin`, `E`, `G`, etc...

Outputs:

- `EI`, `GJ`
- Thickness `tskin` etc
- Weights `Wtank`
- Moments `xWtank`
- Tank Volume `tankVol`
"""
function tankw(gee,Nland,Wfuel, densityFuel,Wpadd,Weng,
                      fstring,fframe,ffadd,deltap,
                      Wppinsul,
                      Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
                      bv,lambdav,nvtail,
                      Rfuse,dRfuse,wfb,nfweb,lambdac,
                      xnose,xshell1,xshell2,
                      xhtail,xvtail,
                      cbox,
                      xeng,
                      hfloor,
                      sigskin,sigbend, rhoskin,rhobend,
                      Eskin,Ebend,Gskin)
# -------------------------------------------------------------
#Tank sizing and weight routine


#--- floor beam material properties
#     (assumed same as stringers, but could be different)
      sigfloor = sigbend
      taufloor = sigbend
      rhofloor = rhobend

      rE = Ebend/Eskin

#--- effective pressure-vessel length

      lshell = xshell2 - xshell1
      lfloor = xshell2 - xshell1 + 2.0*Rfuse

      xshell = 0.5*(xshell1+xshell2)

#--- fuselage cross-section geometric parameters
      wfblim = max( min( wfb , Rfuse ) , 0.0 )
      thetafb = asin(wfblim/Rfuse)
      hfb = sqrt(Rfuse^2 - wfb^2)
      sin2t = 2.0*hfb*wfb/Rfuse^2
      cost  = hfb/Rfuse

      perim = (2.0*pi + 4.0*thetafb)*Rfuse + 2.0*dRfuse

#--------------------------------------------------------------------
#--- fuselage skin and center web thicknesses to withstand pressure load
      tskin =     deltap*Rfuse/sigskin
      tfweb = 2.0*deltap*wfb  /sigskin

#--- cross-sectional areas
      Askin = (2.0*pi+4.0*nfweb*thetafb)*Rfuse*tskin + 2.0*dRfuse*tskin
      #Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      Atank = (pi + nfweb*(2.0*thetafb + sin2t))*Rfuse^2 + 2.0*Rfuse*dRfuse        #####
#          + 2.0*(Rfuse+nfweb*wfb)*dRfuse        #####

#--- nose + rear bulkhead surface areas
      Sbulk = (2.0*pi + 4.0*nfweb*thetafb)*Rfuse^2
      Snose = (2.0*pi + 4.0*nfweb*thetafb)*Rfuse^2* ( 0.333 + 0.667*(lnose/Rfuse)^1.6 )^0.625

#--- component volumes and volume moments
      Vcyl  = Askin*lshell
      Vbulk = Sbulk*tskin
      #Vfweb = Afweb*lshell

      xVcyl  = Vcyl  * 0.5*(xshell1 + xshell2)
      xVnose = Vnose * 0.5*(xnose   + xshell1)
      xVbulk = Vbulk *     (xshell2 + 0.5*Rfuse)
      #xVfweb = Vfweb * 0.5*(xshell1 + xshell2)

#--- weights and weight moments
      Wskin = rhoskin*gee*(Vcyl)
      #Wfweb = rhoskin*gee* Vfweb
      xWskin = rhoskin*gee*(xVcyl)
      #xWfweb = rhoskin*gee* xVfweb

      Wshell = Wskin*(1.0+fstring+fframe+ffadd) #+ Wfweb
      xWshell = xWskin*(1.0+fstring+fframe+ffadd) #+ xWfweb


#--------------------------------------------------------------------
#--- insulation weight
      Winsul = Wppinsul*((1.1*pi+2.0*thetafb)*Rfuse*lshell + 0.55*(Snose+Sbulk))
      xWinsul = Winsul * 0.5*(xshell1 + xshell2)

#--------------------------------------------------------------------
#--- various weight moments
      #xWfix  = Wfix  * xfix
      xWfuel = Wfuel * 0.5*(xshell1 + xshell2)

#--------------------------------------------------------------------
#--- floor structural sizing
      P = (Wfuel+Wseat) * Nland
      wfloor1 = wfb + Rfuse

      if (wfb == 0.0)
#---- full-width floor
       Smax = 0.50 * P
       Mmax = 0.25 * P*wfloor1
      else
#---- floor with center support
       Smax = (5.0/16.0 ) * P
       Mmax = (9.0/256.0) * P*wfloor1
      end

      Afweb = 1.5*Smax/ taufloor
      Afcap = 2.0*Mmax/(sigfloor*hfloor)

      Vfloor = (Afcap + Afweb) * 2.0*wfloor1
      Wfloor = rhofloor*gee*Vfloor + 2.0*wfloor1*lfloor*Wppfloor
      xWfloor = Wfloor * 0.5*(xshell1 + xshell2)

#--- average floor-beam cap thickness ("smeared" over entire floor)
      tfloor = 0.5*Afcap/lfloor

#--------------------------------------------------------------------
#--- torsional stiffnesses
      GJshell = Gskin*4.0*Afuse^2 * tskin / perim
      GJcone  = Gskin*4.0*Afuse^2 * tcone / perim


#--------------------------------------------------------------------
#--- shell bending inertias
      tshell = tskin*(1.0 + rE*fstring*rhoskin/rhobend)

#--- bending stress to be matched
      sigMh = sigbend - rE*0.5*deltap*Rfuse/tshell
      sigMv = sigbend - rE*0.5*deltap*Rfuse/tshell

#--- rear bulkhead where tail structure attaches
      xbulk = xshell2

#--- horizontal-axis bending moment added material
      Ihshell = ( (pi + nfweb*(2.0*thetafb+sin2t))*Rfuse^2 +
                 8.0*nfweb*cost*0.5*dRfuse*Rfuse +
                 (2.0*pi + 4.0*nfweb*thetafb)*(0.5*dRfuse)^2)*Rfuse*tshell+
              0.66667*nfweb*(hfb+0.5*dRfuse)^3*tfweb

      hfuse = Rfuse + 0.5*dRfuse
      A2 = 1.0/(hfuse*sigMh)*
          Nland*(Wpadd+Wfuel+Wshell+Winsul+Wfloor)*
          0.5/lshell
      A1 = 1.0/(hfuse*sigMh)*
          (Nland*Wtail + rMh*Lhmax)
      A0 = -Ihshell/(rE*hfuse^2)
      Abar2 = A2
      Abar1 = 2.0*A2*xbulk + A1
      Abar0 = A2*xbulk^2 + A1*xtail + A0
      desc = max( 0.0 , Abar1^2 - 4.0*Abar0*Abar2 )
      xhbend = (Abar1 - sqrt(desc))*0.5/Abar2


      Ahbendf = Abar2*xf^2 - Abar1*xf + Abar0
      Ahbendb = Abar2*xb^2 - Abar1*xb + Abar0

      Vhbendf = A2*((xbulk-xf)^3 - (xbulk-xhbend)^3)/3.0 +
               A1*((xtail-xf)^2 - (xtail-xhbend)^2)/2.0 +
               A0*(xhbend-xf)
      Vhbendb = A2*((xbulk-xb)^3 - (xbulk-xhbend)^3)/3.0 +
               A1*((xtail-xb)^2 - (xtail-xhbend)^2)/2.0 +
              A0*(xhbend-xb)
      Vhbendc = 0.5*(Ahbendf+Ahbendb)*cbox


      Whbend = rhobend*gee*(Vhbendf + Vhbendb + Vhbendc)

      xWhbend = Whbend * xwing

      EIhshell = Eskin * Ihshell
      EIhbend  = Ebend * 0.5*(Ahbendf+Ahbendb) * 2.0*hfuse^2

#----------------------------------------------------------------
#--- vertical-axis bending moment added material
      nk = Int(round(nfweb/2.0+0.001))
      ik = mod(Int(round(nfweb+0.001))+1,2)
      ksum = 0.
      for k = 1: nk
        ksum = ksum + float(2*k-ik)^2
      end
      Ivshell = ( (pi + nfweb*(2.0*thetafb - sin2t))*Rfuse^2+
                 8.0*cost*nfweb*wfb*Rfuse +
                 (2.0*pi + 4.0*thetafb)*(nfweb*wfb)^2+
                 4.0*thetafb*wfb^2 * ksum )*Rfuse*tshell

      widf = Rfuse + nfweb*wfb
      B1 = 1.0/(widf*sigMv) * (rMv*Lvmax*nvtail)
      B0 = -Ivshell/(rE*widf^2)
      xvbend = xvtail + B0/B1

      Avbendb = B1*(xtail-xb) + B0
      Vvbendb = B1*((xtail-xb)^2 - (xtail-xvbend)^2)/2.0+
               B0*(xvbend-xb)
      Vvbendc = 0.5*Avbendb*cbox
      Wvbend = rhobend*gee*(Vvbendb + Vvbendc)

      xWvbend = Wvbend * (2.0*xwing + xvbend)/3.0


      EIvshell = Eskin * Ivshell
      EIvbend  = Ebend * Avbendb * 2.0*widf^2


#----------------------------------------------------------------
#--- overall tank weight and moment
      Wtank = Wpadd +
             Wshell + Winsul +
             Whbend + Wvbend + Wfuel

      xWtank = xWpadd +
             xWshell + xWinsul +
             xWhbend + xWvbend

#----------------------------------------------------------------
#--- pressurized tank volume
      tankVol = Afuse*(lshell + 0.67*Rfuse)
      #Questions: Is stiffness a desired output?
      #Question: Can we eliminate Wfloor from the code
      #Question: Eliminate Wpadd? What outputs do we need?
      #Eliminate cone parameters?
      #Will we ever use this on a double bubble config? (nfweb etc.)
      #Replaced Wpay with Wfuel

return  tskin, tfloor, xhbend, xvbend,
                       Wshell, Wcone, Winsul,
                       Whbend, Wvbend,
                       Wtank,
                      xWtank,
                       tankVol
end # fusew
