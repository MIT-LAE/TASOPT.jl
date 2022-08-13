"""
fuseW sizes the fuselage and calculates the component weights

Inputs: 

- Geometry `xnose`, `Rfuse`, etc..
- Fixed weights `Wfix`, `Wpay`, `Wseat` etc...
- Material props `sigskin`, `rhoskin`, `E`, `G`, etc...

Outputs:

- `EI`, `GJ` 
- Thickness `tskin` etc
- Weights `Wfuse`
- Moments `xWfuse`
- Cabin Volume `cabVol`
"""
function fusew(gee,Nland,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,Waftfuel,
                      fstring,fframe,ffadd,deltap,
                      Wpwindow,Wppinsul,Wppfloor,
                      Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
                      bv,lambdav,nvtail,
                      Rfuse,dRfuse,wfb,nfweb,lambdac,
                      xnose,xshell1,xshell2,xconend,
                      xhtail,xvtail,
                      xwing,xwbox,cbox,
                      xfix,xapu,xeng,xfuel,
                      hfloor,
                      sigskin,sigbend, rhoskin,rhobend, 
                      Eskin,Ebend,Gskin)
# -------------------------------------------------------------
#Fuselage sizing and weight routine
#-------------------------------------------------------------
#--- cone material properties
#     (assumed same as skin, but could be different)
      taucone = sigskin
      rhocone = rhoskin

#--- floor beam material properties 
#     (assumed same as stringers, but could be different)
      sigfloor = sigbend
      taufloor = sigbend
      rhofloor = rhobend

      rE = Ebend/Eskin

#--- effective nose length and pressure-vessel length
      lnose  = xshell1 - xnose
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
      Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      Afuse = (pi + nfweb*(2.0*thetafb + sin2t))*Rfuse^2 + 2.0*Rfuse*dRfuse        #####
#          + 2.0*(Rfuse+nfweb*wfb)*dRfuse        #####

#--- nose + rear bulkhead surface areas
      Sbulk = (2.0*pi + 4.0*nfweb*thetafb)*Rfuse^2
      Snose = (2.0*pi + 4.0*nfweb*thetafb)*Rfuse^2* ( 0.333 + 0.667*(lnose/Rfuse)^1.6 )^0.625

#--- component volumes and volume moments
      Vcyl  = Askin*lshell
      Vnose = Snose*tskin
      Vbulk = Sbulk*tskin
      Vfweb = Afweb*lshell

      xVcyl  = Vcyl  * 0.5*(xshell1 + xshell2)
      xVnose = Vnose * 0.5*(xnose   + xshell1)
      xVbulk = Vbulk *     (xshell2 + 0.5*Rfuse)
      xVfweb = Vfweb * 0.5*(xshell1 + xshell2)

#--- weights and weight moments
      Wskin = rhoskin*gee*(Vcyl + Vnose + Vbulk)
      Wfweb = rhoskin*gee* Vfweb
      xWskin = rhoskin*gee*(xVcyl + xVnose + xVbulk)
      xWfweb = rhoskin*gee* xVfweb

      Wshell = Wskin*(1.0+fstring+fframe+ffadd) + Wfweb
      xWshell = xWskin*(1.0+fstring+fframe+ffadd) + xWfweb

#--------------------------------------------------------------------
#--- window weight
      Wwindow = Wpwindow * lshell
      xWwindow = Wwindow * 0.5*(xshell1 + xshell2)
#--------------------------------------------------------------------
#--- insulation weight
      Winsul = Wppinsul*((1.1*pi+2.0*thetafb)*Rfuse*lshell + 0.55*(Snose+Sbulk))
      xWinsul = Winsul * 0.5*(xshell1 + xshell2)

#--------------------------------------------------------------------
#--- various weight moments
      xWfix  = Wfix  * xfix
      xWapu  = Wapu  * xapu
      xWseat = Wseat * 0.5*(xshell1 + xshell2)
      xWpadd = Wpadd * 0.5*(xshell1 + xshell2)

#--------------------------------------------------------------------
#--- floor structural sizing
      P = (Wpay+Wseat) * Nland
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
#--- size tailcone to withstand vertical-tail torsion Qv
      Qv = (nvtail*Lvmax*bv/3.0)*(1.0+2.0*lambdav)/(1.0+lambdav)
      Vcone = (Qv/taucone)* 
             (pi + nfweb* 2.0*thetafb)/
            (pi + nfweb*(2.0*thetafb+sin2t))*
            (xconend-xshell2)/Rfuse *
            2.0/(1.0+lambdac)

      Wcone = rhocone*gee*Vcone*(1.0+fstring+fframe+ffadd)
      xWcone = Wcone * 0.5*(xshell2 + xconend)

      tcone = Qv / (2.0*taucone*Afuse)

#--------------------------------------------------------------------
#--- torsional stiffnesses
      GJshell = Gskin*4.0*Afuse^2 * tskin / perim
      GJcone  = Gskin*4.0*Afuse^2 * tcone / perim

#--------------------------------------------------------------------
#--- lumped tail weight and location  
#      (Weng=0 if there are no tail-mounted engines)
      Wtail = Whtail + Wvtail + Wcone + Wapu + Weng + Waftfuel
      xtail = (  xhtail*Whtail +
               xvtail*Wvtail +
               xWcone +
               xapu*Wapu +
               xeng*Weng +
               xfuel*Waftfuel) / Wtail

#--------------------------------------------------------------------
#--- shell bending inertias
      tshell = tskin*(1.0 + rE*fstring*rhoskin/rhobend)

#--- bending stress to be matched
      sigMh = sigbend - rE*0.5*deltap*Rfuse/tshell
      sigMv = sigbend - rE*0.5*deltap*Rfuse/tshell

#--- rear bulkhead where tail structure attaches
      xbulk = xshell2

#----------------------------------------------------------------
#--- horizontal-axis bending moment added material
      Ihshell = ( (pi + nfweb*(2.0*thetafb+sin2t))*Rfuse^2 +
                 8.0*nfweb*cost*0.5*dRfuse*Rfuse +
                 (2.0*pi + 4.0*nfweb*thetafb)*(0.5*dRfuse)^2)*Rfuse*tshell+
              0.66667*nfweb*(hfb+0.5*dRfuse)^3*tfweb

      hfuse = Rfuse + 0.5*dRfuse
      A2 = 1.0/(hfuse*sigMh)*
          Nland*(Wpay+Wpadd+Wshell+Wwindow+Winsul+Wfloor+Wseat)*
          0.5/lshell
      A1 = 1.0/(hfuse*sigMh)*
          (Nland*Wtail + rMh*Lhmax)
      A0 = -Ihshell/(rE*hfuse^2)
      Abar2 = A2
      Abar1 = 2.0*A2*xbulk + A1
      Abar0 = A2*xbulk^2 + A1*xtail + A0
      desc = max( 0.0 , Abar1^2 - 4.0*Abar0*Abar2 )
      xhbend = (Abar1 - sqrt(desc))*0.5/Abar2

      dxwing = xwing - xwbox
      xf = xwing + dxwing + 0.5*cbox
      xb = xwing - dxwing + 0.5*cbox

      Ahbendf = max(Abar2*xf^2 - Abar1*xf + Abar0, 0)
      Ahbendb = max(Abar2*xb^2 - Abar1*xb + Abar0, 0)

      Vhbendf = A2*((xbulk-xf)^3 - (xbulk-xhbend)^3)/3.0 +
               A1*((xtail-xf)^2 - (xtail-xhbend)^2)/2.0 +
               A0*(xhbend-xf)
      Vhbendb = A2*((xbulk-xb)^3 - (xbulk-xhbend)^3)/3.0 +
               A1*((xtail-xb)^2 - (xtail-xhbend)^2)/2.0 +
              A0*(xhbend-xb)
      Vhbendc = 0.5*(Ahbendf+Ahbendb)*cbox


      Whbend = rhobend*gee*(Vhbendf + Vhbendb + Vhbendc)

      xWhbend = Whbend *      xwing

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
                 # ^This simplifies to π×Rfuse³×tshell 
                 # which is the area moment of inertia for a thin walled tube

      widf = Rfuse + nfweb*wfb
      B1 = 1.0/(widf*sigMv) * (rMv*Lvmax*nvtail)
      B0 = -Ivshell/(rE*widf^2)
      xvbend = xvtail + B0/B1 # point where Avbend = 0 

      Avbendb = max(B1*(xtail-xb) + B0, 0)
      Vvbendb = max(B1*((xtail-xb)^2 - (xtail-xvbend)^2)/2.0+
               B0*(xvbend-xb), 0)
      Vvbendc = 0.5*Avbendb*cbox
      Wvbend = rhobend*gee*(Vvbendb + Vvbendc)
      # println("W, Vvb, Vvc Avbend = $Wvbend, $Vvbendb, $Vvbendc, $Avbendb")

      xWvbend = Wvbend * (2.0*xwing + xvbend)/3.0


      EIvshell = Eskin * Ivshell
      EIvbend  = Ebend * Avbendb * 2.0*widf^2

#c      write(*,*) B0*(xvbend-xwbox)
#c      write(*,*) B1*((xvtail-xwbox)^2 - (xvtail-xvbend)^2)/2.0
#c      write(*,*) 'I', Ihshell, Ivshell, tshell
#c      write(*,*) 'A', A0, A1, A2
#c      write(*,*) 'B', B0, B1
#c      write(*,*) 'Lhmax Lvmax', Lhmax/4.45, Lvmax/4.45
#
#c      write(*,*) Ab1^2 - 4.0*Ab0*Ab2
#c      write(*,*) A2 *(xbulk-xhbend)^2 + A1*(xhtail-xhbend) + A0
#c      write(*,*) Ab2*xhbend^2 + Ab1*xhbend + Ab0
#
#c      write(*,*) 'R dR wfb', Rfuse, dRfuse, wfb
#c      write(*,*) 'xt', xwbox, xbulk, xtail
#c      write(*,*) 'xhbend', xhbend, xvbend
#c      pause
#
#c      Ishell = (pi + 2.0*thetafb + sin2t)*Rfuse^3*tshell
#c            + 0.66667*hfb^3*tfweb

#       write(*,*) Ishell, Ihshell, Ivshell

#----------------------------------------------------------------
#--- overall fuse weight and moment
      Wfuse = Wfix + Wapu + Wpadd + Wseat +
             Wshell + Wcone + Wwindow + Winsul + Wfloor+
             Whbend + Wvbend

      xWfuse = xWfix + xWapu + xWpadd + xWseat+
             xWshell + xWcone + xWwindow + xWinsul + xWfloor+
             xWhbend + xWvbend 

#----------------------------------------------------------------
#--- pressurized cabin volume
      cabVol = Afuse*(lshell + 0.67*lnose + 0.67*Rfuse)
      

return  tskin, tcone, tfweb, tfloor, xhbend, xvbend,
                       EIhshell,EIhbend,
                       EIvshell,EIvbend,
                       GJshell ,GJcone,
                       Wshell, Wcone, Wwindow, Winsul, Wfloor,
                       Whbend, Wvbend,
                       Wfuse,
                      xWfuse,
                       cabVol
end # fusew
