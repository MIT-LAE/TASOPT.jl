"""
      fusew(gee, Nland, Wfix, Wpay, Wpadd, Wseat, Wapu, Weng, Waftfuel,
            fuse.weight_frac_string, fuse.weight_frac_frame, ffadd, deltap,
            Wpwindow, Wppinsul, Wppfloor, fuse.n_decks,
            Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
            bv, lambdav, nvtail,
            fuse.layout.radius, fuse.layout.bubble_lower_downward_shift, wfb, fuse.layout.n_webs, fuse.layout.tailcone_taper_ratio,
            fuse.layout.x_nose, fuse.layout.x_pressure_shell_fwd, fuse.layout.x_pressure_shell_aft, fuse.layout.x_cone_end,
            xhtail, xvtail,
            xwing, xwbox, cbox,
            xfix, xapu, xeng, xfuel,
            fuse.floor.thickness,
            fuse.skin.σ, fuse.bending_h.σ, fuse.skin.ρ, fuse.bending_h.ρ, 
            Eskin, Ebend, Gskin)

`fusew` sizes the fuselage and calculates the component weights and structural properties.
It takes inputs related to geometry, fixed weights, material properties, and more to compute the fuselage dimensions, weights, and other parameters.
       
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `Nland::Float64`: load factor for floor beam structural sizing.
      Fixed weights of various components:
      - `Wfix::Float64`: Fixed weight of the structure.
      - `Wpay::Float64`: Fixed weight of payload.
      - `Wpadd::Float64`: Fixed weight of additional equipment.
      - `Wseat::Float64`: Fixed weight of seats.
      - `Wapu::Float64`: Fixed weight of auxiliary power unit.
      - `Weng::Float64`: Fixed weight of engines.
      - `Waftfuel::Float64`: Fixed weight of aft fuel storage.

      Factors for stringers, frames, and additional structural components:
      - `fuse.weight_frac_string::Float64`: Factor for stringers.
      - `fuse.weight_frac_frame::Float64`: Factor for frames.
      - `ffadd::Float64`: Factor for additional structural components.
      
      Pressure differential:
      - `deltap::Float64`: Pressure differential.

      Weights of window, insulation, and floor:
      - `Wpwindow::Float64`: Weight of windows.
      - `Wppinsul::Float64`: Weight of insulation.
      - `Wppfloor::Float64`: Weight of floor per unit area.

      Vertical tail parameters:
      - `Whtail::Float64`: Weight of horizontal tail components.
      - `Wvtail::Float64`: Weight of vertical tail components.
      - `rMh::Float64`: Horizontal tail moment arm.
      - `rMv::Float64`: Vertical tail moment arm.
      - `Lhmax::Float64`: Maximum horizontal tail length.
      - `Lvmax::Float64`: Maximum vertical tail length.
      - `bv::Float64`: Vertical tail span.
      - `lambdav::Float64`: Vertical tail taper ratio.
      - `nvtail::Integer`: Number of vertical tail units.

      Fuselage parameters:
      - `fuse.layout.radius::Float64`: Fuselage radius.
      - `fuse.layout.bubble_lower_downward_shift::Float64`: Fuselage thickness.
      - `wfb::Float64`: Fuselage width.
      - `fuse.layout.n_webs::Integer`: Number of fuselage webs.
      - `fuse.layout.tailcone_taper_ratio::Float64`: Fuselage taper ratio.

      Geometric parameters and locations:
      - `fuse.layout.x_nose::Float64`: X location of the nose.
      - `fuse.layout.x_pressure_shell_fwd::Float64`: X location of the first shell point.
      - `fuse.layout.x_pressure_shell_aft::Float64`: X location of the second shell point.
      - `fuse.layout.x_cone_end::Float64`: X location of the cone end.
      - `xhtail::Float64`: X location of horizontal tail components.
      - `xvtail::Float64`: X location of vertical tail components.
      - `xwing::Float64`: X location of the wing.
      - `xwbox::Float64`: X location of the wing box.
      - `cbox::Float64`: Wing box width.
      - `xfix::Float64`: X location of fixed components.
      - `xapu::Float64`: X location of auxiliary power unit.
      - `xeng::Float64`: X location of engines.
      - `xfuel::Float64`: X location of fuel storage.
      - `fuse.floor.thickness::Float64`: Height of the floor.

      Material properties:
      - `fuse.skin.σ::Float64`: Skin material stress.
      - `fuse.bending_h.σ::Float64`: Bending material stress.
      - `fuse.skin.ρ::Float64`: Skin material density.
      - `fuse.bending_h.ρ::Float64`: Bending material density.
      - `Eskin::Float64`: Skin material Young's modulus.
      - `Ebend::Float64`: Bending material Young's modulus.
      - `Gskin::Float64`: Skin material shear modulus.

      **Outputs:**

      Thicknesses and locations:
      - `tskin::Float64`: Fuselage skin thickness.
      - `tcone::Float64`: Thickness of the tail cone.
      - `tfweb::Float64`: Thickness of fuselage webs.
      - `tfloor::Float64`: Floor beam thickness.
      - `xhbend::Float64`: X location of added material for horizontal-axis bending.
      - `xvbend::Float64`: X location of added material for vertical-axis bending.

      Bending and torsion inertias:
      - `EIhshell::Float64`: Bending inertia for horizontal shell.
      - `EIhbend::Float64`: Bending inertia for horizontal axis bending.
      - `EIvshell::Float64`: Bending inertia for vertical shell.
      - `EIvbend::Float64`: Bending inertia for vertical axis bending.
      - `GJshell::Float64`: Torsional stiffness for horizontal shell.
      - `GJcone::Float64`: Torsional stiffness for tail cone.

      Weights of components and total fuselage weight:
      - `Wshell::Float64`: Weight of fuselage shell components.
      - `Wcone::Float64`: Weight of tail cone.
      - `Wwindow::Float64`: Weight of windows.
      - `Winsul::Float64`: Weight of insulation.
      - `Wfloor::Float64`: Weight of floor.
      - `Whbend::Float64`: Weight of horizontal-axis bending material.
      - `Wvbend::Float64`: Weight of vertical-axis bending material.
      - `Wfuse::Float64`: Total weight of the fuselage.

      Moments
      - `xWfuse::Float64`: Moments.

      Pressurized cabin volume:
      - `cabVol::Float64`: Pressurized cabin volume.

See [here](@ref fuselage) or Section 2.2 of the [TASOPT Technical Description](@ref dreladocs).
"""
function fusew!(fuse,Nland,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
      ifwing, nftanks, 
      Waftfuel, Wftank, ltank, xftankaft, tank_placement,
      ffadd,deltap,
      Wpwindow,Wppinsul,Wppfloor, 
      Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
      bv,lambdav,nvtail,
      xhtail,xvtail,
      xwing,xwbox,cbox,
      xfix,xapu,xeng,xfuel)

      Eskin = fuse.material.E #parg[igEcap]
      Ebend = Eskin * fuse.fuse_shell_modulus_ratio
      Gskin = Eskin * 0.5 / (1.0 + 0.3)

#--- cone material properties
#     (assumed same as skin, but could be different)
      taucone = fuse.skin.σ
      rhocone = fuse.skin.ρ

#--- floor beam material properties 
#     (assumed same as stringers, but could be different)
      sigfloor = fuse.bending_h.σ
      taufloor = fuse.bending_h.σ
      rhofloor = fuse.bending_h.ρ

      rE = Ebend/Eskin

#--- effective nose length and pressure-vessel length
      lnose  = fuse.layout.x_pressure_shell_fwd - fuse.layout.x_nose
      lshell = fuse.layout.x_pressure_shell_aft - fuse.layout.x_pressure_shell_fwd
      lfloor = fuse.layout.x_pressure_shell_aft - fuse.layout.x_pressure_shell_fwd + 2.0*fuse.layout.radius

      xshell = 0.5*(fuse.layout.x_pressure_shell_fwd+fuse.layout.x_pressure_shell_aft)

#--- fuselage cross-section geometric parameters
      wfblim = max( min( fuse.layout.bubble_center_y_offset , fuse.layout.radius ) , 0.0 )
      thetafb = asin(wfblim/fuse.layout.radius) #θfb fuselage bubble subtended half-angle
      hfb = sqrt(fuse.layout.radius^2 - fuse.layout.bubble_center_y_offset^2)
      sin2t = 2.0*hfb*fuse.layout.bubble_center_y_offset/fuse.layout.radius^2 #sin(2θ) = 2sinθcosθ
      cost  = hfb/fuse.layout.radius

      perim = (2.0*pi + 4.0*thetafb)*fuse.layout.radius + 2.0*fuse.layout.bubble_lower_downward_shift

#--------------------------------------------------------------------
#--- fuselage skin and center web thicknesses to withstand pressure load
      fuse.skin.thickness =     deltap*fuse.layout.radius/fuse.skin.σ
      fuse.layout.thickness_webs = 2.0*deltap*fuse.layout.bubble_center_y_offset  /fuse.skin.σ

#--- cross-sectional areas
      Askin = (2.0*pi+4.0*fuse.layout.n_webs*thetafb)*fuse.layout.radius*fuse.skin.thickness + 2.0*fuse.layout.bubble_lower_downward_shift*fuse.skin.thickness
      Afweb = fuse.layout.n_webs*(2.0*hfb+fuse.layout.bubble_lower_downward_shift)*fuse.layout.thickness_webs
      Afuse = (pi + fuse.layout.n_webs*(2.0*thetafb + sin2t))*fuse.layout.radius^2 + 2.0*fuse.layout.radius*fuse.layout.bubble_lower_downward_shift        #####
#          + 2.0*(fuse.layout.radius+fuse.layout.n_webs*fuse.layout.bubble_center_y_offset)*fuse.layout.bubble_lower_downward_shift        #####

#--- nose + rear bulkhead surface areas
      Sbulk = (2.0*pi + 4.0*fuse.layout.n_webs*thetafb)*fuse.layout.radius^2
      Snose = (2.0*pi + 4.0*fuse.layout.n_webs*thetafb)*fuse.layout.radius^2* ( 0.333 + 0.667*(lnose/fuse.layout.radius)^1.6 )^0.625

#--- component volumes and volume moments
      Vcyl  = Askin*lshell
      Vnose = Snose*fuse.skin.thickness
      Vbulk = Sbulk*fuse.skin.thickness
      Vfweb = Afweb*lshell

      xVcyl  = Vcyl  * 0.5*(fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft)
      xVnose = Vnose * 0.5*(fuse.layout.x_nose   + fuse.layout.x_pressure_shell_fwd)
      xVbulk = Vbulk *     (fuse.layout.x_pressure_shell_aft + 0.5*fuse.layout.radius)
      xVfweb = Vfweb * 0.5*(fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft)

#--- weights and weight moments
      Wskin = fuse.skin.ρ*gee*(Vcyl + Vnose + Vbulk)
      Wfweb = fuse.skin.ρ*gee* Vfweb
      xWskin = fuse.skin.ρ*gee*(xVcyl + xVnose + xVbulk)
      xWfweb = fuse.skin.ρ*gee* xVfweb

      fuse.shell.weight = Wskin*(1.0+fuse.weight_frac_string+fuse.weight_frac_frame+ffadd) + Wfweb
      xWshell = xWskin*(1.0+fuse.weight_frac_string+fuse.weight_frac_frame+ffadd) + xWfweb

#--------------------------------------------------------------------
#--- window weight

      if nftanks == 1
            if tank_placement == "front" #If tank is at the front
                  xcabin = 0.5 * (fuse.layout.x_pressure_shell_fwd + ltank + 2.0*ft_to_m + fuse.layout.x_pressure_shell_aft)
            else #tank is at rear
                  xcabin = 0.5 * (fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft - (ltank + 2.0*ft_to_m))
            end
      else
            xcabin = 0.5 * (fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft) #two or zero tanks
      end
      lcabin = fuse.layout.x_pressure_shell_aft - fuse.layout.x_pressure_shell_fwd - nftanks * (ltank + 2.0*ft_to_m) #cabin length is smaller if there are fuel tanks

      fuse.window.weight = fuse.n_decks * Wpwindow * lcabin
      xWwindow = fuse.window.weight * xcabin
      
#--------------------------------------------------------------------
#--- insulation weight
      fuse.insulation.weight = Wppinsul*((1.1*pi+2.0*thetafb)*fuse.layout.radius*lshell + 0.55*(Snose+Sbulk))
      xWinsul = fuse.insulation.weight * 0.5*(fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft)

#--------------------------------------------------------------------
#--- various weight moments
      xWfix  = Wfix  * xfix
      xWapu  = Wapu  * xapu
      xWseat = Wseat * xcabin
      xWpadd = Wpadd * xcabin

#--------------------------------------------------------------------
#--- floor structural sizing
      P = (Wpay+Wseat) * Nland / fuse.n_decks #Total load is distributed across all decks
      wfloor1 = fuse.layout.bubble_center_y_offset + fuse.layout.radius

      if (fuse.layout.bubble_center_y_offset == 0.0) 
#---- full-width floor
       Smax = 0.50 * P
       Mmax = 0.25 * P*wfloor1
      else
#---- floor with center support
       Smax = (5.0/16.0 ) * P
       Mmax = (9.0/256.0) * P*wfloor1
      end

      Afweb = 1.5*Smax/ taufloor
      Afcap = 2.0*Mmax/(sigfloor*fuse.layout.floor_depth)

      Vfloor = (Afcap + Afweb) * 2.0*wfloor1
      fuse.floor.weight = fuse.n_decks * (rhofloor*gee*Vfloor + 2.0*wfloor1*lfloor*Wppfloor)
      xWfloor = fuse.floor.weight * 0.5*(fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft)

#--- average floor-beam cap thickness ("smeared" over entire floor)
      fuse.floor.thickness = 0.5*Afcap/lfloor
      
#--------------------------------------------------------------------
#--- size tailcone to withstand vertical-tail torsion Qv
      Qv = (nvtail*Lvmax*bv/3.0)*(1.0+2.0*lambdav)/(1.0+lambdav)
      Vcone = (Qv/taucone)* 
             (pi + fuse.layout.n_webs* 2.0*thetafb)/
            (pi + fuse.layout.n_webs*(2.0*thetafb+sin2t))*
            (fuse.layout.x_cone_end-fuse.layout.x_pressure_shell_aft)/fuse.layout.radius *
            2.0/(1.0+fuse.layout.tailcone_taper_ratio)

      fuse.cone.weight = rhocone*gee*Vcone*(1.0+fuse.weight_frac_string+fuse.weight_frac_frame+ffadd)
      xWcone = fuse.cone.weight * 0.5*(fuse.layout.x_pressure_shell_aft + fuse.layout.x_cone_end)

      fuse.cone.thickness = Qv / (2.0*taucone*Afuse)

#--------------------------------------------------------------------
#--- torsional stiffnesses
      fuse.shell.GJ = Gskin*4.0*Afuse^2 * fuse.skin.thickness / perim
      fuse.cone.GJ  = Gskin*4.0*Afuse^2 * fuse.cone.thickness / perim

#--------------------------------------------------------------------
#--- lumped tail weight and location  
#      (Weng=0 if there are no tail-mounted engines)
      Wtail = Whtail + Wvtail + fuse.cone.weight + Wapu + Waftfuel + Wftank + Weng

      xtail = (  xhtail*Whtail +
               xvtail*Wvtail +
               xWcone +
               xapu*Wapu + 
               xeng*Weng +
               xftankaft*(Waftfuel + Wftank)) / Wtail

#--------------------------------------------------------------------
#--- shell bending inertias
      tshell = fuse.skin.thickness*(1.0 + rE*fuse.weight_frac_string*fuse.skin.ρ/fuse.bending_h.ρ)

#--- bending stress to be matched
      sigMh = fuse.bending_h.σ - rE*0.5*deltap*fuse.layout.radius/tshell
      sigMv = fuse.bending_h.σ - rE*0.5*deltap*fuse.layout.radius/tshell

#--- rear bulkhead where tail structure attaches
      xbulk = fuse.layout.x_pressure_shell_aft

#----------------------------------------------------------------
#--- horizontal-axis bending moment added material
      Ihshell = ( (pi + fuse.layout.n_webs*(2.0*thetafb+sin2t))*fuse.layout.radius^2 +
                 8.0*fuse.layout.n_webs*cost*0.5*fuse.layout.bubble_lower_downward_shift*fuse.layout.radius +
                 (2.0*pi + 4.0*fuse.layout.n_webs*thetafb)*(0.5*fuse.layout.bubble_lower_downward_shift)^2)*fuse.layout.radius*tshell+
              0.66667*fuse.layout.n_webs*(hfb+0.5*fuse.layout.bubble_lower_downward_shift)^3*fuse.layout.thickness_webs

      hfuse = fuse.layout.radius + 0.5*fuse.layout.bubble_lower_downward_shift
      A2 = 1.0/(hfuse*sigMh)*
          Nland*(Wpay+Wpadd+fuse.shell.weight+fuse.window.weight+fuse.insulation.weight+fuse.floor.weight+Wseat)*
          0.5/lshell
      A1 = 1.0/(hfuse*sigMh)*
          (Nland*Wtail + rMh*Lhmax)
      A0 = -Ihshell/(rE*hfuse^2)
      Abar2 = A2
      Abar1 = 2.0*A2*xbulk + A1
      Abar0 = A2*xbulk^2 + A1*xtail + A0
      desc = max( 0.0 , Abar1^2 - 4.0*Abar0*Abar2 )
      fuse.bending_h.x = (Abar1 - sqrt(desc))*0.5/Abar2

      dxwing = xwing - xwbox
      xf = xwing + dxwing + 0.5*cbox
      xb = xwing - dxwing + 0.5*cbox

      Ahbendf = max(Abar2*xf^2 - Abar1*xf + Abar0, 0)
      Ahbendb = max(Abar2*xb^2 - Abar1*xb + Abar0, 0)

      Vhbendf = A2*((xbulk-xf)^3 - (xbulk-fuse.bending_h.x)^3)/3.0 +
               A1*((xtail-xf)^2 - (xtail-fuse.bending_h.x)^2)/2.0 +
               A0*(fuse.bending_h.x-xf)
      Vhbendb = A2*((xbulk-xb)^3 - (xbulk-fuse.bending_h.x)^3)/3.0 +
               A1*((xtail-xb)^2 - (xtail-fuse.bending_h.x)^2)/2.0 +
              A0*(fuse.bending_h.x-xb)
      Vhbendc = 0.5*(Ahbendf+Ahbendb)*cbox


      fuse.bending_h.weight = fuse.bending_h.ρ*gee*(Vhbendf + Vhbendb + Vhbendc)

      xWhbend = fuse.bending_h.weight *      xwing

      fuse.shell.EIh = Eskin * Ihshell
      fuse.bending_h.EIh  = Ebend * 0.5*(Ahbendf+Ahbendb) * 2.0*hfuse^2
      fuse.bending_v.EIh = fuse.bending_h.EIh

#----------------------------------------------------------------
#--- vertical-axis bending moment added material
      nk = Int(round(fuse.layout.n_webs/2.0+0.001))
      ik = mod(Int(round(fuse.layout.n_webs+0.001))+1,2)
      ksum = 0.
      for k = 1: nk
        ksum = ksum + float(2*k-ik)^2
      end
      Ivshell = ( (pi + fuse.layout.n_webs*(2.0*thetafb - sin2t))*fuse.layout.radius^2+
                 8.0*cost*fuse.layout.n_webs*fuse.layout.bubble_center_y_offset*fuse.layout.radius +
                 (2.0*pi + 4.0*thetafb)*(fuse.layout.n_webs*fuse.layout.bubble_center_y_offset)^2+
                 4.0*thetafb*fuse.layout.bubble_center_y_offset^2 * ksum )*fuse.layout.radius*tshell
                 # ^This simplifies to π×fuse.layout.radius³×tshell 
                 # which is the area moment of inertia for a thin walled tube

      widf = fuse.layout.radius + fuse.layout.n_webs*fuse.layout.bubble_center_y_offset
      B1 = 1.0/(widf*sigMv) * (rMv*Lvmax*nvtail)
      B0 = -Ivshell/(rE*widf^2)
      fuse.bending_v.x = xvtail + B0/B1 # point where Avbend = 0 

      Avbendb = max(B1*(xtail-xb) + B0, 0)
      Vvbendb = max(B1*((xtail-xb)^2 - (xtail-fuse.bending_v.x)^2)/2.0+
               B0*(fuse.bending_v.x-xb), 0)
      Vvbendc = 0.5*Avbendb*cbox
      fuse.bending_v.weight = fuse.bending_h.ρ*gee*(Vvbendb + Vvbendc)
      # println("W, Vvb, Vvc Avbend = $Wvbend, $Vvbendb, $Vvbendc, $Avbendb")

      xWvbend = fuse.bending_v.weight * (2.0*xwing + fuse.bending_v.x)/3.0


      fuse.shell.EIv = Eskin * Ivshell
      fuse.bending_h.EIv  = Ebend * Avbendb * 2.0*widf^2
      fuse.bending_v.EIv = fuse.bending_h.EIv

#----------------------------------------------------------------
#--- overall fuse weight and moment
      fuse.weight = Wfix + Wapu + Wpadd + Wseat +
             fuse.shell.weight + fuse.cone.weight + fuse.window.weight + fuse.insulation.weight + fuse.floor.weight+
             fuse.bending_h.weight + fuse.bending_v.weight

      fuse.moment = xWfix + xWapu + xWpadd + xWseat+
             xWshell + xWcone + xWwindow + xWinsul + xWfloor+
             xWhbend + xWvbend 

#----------------------------------------------------------------
#--- pressurized cabin volume

      if nftanks == 0
            cabVol = Afuse*(lshell + 0.67*lnose + 0.67*fuse.layout.radius)
      else #If there is a fuel tank in the fuselage, the pressure vessel has a smaller air volume
            cabVol = Afuse*(lcabin + 0.67*lnose + 0.67*fuse.layout.radius)
      end

return  cabVol
end # fusew