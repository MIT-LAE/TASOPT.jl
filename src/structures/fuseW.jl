"""
      fusew(gee, fuselage.Nland, Wfix, Wpay, Wpadd, Wseat, Wapu, Weng, Waftfuel,
            fuselage.fstring, fuselage.fframe, fuselage.ffadd, deltap,
            Wpwindow, Wppinsul, Wppfloor, fuselage.n_decks,
            Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
            bv, lambdav, nvtail,
            fuselage.layout.fuse_radius, fuselage.layout.bubble_lower_downward_shift, wfb, fuselage.layout.n_webs, fuselage.layout.tailcone_taper_ratio,
            fuselage.layout.x_nose, fuselage.layout.x_pressure_shell_fwd, fuselage.layout.x_pressure_shell_aft, fuselage.layout.x_cone_end,
            xhtail, xvtail,
            xwing, xwbox, cbox,
            xfix, xapu, xeng, xfuel,
            fuselage.floor.thickness,
            fuselage.skin.σ, fuselage.bending_h.σ, fuselage.skin.ρ, fuselage.bending_h.ρ, 
            Eskin, Ebend, Gskin)

`fusew` sizes the fuselage and calculates the component weights and structural properties.
It takes inputs related to geometry, fixed weights, material properties, and more to compute the fuselage dimensions, weights, and other parameters.
       
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `fuselage.Nland::Float64`: load factor for floor beam structural sizing.
      Fixed weights of various components:
      - `Wfix::Float64`: Fixed weight of the structure.
      - `Wpay::Float64`: Fixed weight of payload.
      - `Wpadd::Float64`: Fixed weight of additional equipment.
      - `Wseat::Float64`: Fixed weight of seats.
      - `Wapu::Float64`: Fixed weight of auxiliary power unit.
      - `Weng::Float64`: Fixed weight of engines.
      - `Waftfuel::Float64`: Fixed weight of aft fuel storage.

      Factors for stringers, frames, and additional structural components:
      - `fuselage.fstring::Float64`: Factor for stringers.
      - `fuselage.fframe::Float64`: Factor for frames.
      - `fuselage.ffadd::Float64`: Factor for additional structural components.
      
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
      - `fuselage.layout.fuse_radius::Float64`: Fuselage radius.
      - `fuselage.layout.bubble_lower_downward_shift::Float64`: Fuselage thickness.
      - `wfb::Float64`: Fuselage width.
      - `fuselage.layout.n_webs::Integer`: Number of fuselage webs.
      - `fuselage.layout.tailcone_taper_ratio::Float64`: Fuselage taper ratio.

      Geometric parameters and locations:
      - `fuselage.layout.x_nose::Float64`: X location of the nose.
      - `fuselage.layout.x_pressure_shell_fwd::Float64`: X location of the first shell point.
      - `fuselage.layout.x_pressure_shell_aft::Float64`: X location of the second shell point.
      - `fuselage.layout.x_cone_end::Float64`: X location of the cone end.
      - `xhtail::Float64`: X location of horizontal tail components.
      - `xvtail::Float64`: X location of vertical tail components.
      - `xwing::Float64`: X location of the wing.
      - `xwbox::Float64`: X location of the wing box.
      - `cbox::Float64`: Wing box width.
      - `xfix::Float64`: X location of fixed components.
      - `xapu::Float64`: X location of auxiliary power unit.
      - `xeng::Float64`: X location of engines.
      - `xfuel::Float64`: X location of fuel storage.
      - `fuselage.floor.thickness::Float64`: Height of the floor.

      Material properties:
      - `fuselage.skin.σ::Float64`: Skin material stress.
      - `fuselage.bending_h.σ::Float64`: Bending material stress.
      - `fuselage.skin.ρ::Float64`: Skin material density.
      - `fuselage.bending_h.ρ::Float64`: Bending material density.
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
function fusew(fuselage,Wfix,Wpay,Wpadd,Wseat,Wapu,Weng,
      
      Waftfuel, Wftank, ltank, xftankaft, tank_placement,
      deltap,
      Wpwindow,Wppinsul,Wppfloor, 
      Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
      bv,lambdav,nvtail,xhtail,xvtail,
      xwing,xwbox,cbox,
      xfix,xapu,xeng,xfuel,
      
      
      
     
      Eskin,Ebend,Gskin)

      

      # fusew(Wfix, Wpaymax, Wpadd, Wseat, Wapu, Wengtail, 
            
#             Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
#              Δp,
#             Wpwindow, Wppinsul, Wppfloor, 
#             Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
#             bv, λv, nvtail,
            
            
#             xhtail, xvtail,
#             xwing, xwbox, cbox,
#             xfix, xapu, xeng, xfuel,
#             )

#--- cone material properties
#     (assumed same as skin, but could be different)
#FIX WHY IS THIS SAME AS SKIN
      taucone = fuselage.skin.σ
      rhocone = fuselage.skin.ρ

#--- floor beam material properties 
#     (assumed same as stringers, but could be different)
      sigfloor = fuselage.bending_h.σ
      taufloor = fuselage.bending_h.σ
      rhofloor = fuselage.bending_h.ρ

      rE = Ebend/Eskin

#--- effective nose length and pressure-vessel length
      lnose  = fuselage.layout.x_pressure_shell_fwd - fuselage.layout.x_nose
      lshell = fuselage.layout.x_pressure_shell_aft - fuselage.layout.x_pressure_shell_fwd
      lfloor = fuselage.layout.x_pressure_shell_aft - fuselage.layout.x_pressure_shell_fwd + 2.0*fuselage.layout.fuse_radius

      xshell = 0.5*(fuselage.layout.x_pressure_shell_fwd+fuselage.layout.x_pressure_shell_aft)

#--- fuselage cross-section geometric parameters
      wfblim = max( min( fuselage.layout.bubble_center_y_offset , fuselage.layout.fuse_radius ) , 0.0 )
      thetafb = asin(wfblim/fuselage.layout.fuse_radius) #θfb fuselage bubble subtended half-angle
      hfb = sqrt(fuselage.layout.fuse_radius^2 - fuselage.layout.bubble_center_y_offset^2)
      sin2t = 2.0*hfb*fuselage.layout.bubble_center_y_offset/fuselage.layout.fuse_radius^2 #sin(2θ) = 2sinθcosθ
      cost  = hfb/fuselage.layout.fuse_radius

      perim = (2.0*pi + 4.0*thetafb)*fuselage.layout.fuse_radius + 2.0*fuselage.layout.bubble_lower_downward_shift

#--------------------------------------------------------------------
#--- fuselage skin and center web thicknesses to withstand pressure load
      tskin =     deltap*fuselage.layout.fuse_radius/fuselage.skin.σ
      tfweb = 2.0*deltap*fuselage.layout.bubble_center_y_offset  /fuselage.skin.σ

#--- cross-sectional areas
      Askin = (2.0*pi+4.0*fuselage.layout.n_webs*thetafb)*fuselage.layout.fuse_radius*tskin + 2.0*fuselage.layout.bubble_lower_downward_shift*tskin
      Afweb = fuselage.layout.n_webs*(2.0*hfb+fuselage.layout.bubble_lower_downward_shift)*tfweb
      Afuse = (pi + fuselage.layout.n_webs*(2.0*thetafb + sin2t))*fuselage.layout.fuse_radius^2 + 2.0*fuselage.layout.fuse_radius*fuselage.layout.bubble_lower_downward_shift        #####
#          + 2.0*(fuselage.layout.fuse_radius+fuselage.layout.n_webs*fuselage.layout.bubble_center_y_offset)*fuselage.layout.bubble_lower_downward_shift        #####

#--- nose + rear bulkhead surface areas
      Sbulk = (2.0*pi + 4.0*fuselage.layout.n_webs*thetafb)*fuselage.layout.fuse_radius^2
      Snose = (2.0*pi + 4.0*fuselage.layout.n_webs*thetafb)*fuselage.layout.fuse_radius^2* ( 0.333 + 0.667*(lnose/fuselage.layout.fuse_radius)^1.6 )^0.625

#--- component volumes and volume moments
      Vcyl  = Askin*lshell
      Vnose = Snose*tskin
      Vbulk = Sbulk*tskin
      Vfweb = Afweb*lshell

      xVcyl  = Vcyl  * 0.5*(fuselage.layout.x_pressure_shell_fwd + fuselage.layout.x_pressure_shell_aft)
      xVnose = Vnose * 0.5*(fuselage.layout.x_nose   + fuselage.layout.x_pressure_shell_fwd)
      xVbulk = Vbulk *     (fuselage.layout.x_pressure_shell_aft + 0.5*fuselage.layout.fuse_radius)
      xVfweb = Vfweb * 0.5*(fuselage.layout.x_pressure_shell_fwd + fuselage.layout.x_pressure_shell_aft)

#--- weights and weight moments
      Wskin = fuselage.skin.ρ*gee*(Vcyl + Vnose + Vbulk)
      Wfweb = fuselage.skin.ρ*gee* Vfweb
      xWskin = fuselage.skin.ρ*gee*(xVcyl + xVnose + xVbulk)
      xWfweb = fuselage.skin.ρ*gee* xVfweb

      Wshell = Wskin*(1.0+fuselage.fstring+fuselage.fframe+fuselage.ffadd) + Wfweb
      xWshell = xWskin*(1.0+fuselage.fstring+fuselage.fframe+fuselage.ffadd) + xWfweb

#--------------------------------------------------------------------
#--- window weight

      if fuselage.nftanks == 1
            if tank_placement == "front" #If tank is at the front
                  xcabin = 0.5 * (fuselage.layout.x_pressure_shell_fwd + ltank + 2.0*ft_to_m + fuselage.layout.x_pressure_shell_aft)
            else #tank is at rear
                  xcabin = 0.5 * (fuselage.layout.x_pressure_shell_fwd + fuselage.layout.x_pressure_shell_aft - (ltank + 2.0*ft_to_m))
            end
      else
            xcabin = 0.5 * (fuselage.layout.x_pressure_shell_fwd + fuselage.layout.x_pressure_shell_aft) #two or zero tanks
      end
      lcabin = fuselage.layout.x_pressure_shell_aft - fuselage.layout.x_pressure_shell_fwd - fuselage.nftanks * (ltank + 2.0*ft_to_m) #cabin length is smaller if there are fuel tanks

      Wwindow = fuselage.n_decks * Wpwindow * lcabin
      xWwindow = Wwindow * xcabin
      
#--------------------------------------------------------------------
#--- insulation weight
      Winsul = Wppinsul*((1.1*pi+2.0*thetafb)*fuselage.layout.fuse_radius*lshell + 0.55*(Snose+Sbulk))
      xWinsul = Winsul * 0.5*(fuselage.layout.x_pressure_shell_fwd + fuselage.layout.x_pressure_shell_aft)

#--------------------------------------------------------------------
#--- various weight moments
      xWfix  = Wfix  * xfix
      xWapu  = Wapu  * xapu
      xWseat = Wseat * xcabin
      xWpadd = Wpadd * xcabin

#--------------------------------------------------------------------
#--- floor structural sizing
      P = (Wpay+Wseat) * fuselage.Nland / fuselage.n_decks #Total load is distributed across all decks
      wfloor1 = fuselage.layout.bubble_center_y_offset + fuselage.layout.fuse_radius

      if (fuselage.layout.bubble_center_y_offset == 0.0) 
#---- full-width floor
       Smax = 0.50 * P
       Mmax = 0.25 * P*wfloor1
      else
#---- floor with center support
       Smax = (5.0/16.0 ) * P
       Mmax = (9.0/256.0) * P*wfloor1
      end

      Afweb = 1.5*Smax/ taufloor
      Afcap = 2.0*Mmax/(sigfloor*fuselage.floor.thickness)

      Vfloor = (Afcap + Afweb) * 2.0*wfloor1
      Wfloor = fuselage.n_decks * (rhofloor*gee*Vfloor + 2.0*wfloor1*lfloor*Wppfloor)
      xWfloor = Wfloor * 0.5*(fuselage.layout.x_pressure_shell_fwd + fuselage.layout.x_pressure_shell_aft)

#--- average floor-beam cap thickness ("smeared" over entire floor)
      tfloor = 0.5*Afcap/lfloor
      
#--------------------------------------------------------------------
#--- size tailcone to withstand vertical-tail torsion Qv
      Qv = (nvtail*Lvmax*bv/3.0)*(1.0+2.0*lambdav)/(1.0+lambdav)
      Vcone = (Qv/taucone)* 
             (pi + fuselage.layout.n_webs* 2.0*thetafb)/
            (pi + fuselage.layout.n_webs*(2.0*thetafb+sin2t))*
            (fuselage.layout.x_cone_end-fuselage.layout.x_pressure_shell_aft)/fuselage.layout.fuse_radius *
            2.0/(1.0+fuselage.layout.tailcone_taper_ratio)

      Wcone = rhocone*gee*Vcone*(1.0+fuselage.fstring+fuselage.fframe+fuselage.ffadd)
      xWcone = Wcone * 0.5*(fuselage.layout.x_pressure_shell_aft + fuselage.layout.x_cone_end)

      tcone = Qv / (2.0*taucone*Afuse)

#--------------------------------------------------------------------
#--- torsional stiffnesses
      GJshell = Gskin*4.0*Afuse^2 * tskin / perim
      GJcone  = Gskin*4.0*Afuse^2 * tcone / perim

#--------------------------------------------------------------------
#--- lumped tail weight and location  
#      (Weng=0 if there are no tail-mounted engines)
      Wtail = Whtail + Wvtail + Wcone + Wapu + Waftfuel + Wftank + Weng

      xtail = (  xhtail*Whtail +
               xvtail*Wvtail +
               xWcone +
               xapu*Wapu + 
               xeng*Weng +
               xftankaft*(Waftfuel + Wftank)) / Wtail

#--------------------------------------------------------------------
#--- shell bending inertias
      tshell = tskin*(1.0 + rE*fuselage.fstring*fuselage.skin.ρ/fuselage.bending_h.ρ)

#--- bending stress to be matched
      sigMh = fuselage.bending_h.σ - rE*0.5*deltap*fuselage.layout.fuse_radius/tshell
      sigMv = fuselage.bending_h.σ - rE*0.5*deltap*fuselage.layout.fuse_radius/tshell

#--- rear bulkhead where tail structure attaches
      xbulk = fuselage.layout.x_pressure_shell_aft

#----------------------------------------------------------------
#--- horizontal-axis bending moment added material
      Ihshell = ( (pi + fuselage.layout.n_webs*(2.0*thetafb+sin2t))*fuselage.layout.fuse_radius^2 +
                 8.0*fuselage.layout.n_webs*cost*0.5*fuselage.layout.bubble_lower_downward_shift*fuselage.layout.fuse_radius +
                 (2.0*pi + 4.0*fuselage.layout.n_webs*thetafb)*(0.5*fuselage.layout.bubble_lower_downward_shift)^2)*fuselage.layout.fuse_radius*tshell+
              0.66667*fuselage.layout.n_webs*(hfb+0.5*fuselage.layout.bubble_lower_downward_shift)^3*tfweb

      hfuse = fuselage.layout.fuse_radius + 0.5*fuselage.layout.bubble_lower_downward_shift
      A2 = 1.0/(hfuse*sigMh)*
          fuselage.Nland*(Wpay+Wpadd+Wshell+Wwindow+Winsul+Wfloor+Wseat)*
          0.5/lshell
      A1 = 1.0/(hfuse*sigMh)*
          (fuselage.Nland*Wtail + rMh*Lhmax)
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


      Whbend = fuselage.bending_h.ρ*gee*(Vhbendf + Vhbendb + Vhbendc)

      xWhbend = Whbend *      xwing

      EIhshell = Eskin * Ihshell
      EIhbend  = Ebend * 0.5*(Ahbendf+Ahbendb) * 2.0*hfuse^2

#----------------------------------------------------------------
#--- vertical-axis bending moment added material
      nk = Int(round(fuselage.layout.n_webs/2.0+0.001))
      ik = mod(Int(round(fuselage.layout.n_webs+0.001))+1,2)
      ksum = 0.
      for k = 1: nk
        ksum = ksum + float(2*k-ik)^2
      end
      Ivshell = ( (pi + fuselage.layout.n_webs*(2.0*thetafb - sin2t))*fuselage.layout.fuse_radius^2+
                 8.0*cost*fuselage.layout.n_webs*fuselage.layout.bubble_center_y_offset*fuselage.layout.fuse_radius +
                 (2.0*pi + 4.0*thetafb)*(fuselage.layout.n_webs*fuselage.layout.bubble_center_y_offset)^2+
                 4.0*thetafb*fuselage.layout.bubble_center_y_offset^2 * ksum )*fuselage.layout.fuse_radius*tshell
                 # ^This simplifies to π×fuselage.layout.fuse_radius³×tshell 
                 # which is the area moment of inertia for a thin walled tube

      widf = fuselage.layout.fuse_radius + fuselage.layout.n_webs*fuselage.layout.bubble_center_y_offset
      B1 = 1.0/(widf*sigMv) * (rMv*Lvmax*nvtail)
      B0 = -Ivshell/(rE*widf^2)
      xvbend = xvtail + B0/B1 # point where Avbend = 0 

      Avbendb = max(B1*(xtail-xb) + B0, 0)
      Vvbendb = max(B1*((xtail-xb)^2 - (xtail-xvbend)^2)/2.0+
               B0*(xvbend-xb), 0)
      Vvbendc = 0.5*Avbendb*cbox
      Wvbend = fuselage.bending_h.ρ*gee*(Vvbendb + Vvbendc)
      # println("W, Vvb, Vvc Avbend = $Wvbend, $Vvbendb, $Vvbendc, $Avbendb")

      xWvbend = Wvbend * (2.0*xwing + xvbend)/3.0


      EIvshell = Eskin * Ivshell
      EIvbend  = Ebend * Avbendb * 2.0*widf^2

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

      if fuselage.nftanks == 0
            cabVol = Afuse*(lshell + 0.67*lnose + 0.67*fuselage.layout.fuse_radius)
      else #If there is a fuel tank in the fuselage, the pressure vessel has a smaller air volume
            cabVol = Afuse*(lcabin + 0.67*lnose + 0.67*fuselage.layout.fuse_radius)
      end


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