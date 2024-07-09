"""
      fusew!(fuse,Nland,Wpay,Weng, nftanks, 
      Waftfuel, Wftank, ltank, xftankaft, tank_placement,deltap,
      Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
      bv,lambdav,nvtail,
      xhtail,xvtail,
      xwing,xwbox,cbox,
      xeng)

`fusew` sizes the fuselage and calculates the component weights and structural properties.
It takes inputs related to geometry, fixed weights, material properties, and more to compute the fuselage dimensions, weights, and other parameters.
       
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `fuselage::struct`: TASOPT.fuselage structure containing layout, internal and external load information
      - `Nland::Float64`: load factor for floor beam structural sizing.
      Cabin and Payload Information:
      - `Wpay::Float64`: Fixed weight of payload.
      - `deltap::Float64`: Pressure differential.
      Engine Parameters:
      - `Weng::Float64`: Fixed weight of engines.
      - `xeng::Float64`: X location of engines.
      Fuel Tank Parameters
      - `nftanks::Int64`: Number of fuel tanks
      - `Waftfuel::Float64`: Fixed weight of aft fuel storage.
      - `Wftank::Float64`: Fixed weight of fuel tank.
      - `ltank::Float64`: Length of fuel tank.
      - `xtankaft::Float64`: X location of aft fuel storage.
      - `tank_placement::String`: Location of tank in fuselage
      Tail parameters:
      - `Whtail::Float64`: Weight of horizontal tail components.
      - `Wvtail::Float64`: Weight of vertical tail components.
      - `rMh::Float64`: Horizontal tail moment arm.
      - `rMv::Float64`: Vertical tail moment arm.
      - `Lhmax::Float64`: Maximum horizontal tail length.
      - `Lvmax::Float64`: Maximum vertical tail length.
      - `bv::Float64`: Vertical tail span.
      - `lambdav::Float64`: Vertical tail taper ratio.
      - `nvtail::Integer`: Number of vertical tail units.
      - `xhtail::Float64`: X location of horizontal tail components.
      - `xvtail::Float64`: X location of vertical tail components.
      Wing Parameters
      - `xwing::Float64`: X location of the wing.
      - `xwbox::Float64`: X location of the wing box.
      - `cbox::Float64`: Wing box width.

      **Outputs:**
      Pressurized cabin volume:
      - `cabVol::Float64`: Pressurized cabin volume.

See [here](@ref fuselage) or Section 2.2 of the [TASOPT Technical Description](@ref dreladocs).
"""
function fusew!(fuse,Nland,Wpay,Weng, nftanks, 
      Waftfuel, Wftank, ltank, xftankaft, tank_placement,deltap,
      Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
      bv,lambdav,nvtail,
      xhtail,xvtail,
      xwing,xwbox,cbox,
      xeng)

      layout = fuse.layout

      Eskin = fuse.material.E #parg[igEcap]
      Ebend = Eskin * fuse.ratio_young_mod_fuse_bending
      Gskin = Eskin * 0.5 / (1.0 + 0.3)

#--- cone material properties
#     (assumed same as skin, but could be different)
      taucone = fuse.skin.σ
      rhocone = fuse.skin.ρ

#--- floor beam material properties 
#     (assumed same as stringers, but could be different)
      sigfloor = fuse.bendingmaterial_h.σ
      taufloor = fuse.bendingmaterial_h.σ
      rhofloor = fuse.bendingmaterial_h.ρ

      rE = Ebend/Eskin

#--- effective nose length and pressure-vessel length
      lnose  = layout.x_pressure_shell_fwd - layout.x_nose
      lshell = layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd
      lfloor = layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd + 2.0*layout.radius

      xshell = 0.5*(layout.x_pressure_shell_fwd+layout.x_pressure_shell_aft)

#--- fuselage cross-section geometric parameters
      wfblim = max( min( layout.bubble_center_y_offset , layout.radius ) , 0.0 )
      thetafb = asin(wfblim/layout.radius) #θfb fuselage bubble subtended half-angle
      hfb = sqrt(layout.radius^2 - layout.bubble_center_y_offset^2)
      sin2t = 2.0*hfb*layout.bubble_center_y_offset/layout.radius^2 #sin(2θ) = 2sinθcosθ
      cost  = hfb/layout.radius

      perim = (2.0*pi + 4.0*thetafb)*layout.radius + 2.0*layout.bubble_lower_downward_shift

#--------------------------------------------------------------------
#--- fuselage skin and center web thicknesses to withstand pressure load
      fuse.skin.thickness =     deltap*layout.radius/fuse.skin.σ
      layout.thickness_webs = 2.0*deltap*layout.bubble_center_y_offset  /fuse.skin.σ

#--- cross-sectional areas
      Askin = (2.0*pi+4.0*layout.n_webs*thetafb)*layout.radius*fuse.skin.thickness + 2.0*layout.bubble_lower_downward_shift*fuse.skin.thickness
      Afweb = layout.n_webs*(2.0*hfb+layout.bubble_lower_downward_shift)*layout.thickness_webs
      Afuse = (pi + layout.n_webs*(2.0*thetafb + sin2t))*layout.radius^2 + 2.0*layout.radius*layout.bubble_lower_downward_shift        #####
#          + 2.0*(layout.radius+layout.n_webs*layout.bubble_center_y_offset)*layout.bubble_lower_downward_shift        #####

#--- nose + rear bulkhead surface areas
      Sbulk = (2.0*pi + 4.0*layout.n_webs*thetafb)*layout.radius^2
      Snose = (2.0*pi + 4.0*layout.n_webs*thetafb)*layout.radius^2* ( 0.333 + 0.667*(lnose/layout.radius)^1.6 )^0.625

#--- component volumes and volume moments
      Vcyl  = Askin*lshell
      Vnose = Snose*fuse.skin.thickness
      Vbulk = Sbulk*fuse.skin.thickness
      Vfweb = Afweb*lshell

      xVcyl  = Vcyl  * 0.5*(layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)
      xVnose = Vnose * 0.5*(layout.x_nose   + layout.x_pressure_shell_fwd)
      xVbulk = Vbulk *     (layout.x_pressure_shell_aft + 0.5*layout.radius)
      xVfweb = Vfweb * 0.5*(layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)

#--- weights and weight moments
      Wskin = fuse.skin.ρ*gee*(Vcyl + Vnose + Vbulk)
      Wfweb = fuse.skin.ρ*gee* Vfweb
      xWskin = fuse.skin.ρ*gee*(xVcyl + xVnose + xVbulk)
      xWfweb = fuse.skin.ρ*gee* xVfweb

      fuse.shell.weight = Wskin*(1.0+fuse.weight_frac_stringers+fuse.weight_frac_frame+fuse.weight_frac_skin_addl) + Wfweb
      xWshell = xWskin*(1.0+fuse.weight_frac_stringers+fuse.weight_frac_frame+fuse.weight_frac_skin_addl) + xWfweb

#--------------------------------------------------------------------
#--- window weight

      if nftanks == 1
            if tank_placement == "front" #If tank is at the front
                  xcabin = 0.5 * (layout.x_pressure_shell_fwd + ltank + 2.0*ft_to_m + layout.x_pressure_shell_aft)
            else #tank is at rear
                  xcabin = 0.5 * (layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft - (ltank + 2.0*ft_to_m))
            end
      else
            xcabin = 0.5 * (layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft) #two or zero tanks
      end
      lcabin = layout.x_pressure_shell_aft - layout.x_pressure_shell_fwd - nftanks * (ltank + 2.0*ft_to_m) #cabin length is smaller if there are fuel tanks

      fuse.window.weight = fuse.n_decks * fuse.window.W_per_length * lcabin
      xWwindow = fuse.window.weight * xcabin
      
#--------------------------------------------------------------------
#--- insulation weight
      fuse.insulation.weight = fuse.insulation.W_per_area*((1.1*pi+2.0*thetafb)*layout.radius*lshell + 0.55*(Snose+Sbulk))
      xWinsul = fuse.insulation.weight * 0.5*(layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)

#--------------------------------------------------------------------
#--- various weight moments
      xWfix  = y_moment(fuse.fixed)
      xWapu  = y_moment(fuse.APU)
      xWseat = fuse.seat.W * xcabin
      xWpadd = fuse.added_payload.W * xcabin

#--------------------------------------------------------------------
#--- floor structural sizing
      P = (Wpay+fuse.seat.W) * Nland / fuse.n_decks #Total load is distributed across all decks
      wfloor1 = layout.bubble_center_y_offset + layout.radius

      if (layout.bubble_center_y_offset == 0.0) 
#---- full-width floor
       Smax = 0.50 * P
       Mmax = 0.25 * P*wfloor1
      else
#---- floor with center support
       Smax = (5.0/16.0 ) * P
       Mmax = (9.0/256.0) * P*wfloor1
      end

      Afweb = 1.5*Smax/ taufloor
      Afcap = 2.0*Mmax/(sigfloor*layout.floor_depth)

      Vfloor = (Afcap + Afweb) * 2.0*wfloor1
      fuse.floor.weight = fuse.n_decks * (rhofloor*gee*Vfloor + 2.0*wfloor1*lfloor*fuse.floor.W_per_area)
      xWfloor = fuse.floor.weight * 0.5*(layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)

#--- average floor-beam cap thickness ("smeared" over entire floor)
      fuse.floor.thickness = 0.5*Afcap/lfloor
      
#--------------------------------------------------------------------
#--- size tailcone to withstand vertical-tail torsion Qv
      Qv = (nvtail*Lvmax*bv/3.0)*(1.0+2.0*lambdav)/(1.0+lambdav)
      Vcone = (Qv/taucone)* 
             (pi + layout.n_webs* 2.0*thetafb)/
            (pi + layout.n_webs*(2.0*thetafb+sin2t))*
            (layout.x_cone_end-layout.x_pressure_shell_aft)/layout.radius *
            2.0/(1.0+layout.taper_tailcone)

      fuse.cone.weight = rhocone*gee*Vcone*(1.0+fuse.weight_frac_stringers+fuse.weight_frac_frame+fuse.weight_frac_skin_addl)
      xWcone = fuse.cone.weight * 0.5*(layout.x_pressure_shell_aft + layout.x_cone_end)

      fuse.cone.thickness = Qv / (2.0*taucone*Afuse)

#--------------------------------------------------------------------
#--- torsional stiffnesses
      fuse.shell.GJ = Gskin*4.0*Afuse^2 * fuse.skin.thickness / perim
      fuse.cone.GJ  = Gskin*4.0*Afuse^2 * fuse.cone.thickness / perim

#--------------------------------------------------------------------
#--- lumped tail weight and location  
#      (Weng=0 if there are no tail-mounted engines)
      Wtail = Whtail + Wvtail + fuse.cone.weight + fuse.APU.W + Waftfuel + Wftank + Weng

      xtail = (  xhtail*Whtail +
               xvtail*Wvtail +
               xWcone +
               xWapu + 
               xeng*Weng +
               xftankaft*(Waftfuel + Wftank)) / Wtail

#--------------------------------------------------------------------
#--- shell bending inertias
      tshell = fuse.skin.thickness*(1.0 + rE*fuse.weight_frac_stringers*fuse.skin.ρ/fuse.bendingmaterial_h.ρ)

#--- bending stress to be matched
      sigMh = fuse.bendingmaterial_h.σ - rE*0.5*deltap*layout.radius/tshell
      sigMv = fuse.bendingmaterial_h.σ - rE*0.5*deltap*layout.radius/tshell

#--- rear bulkhead where tail structure attaches
      xbulk = layout.x_pressure_shell_aft

#----------------------------------------------------------------
#--- horizontal-axis bending moment added material
# # layout.cross_section.n_webs = 0
# A = ((pi + layout.n_webs*(2.0*thetafb+sin2t))*layout.radius^2 +
# (2.0*pi + 4.0*layout.n_webs*thetafb)*(0.5*layout.bubble_lower_downward_shift)^2)*layout.radius*tshell
# B = 8.0*layout.n_webs*cost*0.5*layout.bubble_lower_downward_shift*layout.radius *layout.radius*tshell

      Ihshell = ( (pi + layout.n_webs*(2.0*thetafb+sin2t))*layout.radius^2 +
                 8.0*layout.n_webs*cost*0.5*layout.bubble_lower_downward_shift*layout.radius +
                 (2.0*pi + 4.0*layout.n_webs*thetafb)*(0.5*layout.bubble_lower_downward_shift)^2)*layout.radius*tshell+
              0.66667*layout.n_webs*(hfb+0.5*layout.bubble_lower_downward_shift)^3*layout.thickness_webs
      # println("Ihshell = $Ihshell, $A, $B, $(A+B)")

      hfuse = layout.radius + 0.5*layout.bubble_lower_downward_shift
      A2 = 1.0/(hfuse*sigMh)*
          Nland*(Wpay+fuse.added_payload.W+fuse.shell.weight+fuse.window.weight+fuse.insulation.weight+fuse.floor.weight+fuse.seat.W)*
          0.5/lshell
      A1 = 1.0/(hfuse*sigMh)*
          (Nland*Wtail + rMh*Lhmax)
      A0 = -Ihshell/(rE*hfuse^2)
      Abar2 = A2
      Abar1 = 2.0*A2*xbulk + A1
      Abar0 = A2*xbulk^2 + A1*xtail + A0
      desc = max( 0.0 , Abar1^2 - 4.0*Abar0*Abar2 )
      fuse.bendingmaterial_h.x = (Abar1 - sqrt(desc))*0.5/Abar2

      dxwing = xwing - xwbox
      xf = xwing + dxwing + 0.5*cbox
      xb = xwing - dxwing + 0.5*cbox

      Ahbendf = max(Abar2*xf^2 - Abar1*xf + Abar0, 0)
      Ahbendb = max(Abar2*xb^2 - Abar1*xb + Abar0, 0)

      Vhbendf = A2*((xbulk-xf)^3 - (xbulk-fuse.bendingmaterial_h.x)^3)/3.0 +
               A1*((xtail-xf)^2 - (xtail-fuse.bendingmaterial_h.x)^2)/2.0 +
               A0*(fuse.bendingmaterial_h.x-xf)
      Vhbendb = A2*((xbulk-xb)^3 - (xbulk-fuse.bendingmaterial_h.x)^3)/3.0 +
               A1*((xtail-xb)^2 - (xtail-fuse.bendingmaterial_h.x)^2)/2.0 +
              A0*(fuse.bendingmaterial_h.x-xb)
      Vhbendc = 0.5*(Ahbendf+Ahbendb)*cbox


      fuse.bendingmaterial_h.weight = fuse.bendingmaterial_h.ρ*gee*(Vhbendf + Vhbendb + Vhbendc)

      xWhbend = fuse.bendingmaterial_h.weight *      xwing

      fuse.shell.EIh = Eskin * Ihshell
      fuse.bendingmaterial_h.EIh  = Ebend * 0.5*(Ahbendf+Ahbendb) * 2.0*hfuse^2
      fuse.bendingmaterial_v.EIh = fuse.bendingmaterial_h.EIh

#----------------------------------------------------------------
#--- vertical-axis bending moment added material
      nk = Int(round(layout.n_webs/2.0+0.001))
      ik = mod(Int(round(layout.n_webs+0.001))+1,2)
      ksum = 0.
      for k = 1: nk
        ksum = ksum + float(2*k-ik)^2
      end
      Ivshell = ( (pi + layout.n_webs*(2.0*thetafb - sin2t))*layout.radius^2+
                 8.0*cost*layout.n_webs*layout.bubble_center_y_offset*layout.radius +
                 (2.0*pi + 4.0*thetafb)*(layout.n_webs*layout.bubble_center_y_offset)^2+
                 4.0*thetafb*layout.bubble_center_y_offset^2 * ksum )*layout.radius*tshell
                 # ^This simplifies to π×layout.radius³×tshell 
                 # which is the area moment of inertia for a thin walled tube

      widf = layout.radius + layout.n_webs*layout.bubble_center_y_offset
      B1 = 1.0/(widf*sigMv) * (rMv*Lvmax*nvtail)
      B0 = -Ivshell/(rE*widf^2)
      fuse.bendingmaterial_v.x = xvtail + B0/B1 # point where Avbend = 0 

      Avbendb = max(B1*(xtail-xb) + B0, 0)
      Vvbendb = max(B1*((xtail-xb)^2 - (xtail-fuse.bendingmaterial_v.x)^2)/2.0+
               B0*(fuse.bendingmaterial_v.x-xb), 0)
      Vvbendc = 0.5*Avbendb*cbox
      fuse.bendingmaterial_v.weight = fuse.bendingmaterial_h.ρ*gee*(Vvbendb + Vvbendc)
      # println("W, Vvb, Vvc Avbend = $Wvbend, $Vvbendb, $Vvbendc, $Avbendb")

      xWvbend = fuse.bendingmaterial_v.weight * (2.0*xwing + fuse.bendingmaterial_v.x)/3.0


      fuse.shell.EIv = Eskin * Ivshell
      fuse.bendingmaterial_h.EIv  = Ebend * Avbendb * 2.0*widf^2
      fuse.bendingmaterial_v.EIv = fuse.bendingmaterial_h.EIv

#----------------------------------------------------------------
#--- overall fuse weight and moment
      fuse.weight = fuse.fixed.W + fuse.APU.W + fuse.added_payload.W + fuse.seat.W +
             fuse.shell.weight + fuse.cone.weight + fuse.window.weight + fuse.insulation.weight + fuse.floor.weight+
             fuse.bendingmaterial_h.weight + fuse.bendingmaterial_v.weight

      fuse.moment = xWfix + xWapu + xWpadd + xWseat+
             xWshell + xWcone + xWwindow + xWinsul + xWfloor+
             xWhbend + xWvbend 

#----------------------------------------------------------------
#--- pressurized cabin volume

      if nftanks == 0
            cabVol = Afuse*(lshell + 0.67*lnose + 0.67*layout.radius)
      else #If there is a fuel tank in the fuselage, the pressure vessel has a smaller air volume
            cabVol = Afuse*(lcabin + 0.67*lnose + 0.67*layout.radius)
      end

return  cabVol
end # fusew