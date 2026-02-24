"""
      fusew!(fuse,Nland,Wpay,Weng, nftanks, 
      Waftfuel, Wftank, ltank, xftankaft, tank_placement,deltap,
      Whtail,Wvtail,rMh,rMv,Lhmax,Lvmax,
      bv, λv,nvtail,
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
      - ` λv::Float64`: Vertical tail taper ratio.
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
      bv, λv,nvtail,
      xhtail,xvtail,
      xwing,xwbox,cbox,
      xeng)

      layout = fuse.layout

      Eskin = fuse.material.E #parg[igEcap]
      Ebend = Eskin * fuse.ratio_young_mod_fuse_bending
      Gskin = Eskin * 0.5 / (1.0 + 0.3)

#--- cone material properties
#     (assumed same as skin, but could be different)

      rE = Ebend/Eskin

      ## Size primary shell
      A_fuse, fuse.insulation =  calculate_shell_geometry!(fuse, deltap)
      x_cabin, l_cabin = get_cabin_dimensions(layout,tank_placement, nftanks, ltank)
      fuse.window = size_windows(x_cabin, l_cabin, fuse.window_W_per_length, fuse.n_decks)

#--------------------------------------------------------------------
#--- various weight moments
      xWfix  = y_moment(fuse.fixed)
      xWapu  = y_moment(fuse.APU)
      xWseat = fuse.seat.W * x_cabin
      xWpadd = fuse.added_payload.W * x_cabin

#--------------------------------------------------------------------
#--- floor structural sizing
      P = (Wpay+fuse.seat.W) * Nland / fuse.n_decks #Total load is distributed across all decks
      floor_half_width = layout.bubble_center_y_offset + layout.radius
      xWfloor = size_floor!(fuse, P, layout.l_floor, floor_half_width, mid_span_support = false)

#--------------------------------------------------------------------
#--- size tailcone to withstand vertical-tail torsion Qv
      size_tailcone(fuse, nvtail, Lvmax, bv, λv)

      thetafb, hfb, sin2t, cost, web_length = web_geometry(layout.cross_section)

      xWcone = y_moment(fuse.cone.weight)

#--------------------------------------------------------------------
#--- torsional stiffnesses
      fuse.shell.GJ = Gskin*4.0*A_fuse^2 * fuse.skin.thickness / get_perimeter(layout.cross_section)
#--------------------------------------------------------------------
#--- lumped tail weight and location  
#      (Weng=0 if there are no tail-mounted engines)
      Wtail = Whtail + Wvtail + fuse.cone.weight.W + fuse.APU.W + Waftfuel + Wftank + Weng

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
          Nland*(Wpay+fuse.added_payload.W+fuse.shell.weight.W+fuse.window.W+fuse.insulation.W+fuse.floor.weight.W+fuse.seat.W)*
          0.5/layout.l_shell
      A1 = 1.0/(hfuse*sigMh)*
          (Nland*Wtail + rMh*Lhmax)
      A0 = -Ihshell/(rE*hfuse^2)
      Abar2 = A2
      Abar1 = 2.0*A2*xbulk + A1
      Abar0 = A2*xbulk^2 + A1*xtail + A0
      desc = max( 0.0 , Abar1^2 - 4.0*Abar0*Abar2 )
      bending_h_x = (Abar1 - sqrt(desc))*0.5/Abar2

      dxwing = xwing - xwbox
      xf = xwing + dxwing + 0.5*cbox
      xb = xwing - dxwing + 0.5*cbox

      Ahbendf = max(Abar2*xf^2 - Abar1*xf + Abar0, 0)
      Ahbendb = max(Abar2*xb^2 - Abar1*xb + Abar0, 0)

      Vhbendf = A2*((xbulk-xf)^3 - (xbulk-fuse.bendingmaterial_h.weight.x)^3)/3.0 +
               A1*((xtail-xf)^2 - (xtail-fuse.bendingmaterial_h.weight.x)^2)/2.0 +
               A0*(fuse.bendingmaterial_h.weight.x-xf)
      Vhbendb = A2*((xbulk-xb)^3 - (xbulk-fuse.bendingmaterial_h.weight.x)^3)/3.0 +
               A1*((xtail-xb)^2 - (xtail-fuse.bendingmaterial_h.weight.x)^2)/2.0 +
              A0*(fuse.bendingmaterial_h.weight.x-xb)
      Vhbendc = 0.5*(Ahbendf+Ahbendb)*cbox


      fuse.bendingmaterial_h.weight = Weight(W = fuse.bendingmaterial_h.ρ*gee*(Vhbendf + Vhbendb + Vhbendc), 
                                                x = bending_h_x)

      xWhbend = fuse.bendingmaterial_h.weight.W *      xwing

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
      bending_v_x = xvtail + B0/B1 # point where Avbend = 0 

      Avbendb = max(B1*(xtail-xb) + B0, 0)
      Vvbendb = max(B1*((xtail-xb)^2 - (xtail-fuse.bendingmaterial_v.weight.x)^2)/2.0+
               B0*(fuse.bendingmaterial_v.weight.x-xb), 0)
      Vvbendc = 0.5*Avbendb*cbox
      fuse.bendingmaterial_v.weight = Weight(W = fuse.bendingmaterial_h.ρ*gee*(Vvbendb + Vvbendc), 
      x = bending_v_x)
      # println("W, Vvb, Vvc Avbend = $Wvbend, $Vvbendb, $Vvbendc, $Avbendb")

      xWvbend = fuse.bendingmaterial_v.weight.W * (2.0*xwing + fuse.bendingmaterial_v.weight.x)/3.0


      fuse.shell.EIv = Eskin * Ivshell
      fuse.bendingmaterial_h.EIv  = Ebend * Avbendb * 2.0*widf^2
      fuse.bendingmaterial_v.EIv = fuse.bendingmaterial_h.EIv

#----------------------------------------------------------------
#--- overall fuse weight and moment
      fuse.weight = fuse.fixed.W + fuse.APU.W + fuse.added_payload.W + fuse.seat.W +
             fuse.shell.weight.W + fuse.cone.weight.W + fuse.window.W + fuse.insulation.W + fuse.floor.weight.W+
             fuse.bendingmaterial_h.weight.W + fuse.bendingmaterial_v.weight.W

      fuse.moment = xWfix + xWapu + xWpadd + xWseat+
      y_moment(fuse.shell.weight) + xWcone + y_moment(fuse.window) + y_moment(fuse.insulation) + xWfloor+
             xWhbend + xWvbend 

#----------------------------------------------------------------
#--- pressurized cabin volume

      if nftanks == 0
            cabVol = A_fuse*(layout.l_shell + 0.67*layout.l_nose + 0.67*layout.radius)
      else #If there is a fuel tank in the fuselage, the pressure vessel has a smaller air volume
            cabVol = A_fuse*(l_cabin + 0.67*layout.l_nose + 0.67*layout.radius)
      end

return  cabVol
end # fusew


"""
    size_floor(floor::StructuralMember, load, floor_depth, floor_half_width;
      mid_span_support::Bool)

"""
function size_floor!(fuse, load, floor_length, floor_half_width;
      mid_span_support::Bool)

      floor = fuse.floor
      layout = fuse.layout
      # P = (Wpay + fuse.seat.W) * Nland / fuse.n_decks #Total load is distributed across all decks

      # Assume pinned ends (statically indeterminate problem)
      if mid_span_support # floor with center support
            Smax = (5.0 / 16.0) * load
            Mmax = (9.0 / 256.0) * load * floor_half_width
      else
            # full-width floor
            Smax = 0.50 * load
            Mmax = 0.25 * load * floor_half_width
      end

      # Parabolic loading of the web gives us the 1.5 coeff
      Afweb = 1.5 * Smax / floor.material.τmax
      # σ = Mmax*y/I = Mmax*h/2/I
      # I ≈ Acap*h²/4
      # σ ≈ 2*Mmax/Acap/h
      Afcap = 2.0 * Mmax / (floor.material.σmax * layout.floor_depth)

      Vfloor = (Afcap + Afweb) * 2.0 * floor_half_width
      W_floor_beam = floor.material.ρ * gee * Vfloor
      W_floor_panels = 2 * floor_half_width * floor_length * fuse.floor_W_per_area

      floor_weight = fuse.n_decks * (W_floor_beam + W_floor_panels)
      xWfloor = fuse.floor.weight.W * 0.5 * (layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft)

      # average floor-beam cap thickness ("smeared" over entire floor)
      fuse.floor.thickness = 0.5*Afcap/layout.l_floor # half since two flanges in an I beam
      fuse.floor.weight.W = floor_weight

      return xWfloor
end  # function size_floor!

# size_floor(cs::SingleBubble) = size_floor( mid_span_support = false)

"""
    size_tailcone(fuse::Fuselage, n_vertical_tails, L_vmax, b_v, λv)

TBW
"""
function size_tailcone(fuse::Fuselage, n_vertical_tails, L_vmax, b_v, λv)

      layout = fuse.layout
      cone = fuse.cone
      #Calculate torsional moment from vertical tail
      Qv = n_vertical_tails * (L_vmax * b_v / 3.0) * (1.0 + 2.0 * λv) / (1.0 + λv)
  
      θ_web, h_web, sin2θ, _, web_length = web_geometry(layout.cross_section)
      n_webs = layout.n_webs
  
      # Get cone volume
      V_cone = (Qv / cone.τ) *
               (π + 2 * n_webs * θ_web) / (π + n_webs * (2θ_web + sin2θ)) *
               (layout.x_cone_end - layout.x_pressure_shell_aft) / layout.radius *
               (2.0 / (1.0 + layout.taper_tailcone))

      x_cone = 0.5 * (layout.x_pressure_shell_aft + layout.x_cone_end)
  
      W_frac_add =
          (fuse.weight_frac_stringers + fuse.weight_frac_frame + fuse.weight_frac_skin_addl)
  
      cone.weight = Weight(W = cone.ρ * gee * V_cone * (1 + W_frac_add), x = x_cone)
  
      A_fuse = area(layout.cross_section)
      perimeter = get_perimeter(layout.cross_section)
  
      cone.thickness = Qv / (2.0 * cone.τ * A_fuse)
      cone.GJ = cone.material.G * 4.0 * A_fuse^2 * cone.thickness / perimeter
  
  end  # function size_tailcone