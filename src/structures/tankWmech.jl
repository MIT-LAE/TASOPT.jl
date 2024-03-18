"""
      size_inner_tank(fuse_tank, t_cond::Vector{Float64}, ρfuel::Float64,
                        Rfuse::Float64, dRfuse::Float64, wfb::Float64, nfweb::Float64,
                        Wfuel::Float64)

`size_inner_tank` calculates the weight of the cryogenic fuel tank for a LH-fueled aircraft.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `fuse_tank::Struct`: structure with tank parameters.
      - `t_cond::Float64`: Vector with tank isulation layer thickness. Provided separately from fuse_tank as it changes during 
      non-linear solve process.

      **Outputs:**
      - `Wtank_total::Float64`: Total tank weight including fuel (N).
      - `l_cyl::Float64`: Length of the cylindrical portion (m).
      - `tskin::Float64`: Thickness of the tank's skin (m).
      - `Rtank_outer::Float64`: Outer radius of the tank (m).
      - `Vfuel::Float64`: Volume of fuel (m^3).
      - `Wtank::Float64`: Weight of the empty tank (N).
      - `Wfuel_tot::Float64`: Total weight of fuel (N).
      - `Winsul_sum::Float64`: Sum of insulation weight (N).
      - `t_head::Float64`: Thickness of the tank's head (m).
      - `Whead::Float64`: Weight of the tank's head (N).
      - `Wcyl::Float64`: Weight of the tank's cylinder (N).
      - `Winsul::Float64`: Weight of insulation (N).
      - `Shead_insul::Float64`: Insulated surface area of the head (m^2).
      - `l_tank::Float64`: Total length of the tank (m).

NOTE: Al alloy 2219 has been recommended as tank material (from H2 tank paper in OneNote)

See [here](@ref fueltanks).
"""

function size_inner_tank(fuse_tank, t_cond::Vector{Float64})

      #Unpack parameters in fuse_tank
      Rfuse = fuse_tank.Rfuse
      dRfuse = fuse_tank.dRfuse
      wfb = fuse_tank.wfb
      nfweb = fuse_tank.nfweb

      Wfuel = fuse_tank.Wfuelintank

      ρfuel = fuse_tank.rhofuel
      ftankstiff = fuse_tank.ftankstiff
      ftankadd = fuse_tank.ftankadd
      Δp = fuse_tank.ptank
      sigskin = fuse_tank.UTSinner
      material_insul = fuse_tank.material_insul
      rhoskin = fuse_tank.rhoinner
      clearance_fuse = fuse_tank.clearance_fuse
      AR = fuse_tank.ARtank
      ullage_frac = fuse_tank.ullage_frac
      weld_eff = fuse_tank.ew

# Total thickness:
      thickness_insul = sum(t_cond)

      sa = sigskin / 4 #Maximum allowable stress is 1/4 Ultimate tensile strength (Barron 1985, p. 359)
      
      Rtank_outer = Rfuse - thickness_insul - clearance_fuse

      tskin = Δp * (2 * Rtank_outer) / (2 * sa * weld_eff + 0.8 * Δp) #(7.1) in Barron (1985)

      Rtank = Rtank_outer - tskin
      #tfweb = 2.0 * Δp * wfb  / ew
      Lhead = Rtank / AR       # eg. for a 2:1 ellipsoid majorax/minorax = 2/1 ⟹ R/Lhead = 2/1 
      
      K = (1/6) * (AR^2 + 2) # Aspect ratio of 2:1 for the head (# Barron pg 359) 
      t_head = Δp* (2*Rtank_outer) * K/ (2 * sa * weld_eff + 2 * Δp * (K - 0.1)) #(7.2) in Barron (1985)

#--- Calculate length of cylindrical portion
      Wfuel_tot = Wfuel #Wfuel already includes the amount that boils off
      Vfuel = Wfuel_tot / (gee * ρfuel)
      Vinternal = (1 + ullage_frac)*Vfuel  # required interal volume
      V_ellipsoid = 2π * (Rtank^3 / AR) / 3  # Half the vol of std ellipsoid = 1/2×(4π/3 ×(abc)) where a,b,c are the semi-axes length. Here a = R/AR, b=c=R
                                       # Also see: https://neutrium.net/equipment/volume-and-wetted-area-of-partially-filled-horizontal-vessels/
      V_cylinder = Vinternal - 2*V_ellipsoid
      l_cyl = V_cylinder / (π * (Rtank^2)) #required length of cylindrical portion

#--- tank cross-section geometric parameters
      wfblim = max( min( wfb , Rtank) , 0.0 )
      thetafb = asin(wfblim / Rtank)

#--- areas
      Scyl = (2.0*π+4.0*nfweb*thetafb)*Rtank*l_cyl + 2.0*dRfuse*l_cyl # Surface area of cylindrical part
      #Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      #Atank = (π + nfweb*(2.0*thetafb + sin2t))*Rtank^2 + 2.0*Rtank*dRfuse #+ 2.0*(Rtank+nfweb*wfb)*dRfuse
      Shead = (2.0*π + 4.0*nfweb*thetafb)*Rtank^2* ( 0.333 + 0.667*(Lhead/Rtank)^1.6 )^0.625 # This form is better for insul thickness 
                                                                                          # but just as a note to reader this comes from  semi- oblate spheroid surf area is ≈ 2π×R²[1/3 + 2/3×(1/AR)^1.6]^(1/1.6)
#--- component volumes
      Vcyl  = Scyl*tskin    # volume of the metal in the cylindrical part
      #Vhead = Shead*tskin
      Vhead = Shead * t_head # volume of head

      Sinternal = Scyl + Shead

#--- weights and weight moments
      Whead = rhoskin*gee*Vhead
      Wcyl  = rhoskin*gee*Vcyl
      Wtank = (Wcyl + 2*Whead) *(1.0 + ftankstiff + ftankadd) # What is an appropriate mass addtion from fasteners/ supports

#--- insulation weight!
      N = length(t_cond)
      Vcyl_insul = zeros(N)
      Winsul = zeros(N)
      Shead_insul = zeros(N + 1) #add one for first (tank wall) surface 
      Vhead_insul = zeros(N)
      rho_insul = zeros(N)
      L = Lhead + tskin

      #Assemble array with layer densities
      for i = 1:N
            rho_insul[i] = insulation_density_calc(material_insul[i])
      end

      Ro = Ri = Rtank_outer # Start calculating insulation from the outer wall of the metal tank ∴Ri of insul = outer R of tank
      Shead_insul[1] = (2.0*π + 4.0*nfweb*thetafb)*(Ro)^2* ( 0.333 + 0.667*(L/Ro)^1.6 )^0.625
      
      for n in 1:N
            
            Ro = Ro + t_cond[n]
            L  = L  + t_cond[n]
            # println("AR ≈ $(Ro/L)")
            Vcyl_insul[n]  = (π * ( Ro^2 - Ri^2 ) * l_cyl)
            Shead_insul[n+1] = (2.0*π + 4.0*nfweb*thetafb)*(Ro)^2* ( 0.333 + 0.667*(L/Ro)^1.6 )^0.625

            Area_coeff = Shead_insul[n+1] / Ro^2 #coefficient that relates area and radius squared
            Vhead_insul[n] = ((Shead_insul[n] + Shead_insul[n+1])/2 - Area_coeff/(6) * t_cond[n]^2) * t_cond[n] #Closed-form solution
            
            Winsul[n] = (Vcyl_insul[n] + 2*Vhead_insul[n]) * rho_insul[n] * gee
            # println("AR = $(Ro/L)")
            Ri = Ro
      end

      Winsul_sum = sum(Winsul)
      Wtank = (Wtank + Winsul_sum)
      l_tank = l_cyl + 2*Lhead
#--- overall tank weight
      Wtank_total = Wtank + Wfuel_tot

return  Wtank_total, l_cyl, tskin, Rtank_outer, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul, Sinternal, Shead_insul, l_tank
end

function size_outer_tank(fuse_tank, l_cyl)
      #Unpack parameters in fuse_tank
      poiss = fuse_tank.poissouter
      Eouter = fuse_tank.Eouter
      ρouter = fuse_tank.rhoouter
      ftankstiff = fuse_tank.ftankstiff
      ftankadd = fuse_tank.ftankadd
      wfb = fuse_tank.wfb
      nfweb = fuse_tank.nfweb
      ARtank = fuse_tank.ARtank

      pc = 4 * pref #4*p_atm; Collapsing pressure, Eq. (7.11) in Barron (1985)

      #Assumed parameters
      L_Do = 0.5 #TODO get rid of this hack and size stiffeners properly

      #Calculate outer tank geometry
      Rtank_outer = fuse_tank.Rfuse - fuse_tank.clearance_fuse
      Do = 2 * Rtank_outer #outside diameter

      #Find cylinder wall thickness. This applies to a short cylinder.
      pressure_res(t_D) = 2.42*Eouter*(t_D)^(5/2) / ( (1 - poiss^2)^(3/4) * (L_Do - 0.45*sqrt(t_D)) ) - pc
      t_Do = find_zero(pressure_res, 1e-3) #Find root with Roots.jl
      t_cyl = t_Do * Do

      #Find head wall thickness
      if ARtank == 2.0
            K1 = 0.90# See table 7.6 for D/D1=2.0 in Barron p. 367
      elseif ARtank == 1.0
            K1 = 0.50
      else  
            println("ARtank of heads not supported, see size_outer_tank()")
            K1=1.0
      end
      t_head = K1 * Do * sqrt(pc * sqrt(3*(1 - poiss^2))/ (0.5*Eouter))

      ## Areas
      wfblim = max( min( wfb , Rtank_outer) , 0.0 )
      thetafb = asin(wfblim / Rtank_outer)

      Shead = (2.0*π + 4.0*nfweb*thetafb)*(Rtank_outer)^2* (0.333 + 0.667*(1/ARtank)^1.6 )^0.625
      Scyl  = 2π*Rtank_outer*l_cyl  # Cross-sectional area

      S_outer = Shead + 2 * Scyl

      ## Volume and Weight
      Vcyl  = Scyl*t_cyl
      Vhead = Shead*t_head

      Wcyl  = Vcyl *ρouter*gee
      Whead =  Vhead*ρouter*gee

      Wtank =(Wcyl + 2 * Whead) * (1.0 + ftankstiff + ftankadd)

      return Wtank, Wcyl, Whead, S_outer, Shead, Scyl, t_cyl, t_head
end

"""
      insulation_density_calc(material)

This function calculates the density of different insulation materials.
      
!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `material::String`: material name.

      **Outputs:**
      - `ρ::Float64`: mass density (kg/m^3).
"""
function insulation_density_calc(material::String)
      if lowercase(material) == "rohacell41s"
            ρ = 35.0 #kg/m^3. From Brewer (1991)
      elseif lowercase(material) == "polyurethane27"
            ρ = 27.0 #kg/m^3
      elseif lowercase(material) == "polyurethane32"
            ρ = 32.0 #kg/m^3
      elseif lowercase(material) == "polyurethane35"
            ρ = 35.0 #kg/m^3
      elseif lowercase(material) == "vacuum"
            ρ = 0 #kg/m^3
      else
            error("Insulation materials currently supported are
                  [rohacell41S, polyurethane27, polyurethane32, polyurethane35, vacuum],
                  but you supplied $material")
      end
      return ρ
end
