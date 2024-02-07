"""
      tankWmech(gee, ρfuel,
            ftankstiff, ftankadd, Δp,
            Rfuse, dRfuse, wfb, nfweb,
            sigskin, rho_insul, rhoskin,
            Wfuel, m_boiloff, t_cond, clearance_fuse, AR)

`tankWmech` calculates the weight of the cryogenic fuel tank for a LH-fueled aircraft.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `gee::Float64`: Gravitational acceleration (m/s^2).
      - `ρfuel::Float64`: Fuel density (kg/m^3).
      - `ftankstiff::Float64`: Tank stiffness factor.
      - `ftankadd::Float64`: Additional mass factor for the tank.
      - `Δp::Float64`: Pressure difference (Pa).
      - `Rfuse::Float64`: Fuselage radius (m).
      - `dRfuse::Float64`: Subtraction factor accounting for fuselage flatness (m).
      - `wfb`: Parameters for multiple-bubble configuration.
      - `nfweb`: Number of bubbles.
      - `sigskin::Float64`: Material property.
      - `material_insul::Array{String,1}`: material name for each MLI layer.
      - `rhoskin::Float64`: Material property.
      - `Wfuel::Float64`: Weight of fuel (N).
      - `m_boiloff::Float64`: Mass boiled off during the mission flight (kg).
      - `t_cond::Array{Float64,1}`: Thickness of insulation layers (m).
      - `clearance_fuse::Float64`: Clearance for the fuselage (m).
      - `AR::Float64`: Aspect ratio.

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
function tankWmech(gee::Float64, ρfuel::Float64,
                  ftankstiff::Float64, ftankadd::Float64, Δp::Float64,
                  Rfuse::Float64, dRfuse::Float64, wfb, nfweb,
                  sigskin, material_insul, rhoskin,
                  Wfuel, m_boiloff, t_cond::Array{Float64,1}, clearance_fuse, AR)
      
# Total thickness:
      thickness_insul = sum(t_cond)

# Input paramters:
      weld_eff = 0.9 #lower strength due to welding?
      ullage_frac = 0.1 # V. thesis says ~3% but Barron recommends 10% #TODO: figure out who these people are
      β = 2.0 #TODO: what is this? a safety factor?
# Add an additional pressure factor
      Δpdes = Δp * β

      Rtank_outer = Rfuse - thickness_insul - clearance_fuse

      tskin = Δpdes * (2 * Rtank_outer) / (2 * sigskin * weld_eff + 0.8 * Δpdes) #TODO: this sizes the skin thickness
                                                                              #based on stress; figure out where it comes from

      Rtank = Rtank_outer - tskin
      tfweb = 2.0 * Δpdes * wfb  / sigskin
      Lhead = Rtank / AR       # eg. for a 2:1 ellipsoid majorax/minorax = 2/1 ⟹ R/Lhead = 2/1 
      
      K = (1/6) * (AR^2 + 2) # Aspect ratio of 2:1 for the head (# Barron pg 359) 
      t_head = Δpdes* (2*Rtank_outer) * K/ (2 * sigskin * weld_eff + 2 * Δpdes * (K-.1)) #Verstraete

#--- Calculate length of cylindrical portion
      Wfuel_tot = Wfuel + (m_boiloff * gee)
      Vfuel = Wfuel_tot / (gee * ρfuel)
      Vinternal = (1 + ullage_frac)*Vfuel  # required interal volume
      V_ellipsoid = 2π * (Rtank^3 / AR) / 3  # Half the vol of std ellipsoid = 1/2×(4π/3 ×(abc)) where a,b,c are the semi-axes length. Here a = R/AR, b=c=R
                                       # Also see: https://neutrium.net/equipment/volume-and-wetted-area-of-partially-filled-horizontal-vessels/
      V_cylinder = Vinternal - 2*V_ellipsoid
      l_cyl = V_cylinder / (π * (Rtank^2)) #required length of cylindrical portion

#--- tank cross-section geometric parameters
      wfblim = max( min( wfb , Rtank) , 0.0 )
      thetafb = asin(wfblim / Rtank)
      hfb = sqrt(Rtank^2 - wfb^2)
      sin2t = 2.0 * hfb * wfb / Rtank^2

#--- areas
      Askin = (2.0*π+4.0*nfweb*thetafb)*Rtank*tskin + 2.0*dRfuse*tskin # Cross-sectional area of the cylindrical part
      #Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      #Atank = (π + nfweb*(2.0*thetafb + sin2t))*Rtank^2 + 2.0*Rtank*dRfuse #+ 2.0*(Rtank+nfweb*wfb)*dRfuse
      Shead = (2.0*π + 4.0*nfweb*thetafb)*Rtank^2* ( 0.333 + 0.667*(Lhead/Rtank)^1.6 )^0.625 # This form is better for insul thickness 
                                                                                          # but just as a note to reader this comes from  semi- oblate spheroid surf area is ≈ 2π×R²[1/3 + 2/3×(1/AR)^1.6]^(1/1.6)
#--- component volumes
      Vcyl  = Askin*l_cyl    # volume of the metal in the cylindrical part
      #Vhead = Shead*tskin
      Vhead = Shead * t_head # volume of head

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
            Vhead_insul[n] = (Shead_insul[n] + Shead_insul[n+1])/2  * t_cond[n]
            
            Winsul[n] = (Vcyl_insul[n] + 2*Vhead_insul[n]) * rho_insul[n] * gee
            # println("AR = $(Ro/L)")
            Ri = Ro
      end
      
      Winsul_sum = sum(Winsul)
      Wtank = (Wtank + Winsul_sum)
#--- overall tank weight
      Wtank_total = Wtank + Wfuel_tot
      l_tank = l_cyl + 2*Lhead

return  Wtank_total, l_cyl, tskin, Rtank_outer, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead_insul, l_tank
end

function insulation_density_calc(material)
      if material == "rohacell31"
            ρ = 32.0 #kg/m^3. From manufacturer sheet
      elseif material == "polyurethane"
            ρ = 27.0 #kg/m^3. From Brewer (1991)
      end
      return ρ
end
