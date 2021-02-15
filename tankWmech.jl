"""
tankWmech calculates weight of the tank that holds the LH2
Inputs:
NOTE: Everything is in SI units.
NOTE: Al alloy 2219 has been recommended as tank material (from H2 tank paper in OneNote)
-Rfuse is fuselage radius (m), dRfuse (m) is the subtraction factor that accounts for flatness of fuselage at bottom
-wfb, nfb are parameters for multiple-bubble configuration
-lcv (J/kg) is lower calorific value of fuel
-eta is overall efficiency of gas turbine/turboelectric powertrain etc.
-sigskin, rhoskin are material properties
-thickness_insul is insulation thickness of all layers combined (m)

Outputs:

- Wtank_total: tank weight including fuel (N)
- Wtank: empty tank weight (N)
"""
function tankWmech(gee, ρfuel,
                      fstring, ffadd, Δp,
                      Rfuse, dRfuse, wfb, nfweb,
                      sigskin, rho_insul, rhoskin,
                      Wfuel, m_boiloff, thickness_insul, t_cond, clearance_fuse, AR)

#--- fuselage skin and center web thicknesses to withstand pressure load
      weld_eff = 0.9
      ullage_frac = 0.1
      β = 2.0

      # Add an additional 
      Δpdes = Δp*β

      Rtank_outer = Rfuse - thickness_insul - clearance_fuse

      tskin = Δpdes * (2 * Rtank_outer) / (2 * sigskin * weld_eff + 0.8 * Δpdes)
      Rtank = Rtank_outer - tskin
      tfweb = 2.0 * Δpdes * wfb  / sigskin
      Lhead = Rtank/AR       # eg. for a 2:1 ellipsoid majorax/minorax = 2/1 ⟹ R/Lhead = 2/1 
      K = (1/6) * (AR^2 + 2) # Aspect ratio of 2:1 for the head (# Barron pg 359) 
      t_head = Δpdes* (2*Rtank_outer) * K/ (2 * sigskin * weld_eff + 2 * Δpdes * (K-.1)) #Verstraete

#--- Calculate length of shell
      Wfuel_tot = Wfuel + m_boiloff * gee
      Vfuel = Wfuel_tot / (gee * ρfuel)
      Vinternal = (1 + ullage_frac)*Vfuel # V. thesis says ~3% but Barron recommends 10%
      V_ellipsoid = 2π*(Rtank^3/AR)/3  # Half the vol of std ellipsoid = 1/2×(4π/3 ×(abc)) where a,b,c are the semi-axes length. Here a = R/AR, b=c=R
                                       # Also see: https://neutrium.net/equipment/volume-and-wetted-area-of-partially-filled-horizontal-vessels/
      V_cylinder = Vinternal - 2*V_ellipsoid
      lshell = V_cylinder / (π * (Rtank^2))

#--- tank cross-section geometric parameters
      wfblim = max( min( wfb , Rtank) , 0.0 )
      thetafb = asin(wfblim / Rtank)
      hfb = sqrt(Rtank^2 - wfb^2)
      sin2t = 2.0*hfb*wfb/Rtank^2

#--- areas
      Askin = (2.0*π+4.0*nfweb*thetafb)*Rtank*tskin + 2.0*dRfuse*tskin
      #Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      #Atank = (π + nfweb*(2.0*thetafb + sin2t))*Rtank^2 + 2.0*Rtank*dRfuse + 2.0*(Rtank+nfweb*wfb)*dRfuse
      Shead = (2.0*π + 4.0*nfweb*thetafb)*Rtank^2* ( 0.333 + 0.667*(Lhead/Rtank)^1.6 )^0.625 # This form is better for insul thickness 
                                                                                          # but just as a note to reader this comes from  semi- oblate spheroid surf area is ≈ 2π/3×R²[1/3 + 2/3×(1/AR)^1.6]^(1/1.6)
#--- component volumes
      Vcyl  = Askin*lshell #volume of the metal
      #Vhead = Shead*tskin
      Vhead = Shead * t_head #volume of head

#--- weights and weight moments
      Whead = rhoskin*gee*Vhead
      Wcyl = rhoskin*gee*Vcyl
      Wtank = (Wcyl + 2*Whead) *(1.0 + fstring + ffadd) # What is an appropriate mass addtion from fasteners/ supports

#--- insulation weight!
      N = length(t_cond)
      Vcyl_insul = zeros(N)
      Winsul = zeros(N)
      Shead_insul = zeros(N)
      Vhead_insul = zeros(N)

      Ri = Rtank_outer # Start calculating insulation from the outer wall of the metal tank ∴Ri of insul = outer R of tank
      for n in 1:N
            
            Ro = Rtank_outer + sum(t_cond[1:n])
            L  = Lhead + sum(t_cond[1:n])

            Vcyl_insul[n]  = (π * ( Ro^2 - Ri^2 ) * lshell)
            Shead_insul[n] = (2.0*π + 4.0*nfweb*thetafb)*(Ro)^2* ( 0.333 + 0.667*(L/Ro)^1.6 )^0.625
            Vhead_insul[n] = Shead_insul[n] * t_cond[n]
            Winsul[n] = (Vcyl_insul[n] + 2*Vhead_insul[n]) * rho_insul[n] * gee
            
            Ri = Ro
      end
      Winsul_sum = sum(Winsul)

#--- overall tank weight
      Wtank_total = Wtank + Wfuel_tot + Winsul_sum
      l_tank = lshell + 2*Lhead

return  Wtank_total, lshell, tskin, Rtank_outer, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead_insul, l_tank
end
