"""
tankWmech calculates weight of the tank that holds the LH2
Inputs:
NOTE: Everything is in SI units.
NOTE: Al alloy 2219 has been recommended as tank material (from H2 tank paper in OneNote)
-LD is L/D of airplane
-xshell1, xshell2 are start and end coordinates of tank (in meters)
-gee is gravitational acceleration (m2/s)
-Rfuse is fuselage radius (m), dRfuse (m) is the subtraction factor that accounts for flatness of fuselage at bottom
-wfb, nfb are parameters for multiple-bubble configuration
-m_airplane is airplane mass (kg)
-R (m) is specified range for a given mission
-lcv (J/kg) is lower calorific value of fuel
-eta is overall efficiency of gas turbine/turboelectric powertrain etc.
-sigskin, rhoskin are material properties
-thickness_insul is insulation thickness of all layers combined (m)

Outputs:

- Wtank_total: tank weight including fuel (N)
- Wtank: empty tank weight (N)
"""
function tankWmech(gee, rhoFuel,
                      fstring, ffadd, deltap,
                      Rfuse, dRfuse, wfb, nfweb,
                      sigskin, rho_insul, rhoskin,
                      Wfuel, m_boiloff, thickness_insul, t_cond, clearance_fuse, AR)

#--- fuselage skin and center web thicknesses to withstand pressure load
      Rtank_outer = Rfuse - thickness_insul - clearance_fuse
      #tskin = 2.35 * Δp * Rtank_outer / sigskin
      tskin = 1.1 * Δp * (2 * Rtank_outer) / (2 * sigskin * weld_eff + 0.8 * 1.1 * Δp) #2.35 FOS https://www.nrel.gov/docs/fy02osti/32405b27.pdf
      Rtank = Rtank_outer - tskin
      tfweb = 2.0 * Δp * wfb  / sigskin
      #Rhead = 0.8 * Rtank * 2
      Lhead = Rtank/AR       # eg. for a 2:1 ellipsoid majorax/minorax = 2/1 ⟹ R/Lhead = 2/1 
      K = (1/6) * (AR^2 + 2) #Verstraete aspect ratio of 2:1 for the head (# Barron pg 359)
      t_head = 1.1 * Δp * (2*Rtank_outer) * K/ (2 * sigskin * weld_eff + 2 * 1.1 * Δp * (K-.1)) #Verstraete

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
      perim = (2.0*pi + 4.0*thetafb)*Rtank + 2.0*dRfuse

#--- areas
      Askin = (2.0*pi+4.0*nfweb*thetafb)*Rtank*tskin + 2.0*dRfuse*tskin
      #Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      #Atank = (pi + nfweb*(2.0*thetafb + sin2t))*Rtank^2 + 2.0*Rtank*dRfuse + 2.0*(Rtank+nfweb*wfb)*dRfuse
      Shead = (2.0*pi + 4.0*nfweb*thetafb)*Rtank^2* ( 0.333 + 0.667*(Lhead/Rtank)^1.6 )^0.625
#--- component volumes
      Vcyl  = Askin*lshell #volume of the metal
      #Vhead = Shead*tskin
      Vhead = Shead * t_head #volume of head

#--- weights and weight moments
      Whead = rhoskin*gee*Vhead
      Wcyl = rhoskin*gee*Vcyl
      Wtank = (Wcyl + 2*Whead) #*(1.0+fstring+ffadd)

#--- insulation weight!
      N = length(t_cond)
      Vcyl_insul = zeros(N)
      Winsul = zeros(N)
      Shead_insul = zeros(N)
      Vhead_insul = zeros(N)
      s=0 #thickness of previous layer
      for n in 1:N
            Vinsul[n] = (pi * (((Rtank_outer+sum(t_cond[1:n]))^2)-((Rtank_outer+s)^2)) * lshell)# + (pi * t_cond[n] * (Rtank_outer+sum(t_cond[1:n]))^2)
            Shead_insul[n] = (2.0*pi + 4.0*nfweb*thetafb)*(Rtank+sum(t_cond[1:n]))^2* ( 0.333 + 0.667*((Lhead+sum(t_cond[1:n]))/(Rtank+sum(t_cond[1:n])))^1.6 )^0.625
            Vhead_insul[n] = 2 * Shead_insul[n] * t_cond[n]
            Winsul[n] = (Vinsul[n] + Vhead_insul[n]) * rho_insul[n] * gee
            s = sum(t_cond[1:n])
      end
      Winsul_sum = sum(Winsul)

#--- overall tank weight
      Wtank_total = Wtank + Wfuel + Winsul_sum
      l_tank = lshell + 2*Lhead

return  Wtank_total, lshell, tskin, Rtank_outer, Vfuel, Wtank, Wfuel_tot, Winsul_sum, t_head, Whead, Wcyl, Winsul, Shead_insul, l_tank
end
