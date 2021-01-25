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
                      Wfuel, m_boiloff, thickness_insul, t_cond, clearance_fuse)

#--- fuselage skin and center web thicknesses to withstand pressure load
      Rtank_outer = Rfuse - thickness_insul - clearance_fuse
      #tskin = 2.35 * deltap * Rtank_outer / sigskin
      tskin = 2.35 * 1.1 * deltap * 2 * Rtank_outer / (2 * sigskin * 0.8 + 0.8 * 1.1 * deltap) #2.35 FOS https://www.nrel.gov/docs/fy02osti/32405b27.pdf
      Rtank = Rtank_outer - tskin
      tfweb = 2.0 * deltap * wfb  / sigskin
      Rhead = 0.8 * Rtank
      K = (1/6) * (2+1.6)
      t_head = 2.35 * 1.1 * deltap * 2 * Rtank_outer * K/ (2 * sigskin * 0.8 + 2 * 1.1 * deltap * (K-1))

#--- Calculate length of shell

      Vfuel = Wfuel / (gee * rhoFuel)
      Vfuel_allowance = (0.031 * Vfuel) + Vfuel #Recommended by Verstraete
      V_ellipsoid = 8 * (Rtank^3) * 0.5 * pi / 12  #https://neutrium.net/equipment/volume-and-wetted-area-of-partially-filled-horizontal-vessels/
      V_remaining = Vfuel - V_ellipsoid
      lshell = Vfuel_allowance / (pi * (Rtank^2))

#--- tank cross-section geometric parameters
      wfblim = max( min( wfb , Rtank) , 0.0 )
      thetafb = asin(wfblim / Rtank)
      hfb = sqrt(Rtank^2 - wfb^2)
      sin2t = 2.0*hfb*wfb/Rtank^2
      perim = (2.0*pi + 4.0*thetafb)*Rtank + 2.0*dRfuse

#--- areas
      Askin = (2.0*pi+4.0*nfweb*thetafb)*Rtank*tskin + 2.0*dRfuse*tskin
      Afweb = nfweb*(2.0*hfb+dRfuse)*tfweb
      Atank = (pi + nfweb*(2.0*thetafb + sin2t))*Rtank^2 + 2.0*Rtank*dRfuse + 2.0*(Rtank+nfweb*wfb)*dRfuse
      Shead = (2.0*pi + 4.0*nfweb*thetafb)*Rtank^2* ( 0.333 + 0.667*(Rhead/Rtank)^1.6 )^0.625
#--- component volumes
      Vcyl  = Askin*lshell
      #Vhead = Shead*tskin
      Vhead = Shead * t_head

#--- weights and weight moments
      Wtank = rhoskin*gee*(Vcyl+Vhead)*(1.0+fstring+ffadd)
      #Wtank = Wtank*(1.0+fstring+ffadd)

#--- insulation weight!
      N = length(t_cond)
      Vinsul = zeros(N)
      Winsul = zeros(N)
      s=0 #thickness of previous layer
      for n in 1:N
            Vinsul[n] = pi * (((Rtank_outer+sum(t_cond[1:n]))^2)-((Rtank_outer+s)^2)) * lshell
            Winsul[n] = Vinsul[n] * rho_insul[n] * gee
            s = sum(t_cond[1:n])
      end
      Winsul_sum = sum(Winsul)
      #Winsul = Wppinsul*(1.1*pi+2.0*thetafb)*Rtank*lshell

#--- overall tank weight
      Wtank_total = Wtank + Wfuel + Winsul_sum + 20 * gee #20kg allowance according to Verstraete

#--- pressurized tank volume
      #tankVol = Atank*(lshell + 0.67*Rfuse)

return  Wtank_total, lshell, tskin, Rtank, Vfuel, Wtank
end
