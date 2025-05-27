using StaticArrays

"""
    fuselage_drag!(fuse, parg, para, ip)

Calculates surface velocities, boundary layer, wake 
for a quasi-axisymmetric body in compressible flow. Formerly, `fusebl!()`.

A compressible source line represents the potential flow.
An integral BL formulation with lateral divergence represents the surface BL and wake.
An added-source distribution represents the viscous displacement influence on the potential flow.
The body shape is defined by its area and perimeter distributions `A(x)`,  `b0(x)`,
which are defined by the various geometric parameters in `parg`. 

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `fuse::Fuselage`: Fuselage model object.
      - `parg::AbstractArray{Float64}`: Vector of `aircraft` model geometry parameters.
      - `para::AbstractArray{Float64}`: Vector of `aircraft` model aerodynamic parameters.
      - `parm::AbstractArray{Float64}`: Vector of `aircraft` model mission parameters.
      - `ip::Integer`: Index of flight point in `par` arrays.
      
      **Outputs:**
      - No explicit outputs. Computed drag values are saved to `para` of `aircraft` model.

See [Simplified Viscous/Inviscid Analysis for Nearly-Axisymmetric Bodies](../assets/drela_TASOPT_2p16/axibl.pdf).
See also [`_axisymm_BL`](@ref) and [`_axisymm_flow`](@ref).

!!! compat "Future Changes"
      In an upcoming revision, an `aircraft` `struct` and auxiliary indices will be passed in lieu of pre-sliced `par` arrays.

"""
function fuselage_drag!(fuse, parm, para, ip)
      
#     nc,     # number of control points for fuselage potential-flow problem
#     nbl,    # number of BL+wake points
#     iblte   # index of TE point

      nbldim = 60
      xbl  = @MVector zeros(nbldim)     # body x coordinates
      zbl  = @MVector zeros(nbldim)     # body z coordinates
      sbl  = @MVector zeros(nbldim)     # body + wake arc length  (BL coordinate)
      dybl = @MVector zeros(nbldim)     # body y-offset of edge-type tail
      bbl  = @MVector zeros(nbldim)     # transverse width (body perimeter)
      rnbl = @MVector zeros(nbldim)     # dr/dn  (cosine of body contour angle from axis)
      uinv = @MVector zeros(nbldim)     # inviscid edge velocity (in absence of delta*)
      uebl = @MVector zeros(nbldim)     # actual edge velocity
      dsbl = @MVector zeros(nbldim)     # displacement thickness delta*
      thbl = @MVector zeros(nbldim)     # momentum thickness theta
      tsbl = @MVector zeros(nbldim)     # K.E. thickness  theta*
      dcbl = @MVector zeros(nbldim)     # density flux thickness  delta**
      ctbl = @MVector zeros(nbldim)     # max shear stress coefficient
      cqbl = @MVector zeros(nbldim)     # equilib.max shear stress coefficient
      cfbl = @MVector zeros(nbldim)     # skin friction coefficient
      cdbl = @MVector zeros(nbldim)     # dissipation coefficient
      hkbl = @MVector zeros(nbldim)     # kinematic shape parameter
      phbl = @MVector zeros(nbldim)     # running integral of dissipation

#       include 'constants.inc'

      Pend, Pinf, KTE, Kend, Kinf = zeros(5)

      Vol = 0.0

      tapers_to_edge = compare_strings(fuse.layout.opt_tapers_to, "edge")
      xnose = fuse.layout.x_nose
      xend = fuse.layout.x_end
      xblend1 = fuse.layout.x_start_cylinder
      xblend2 = fuse.layout.x_end_cylinder

      Mach  = para[iaMach, ip]
      altkm = para[iaalt, ip]/1000.0
      ΔTatmos = parm[imDeltaTatm] #atmosphere temperature difference
      T0,p0,rho0,a0,mu0 = atmos(altkm, ΔTatmos) #get atmospheric parameters
    
      Reunit = Mach*a0 * rho0/mu0

      wfb = fuse.layout.bubble_center_y_offset
      Rfuse = fuse.layout.radius
      dRfuse = fuse.layout.bubble_lower_downward_shift

#---- fuselage cross-section geometric parameters
      wfblim = max( min( wfb , Rfuse ) , 0.0 )
      thetafb = asin(wfblim/Rfuse)
      hfb = sqrt(Rfuse^2 - wfb^2)
      sin2t = 2.0*hfb*wfb/Rfuse^2
      Sfuse = (pi + 2.0*thetafb + sin2t)*Rfuse^2 + 2.0*Rfuse*dRfuse

      anose = fuse.layout.nose_radius
      btail = fuse.layout.tail_radius

#---- calculate potential-flow surface velocity uinv(.) using PG source line
      nc = 30

      nbl, iblte =  _axisymm_flow(xnose,xend,xblend1,xblend2,Sfuse, 
                            anose,btail,Int(tapers_to_edge),
                            Mach, nc, nbldim,  xbl,zbl,sbl,dybl,uinv)
     
#---- fuselage volume and perimeter
      @inbounds for i = 1:iblte-1
            ra = 0.5*(zbl[i+1]+zbl[i])
            Vol += π*ra^2*(xbl[i+1] - xbl[i])
      end
      fuse.volume = Vol
      
      if !(tapers_to_edge)  #i.e., tapers to point
       @inbounds for  ibl = 1: nbl
         bbl[ibl] = 2.0*pi*zbl[ibl]
       end
      else #tapers to edge
       @inbounds for  ibl = 1: nbl
         bbl[ibl] = 2.0*pi*zbl[ibl] + 4.0*dybl[ibl]
       end
      end

#---- dr/dn  (cosine of body contour angle from axis)
      @inbounds for  ibl = 2: nbl-1
        dxm = xbl[ibl] - xbl[ibl-1]
        dzm = zbl[ibl] - zbl[ibl-1]
        dsm = sbl[ibl] - sbl[ibl-1]

        dxp = xbl[ibl+1] - xbl[ibl]
        dzp = zbl[ibl+1] - zbl[ibl]
        dsp = sbl[ibl+1] - sbl[ibl]

        dxo = (dxm/dsm)*dsp + (dxp/dsp)*dsm
        dzo = (dzm/dsm)*dsp + (dzp/dsp)*dsm

        rnbl[ibl] = dxo / sqrt(dxo^2 + dzo^2)
      end
#---- extrapolate to endpoints
      rnbl[1] = rnbl[2] - (0.5*(sbl[3]+sbl[2]) - sbl[1]) * (rnbl[3]-rnbl[2])/(sbl[3]-sbl[2])
      rnbl[nbl] = rnbl[nbl-1] + (sbl[nbl] - 0.5*(sbl[nbl-1]+sbl[nbl-2])) * (rnbl[nbl-1]-rnbl[nbl-2])/(sbl[nbl-1]-sbl[nbl-2])
      rnbl[1] = max( rnbl[1] , 0.0 )
      rnbl[nbl] = max( rnbl[nbl] , 0.0 )

#---- perform viscous/inviscid BL calculation driven by uinv(.)
      fex = para[iafexcdf, ip]
     
      uebl, dsbl, thbl, tsbl, dcbl,
      cfbl, cdbl, ctbl, hkbl, phbl  = _axisymm_BL(nbldim, nbl,iblte, 
					sbl, bbl, rnbl, uinv, Reunit, Mach, fex)

      gam = 1.4 #gamSL
      gmi = gam - 1.0

#---- KE defect at TE and surface dissipation
      i = iblte
      trbl = 1.0 + 0.5*gmi*Mach^2*(1.0 - uebl[i]^2) 
      rhbl = trbl^(1.0/gmi)
      KTE = 0.5*rhbl*uebl[i]^3 * tsbl[i] * (bbl[i] + 2.0*pi*dsbl[i])
      Difsurf = phbl[i]

#---- momentum and KE defects and accumulated dissipation at end of wake
      i = nbl
      trbl = 1.0 + 0.5*gmi*Mach^2*(1.0 - uebl[i]^2) 
      rhbl = trbl^(1.0/gmi)
      Pend =     rhbl*uebl[i]^2 * thbl[i] * (bbl[i] + 2.0*pi*dsbl[i])
      Kend = 0.5*rhbl*uebl[i]^3 * tsbl[i] * (bbl[i] + 2.0*pi*dsbl[i])
      Difend = phbl[i]

#---- far-downstream momentum defect via Squire-Young (note that Vinf = 1 here)
      Hend = dsbl[i]/thbl[i]
      Hinf = 1.0 + gmi*Mach^2
      Havg = 0.5*(Hend+Hinf)
      Pinf = Pend * uebl[i]^Havg

#---- far-downstream KE defect  0.5 rho V^3 Theta*  (note that Vinf = 1 here)
      Kinf = Pinf

#---- additional dissipation downstream of last wake point
      tsinf = 2.0*Kinf
      dcinf = 0.5*gmi*Mach^2 * tsinf
      dDif = Kinf - Kend + 0.5*(dcinf + dcbl[i]) * (1.0 - uebl[i])
      Difinf = Difend + dDif

#---- wake dissipation
      Difwake = Difinf - Difsurf

#---- store dissipation, KE, and drag areas
      Vinf = 1.0
      qinf = 0.5
      para[iaDAfsurf, ip] = Difsurf/(qinf*Vinf)
      para[iaDAfwake, ip] = Difwake/(qinf*Vinf)
      para[iaKAfTE  , ip] = KTE/(qinf*Vinf)
      para[iaPAfinf , ip] = Pinf/qinf
    
end # fuselage_drag!

"""
    broadcast_fuselage_drag!(para, ip)

Broadcasts fuselage drag calculations across all flight points in `para` arrays.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `para::AbstractArray{Float64}`: Vector of `aircraft` model aerodynamic parameters.
      - `ip::Integer`: Index of flight point in `par` arrays.
      
      **Outputs:**
      - No explicit outputs. Values are saved to `para` of `aircraft` model.
"""
function broadcast_fuselage_drag!(para, ip)
    #---- assume K.E., dissipation, drag areas will be the same for all points
    KAfTE   = para[iaKAfTE  , ip] # Kinetic energy area at T.E.
    DAfsurf = para[iaDAfsurf, ip] # Surface dissapation area 
    DAfwake = para[iaDAfwake, ip] # Wake dissapation area
    PAfinf  = para[iaPAfinf , ip] # Momentum area at ∞

    # Assume K.E., Disspation and momentum areas are const. for all mission points:
    para[iaKAfTE  , :] .= KAfTE
    para[iaDAfsurf, :] .= DAfsurf
    para[iaDAfwake, :] .= DAfwake
    para[iaPAfinf , :] .= PAfinf
end