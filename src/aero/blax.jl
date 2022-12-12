"""
    blax(ndim, n,ite, xi, bi, rni, uinv, Reyn, Mach, fexcr)
     
Axisymmetric boundary layer + wake calculation routine.
Uses specified inviscid velocity, corrects for viscous
displacement to allow calculation of separated flow.

# Inputs

- `ndim::Integer`: physical array dimension.
- `n::Integer`: number of BL+wake points.
- `ite::Integer`: index of trailing edge point, start of wake.
- `xi::Array{Float64}`: arc length array (BL coordinate).
- `bi::Array{Float64}`: lateral width of BL (body perimeter, zero for wake) = 1 for 2D.
- `rni::Array{Float64}`: ``dr/dn`` ``= 0`` for 2D.
- `uinv::Array{Float64}`: inviscid velocity.
- `Reyn::Float64`: Reynolds number,  ``rho_ref u_ref l_ref / mu_ref``.
- `Mach::Float64`: Mach number,  ``u_ref / a_ref``.
- `fexcr::Float64`: excrescence multiplier, applied to wall `Cf` ``= 1`` for smooth wall.

Assumed units for all quantities:
- `l_ref`: same unit as used for input `xi`,`bi`.
- `u_ref`: freestream velocity.
- `a_ref`: freestream speed of sound.
- `rho_ref`: freestream density.
- `mu_ref`: freestream viscosity.
    
# Outputs
- `uei::Array{Float64}`:  edge velocity, ( ``= uinv[i] + displacement correction`` ).
- `dsi::Array{Float64}`:  displacement thickness (``\\delta^*``).
- `thi::Array{Float64}`:  momentum thickness     (``\\theta``).
- `tsi::Array{Float64}`:  kinetic energy thickness (``\\theta^*``).
- `dci::Array{Float64}`:  density flux thickness   (``\\delta^{**}``).
- `cfi::Array{Float64}`:  skin friction coefficient, normalized with local `rho`, `u`.
- `cdi::Array{Float64}`:  dissipation coefficient  , normalized with local `rho`,`u`.
- `cti::Array{Float64}`:  max shear-stress coefficient, normalized with local `rho`,`u`.
- `hki::Array{Float64}`:  kinematic shape parameter.
- `phi::Array{Float64}`:  integrated dissipation.


Other outputs of interest can be computed as follows.
  These are in units of `l_ref`, `rho_ref`, `u_ref`

- `beff`: Effective perimeter,  ``bi  +  2 pi dsi rni``.
- `rhi`: Edge density, ``(1 + 0.5*(gam-1)*Mach^2*(1.0-uei^2))^(1/(gam-1))``.
- `mdef`: Total mass defect, ``rhi uei   dsi beff``.
- `Pdef`: Total mom. defect, ``rhi uei^2 thi beff``.
- `Edef`: Total KE defect, ``rhi uei^3 tsi beff / 2``.
- `rhi`: Wall shear force/span, ``tw b =   uei^2 cfi beff / 2``.
- `Diss`: Dissipation integral, ``rhi uei^3 cdi beff``.

Body profile drag `Dp` is the far-downstream momentum defect `Pinf`,
best obtained by applying Squire-Young to the last wake point ``i = n`` :

- ``Pend = rhi*uei^2 * thi * beff``
- ``Hend = dsi/thi``
- ``Hinf = 1.0 + (gam-1)*Mach^2``
- ``Havg = 0.5*(Hend+Hinf)``
- ``Pinf = Pend * uei^Havg  =  Dp``

See Appendix E.4 of TASOPT docs.
See also [`blsys`](@ref) and [`blvar`](@ref).
"""
function blax(ndim, n,ite, xi, bi, rni, uinv, Reyn, Mach, fexcr)


      uei  = @MVector zeros(ndim) 
      rhi  = @MVector zeros(ndim)
      dsi  = @MVector zeros(ndim)
      thi  = @MVector zeros(ndim)
      tsi  = @MVector zeros(ndim)
      dci  = @MVector zeros(ndim)
      cfi  = @MVector zeros(ndim)
      cdi  = @MVector zeros(ndim)
      cti  = @MVector zeros(ndim)
      hki  = @MVector zeros(ndim)
      phi  = @MVector zeros(ndim)

      # aa = @MMatrix zeros(3,3)
      # bb = @MMatrix zeros(3,3)
      # rr = @MVector zeros(3)
        
      hm,  hm_thm, hm_dsm = 1e20*ones(3)
      hkm, hkm_thm, hkm_dsm, hkm_uem = 0*1e20*ones(4)
      hcm, hcm_thm, hcm_dsm, hcm_uem = 1e20*ones(4)
      hsm, hsm_thm, hsm_dsm, hsm_uem = 1e20*ones(4)
      cfm, cfm_thm, cfm_dsm, cfm_uem = 1e20*ones(4)
      dim, dim_thm, dim_dsm, dim_uem = 1e20*ones(4)
      h , h_th, h_ds,
      hk, hk_th, hk_ds, hk_ue,
      hc, hc_th, hc_ds, hc_ue,
      hs, hs_th, hs_ds, hs_ue,
      cf, cf_th, cf_ds, cf_ue,
      di, di_th, di_ds, di_ue = 1e20*ones(23)
      

      idim = 60
      mdi      = zeros(idim)
      uvis     = zeros(idim)
      uvis_mdi = zeros(idim,idim)
      dthi     = zeros(idim)
      dmdi     = zeros(idim)
      duei     = zeros(idim)
      ddsi     = zeros(idim)

      kdim = 3*n # instead of idim
      asys = zeros(kdim,kdim)
      rsys = zeros(kdim)
    
    simi, lami, wake, direct = true, true, true, true

#---- pi  and  1/(4 pi)
      qopi = 1/(4*pi)
      gam  = 1.4 
      gmi = gam - 1.0

      hksep =  2.9 
      eps =  1.0e-6 # Tolerance for Newton solver

      if(n > idim)
       println("BLAX: Local array overflow.  Increase idim to", n)
       quit()
      end
    
#---- initialize ue to inviscid uinv, and initialize rhoe
      uei[1:n] = uinv[1:n]
      rhi = @. (1 + 0.5*gmi*Mach^2 * (1.0 - uei^2))^(1.0/gmi)

    #  @inbounds for i = 1:n
    #   uei[i] = uinv[i]
    #   trat = 1.0 + 0.5*gmi*Mach^2 * (1.0-uei[i]^2)
    #   rhi[i] = trat^(1.0/gmi)
    #  end

#---- first point is not calculated if xi=0 there
      if(xi[1] == 0.0)
       thi[1] = 0.
       dsi[1] = 0.
       mdi[1] = 0.
       phi[1] = 0.
      end

# =============================================================================
#---- sweep downstream to initialize BL via standard BL march, 
#-     fudging by specifying plausible Hk if separation occurs (Goldstein singularity)
#-    (this is to provide an intial guess only, later the global system is solved simultaneously via Newton iteration)

   direct = true
 @inbounds for i = 2:n

        simi = xi[i-1] == 0.0
#c      lami = simi    # laminar   similarity station
        lami = false   # turbulent similarity station
        wake = i  > ite

        x = xi[i] #xlocation
        b = bi[i] #perimeter
        rn = rni[i] #dr/dn
#[prash] previous points x minus
        xm = xi[i-1]
        bm = bi[i-1]
        rnm = rni[i-1]

#------ set i-1 station, and initialize station i for iteration
        if(simi)
         uem = 0.
         thm = 0.
         dsm = 0.

         ue = uei[i]
         rex = ue*x * Reyn
         th = 0.4*x / sqrt(rex)
         ds = th * 2.0

        else
         uem = uei[i-1]
         thm = thi[i-1]
         dsm = dsi[i-1]

         th = thi[i-1]
         ds = dsi[i-1]
         if(direct)
#-------- previous station was direct... use specified ue as initial guess
          ue = uei[i]
         else
#-------- previous station was inverse... use previous-station ue as initial guess
          ue = uei[i-1]
         end

        end
#------ always try direct mode first
        direct = true
        hkprev = 0.

      @inbounds for iter = 1:20 
          (h , h_th, h_ds,
          hk, hk_th, hk_ds, hk_ue,
          hc, hc_th, hc_ds, hc_ue,
          hs, hs_th, hs_ds, hs_ue,
          cf, cf_th, cf_ds, cf_ue,
          di, di_th, di_ds, di_ue ) = blvar(simi,lami,wake, Reyn,Mach, fexcr, x, th ,ds ,ue) 
            
          if(iter == 1)
           hkprev = hk
          else
           if(hk > hksep) # & hk > hkprev)
#---------- Hk limit exceeded... switch to inverse mode for this point
            direct = false
           end
          end

#-------- set up 2-point differenced BL equation system
          aa, bb, rr = blsys(simi,lami,wake,direct, Mach, uinv[i],hksep,
                      x,b,rn,th,ds,ue, 
                      h , h_th, h_ds,
                      hk, hk_th, hk_ds, hk_ue,
                      hc, hc_th, hc_ds, hc_ue,
                      hs, hs_th, hs_ds, hs_ue,
                      cf, cf_th, cf_ds, cf_ue,
                      di, di_th, di_ds, di_ue,
                      xm,bm,rnm,thm,dsm,uem, 
                      hm , hm_thm, hm_dsm,
                      hkm, hkm_thm, hkm_dsm, hkm_uem,
                      hcm, hcm_thm, hcm_dsm, hcm_uem,
                      hsm, hsm_thm, hsm_dsm, hsm_uem,
                      cfm, cfm_thm, cfm_dsm, cfm_uem,
                      dim, dim_thm, dim_dsm, dim_uem)

          rr = aa\rr #solve linear system
          dth = -rr[1]
          dds = -rr[2]
          due = -rr[3]

          rlx = 1.0

          if(rlx*dth > 1.6*th)  rlx =  1.6*th/dth;  end
          if(rlx*dth < -0.6*th) rlx = -0.6*th/dth;  end
          if(rlx*dds >  2.5*ds) rlx =  2.5*ds/dds;  end
          if(rlx*dds < -0.4*ds) rlx = -0.4*ds/dds;  end
          if(rlx*due >  0.2*ue) rlx =  0.2*ue/due;  end
          if(rlx*due < -0.1*ue) rlx = -0.1*ue/due;  end

          dmax = max( abs(dth)/th , abs(dds)/ds , abs(due)/ue )
            

          th = th + rlx*dth
          ds = ds + rlx*dds
          ue = ue + rlx*due
            
          if(dmax < eps) break end
            
        end # end iter loop

      uei[i] = ue
      dsi[i] = ds
      thi[i] = th
      tsi[i] = hs*th
      dci[i] = hc*th
      cfi[i] = cf
      cdi[i] = di*hs/2.0
      cti[i] = 0.03 * 0.5*hs*((hk-1.0)/hk)^2
      hki[i] = hk

      dib  = cdi[i]  *rhi[i]  *uei[i]^3   * (b  + 2.0*pi*ds *rn )
      dibm = cdi[i-1]*rhi[i-1]*uei[i-1]^3 * (bm + 2.0*pi*dsm*rnm)
      phi[i] = phi[i-1] + 0.5*(dib + dibm) * (x - xm)
      mdi[i] = ue*ds*(b + 2.0*pi*ds)

      thm = th
      dsm = ds
      hsm = hs
      uem = ue

      hm = h
      hm_thm = h_th
      hm_dsm = h_ds

      hkm = hk
      hkm_thm = hk_th
      hkm_dsm = hk_ds
      hkm_uem = hk_ue

      hsm = hs
      hsm_thm = hs_th
      hsm_dsm = hs_ds
      hsm_uem = hs_ue

      hcm = hc
      hcm_thm = hc_th
      hcm_dsm = hc_ds
      hcm_uem = hc_ue

      dim     = di
      dim_thm = di_th
      dim_dsm = di_ds
      dim_uem = di_ue

      if(i >= ite)
       cfm     = 0.
       cfm_thm = 0.
       cfm_dsm = 0.
       cfm_uem = 0.
      else
       cfm     = cf
       cfm_thm = cf_th
       cfm_dsm = cf_ds
       cfm_uem = cf_ue
      end

    end #end i loop


# =============================================================================
#---- perform global viscous/inviscid iteration with source model 
#-     for viscous displacement

#---- size of system
      nsys = 3*n

#---- clear and accumulate source-influence matrix
      uvis_mdi=zeros(n,n)
    @inbounds for i = 1:n
      @inbounds for j = 1: n-1
          dx = xi[i] - 0.5*(xi[j+1]+xi[j])
          uvis_mdi[i,j+1] = uvis_mdi[i,j+1] + qopi/(dx*abs(dx))
          uvis_mdi[i,j  ] = uvis_mdi[i,j  ] - qopi/(dx*abs(dx))
        end
      end

#    print("15, :", uvis_mdi[15,:])
#------------------------------------------------------------
#---- global Newton iteration
     npass = 20
   @inbounds for ipass = 1:npass 

#---- clear system matrix and righthand side
      rsys = zeros(nsys)
      asys = zeros(nsys,nsys)
      

#---- first point variables are held frozen... put 1's on diagonal
      asys[1,1] = 1.0
      asys[2,2] = 1.0
      asys[3,3] = 1.0


#---- set current uvis corresponding to current mdi
      dsi[1] = 0.
    
    @inbounds for i = 2: n
        uvis[i] = uinv[i]
      @inbounds for j = 1: n
          uvis[i] = uvis[i] + uvis_mdi[i,j]*mdi[j]
        end
      end

#---- sweep downstream to set up BL equations
for i = 2:n
        
        simi = xi[i-1] == 0.0
#c      lami = simi
        lami = false
        wake = i  > ite

        x  = xi[i]
        b  = bi[i]
        rn = rni[i]
        th = thi[i]
        md = mdi[i]
        ue = uei[i]

        xm  = xi[i-1]
        bm  = bi[i-1]
        rnm = rni[i-1]
        thm = thi[i-1]
        mdm = mdi[i-1]
        uem = uei[i-1]


        bp = b * 0.25/pi
        mpu = 0.5*md/(pi*ue)
        if(rn <= 1.0e-6)
         ds    =  md/(ue*b)
         ds_md = 1.0/(ue*b)
         ds_ue = -ds/ue
        else
         ds    =    (sqrt(bp^2 + mpu) - bp) / rn
         ds_md = 0.5/sqrt(bp^2 + mpu) * ( mpu/md) / rn
         ds_ue = 0.5/sqrt(bp^2 + mpu) * (-mpu/ue) / rn
        end

        if(simi)
         dsm = 0.
         dsm_mdm = 0.
         dsm_uem = 0.
        else
         if(rnm <= 1.0e-6)
          dsm     =  mdm/(uem*bm)
          dsm_mdm = 1.0/(uem*bm)
          dsm_uem = -dsm/uem
         else
          bpm = bm * 0.25/pi
          mpum = 0.5*mdm/(pi*uem)
          dsm     =    (sqrt(bpm^2 + mpum) - bpm) / rnm
          dsm_mdm = 0.5/sqrt(bpm^2 + mpum) * ( mpum/mdm) / rnm
          dsm_uem = 0.5/sqrt(bpm^2 + mpum) * (-mpum/uem) / rnm
         end
        end
            
        (h , h_th, h_ds,
        hk, hk_th, hk_ds, hk_ue,
        hc, hc_th, hc_ds, hc_ue,
        hs, hs_th, hs_ds, hs_ue,
        cf, cf_th, cf_ds, cf_ue,
        di, di_th, di_ds, di_ue ) = blvar(simi,lami,wake, Reyn,Mach, fexcr,x, th ,ds ,ue)

        direct = true
        aa,bb,rr = blsys(simi,lami,wake,direct, Mach, uinv[i],hksep,
                      x,b,rn,th,ds,ue, 
                      h , h_th, h_ds,
                      hk, hk_th, hk_ds, hk_ue,
                      hc, hc_th, hc_ds, hc_ue,
                      hs, hs_th, hs_ds, hs_ue,
                      cf, cf_th, cf_ds, cf_ue,
                      di, di_th, di_ds, di_ue,
                      xm,bm,rnm,thm,dsm,uem, 
                      hm , hm_thm, hm_dsm,
                      hkm, hkm_thm, hkm_dsm, hkm_uem,
                      hcm, hcm_thm, hcm_dsm, hcm_uem,
                      hsm, hsm_thm, hsm_dsm, hsm_uem,
                      cfm, cfm_thm, cfm_dsm, cfm_uem,
                      dim, dim_thm, dim_dsm, dim_uem)

#------ put BL equations of small 3x3 system into big system
      @inbounds for k = 1:2
          ksys = 3*(i-1) + k
          rsys[ksys] = rr[k]

          r_th = aa[k,1]
          r_ds = aa[k,2]
          r_ue = aa[k,3]

          r_thm = bb[k,1]
          r_dsm = bb[k,2]
          r_uem = bb[k,3]

          lsys = 3*(i-1) + 1
          asys[ksys,lsys  ] = r_th
          asys[ksys,lsys-3] = r_thm

          lsys = 3*(i-1) + 2
          asys[ksys,lsys  ] = r_ds *ds_md
          asys[ksys,lsys-3] = r_dsm*dsm_mdm

          lsys = 3*(i-1) + 3
          asys[ksys,lsys  ] = r_ds *ds_ue   + r_ue
          asys[ksys,lsys-3] = r_dsm*dsm_uem + r_uem
        end

#------ set up 3rd ue equation, ue = uvis
        k = 3
        ksys = 3*(i-1) + k
        rsys[ksys] = ue - uvis[i]

      @inbounds for j = 1:n
          lsys = 3*(j-1) + 2
          asys[ksys,lsys] = -uvis_mdi[i,j]
        end

        lsys = 3*(i-1) + 3
        asys[ksys,lsys] = asys[ksys,lsys] + 1.0

#------ also store dependent variables for returning
        tsi[i] = hs*th
        dci[i] = hc*th
        cfi[i] = cf
        cdi[i] = di*hs/2.0
        cti[i] = 0.03 * 0.5*hs*((hk-1.0)/hk)^2
        hki[i] = hk

#------ set dependent i-1 variables for next interval
        hm = h
        hm_thm = h_th
        hm_dsm = h_ds
  
        hkm = hk
        hkm_thm = hk_th
        hkm_dsm = hk_ds
        hkm_uem = hk_ue
  
        hsm = hs
        hsm_thm = hs_th
        hsm_dsm = hs_ds
        hsm_uem = hs_ue
  
        hcm = hc
        hcm_thm = hc_th
        hcm_dsm = hc_ds
        hcm_uem = hc_ue
  
        if(i >= ite)
         cfm     = 0.
         cfm_thm = 0.
         cfm_dsm = 0.
         cfm_uem = 0.
        else
         cfm     = cf
         cfm_thm = cf_th
         cfm_dsm = cf_ds
         cfm_uem = cf_ue
        end
  
        dim     = di
        dim_thm = di_th
        dim_dsm = di_ds
        dim_uem = di_ue

end #end i = 2:n loop


#      do i = 1, n
#        k1 = 3*(i-1) + 1
#        k2 = 3*(i-1) + 2
#        k3 = 3*(i-1) + 3
#        write(*,'(1x,i4,3e12.4,3x,3g15.7)')
#     &    i, rsys(k1),rsys(k2),rsys(k3), uei[i], uvis[i], uei[i]-uvis[i]
#      end

#---- solve Newton system
      rsys = asys\rsys #call gaussn(kdim,nsys,asys,rsys,1)

#---- set Newton changes, set max change >dmax, and set limiter rlx
        dmax = 0.
        rlx = 1.0
      @inbounds for i = 2:n
          b  = bi[i]
          th = thi[i]
          md = mdi[i]
          ue = uei[i]

          kth = 3*(i-1) + 1
          kmd = 3*(i-1) + 2
          kue = 3*(i-1) + 3 
          dthi[i] = -rsys[kth]
          dmdi[i] = -rsys[kmd]
          duei[i] = -rsys[kue]

          bp = b * 0.25/pi
          mpu = 0.5*md/(pi*ue)
          ds    =     sqrt(bp^2 + mpu) - bp
          ds_md = 0.5/sqrt(bp^2 + mpu) * ( mpu/md)
          ds_ue = 0.5/sqrt(bp^2 + mpu) * (-mpu/ue)

          ddsi[i] = ds_md*dmdi[i] + ds_ue*duei[i]

          if(rlx*dthi[i] >  1.6*th) rlx =  1.6*th/dthi[i] end
          if(rlx*dthi[i] < -0.6*th) rlx = -0.6*th/dthi[i] end
          if(rlx*ddsi[i] >  2.5*ds) rlx =  2.5*ds/ddsi[i] end
          if(rlx*ddsi[i] < -0.4*ds) rlx = -0.4*ds/ddsi[i] end
          if(rlx*duei[i] >  0.2*ue) rlx =  0.2*ue/duei[i] end
          if(rlx*duei[i] < -0.1*ue) rlx = -0.1*ue/duei[i] end
            
#           println(i, ", ", rlx)
          dmax = max( dmax , abs(dthi[i])/th , abs(ddsi[i])/ds , abs(duei[i])/ue)
#           println("dmax = $dmax, ue = $ue")
        end

#------ perform Newton updated, with limiter if rlx < 1
      @inbounds for i = 2:n
          thi[i] = thi[i] + rlx*dthi[i]
          dsi[i] = dsi[i] + rlx*ddsi[i]
          uei[i] = uei[i] + rlx*duei[i]
          mdi[i] = uei[i]*dsi[i]*(bi[i] + 2.0*pi*dsi[i]*rni[i])
        end

        if(dmax <  eps) break end

end


#---- integrate running dissipation
      phi[1] = 0.
    @inbounds for i = 2: n
        x  = xi[i]
        xm = xi[i-1]

        dib  = cdi[i]  *rhi[i]  *uei[i]^3 * (bi[i]   + 2.0*pi*dsi[i]  *rni[i])
        dibm = cdi[i-1]*rhi[i-1]*uei[i-1]^3 * (bi[i-1] + 2.0*pi*dsi[i-1]*rni[i-1])

        phi[i] = phi[i-1] + 0.5*(dib + dibm) * (x - xm)
      end

      return uei, dsi, thi, tsi, dci, cfi, cdi, cti, hki, phi 

end # blax

