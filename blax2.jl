
function blax2(ndim, n,ite, xᵢ,bi,rni, uinv, Reyn, Mach, fexcr )
#-----------------------------------------------------------------
#     Axᵢsymmetric boundary layer + wake calculation routine.
#     Uses specified inviscid velocity, corrects for viscous
#     displacement to allow calculation of separated flow.
#
#  Inputs
#  ------
#    ndim    physical array dimension
#    n       number of BL+wake points
#    ite     index of trailing edge point, start of wake
#    xᵢ(.)   arc length array (BL coordinate)
#    bi(.)   lateral width of BL (body perimeter, zero for wake)
#             = 1 for 2D
#    rni(.)  dr/dn
#             = 0 for 2D
#    uinv(.) inviscid velocity
#    Reyn    Reynolds number,  ρ_ref u_ref l_ref / mu_ref
#    Mach    Mach number    ,  u_ref / a_ref
#    fexcr   excrescence multiplier, applied to wall Cf 
#             = 1 for smooth wall
#
#  Assumed units for all quantities:
#    l_ref   same unit as used for input xᵢ,bi
#    u_ref   freestream velocity
#    a_ref   freestream speed of sound
#    ρ_ref freestream density
#    mu_ref  freestream viscosity
#    
#  Outputs
#  -------
#    uₑᵢ[i]  edge velocity, ( = uinv[i] + displacement correction )
#    δᵢ[i]  displacement θᵢckness
#    θᵢ[i]  momentum θᵢckness
#    θ∗ᵢ[i]  kinetic energy θᵢckness
#    δᵢ⁺⁺[i]  density flux θᵢckness
#    cfi[i]  skin friction coefficient, normalized with local ρ,u
#    cdi[i]  dissipation coefficient  , normalized with local ρ,u
#    cti[i]  max shear-stress coefficient, normalized with local ρ,u
#    hki[i]  kinematic shape parameter
#    ϕ[i]  integrated dissipation
#
#
#  Other outputs of interest can be computed as follows.
#  These are in units of l_ref, ρ_ref, u_ref
#
#  Effective perimeter  : beff  =  bi  +  2 π δᵢ rni
#  Edge density         : ρᵢ = (1 + 0.5*(ɣ-1)*Mach^2*(1.0-uₑᵢ^2))^(1/(ɣ-1))
#  Total mass defect    : mdef =  ρᵢ uₑᵢ   δᵢ beff
#  Total mom. defect    : Pdef =  ρᵢ uₑᵢ^2 θᵢ beff
#  Total KE defect      : Edef =  ρᵢ uₑᵢ^3 θ∗ᵢ beff / 2
#  Wall shear force/span: tw b =  ρᵢ uₑᵢ^2 cfi beff / 2
#  Dissipation integral : Diss =  ρᵢ uₑᵢ^3 cdi beff
#
#  Body profile drag Dp is the far-downstream momentum defect P∞,
#  best obtained by applying Squire-Young to the last wake point i = n :
#
#   Pend = ρᵢ*uₑᵢ^2 * θᵢ * beff
#   Hend = δᵢ/θᵢ
#   H∞ = 1.0 + (ɣ-1)*Mach^2
#   Havg = 0.5*(Hend+H∞)
#   P∞ = Pend * uₑᵢ^Havg  =  Dp
#
#-----------------------------------------------------------------

#      real xᵢ(ndim), bi(ndim), rni(ndim), uinv(ndim),

      uₑᵢ  = zeros(ndim)
      ρᵢ  = zeros(ndim)
      δᵢ  = zeros(ndim)
      θᵢ  = zeros(ndim)
      θ∗ᵢ  = zeros(ndim)
      δᵢ⁺⁺  = zeros(ndim)
      cfi  = zeros(ndim)
      cdi  = zeros(ndim)
      cti  = zeros(ndim)
      hki  = zeros(ndim)
      ϕ  = zeros(ndim)

      aa = zeros(3,3)
      bb = zeros(3,3)
      rr = zeros(3)
    
      hm,  hm_thm, hm_dsm = 1e20*ones(3)
      hkm, hkm_thm, hkm_dsm, hkm_uem = 0*1e20*ones(4)
      hcm, hcm_thm, hcm_dsm, hcm_uem = 1e20*ones(4)
      hsm, hsm_thm, hsm_dsm, hsm_uem = 1e20*ones(4)
      cfm, cfm_thm, cfm_dsm, cfm_uem = 1e20*ones(4)
      dim, dim_thm, dim_dsm, dim_uem = 1e20*ones(4)
      h , ∂h∂θ, ∂h∂δ,
      hk, hk_th, hk_ds, hk_ue,
      hc, hc_th, hc_ds, hc_ue,
      hs, hs_th, hs_ds, hs_ue,
      cf, cf_th, cf_ds, cf_ue,
      di, di_th, di_ds, di_ue = 1e20*ones(23)


      idim = 60
      mdi      = zeros(idim)
      uvis     = zeros(idim)
      uvis_mdi = zeros(idim,idim)
      dθᵢ     = zeros(idim)
      dmdi     = zeros(idim)
      duₑᵢ     = zeros(idim)
      dδᵢ     = zeros(idim)

      kdim=3*n #idim
      asys = zeros(kdim,kdim)
      rsys = zeros(kdim)

      simi, lami, wake, direct = true, true, true, true

#---- π  and  1/(4 π)
      π¼  = 0.07957747154594766788444188168625718
      ɣ   = 1.4

      hksep = 2.9

      eps = 1.0e-6

      if(n > idim) 
       println("BLAX: Local array overflow.  Increase idim to", n)
       quit()
      end

      gmi = ɣ - 1.0

#---- initialize ue to inviscid uinv, and initialize ρe
      for i = 1: n
        uₑᵢ[i] = uinv[i]

        trat = 1.0 + 0.5*gmi*Mach^2 * (1.0-uₑᵢ[i]^2)
        ρᵢ[i] = trat^(1.0/gmi)
      end

#---- first point is not calculated if xᵢ=0 there
      if(xᵢ[1] == 0.0) 
       θᵢ[1] = 0.
       δᵢ[1] = 0.
       mdi[1] = 0.
       ϕ[1] = 0.
      end

# =============================================================================
#---- sweep downstream to initialize BL via standard BL march, 
#-     fudging by specifying plausible Hk if separation occurs

direct = true
for i = 2: n #BL march loop

        simi = xᵢ[i-1] == 0.0

#c      lami = simi     # laminar   similarity station
        lami = false  # turbulent similarity station

        wake = i  > ite

        x = xᵢ[i]
        b = bi[i]
        rn = rni[i]

        xm = xᵢ[i-1]
        bm = bi[i-1]
        rnm = rni[i-1]

#------ set i-1 station, and initialize station i for iteration
        if(simi) 
         uem = 0.
         thm = 0.
         dsm = 0.

         ue = uₑᵢ[i]
         rex = ue*x * Reyn
         th = 0.4*x / sqrt(rex)
         ds = th * 2.0

        else
         uem = uₑᵢ[i-1]
         thm = θᵢ[i-1]
         dsm = δᵢ[i-1]

         th = θᵢ[i-1]
         ds = δᵢ[i-1]
         if(direct) 
#-------- previous station was direct... use specified ue as initial guess
          ue = uₑᵢ[i]
         else
#-------- previous station was inverse... use previous-station ue as initial guess
          ue = uₑᵢ[i-1]
         end

        end
 
#------ always try direct mode first
        direct = true
        hkprev = 0.

        for iter = 1: 20 #Newton iteration
          
          (h , ∂h∂θ, ∂h∂δ,
          hk, hk_th, hk_ds, hk_ue,
          hc, hc_th, hc_ds, hc_ue,
          hs, hs_th, hs_ds, hs_ue,
          cf, cf_th, cf_ds, cf_ue,
          di, di_th, di_ds, di_ue ) = blvar(simi,lami,wake, Reyn,Mach, fexcr, x, th ,ds ,ue) 
            
          if(iter==1) 
           hkprev = hk
          else
           if(hk > hksep)  # .and. hk > hkprev) 
#---------- Hk limit exceeded... switch to inverse mode for θᵢs point
            direct = false
           end
          end

#-------- set up 2-point differenced BL equation system
          aa, bb, rr = blsys(simi,lami,wake,direct, Mach, uinv[i],hksep,
                      x,b,rn,th,ds,ue, 
                      h , ∂h∂θ, ∂h∂δ,
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
  
          rr = aa\rr
          dth = -rr[1]
          dds = -rr[2]
          due = -rr[3]

          rlx = 1.0
            
            if(rlx*dth >  1.6*th) rlx =  1.6*th/dth ; end
            if(rlx*dth < -0.6*th) rlx = -0.6*th/dth ; end
            if(rlx*dds >  2.5*ds) rlx =  2.5*ds/dds ; end
            if(rlx*dds < -0.4*ds) rlx = -0.4*ds/dds ; end
            if(rlx*due >  0.2*ue) rlx =  0.2*ue/due ; end
            if(rlx*due < -0.1*ue) rlx = -0.1*ue/due ; end

          dmax = max( abs(dth)/th , abs(dds)/ds , abs(due)/ue )

#          if(iter==1) write(*,*)
#          rt = ue*th*Reyn
#          write(*,1200) iter, direct, dmax, rlx, ue, hk, rt, cf, di
# 1200     format(1x,i4, l2, e12.4, 9g13.5)

          th = th + rlx*dth
          ds = ds + rlx*dds
          ue = ue + rlx*due

          if(dmax < eps) 
            break; 
          end
       end #newton iteration


      uₑᵢ[i] = ue
      δᵢ[i] = ds
      θᵢ[i] = th
      θ∗ᵢ[i] = hs*th
      δᵢ⁺⁺[i] = hc*th
      cfi[i] = cf
      cdi[i] = di*hs/2.0
      cti[i] = 0.03 * 0.5*hs*((hk-1.0)/hk)^2
      hki[i] = hk

      dib  = cdi[i]  *ρᵢ[i]  *uₑᵢ[i]^3   * (b  + 2.0*π*ds *rn )
      dibm = cdi[i-1]*ρᵢ[i-1]*uₑᵢ[i-1]^3 * (bm + 2.0*π*dsm*rnm)
      ϕ[i] = ϕ[i-1] + 0.5*(dib + dibm) * (x - xm)
      mdi[i] = ue*ds*(b + 2.0*π*ds)


      thm = th
      dsm = ds
      hsm = hs
      uem = ue

      hm = h
      hm_thm = ∂h∂θ
      hm_dsm = ∂h∂δ

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

      if(i>=ite) 
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

end # BL march loop


# =============================================================================
#---- perform global viscous/inviscid iteration with source model 
#-     for viscous displacement

#---- size of system
      nsys = 3*n

#---- clear and accumulate source-influence matrix
      for i = 1: n
        for j = 1: n
          uvis_mdi[i,j] = 0.
        end
      end
    
      for i = 1: n
        for j = 1: n-1
          dx = xᵢ[i] - 0.5*(xᵢ[j+1]+xᵢ[j])
#         dm = mdi(j+1) - mdi[j]
#         du = dm * π¼ / (dx*abs(dx))
#         duₑᵢ = duₑᵢ + du
          uvis_mdi[i,j+1] = uvis_mdi[i,j+1] + π¼/(dx*abs(dx))
          uvis_mdi[i,j  ] = uvis_mdi[i,j  ] - π¼/(dx*abs(dx))
        end
      end

#      i = 15
#      for j = 1: n
#         write(*,*) i, j, uvis_mdi[i,j]
#      end

#------------------------------------------------------------
#---- global Newton iteration
      npass = 20
      for ipass = 1: npass #Newton iteration

#---- clear system matrix and righthand side
      for k = 1: nsys
        rsys[k] = 0.
        for l = 1: nsys
          asys[k,l] = 0.
        end
      end

#---- first point variables are held frozen... put 1's on diagonal
      asys[1,1] = 1.0
      asys[2,2] = 1.0
      asys[3,3] = 1.0


#---- set current uvis corresponding to current mdi
      δᵢ[1] = 0.
      for i = 2: n
        uvis[i] = uinv[i]
        for j = 1: n
          uvis[i] = uvis[i] + uvis_mdi[i,j]*mdi[j]
        end
      end

#---- sweep downstream to set up BL equations
      for i = 2: n  # Set up BL loop
        simi = xᵢ[i-1] == 0.0

#c      lami = simi
        lami = false

        wake = i  > ite

        x  = xᵢ[i]
        b  = bi[i]
        rn = rni[i]
        th = θᵢ[i]
        md = mdi[i]
        ue = uₑᵢ[i]

        xm  = xᵢ[i-1]
        bm  = bi[i-1]
        rnm = rni[i-1]
        thm = θᵢ[i-1]
        mdm = mdi[i-1]
        uem = uₑᵢ[i-1]

#c      md = ue*ds*(b + 2.0*π*ds*rn)

#c      md/(ue*2*π) = ds*(b/2π + ds*rn)
#c      0.5 rn ds^2 + (b/4π)*ds - md/(4π ue) = 0
#c      ds = ( sqrt(bp^2 + 0.5*md/(π*ue)) - bp ) / rn
        bp = b * 0.25/π
        mpu = 0.5*md/(π*ue)
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
          bpm = bm * 0.25/π
          mpum = 0.5*mdm/(π*uem)
          dsm     =    (sqrt(bpm^2 + mpum) - bpm) / rnm
          dsm_mdm = 0.5/sqrt(bpm^2 + mpum) * ( mpum/mdm) / rnm
          dsm_uem = 0.5/sqrt(bpm^2 + mpum) * (-mpum/uem) / rnm
         end
        end

#        if(.not.simi) 
#        call blvar(simi,lami,wake, Reyn,Mach, fexcr,
#     &                 xm, thm ,dsm ,uem , 
#     &                 hm, hm_thm, hm_dsm,
#     &                 hkm, hkm_thm, hkm_dsm, hkm_uem,
#     &                 hcm, hcm_thm, hcm_dsm, hcm_uem,
#     &                 hsm, hsm_thm, hsm_dsm, hsm_uem,
#     &                 cfm, cfm_thm, cfm_dsm, cfm_uem,
#     &                 dim, dim_thm, dim_dsm, dim_uem )
#        end

        (h , ∂h∂θ, ∂h∂δ,
        hk, hk_th, hk_ds, hk_ue,
        hc, hc_th, hc_ds, hc_ue,
        hs, hs_th, hs_ds, hs_ue,
        cf, cf_th, cf_ds, cf_ue,
        di, di_th, di_ds, di_ue ) = blvar(simi,lami,wake, Reyn,Mach, fexcr,x, th ,ds ,ue)
        direct = true
            
        aa,bb,rr = blsys(simi,lami,wake,direct, Mach, uinv[i],hksep,
                      x,b,rn,th,ds,ue, 
                      h , ∂h∂θ, ∂h∂δ,
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
        for k = 1: 2
          ksys = 3*(i-1)+ k
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

        for j = 1: n
          lsys = 3*(j-1) + 2
          asys[ksys,lsys] = -uvis_mdi[i,j]
        end

        lsys = 3*(i-1) + 3
        asys[ksys,lsys] = asys[ksys,lsys] + 1.0

#------ also store dependent variables for returning
        θ∗ᵢ[i] = hs*th
        δᵢ⁺⁺[i] = hc*th
        cfi[i] = cf
        cdi[i] = di*hs/2.0
        cti[i] = 0.03 * 0.5*hs*((hk-1.0)/hk)^2
        hki[i] = hk

#------ set dependent i-1 variables for next interval
        hm = h
        hm_thm = ∂h∂θ
        hm_dsm = ∂h∂δ
  
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
  
        if(i>=ite) 
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

        end #end BL loop


#      for i = 1: n
#        k1 = 3*[i-1] + 1
#        k2 = 3*[i-1] + 2
#        k3 = 3*[i-1] + 3
#        write(*,'(1x,i4,3e12.4,3x,3g15.7)')
#     &    i, rsys(k1),rsys(k2),rsys(k3), uₑᵢ[i], uvis[i], uₑᵢ[i]-uvis[i]
#      end

#---- solve Newton system
#      call gaussn(kdim,nsys,asys,rsys,1)

rsys = asys\rsys

#---- set Newton changes, set max change dmax, and set limiter rlx
        dmax = 0.
        rlx = 1.0
        for i = 2: n
          b  = bi[i]
          th = θᵢ[i]
          md = mdi[i]
          ue = uₑᵢ[i]

          kth = 3*(i-1) + 1
          kmd = 3*(i-1) + 2
          kue = 3*(i-1) + 3 
          dθᵢ[i] = -rsys[kth]
          dmdi[i] = -rsys[kmd]
          duₑᵢ[i] = -rsys[kue]

          bp = b * 0.25/π
          mpu = 0.5*md/(π*ue)
          ds    =     sqrt(bp^2 + mpu) - bp
          ds_md = 0.5/sqrt(bp^2 + mpu) * ( mpu/md)
          ds_ue = 0.5/sqrt(bp^2 + mpu) * (-mpu/ue)

#          ds    =  md/(b*ue)
#          ds_md = 1.0/(b*ue)
#          ds_ue = -ds/ue

          dδᵢ[i] = ds_md*dmdi[i] + ds_ue*duₑᵢ[i]

          if(rlx*dθᵢ[i] >  1.6*th) rlx =  1.6*th/dθᵢ[i]; end 
          if(rlx*dθᵢ[i] < -0.6*th) rlx = -0.6*th/dθᵢ[i]; end
          if(rlx*dδᵢ[i] >  2.5*ds) rlx =  2.5*ds/dδᵢ[i]; end
          if(rlx*dδᵢ[i] < -0.4*ds) rlx = -0.4*ds/dδᵢ[i]; end
          if(rlx*duₑᵢ[i] >  0.2*ue) rlx =  0.2*ue/duₑᵢ[i]; end
          if(rlx*duₑᵢ[i] < -0.1*ue) rlx = -0.1*ue/duₑᵢ[i]; end

          dmax = max( dmax , abs(dθᵢ[i])/th , 
	abs(dδᵢ[i])/ds , 
	abs(duₑᵢ[i])/ue   )
        end

#------ perform Newton updated, with limiter if rlx < 1
        for i = 2: n
          θᵢ[i] = θᵢ[i] + rlx*dθᵢ[i]
          δᵢ[i] = δᵢ[i] + rlx*dδᵢ[i]
          uₑᵢ[i] = uₑᵢ[i] + rlx*duₑᵢ[i]
          mdi[i] = uₑᵢ[i]*δᵢ[i]*(bi[i] + 2.0*π*δᵢ[i]*rni[i])
        end

#        for i = 2: n
#          θᵢ[i] = θᵢ[i] + rlx*dθᵢ[i]
#cc        δᵢ[i] = δᵢ[i] + rlx*dδᵢ[i]
#          uₑᵢ[i] = uₑᵢ[i] + rlx*duₑᵢ[i]
#          mdi[i] = mdi[i] + rlx*dmdi[i]
#          δᵢ[i] = mdi[i]/(bi[i]*uₑᵢ[i])
#        end

#        if(ipass==1) write(*,*)
#        write(*,2200) ipass, dmax, rlx
# 2200   format(1x,i4, e12.4, f8.4)

            if(dmax < eps) break; end

        end #End newton iteration



#---- integrate running dissipation
      ϕ[1] = 0.
      for i = 2: n
        x  = xᵢ[i]
        xm = xᵢ[i-1]

        dib  = cdi[i]  *ρᵢ[i]  *uₑᵢ[i]^3   * (bi[i]   + 2.0*π*δᵢ[i]  *rni[i])
        dibm = cdi[i-1]*ρᵢ[i-1]*uₑᵢ[i-1]^3 * (bi[i-1] + 2.0*π*δᵢ[i-1]*rni[i-1])

        ϕ[i] = ϕ[i-1] + 0.5*(dib + dibm) * (x - xm)
      end

      return uₑᵢ, δᵢ, θᵢ, θ∗ᵢ, δᵢ⁺⁺, cfi, cdi, cti, hki, ϕ
      end # blax
