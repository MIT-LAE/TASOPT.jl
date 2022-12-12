"""
    trefftz1(nsurf, npout, npinn, npimg, Sref, bref,
    b,bs,bo,bop, zcent, 
    po,gammat,gammas, fLo,ktip, Lspec,CLsurfsp,)

Trefftz plane routine for the induced drag computation.

Inputs
- `nsurf::Integer`: number of surfaces (typically wing and horizontal tail)
- `npout`::Integer: number of spanwise intervals (outer panel)
- `npinn`::Integer:  "     "       (inner panel)
- `npimg`::Integer:  "     "       (image inside fuselage)
- `b::Float64`: span
- `bs::Float64`: wing-break span.
- `bo::Float64`: wing-root span.
- `bop::Float64`: span of wing-root streamline in Trefftz plane
- `zcent`:
- `gammat`:
- `gammas`:
- `fLo`: wing root load adjustment factors.
- `ktip`: wing tip load adjustment factors.
- `Lspec`:
- `CLsurfsp`:

See section A 2.13 of TASOPT docs.
"""
function trefftz1(nsurf, npout, npinn, npimg,
	Sref, bref,
	b,bs,bo,bop, zcent,
	po,gammat,gammas, fLo,ktip,
	Lspec,CLsurfsp,)
#	CLsurf,CL,CD,spanef,
#	idim,ifrst,ilast,
#	yc,zc,gc,vc,wc,vnc, ycp,zcp)

#      implicit real(a-h,l-z)
#
#      integer nsurf, npout(nsurf), npinn(nsurf), npimg(nsurf)
#      real b(nsurf), bs(nsurf), bo(nsurf), bop(nsurf)
#      real zcent(nsurf)
#      real po(nsurf)
#      real gammat(nsurf), gammas(nsurf)
#      integer ktip
#      real CLsurfsp(nsurf), CLsurf(nsurf)
#      real CL, CD
#      logical Lspec
#
#      real yc(idim), ycp(idim),
#	zc(idim), zcp(idim)
#      real gc(idim), vc(idim), wc(idim), vnc(idim)
#      integer ifrst(nsurf), ilast(nsurf)
#
#
#      parameter (jdim=360)
#      real t(jdim), 
#	y(jdim), yp(jdim),
#	z(jdim), zp(jdim)
#      real gw(jdim)
#
#      common  /com_lu/ ilu
#
#      data π / 3.1415926535897932384626 /

#---- center resolution increase factor (0..1)
#     bunch = 0.75 
      bunch =  0.5 
#c    bunch =  0. 
      
      idim::Int = 360
      jdim::Int = 360

      ifrst = zeros(Int, nsurf)
      ilast = zeros(Int, nsurf)
      
      t     = zeros(Float64, jdim)
      y     = zeros(Float64, jdim)
      yp    = zeros(Float64, jdim)
      z     = zeros(Float64, jdim)
      zp    = zeros(Float64, jdim)
      gw    = zeros(Float64, jdim)
      
      yc    = zeros(Float64, jdim)
      ycp   = zeros(Float64, jdim)
      zc    = zeros(Float64, jdim)
      zcp   = zeros(Float64, jdim)
      gc    = zeros(Float64, jdim)
      vc    = zeros(Float64, jdim)
      wc    = zeros(Float64, jdim)
      vnc   = zeros(Float64, jdim)
      
      CLsurf= zeros(Float64, nsurf)
      
      isum = 0
      @inbounds for  isurf = 1: nsurf
        isum = isum + npout[isurf] + npinn[isurf] + npimg[isurf] + 1
      end

      if(isum > idim) 
	      println("TREFFTZ: Passed array overflow. Increase idim to ",isum)
       quit()
      end

      if(isum > jdim) 
	      println("TREFFTZ: Local array overflow. Increase jdim to ", isum)
       quit()
      end


      i = 0
@inbounds for  isurf = 1: nsurf

#---- eta at center, side-of-body, wing break, tip
      e0 = 0.0
      eo = bo[isurf]/b[isurf]
      es = bs[isurf]/b[isurf]
      e1 = 1.0

      eop = bop[isurf]/b[isurf]

      t0 = 1.0
      to = acos(eo)/(0.5*π)
      ts = acos(es)/(0.5*π)
      t1 = 0.0

      if(bunch > 0.0) 
       to = (1.0+bunch - sqrt((1.0+bunch)^2 - 4.0*bunch*to))*0.5/bunch
       ts = (1.0+bunch - sqrt((1.0+bunch)^2 - 4.0*bunch*ts))*0.5/bunch
      end

      y0 = 0.
      yo = 0.5*bo[isurf]
      ys = 0.5*bs[isurf]
      y1 = 0.5*b[isurf]

      yop = 0.5*bop[isurf]

      z0 = zcent[isurf]
      zo = zcent[isurf]
      zs = zcent[isurf]
      z1 = zcent[isurf]

      g0 = po[isurf]
      go = po[isurf]
      gs = po[isurf]*gammas[isurf]
      g1 = po[isurf]*gammat[isurf]


#      k0 = 1
#      ko = 1 + int( float(n[isurf])*(to-t0)/(t1-t0) + 0.5 )
#      ks = 1 + int( float(n[isurf])*(ts-t0)/(t1-t0) + 0.5 )
#      k1 = 1 + n[isurf]

      k0 = 1
      ko = 1 + npimg[isurf]
      ks = 1 + npimg[isurf] + npinn[isurf]
      k1 = 1 + npimg[isurf] + npinn[isurf] + npout[isurf]

#      write(*,*) y0, yo, ys, y1
#      write(*,*) t0, to, ts, t1
#      write(*,*) 1, ko, ks, n[isurf]+1

      i = i+1
      ifrst[isurf] = i
      t[i] = t0
      y[i] = y0
      z[i] = z0

      yp[i] = y0
      zp[i] = z0

#---- set points over fuselage
      @inbounds for  k = k0+1: ko
        i = i+1

        fk = float(k-k0)/float(ko-k0)
        t[i] = t0*(1.0-fk) + to*fk
        tc = 0.5*(t[i-1]+t[i])

        e  = cos(0.5*π*(t[i] + bunch*t[i]*(1.0-t[i])))
        ec = cos(0.5*π*(tc   + bunch*tc  *(1.0-tc  )))
        if(eo-e0 == 0.0) 
         fi = 1.0
         fc = 0.5
        else
         fi = (e -e0)/(eo-e0)
         fc = (ec-e0)/(eo-e0)
        end

        y[i] = y0*(1.0-fi) + yo*fi
        z[i] = z0*(1.0-fi) + zo*fi

        yc[i-1] =  y0*(1.0-fc) + yo*fc
        zc[i-1] =  z0*(1.0-fc) + zo*fc
        gc[i-1] = (g0*(1.0-fc) + go*fc) * sqrt(1.0-ec^ktip)

        yexp = (eo/eop)^2

        yp[i] = yop * (y[i]/yo)^yexp
        zp[i] = z[i]


        ycp[i-1] = yop * (yc[i-1]/yo)^yexp
        zcp[i-1] = zc[i-1]
      end

      @inbounds for  k = ko+1: ks
        i = i+1

        fk = float(k-ko)/float(ks-ko)
        t[i] = to*(1.0-fk) + ts*fk
        tc = 0.5*(t[i-1]+t[i])

        e  = cos(0.5*π*(t[i] + bunch*t[i]*(1.0-t[i])))
        ec = cos(0.5*π*(tc   + bunch*tc  *(1.0-tc  )))
        if(es-eo == 0.0) 
         fi = 1.0
         fc = 0.5
        else
         fi = (e -eo)/(es-eo)
         fc = (ec-eo)/(es-eo)
        end

        y[i] = yo*(1.0-fi) + ys*fi
        z[i] = zo*(1.0-fi) + zs*fi

        yc[i-1] =  yo*(1.0-fc) + ys*fc
        zc[i-1] =  zo*(1.0-fc) + zs*fc
        gc[i-1] = (go*(1.0-fc) + gs*fc) * sqrt(1.0-ec^ktip)

        yp[i] = sqrt(y[i]^2 - yo^2 + yop^2)
        zp[i] = z[i]

        ycp[i-1] = sqrt(yc[i-1]^2 - yo^2 + yop^2)
        zcp[i-1] = zc[i-1]
      end

      @inbounds for  k = ks+1: k1
        i = i+1

        fk = float(k-ks)/float(k1-ks)
        t[i] = ts*(1.0-fk) + t1*fk
        tc = 0.5*(t[i-1]+t[i])

        e  = cos(0.5*π*(t[i] + bunch*t[i]*(1.0-t[i])))
        ec = cos(0.5*π*(tc   + bunch*tc  *(1.0-tc  )))
        if(e1-es == 0.0) 
         fi = 1.0
         fc = 0.5
        else
         fi = (e -es)/(e1-es)
         fc = (ec-es)/(e1-es)
        end
        y[i] = ys*(1.0-fi) + y1*fi
        z[i] = zs*(1.0-fi) + z1*fi

        yc[i-1] =  ys*(1.0-fc) + y1*fc
        zc[i-1] =  zs*(1.0-fc) + z1*fc
        gc[i-1] = (gs*(1.0-fc) + g1*fc) * sqrt(1.0-ec^ktip)

        yp[i] = sqrt(y[i]^2 - yo^2 + yop^2)
        zp[i] = z[i]

        ycp[i-1] = sqrt(yc[i-1]^2 - yo^2 + yop^2)
        zcp[i-1] = zc[i-1]
      end
      ilast[isurf] = i

#---- dummy control point between surfaces
      yc[i] = 0.
      zc[i] = 0.
      gc[i] = 0.

 end # nsurf loop

 ii = ilast[nsurf]

      @inbounds for  isurf = 1: nsurf
        i = ifrst[isurf]
        gw[i] = 0.
        @inbounds for  i = ifrst[isurf]+1: ilast[isurf]-1
          gw[i] = gc[i-1] - gc[i]
        end
        i = ilast[isurf]
        gw[i] = gc[i-1]
      end


      if(Lspec) 
#----- scale circulations to get specified lift for each surface
       @inbounds for  isurf = 1: nsurf
         cltest = 0.
         @inbounds for  i = ifrst[isurf]: ilast[isurf]-1
           dy = yp[i+1] - yp[i]
           cltest = cltest + 2.0*gc[i]*dy * bref/(0.5*Sref)
         end

         gfac = CLsurfsp[isurf]/cltest
#c         write(*,*) isurf, gfac

         @inbounds for  i = ifrst[isurf]: ilast[isurf]
           gc[i] = gc[i]*gfac
           gw[i] = gw[i]*gfac
         end
       end
      end
      

      CL = 0.
      CD = 0.
      @inbounds for  isurf = 1: nsurf

      CLsurf[isurf] = 0.

      @inbounds for  i = ifrst[isurf]: ilast[isurf]-1
        dy = yp[i+1] - yp[i]
        dz = zp[i+1] - zp[i]
        ds = sqrt(dy^2 + dz^2)

        vsum = 0.
        wsum = 0.
        @inbounds for  j = 1: ii
#-------- velocities of wake vortex, at control point
          yb = ycp[i] - yp[j]
          zb = zcp[i] - zp[j]
          rsq = yb^2 + zb^2
          vsum = vsum - gw[j] * zb/rsq
          wsum = wsum + gw[j] * yb/rsq

#-------- velocities of wake vortex image, at control point
          yb = ycp[i] + yp[j]
          zb = zcp[i] - zp[j]
          rsq = yb^2 + zb^2
          vsum = vsum + gw[j] * zb/rsq
          wsum = wsum - gw[j] * yb/rsq
        end
        vc[i] = vsum*bref/(2.0*π)
        wc[i] = wsum*bref/(2.0*π)

        vnc[i] = (dy*wc[i] - dz*vc[i])/ds

        dCL = 2.0*gc[i]       *dy * bref/(0.5*Sref)
        dCD =    -gc[i]*vnc[i]*ds * bref/(0.5*Sref)

        CL = CL + dCL
        CD = CD + dCD

        CLsurf[isurf] = CLsurf[isurf] + dCL
      end

 end #end for loop

      AR = bref^2/Sref
      CDell = CL^2 / (π*AR)

      spanef = CDell / CD

      return CLsurf, CL, CD , spanef
      end
