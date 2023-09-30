# Investigating potential to have wrappers to make code cleaner
struct tfp_workarrays
  t::Vector{Float64}
  y::Vector{Float64}
  yp::Vector{Float64}
  z::Vector{Float64}
  zp::Vector{Float64}
  gw::Vector{Float64}
  yc::Vector{Float64}
  ycp::Vector{Float64}
  zc::Vector{Float64}
  zcp::Vector{Float64}
  gc::Vector{Float64}
  vc::Vector{Float64}
  wc::Vector{Float64}
  vnc::Vector{Float64}
end

struct tfp
  nsurf::Int
  npout::Vector{Int}
  npinn::Vector{Int}
  npimg::Vector{Int}
  Sref::Float64
  bref::Float64
  b::Vector{Float64}
  bs::Vector{Float64}
  bo::Vector{Float64}
  bop::Vector{Float64}
  zcent::Vector{Float64}
  po::Vector{Float64}
  γt::Vector{Float64}
  γs::Vector{Float64}
  CLsurfsp::Vector{Float64}
  ktip::Float64
  workarrays::tfp_workarrays
end

function trefftz(tf::tfp)
  
  trefftz1(tf.nsurf, tf.npout, tf.npinn, tf.npimg,
	tf.Sref, tf.bref,
	tf.b, tf.bs, tf.bo, tf.bop, tf.zcent,
	tf.po,tf.γt,tf.γs, 1, tf.ktip,
	true,tf.CLsurfsp,
  tf.workarrays.t, tf.workarrays.y, tf.workarrays.yp, tf.workarrays.z, tf.workarrays.zp, 
  tf.workarrays.gw, tf.workarrays.yc, tf.workarrays.ycp, tf.workarrays.zc, tf.workarrays.zcp,
  tf.workarrays.gc, tf.workarrays.vc, tf.workarrays.wc, tf.workarrays.vnc)
  ## Note: do not do the following. This will slow the code down. Interestingly seems to have more allocations than the above.
  # [getproperty(tf.workarrays,f) for f in fieldnames(typeof(tf.workarrays))]...)

end

"""
    trefftz1(nsurf, npout, npinn, npimg, 
          Sref, bref, b, bs, bo, bop, 
          zcent, po, gammat, gammas, 
          fLo,ktip, Lspec, CLsurfsp)

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
	po,gammat,gammas, fLo, ktip,
	Lspec,CLsurfsp,t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)

#---- center resolution increase factor (0..1)
#     bunch = 0.75 
      bunch =  0.5 
#     bunch =  0. 

      ifrst = zeros(Int, nsurf)
      ilast = zeros(Int, nsurf)
      
      CLsurf= zeros(Float64, nsurf)
      
      # Calcualte total number of points
      # outboard point + inboard + points within fuselage + dummy between surfaces
      isum = 0
      @inbounds for  isurf = 1: nsurf
        isum = isum + npout[isurf] + npinn[isurf] + npimg[isurf] + 1
      end

      if(isum > idim) 
	      println("TREFFTZ: Passed array overflow. Increase idim to ",isum)
        exit()
      end

      if(isum > jdim) 
	      println("TREFFTZ: Local array overflow. Increase jdim to ", isum)
        exit()
      end


      i = 0
@inbounds for  isurf = 1: nsurf

#---- η at center, side-of-body, wing break, tip
      e0 = 0.0
      eo = bo[isurf]/b[isurf]
      es = bs[isurf]/b[isurf]
      e1 = 1.0

      # ηₒ' fuse non-dim y location that is constricted in the Trefftz plane 
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

      k0 = 1
      ko = 1 + npimg[isurf]
      ks = 1 + npimg[isurf] + npinn[isurf]
      k1 = 1 + npimg[isurf] + npinn[isurf] + npout[isurf]

      i = i+1
      ifrst[isurf] = i
      t[i] = t0
      y[i] = y0
      z[i] = z0

      yp[i] = y0
      zp[i] = z0

#---- set points over fuselage
      @inbounds for  k = k0+1 : ko
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
        gc[i-1] = (g0*(1.0-fc) + go*fc) * sqrt(1.0-ec^ktip) #Eq. 389

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

        yp[i] = sqrt(y[i]^2 - yo^2 + yop^2) #Calculate y'[i] from Eq 391. conservtion of mass effectively
        zp[i] = z[i]

        ycp[i-1] = sqrt(yc[i-1]^2 - yo^2 + yop^2)
        zcp[i-1] = zc[i-1]
      end

      @inbounds for  k = ks+1 : k1
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
#        println("$isurf, $gfac")

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
