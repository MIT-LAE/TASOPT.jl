"""
      induced_drag!(para, wing, htail)

Computes the induced drag via the Trefftz plane. Calls [`_trefftz_analysis`](@ref). Formerly, `cditrp!()`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `para::AbstractArray{Float64}`: Array of `aircraft` model aerodynamic parameters.
      - `wing::TASOPT.Wing`: Wing Structure.
      - `htail::TASOPT.Tail`: Htail Structure.

      **Outputs:**
      - No explicit outputs. Computed induced drag value and span efficiency are saved to `para` of `aircraft` model.

!!! compat "Future Changes"
      In an upcoming revision, an `aircraft` `struct` and auxiliary indices will be passed in lieu of pre-sliced `par` arrays.

"""
function induced_drag!(para, wing, htail)

      CL = para[iaCL]

      if (CL == 0.0) 
	    para[iaCDi]  = 0.
	    para[iaspaneff] = 1.0
       return
      end

      CLhtail = para[iaCLh]*htail.layout.S/wing.layout.S
      # println("CLhtail: $(para[iaCLh]) $(htail.layout.S) $(parg[igS])")
      bref = wing.layout.span
      Sref = wing.layout.S

      Mach = para[iaMach]

      specifies_CL = true

      b        = zeros(Float64, 2)
      bs       = zeros(Float64, 2)
      bo       = zeros(Float64, 2)
      bop      = zeros(Float64, 2)
      zcent    = zeros(Float64, 2)
      gammas   = zeros(Float64, 2)
      gammat   = zeros(Float64, 2)
      po       = zeros(Float64, 2)
      CLsurfsp = zeros(Float64, 2)

      npout = zeros(Float64, 2) # outer panel
      npinn = zeros(Float64, 2) # inner panel
      npimg = zeros(Float64, 2) # image inside fuselage

      #Alternatively can define as b  = [parg[igb], parg[igbh]] for both wing and tail simultaneously 
#---- wing wake parameters
      fLo =  wing.fuse_lift_carryover
#      fLo = 0.0

#---- span, wing-break span, wing-root span
      b[1]  = wing.layout.span
      bs[1] = wing.layout.break_span
      bo[1] = wing.layout.root_span

#---- span of wing-root streamline in Trefftz Plane
      bop[1] = wing.layout.root_span * 0.2

      zcent[1]  = wing.layout.z
      gammas[1] = wing.inboard.λ*para[iarcls]
      gammat[1] = wing.outboard.λ*para[iarclt]
      po[1]     = 1.0
      CLsurfsp[1] = CL - CLhtail

#---- velocity-change fractions at wing spanwise locations due to fuselage flow
      fduo = para[iafduo]
      fdus = para[iafdus]
      fdut = para[iafdut]

#---- horizontal tail wake parameters
      b[2]   = htail.layout.span
      bs[2]  = htail.layout.root_span
      bo[2]  = htail.layout.root_span
      bop[2] = htail.layout.root_span

      zcent[2] = htail.layout.z
      gammas[2] = 1.0
      gammat[2] = htail.outboard.λ
      po[2]     = 1.0
      CLsurfsp[2] = CLhtail


#---- number of surfaces  (wing, horizontal tail)
      nsurf = 2

#---- number of spanwise intervals
# -----------------------------------------------------------------------
# The below # of panels (43 total) gives about 1.8% higher CDi relative to the 
# finest one here with ~360 panels
      npout[1] = 20  # outer panel
      npinn[1] = 6   # inner panel
      npimg[1] = 3   # image inside fuselage
  
      npout[2] = 10  # outer panel
      npinn[2] = 0   # inner panel
      if (bo[2] == 0.0) 
       npimg[2] = 0
      else
       npimg[2] = 2   # image inside fuselage  (or inner panel if T-tail)
      end
# -----------------------------------------------------------------------
# The below # of panels (84 total) gives about 0.60% higher CDi relative to the 
# finest one here with ~360 panels
#     npout[1] = 40  # outer panel
#     npinn[1] = 12  # inner panel
#     npimg[1] = 6   # image inside fuselage
# 
#     npout[2] = 20  # outer panel
#     npinn[2] = 0   # inner panel
#     npimg[2] = 4   # image inside fuselage  (or inner panel if T-tail)
# -----------------------------------------------------------------------
#      npout[1] = 160  # outer panel
#      npinn[1] = 48   # inner panel
#      npimg[1] = 24   # image inside fuselage
# 
#      npout[2] = 80  # outer panel
#      npinn[2] = 0   # inner panel
#      npimg[2] = 16  # image inside fuselage  (or inner panel if T-tail)
# -----------------------------------------------------------------------

      ktip = 16
      #CLsurf = zeros(Float64, nsurf)
      # println("$nsurf, $npout, $npinn, $npimg, 
	# $Sref, $bref,
	# $b,$bs,$bo,$bop, $zcent,
	# $po,$gammat,gammas, $fLo, $ktip,
      # $specifies_CL,$CLsurfsp")


      CLsurf, CLtp, CDtp, sefftp = _trefftz_analysis(nsurf, npout, npinn, npimg, 
	Sref, bref,
	b,bs,bo,bop, zcent,
	po,gammat,gammas, fLo, ktip,
      specifies_CL,CLsurfsp, t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)
      
      # println("$CLsurf, $CLtp, $CDtp, $sefftp")

      para[iaCDi] = CDtp
      para[iaspaneff] = sefftp

      return
end # induced_drag!


# Trefftz analysis functions -------------------
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
  
  _trefftz_analysis(tf.nsurf, tf.npout, tf.npinn, tf.npimg,
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
    _trefftz_analysis(nsurf, npout, npinn, npimg, 
          Sref, bref, b, bs, bo, bop, 
          zcent, po, gammat, gammas, 
          fLo,ktip, specifies_CL, CLsurfsp)

Trefftz plane routine for the induced drag computation of `nsurf` number of surfaces. Formerly, `trefftz1()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `nsurf::Integer`: Number of surfaces (typically wing and horizontal tail).
    - `npout::Integer`: Number of spanwise intervals (outer panel).
    - `npinn::Integer`:  "     "       (inner panel).
    - `npimg::Integer`:  "     "       (image inside fuselage).
    - `b::Float64`: Span.
    - `bs::Float64`: Wing-break span.
    - `bo::Float64`: Wing-root span.
    - `bop::Float64`: Span of wing-root streamline in Trefftz plane
    - `zcent::Vector{Float64}`: Vertical position at centerline for each surface.
    - `gammat::Vector{Float64}`, gammas::Vector{Float64}`: Wing lift distribution "taper" ratios for outer and inner wing sections, respectively.
    - `fLo`: wing root load adjustment factors (currently not implemented).
    - `ktip::Float64`: wing tip load adjustment factors.
    - `specifies_CL::Bool`: Flag for specified lift calculation (scales vorticities to achieve `CLsurfsp` before computing induced drag).
    - `CLsurfsp::Vector{Float64}`: Prescribed surface lift coefficient.

    **Outputs:**
    -`CLsurf::Vector{Float64}`: Lift coefficients for each surface. 
    -`CL::Float64`: Sum of lift coefficients of all surfaces.
    -`CD::Float64`: Sum of induced drag coefficients of all surfaces.
    -`spanef::Float64`: Span efficiency of combined surfaces (``= (CL^2 / (π*AR))/CD``).

See [theory above](@ref trefftz) or Sections 2.14.7 and 3.8.1 of the [TASOPT Technical Desc](@ref dreladocs).
"""
function _trefftz_analysis(nsurf, npout, npinn, npimg,
	Sref, bref,
	b,bs,bo,bop, zcent,
	po, gammat, gammas, fLo, ktip,
	specifies_CL,CLsurfsp,t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)

#---- center resolution increase factor (0..1)
#     bunch = 0.75 
      bunch =  0.5 
#     bunch =  0. 

      ifrst = zeros(Int, nsurf)
      ilast = zeros(Int, nsurf)
      
      CLsurf= zeros(Float64, nsurf)
      
      # Calculate total number of points
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


      i::Int64 = 0

@inbounds for  isurf = 1: nsurf

#---- η at center, side-of-body, wing break, tip
      e0 = 0.0
      eo = bo[isurf]/b[isurf]
      es = bs[isurf]/b[isurf]
      e1 = 1.0

      # ηₒ' fuse non-dim y location that is constricted in the Trefftz plane 
      eop = bop[isurf]/b[isurf]

      # Inverting the cosine spacing and 
      # then normalizing by π/2 so θ spans from 1.0 to 0.0
      # You want these specific θs since you know what Γ is at these points from 
      # the piece wise lift distributions
      t0 = 1.0 # which is basically acos(0.0)/(π/2)
      to = acos(eo)/(0.5*π)
      ts = acos(es)/(0.5*π)
      t1 = 0.0 # which is basically acos(1.0)

      # This is to transform the to and ts angles to the right values such
      # that the bunch transform gives us back exactly the right ηo and ηs.
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
      @inbounds for  k = k0+1 : ko #start at k0+1 cause you already set stuff for k0.
        i = i+1

        fk = (k-k0)/(ko-k0)
        t[i] = t0*(1.0-fk) + to*fk #field points at i
        tc = 0.5*(t[i-1]+t[i]) #collocation point at i+1/2 points

        e  = cos(0.5*π*(t[i] + bunch*t[i]*(1.0-t[i])))
        ec = cos(0.5*π*(tc   + bunch*tc  *(1.0-tc  )))
       
        # E.g., eo can be = e0 if we have something like a T-tail
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
      #Inner panel
      @inbounds for  k = ko+1: ks
        i = i+1

        fk = (k-ko)/(ks-ko)
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
      # Outer panel
      @inbounds for  k = ks+1 : k1
        i = i+1

        fk = (k-ks)/(k1-ks)
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
        gw[i] = 0. #Circulation in the wake
        @inbounds for  i = ifrst[isurf]+1: ilast[isurf]-1
          gw[i] = gc[i-1] - gc[i]
        end
        i = ilast[isurf]
        gw[i] = gc[i-1]
      end


      if(specifies_CL) 
#----- scale circulations to get specified lift for each surface
       @inbounds for  isurf = 1: nsurf
         cltest = 0.
         @inbounds for  i = ifrst[isurf]: ilast[isurf]-1
           dy = yp[i+1] - yp[i]
           cltest = cltest + gc[i]*dy
         end
         
         cltest = cltest * 2.0 * bref/(0.5*Sref)

         ## Calculate the scaling factor for the circulations
         gfac = CLsurfsp[isurf]/(cltest)
        # println("$isurf, $gfac, $cltest, $(CLsurfsp[isurf])")

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
      CDell = CL^2 / (π*AR) #CD for elliptical loading here. 
      # For general CD = CL²/(π*AR*spaneff) ∴spaneff calculated as:
      spanef = CDell / CD

      return CLsurf, CL, CD , spanef
      end
