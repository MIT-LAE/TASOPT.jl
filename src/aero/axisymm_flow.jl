"""
    _axisymm_flow(xnose, xend, xblend1, xblend2, Amax, 
	      anose, btail, iclose,
	      Mach, nc, nldim,
            xl, zl, sl, dyl, ql)

Calculates compressible potential flow about a quasi-axisymmetric body, 
using a simple piecewise-constant source line. Formerly, `axisol!()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `xnose::Float64`: X (axial) location of nose point. 
    - `xend::Float64`: X location of tail point.
    - `xblend1::Float64`: X location of nose-section blend point.
    - `xblend2::Float64`: X location of tail-section blend point.
    - `Amax::Float64`: Maximum cross-sectional area.
    - `anose::Float64`: Nose-section shape exponent.
    - `btail::Float64`: Tail-section shape exponent.
    - `iclose::Integer`: If 0, tail tapers to a point, otherwise to an edge.
    - `Mach::Float64`: Freestream Mach number for Prandtl-Glauert.
    - `nc::Integer`: Number of control points to be used.
    - `nldim::Integer`: Max dimension of passed arrays.

    **Outputs:**
    - `nl::Integer`: Number of output surface and wake points.
    - `ilte::Integer`: Index of TE point.
    - `xl::Array{Float64}`: X (axial) locations of surface segment endpoints.
    - `zl::Array{Float64}`: Z (vertical) locations of surface segment endpoints.
    - `sl::Array{Float64}`: Arc lengths along surface and wake.
    - `dyl::Array{Float64}`: Half-width of edge-type tail section.
    - `ql::Array{Float64}`: Velocities V/V_inf along surface and wake.

See [theory above](@ref axi) or Section 3 of [Simplified Viscous/Inviscid Analysis for Nearly-Axisymmetric Bodies](../assets/drela_TASOPT_2p16/axibl.pdf). 
See also [`fuselage_drag!`](@ref).
"""
function _axisymm_flow(xnose, xend, xblend1, xblend2, Amax, 
	anose, btail, iclose,
	Mach, nc, nldim,
      xl, zl, sl, dyl, ql)

      
      idim  = 40
      xc = zeros(idim) 
      zc = zeros(idim)
      ac = zeros(idim)
      dyc = zeros(idim)
      src = zeros(idim)

      zxc = zeros(idim)
      nxc = zeros(idim)
      nyc = zeros(idim)
      nzc = zeros(idim)
      
      aa = @MMatrix zeros(idim,idim)
      rr = @MVector zeros(idim)

#cc   ispace =  0     # uniform x spacing
      ispace =  1     # cosine  x spacing

      if(nc > idim) 
	println("_axisymm_flow: Local-array overflow.  Increase idim to", nc)
        quit()
      end

      ilte  = nc + 1
      nl = nc + Int(floor(nc/2)) + 2

      if(nl > nldim) 
	println("_axisymm_flow: Passed-in array overflow.  Increase nldim to", nl)
        quit()
      end

      # xl  = @MVector zeros(nldim)
      # zl  = @MVector zeros(nldim)
      # sl  = @MVector zeros(nldim)
      # dyl = @MVector zeros(nldim)
      # ql  = @MVector zeros(nldim)


      beta = sqrt(1.0 - Mach^2)

#---- max radius of equivalent round body
      Rcyl = sqrt(Amax/pi)

#---- set body geometry points
    @inbounds for  i = 1: ilte
        tl = float(i-1) / float(ilte-1)
        if(ispace==0) 
         frac = tl
        else
         frac = 0.5*(1.0 - cos(pi*tl))
        end
        xl[i] = xnose*(1.0-frac) + xend*frac

        if(i==1 || i==ilte) 
         zl[i] = 0.
        elseif(xl[i] < xblend1) 
         f = 1.0 - (xl[i]-xnose)/(xblend1-xnose)
         zl[i] = Rcyl*(1.0 - f^anose)^(1.0/anose)
        elseif(xl[i] < xblend2) 
         zl[i] = Rcyl
        else
         f = (xl[i]-xblend2)/(xend-xblend2)
         zl[i] = Rcyl*(1.0 - f^btail)
        end

        if(iclose==0) 
          dyl[i] = 0.
        else
          if(xl[i] < xblend2) 
           dyl[i] = 0.
          else
           dyl[i] = Rcyl - zl[i]
          end
        end

      end
      zl[ilte] = 0.25*zl[ilte-1]
 
      if(iclose==0) 
       dyl[ilte] = 0.
      else
       dyl[ilte] = Rcyl - zl[ilte]
      end

#---- wake geometry points
      @inbounds for  i = ilte+1: nl
        xl[i] = 2.0*xend - xl[2*ilte-i]
        zl[i] = 0.125*zl[ilte-1]

        if(iclose==0) 
         dyl[i] = 0.
        else
         dyl[i] = Rcyl - zl[i]
        end
      end

#---- calculate arc lengths
      sl[1] = 0.
      @inbounds for  i = 1: nl-1
        ds = sqrt((xl[i+1]-xl[i])^2 + (zl[i+1]-zl[i])^2)
        sl[i+1] = sl[i] + ds
      end


#----- calculate wetted area, out of interest
#      Awet = 0.
#      for i = 1: ilte-1
#        ds = sqrt((xl[i+1]-xl[i])^2 + (zl[i+1]-zl[i])^2)
#        ra = 0.5*(zl[i+1]+zl[i])
#        Awet = Awet + 2.0*pi*ra*ds
#      end


#---- set control points xc(),zc(.) on actual surface
#-     dyc() is lateral offset over edge-type tail section
      @inbounds for  i = 1: nc
        tlo = float(i-1) / float(ilte-1)
        tlp = float(i  ) / float(ilte-1)

        tc = 0.5*(tlo+tlp)
        if(ispace==0) 
         frac = tc
        else
         frac = 0.5*(1.0 - cos(pi*tc))
        end
        xc[i] = xnose*(1.0-frac) + xend*frac

#c      xc[i] = 0.5*(xl[i] + xl[i+1])
        if    (xc[i] < xblend1) 
#------- nose section
         f = 1.0 - (xc[i]-xnose)/(xblend1-xnose)
         f_x = -1.0/(xblend1-xnose)
         zc[i] = Rcyl*(1.0 - f^anose)^(1.0/anose)
         z_f = -zc[i]/(1.0 - f^anose) * f^(anose-1.0)
         z_x = z_f*f_x
         dyc[i] = 0.

        elseif(xc[i] < xblend2) 
#------- constant section
         zc[i] = Rcyl
         z_x = 0.
         dyc[i] = 0.

        else
#------- tail section
         f = (xc[i]-xblend2)/(xend-xblend2)
         f_x = 1.0/(xend-xblend2)
         zc[i] = Rcyl*(1.0 - f^btail)
         z_f = -Rcyl*btail*f^(btail-1.0)
         z_x = z_f*f_x

         if(iclose==0) 
#-------- taper to a point
          dyc[i] = 0.
         else
#-------- taper to an edge
          dyc[i] = Rcyl - zc[i]
         end
        end

#       dyc[i] = 0.05*zc[i]   ####

#------ surface-normal vector
        nxc[i] =  z_x / sqrt(1.0 + z_x^2)
        nzc[i] = -1.0 / sqrt(1.0 + z_x^2)
        nyc[i] = 0.

        zxc[i] = z_x

#------ meridonal arc length and surface area of i..i+1 segment
        ds = sqrt((xl[i+1]-xl[i])^2 + (zl[i+1]-zl[i])^2)
        ac[i] = ds * 2.0*pi*(zl[i] + 4.0*zc[i] + zl[i+1])/6.0
      end


#---- clear AIC matrix
	rr = zeros(nc)
	aa = zeros(nc,nc)      

#---- go over control points
      @inbounds for  i = 1: nc
#
#------ go over source segments, setting up V.n=0 contributions at control point i
        yc = 0.
        @inbounds for  j = 1: nc
          if(dyc[j] == 0.0) 
#--------- singularity element is source line over xl[i]...xl[i+1]
           u0, v0, w0 =  vline(xc[i],yc,zc[i], xl[j],xl[j+1],beta)
          else
#--------- singularity element is source panel over xl[i]...xl[i+1], y1..y2
           y1 = -dyc[j]
           y2 =  dyc[j]
           u0, v0, w0 = vsurf(xc[i],yc,zc[i], xl[j],xl[j+1],y1,y2,beta)
          end

#-------- add on contribution to V.n=0 equation
          aa[i,j] = aa[i,j] + u0*nxc[i] + v0*nyc[i] + w0*nzc[i]
        end

#------ unit-freestream contribtuion to V.n=0 equation
        rr[i] = nxc[i]
      end

#---- solve for all source strenghts
      rr = aa\rr  #call gaussn(idim,nc,aa,rr,1)
      @inbounds for  i = 1: nc
        src[i] = -rr[i]
      end
      
#      for i = 1: n
#        write(41,*) xc[i], src[i]
#      end

#---- calculate velocities at all xl,zl points
      ql[1] = 0.0
      @inbounds for  i = 2: nl
#------ clear velocity accumulators
        ul = 0.
        vl = 0.
        wl = 0.
#
#------ accumulate source veleocities
        yl = 0.
        @inbounds for  j = 1: nc
          if(dyc[j] == 0.0) 
           u0,v0,w0 =  vline(xl[i],yl,zl[i], xl[j],xl[j+1],beta)
          else
           y1 = -dyc[j]
           y2 =  dyc[j]
           u0,v0,w0 =  vsurf(xl[i],yl,zl[i], xl[j],xl[j+1],y1,y2,beta)
          end
          ul = ul + u0*src[j]
          vl = vl + v0*src[j]
          wl = wl + w0*src[j]
        end
#
#------ add on freestream
        ul = ul + 1.0

#------ set speed
        ql[i] = sqrt(ul^2 + vl^2 + wl^2)

#        write(42,'(1x,9g13.5)') 
#     &     xl[i],zl[i],sl[i], ql[i]
      end

      #Create static arrays for performance
      # sxl  = SVector{nldim}(xl)
      # szl  = SVector{nldim}(zl)
      # ssl  = SVector{nldim}(sl)
      # sdyl = SVector{nldim}(dyl)
      # sql  = SVector{nldim}(ql)

      return nl, ilte #xl, zl, sl, dyl, ql
      end # _axisymm_flow


"""
----------------------------------------------------------
     Sets velocity u0,v0,w0 at location x,y,z,
     due to unit-strength source line segment located 
     at (x1..x2, 0, 0)
----------------------------------------------------------
"""   
function vline(x,y,z,x1,x2,b)
      qopi= 1/(4*pi)

      bsq = b^2

      yb = y*b
      zb = z*b

      R1 = sqrt((x1-x)^2 + yb^2 + zb^2)
      R2 = sqrt((x2-x)^2 + yb^2 + zb^2)

      hbsq = yb^2 + zb^2

#      p0_zb = -((x-x2)/R2 - (x-x1)/R1)*qopi/zb
#      w0 = p0_zb*zb_r
#      p0 = log( (R1-x+x1) / (R2-x+x2) ) * qopi / bsq
#      f0 = (R2 - R1 - (x2-x1)) * qopi

      u0 = qopi*(   1.0/R2 -    1.0/R1)         / bsq
      v0 = qopi*((x2-x)/R2 - (x1-x)/R1)*yb/hbsq / b
      w0 = qopi*((x2-x)/R2 - (x1-x)/R1)*zb/hbsq / b

# R(x,y,z,b) = sqrt(x^2 + (y*b)^2 + (z*b)^2)
# ul(x,y,z,b) = ( 1/R(x-x2,y,z,b) - 1/R(x-x1,y,z,b) ) / b^2
# vl(x,y,z,b) = ((x2-x)/R(x-x2,y,z,b) - (x1-x)/R(x-x1,y,z,b))*(y*b)/((y*b)^2 + (z*b)^2) / b
# wl(x,y,z,b) = ((x2-x)/R(x-x2,y,z,b) - (x1-x)/R(x-x1,y,z,b))*(z*b)/((y*b)^2 + (z*b)^2) / b

      return u0, v0, w0
end # vline

"""
Sets velocity u0,v0,w0 at location x,y,z,
due to unit-strength source panel segment located
at (x1..x2, y1..y2, 0)
"""
function vsurf(x,y,z, x1,x2,y1,y2,b)
      qopi= 1/(4*pi)

      bsq = b^2

      y1b = y1*b
      y2b = y2*b
      zb = z*b

      R11 = sqrt((x1-x)^2 + y1b^2 + zb^2)
      R12 = sqrt((x1-x)^2 + y2b^2 + zb^2)
      R21 = sqrt((x2-x)^2 + y1b^2 + zb^2)
      R22 = sqrt((x2-x)^2 + y2b^2 + zb^2)

      g11 = log( abs(y1b + R11) )
      g12 = log( abs(y2b + R12) )
      g21 = log( abs(y1b + R21) )
      g22 = log( abs(y2b + R22) )

      h11 = atanh( (x1-x) / R11 )
      h12 = atanh( (x1-x) / R12 )
      h21 = atanh( (x2-x) / R21 )
      h22 = atanh( (x2-x) / R22 )

      t11 = atan( (x1-x)*y1b , zb*R11 )
      t12 = atan( (x1-x)*y2b , zb*R12 )
      t21 = atan( (x2-x)*y1b , zb*R21 )
      t22 = atan( (x2-x)*y2b , zb*R22 )

      u0 = qopi*(g22 - g21 - g12 + g11)/(y2b-y1b) / bsq
      v0 = qopi*(h22 - h21 - h12 + h11)/(y2b-y1b) / b
      w0 = qopi*(t22 - t21 - t12 + t11)/(y2b-y1b) / b

# gr(x,y,z,xo,yo,b) = log( abs( (yo-y)*b + R(x-xo,y-yo,z,b) ) )
# hy(x,y,z,xo,yo,b) = atanh( (xo-x) / R(x-xo,y-yo,z,b) )
# tz(x,y,z,xo,yo,b) = atan2( (xo-x)*(yo-y)*b , (z*b)*R(x-xo,y-yo,z,b) )

# us(x,y,z,y1,y2,b) = (gr(x,y,z,x2,y2,b) - gr(x,y,z,x2,y1,b) - gr(x,y,z,x1,y2,b) + gr(x,y,z,x1,y1,b)) / (y2-y1) / b^3
# vs(x,y,z,y1,y2,b) = (hy(x,y,z,x2,y2,b) - hy(x,y,z,x2,y1,b) - hy(x,y,z,x1,y2,b) + hy(x,y,z,x1,y1,b)) / (y2-y1) / b^2
# ws(x,y,z,y1,y2,b) = (tz(x,y,z,x2,y2,b) - tz(x,y,z,x2,y1,b) - tz(x,y,z,x1,y2,b) + tz(x,y,z,x1,y1,b)) / (y2-y1) / b^2

      return u0,v0,w0
end # vsurf


