"""
      induced_drag!(para, wing, htail, trefftz_config)

Computes the induced drag via the Trefftz plane. Calls [`_trefftz_analysis`](@ref). Formerly, `cditrp!()`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `para::AbstractArray{Float64}`: Array of `aircraft` model aerodynamic parameters.
      - `wing::TASOPT.Wing`: Wing Structure.
      - `htail::TASOPT.Tail`: Htail Structure.
      - `trefftz_config::TrefftzPlaneConfig`: Trefftz plane analysis configuration (panel discretization and analysis parameters).

      **Outputs:**
      - No explicit outputs. Computed induced drag value and span efficiency are saved to `para` of `aircraft` model.

!!! compat "Future Changes"
      In an upcoming revision, an `aircraft` `struct` and auxiliary indices will be passed in lieu of pre-sliced `par` arrays.

"""
function induced_drag!(para, wing, htail, trefftz_config::TrefftzPlaneConfig)

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

      gammas   = zeros(Float64, 2)
      gammat   = zeros(Float64, 2)
      po       = zeros(Float64, 2)
      CLsurfsp = zeros(Float64, 2)

#---- wing wake parameters
      fLo =  wing.fuse_lift_carryover

      gammas[1] = wing.inboard.λ*para[iarcls]
      gammat[1] = wing.outboard.λ*para[iarclt]
      po[1]     = 1.0
      CLsurfsp[1] = CL - CLhtail

#---- velocity-change fractions at wing spanwise locations due to fuselage flow
      fduo = para[iafduo]
      fdus = para[iafdus]
      fdut = para[iafdut]

#---- horizontal tail wake parameters
      gammas[2] = 1.0
      gammat[2] = htail.outboard.λ
      po[2]     = 1.0
      CLsurfsp[2] = CLhtail


#---- number of surfaces  (wing, horizontal tail)
      nsurf = 2

      CLsurf, CLtp, CDtp, sefftp = _trefftz_analysis(nsurf, trefftz_config,
	wing, htail,
	Sref, bref,
	po,gammat,gammas, fLo,
      specifies_CL,CLsurfsp, t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)
      
      # println("$CLsurf, $CLtp, $CDtp, $sefftp")

      para[iaCDi] = CDtp
      para[iaspaneff] = sefftp

      return
end # induced_drag!


"""
    bunch_transform(t, bunch)

Apply bunching transformation to input `t ∈ \\[0,1\\]`.

Transforms `t` to cluster points near the center (t=0.5) when `bunch > 0`.
The transformation is: `t_bunched = t + bunch * t * (1 - t)`

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: Original normalized angle ∈ \\[0,1\\].
    - `bunch::Float64`: Clustering factor ∈ \\[0,1\\]. Higher values cluster more points near center.

    **Outputs:**
    - `t_bunched::Float64`: Transformed parameter value.

See also [`inv_bunch_transform`](@ref).
"""
@inline function bunch_transform(t, bunch)
    return t + bunch * t * (1.0 - t)
end

"""
    inv_bunch_transform(t_bunched, bunch)

Inverse of the bunching transformation. Recovers original `t` from bunched value.

Solves `t_bunched = t + bunch * t * (1 - t)` for `t` using the quadratic formula.
See also [`bunch_transform`](@ref).
"""
@inline function inv_bunch_transform(t_bunched, bunch)
    return (1.0 + bunch - sqrt((1.0 + bunch)^2 - 4.0 * bunch * t_bunched)) * 0.5 / bunch
end

@enum WAKE_CONTRACTION_TYPE begin
    FUSEWAKE
    WINGWAKE
end

function get_wake_contraction(wake_type::WAKE_CONTRACTION_TYPE, y, yo, yop)
      if wake_type === FUSEWAKE
            # Power law to contract the stream tubes in the fuselage region.
            @assert y ≤ yo "y must be within fuselage region (y ≤ yo)"
            yexp = (yo / yop)^2
            return yop * (y / yo)^yexp
      else
            # Mass conservation based contraction outside fuselage region.
            @assert y ≥ yo "y must be outside fuselage region (y ≥ yo)"
            return sqrt(y^2 - yo^2 + yop^2)
      end
end

"""
    generate_panel_points!(
        i, t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
        k_start, k_end,
        t_start, t_end,
        e_start, e_end,
        y_start, y_end, z_start, z_end,
        g_start, g_end,
        yo, yop, bunch, ktip
    ) -> i_end

Generate points for a single panel section (image, inboard, or outboard).
This eliminates code duplication between the three panel loops.

Type-stable and efficient for repeated use in different panel sections.
"""
@inline function generate_panel_points!(
    i::Int,
    t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
    k_start::Int, k_end::Int,
    t_start::Float64, t_end::Float64,
    e_start::Float64, e_end::Float64,
    y_start::Float64, y_end::Float64,
    z_start::Float64, z_end::Float64,
    g_start::Float64, g_end::Float64,
    yo::Float64, yop::Float64,
    bunch::Float64, ktip::Float64, 
    wake_contract::WAKE_CONTRACTION_TYPE)

    @inbounds for k = k_start+1:k_end
        i = i + 1

        # Interpolation fraction
        fk = (k - k_start) / (k_end - k_start)

        # Theta values for field point and control point
        t[i] = t_start * (1.0 - fk) + t_end * fk
        tc = 0.5 * (t[i-1] + t[i])

        # Eta values (spanwise position)
        e  = cos(0.5π * bunch_transform(t[i], bunch))
        ec = cos(0.5π * bunch_transform(tc, bunch))

        # Interpolation fractions in eta space
        if e_end - e_start == 0.0
            fi = 1.0
            fc = 0.5
        else
            fi = (e - e_start) / (e_end - e_start)
            fc = (ec - e_start) / (e_end - e_start)
        end

        # Physical coordinates - field points
        y[i] = y_start * (1.0 - fi) + y_end * fi
        z[i] = z_start * (1.0 - fi) + z_end * fi

        # Physical coordinates - control points
        yc[i-1] = y_start * (1.0 - fc) + y_end * fc
        zc[i-1] = z_start * (1.0 - fc) + z_end * fc

        # Circulation with tip rolloff
        gc[i-1] = (g_start * (1.0 - fc) + g_end * fc) * sqrt(1.0 - ec^ktip)

        yp[i] = get_wake_contraction(wake_contract, y[i], yo, yop)
        ycp[i-1] = get_wake_contraction(wake_contract, yc[i-1], yo, yop)

        zp[i] = z[i]
        zcp[i-1] = zc[i-1]
    end

    return i
end

"""
    generate_trefftz_points_single_surface!(
        i_start::Int,
        t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
        surface,
        po, gammat, gammas,
        panels::SurfaceDiscretization,
        trefftz_config::TrefftzPlaneConfig
    ) -> i_end::Int

Generate field points and control points for a single lifting surface in the Trefftz plane.
Returns the final index after point generation.

This is an incremental refactoring step - extracted from _trefftz_analysis for clarity.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `i_start::Int`: Starting index in the arrays.
    - `t, y, yp, z, zp, yc, ycp, zc, zcp, gc`: Work arrays for point storage.
    - `surface`: Wing or tail structure (TASOPT.Wing or TASOPT.Tail) containing geometry.
    - `po::Float64`: Root circulation scaling factor.
    - `gammat::Float64`: Outer section lift distribution taper ratio.
    - `gammas::Float64`: Inner section lift distribution taper ratio.
    - `panels::SurfaceDiscretization`: Panel discretization for this surface.
    - `trefftz_config::TrefftzPlaneConfig`: Configuration with bunching and tip parameters.

    **Outputs:**
    - `i_end::Int`: Final index after point generation.
"""
@inline function generate_trefftz_points_single_surface!(
    i_start::Int,
    t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
    surface,
    po::Float64, gammat::Float64, gammas::Float64,
    panels::SurfaceDiscretization,
    trefftz_config::TrefftzPlaneConfig,
    root_contraction::Float64
)
    # Extract geometry from surface struct
    b = surface.layout.span
    bs = surface.layout.break_span  # uses the getproperty accessor
    bo = surface.layout.root_span
    bop = bo * root_contraction
    zcent = surface.layout.z

    # Extract configuration
    bunch = trefftz_config.bunch
    ktip = trefftz_config.k_tip

    # Extract panel counts
    n_img = panels.n_image_panels
    n_inn = panels.n_inner_panels
    n_out = panels.n_outer_panels

    # Key span positions (η = spanwise location normalized by b)
    e0 = 0.0
    eo = bo/b
    es = bs/b
    e1 = 1.0
    eop = bop/b

    # Transform to θ space
    t0 = 1.0
    to = acos(eo)/(0.5*π)
    ts = acos(es)/(0.5*π)
    t1 = 0.0

    # Apply inverse bunching
    if bunch > 0.0
        to = inv_bunch_transform(to, bunch)
        ts = inv_bunch_transform(ts, bunch)
    end

    # Physical coordinates at key stations
    y0 = 0.0
    yo = 0.5*bo
    ys = 0.5*bs
    y1 = 0.5*b
    yop = 0.5*bop

    z0 = zcent
    zo = zcent
    zs = zcent
    z1 = zcent

    # Circulation at key stations
    g0 = po
    go = po
    gs = po*gammas
    g1 = po*gammat

    # Panel indices
    k0 = 1
    ko = 1 + n_img
    ks = 1 + n_img + n_inn
    k1 = 1 + n_img + n_inn + n_out

    # Start with center point
    i = i_start
    t[i] = t0
    y[i] = y0
    z[i] = z0
    yp[i] = y0
    zp[i] = z0

    # Image panels (fuselage region)
    i = generate_panel_points!(
        i, t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
        k0, ko, t0, to, e0, eo,
        y0, yo, z0, zo, g0, go,
        yo, yop, bunch, ktip, FUSEWAKE)

    # Inboard panels
    i = generate_panel_points!(
        i, t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
        ko, ks, to, ts, eo, es,
        yo, ys, zo, zs, go, gs,
        yo, yop, bunch, ktip, WINGWAKE)

    # Outboard panels
    i = generate_panel_points!(
        i, t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
        ks, k1, ts, t1, es, e1,
        ys, y1, zs, z1, gs, g1,
        yo, yop, bunch, ktip, WINGWAKE)

    return i
end

"""
    _trefftz_analysis(nsurf, trefftz_config, wing, htail, Sref, bref,
          po, gammat, gammas, fLo, specifies_CL, CLsurfsp,
          t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)

Trefftz plane routine for the induced drag computation of `nsurf` number of surfaces. Formerly, `trefftz1()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `nsurf::Integer`: Number of surfaces (typically wing and horizontal tail).
    - `trefftz_config::TrefftzPlaneConfig`: Configuration containing panel discretization and analysis parameters.
    - `wing`: Wing structure (TASOPT.Wing).
    - `htail`: Horizontal tail structure (TASOPT.Tail).
    - `Sref::Float64`: Reference wing area.
    - `bref::Float64`: Reference wing span.
    - `po::Vector{Float64}`: Root circulation scaling factor for each surface.
    - `gammat::Vector{Float64}`: Wing lift distribution "taper" ratios for outer sections.
    - `gammas::Vector{Float64}`: Wing lift distribution "taper" ratios for inner sections.
    - `fLo::Float64`: Wing root load adjustment factor (currently not implemented).
    - `specifies_CL::Bool`: Flag for specified lift calculation (scales vorticities to achieve `CLsurfsp`).
    - `CLsurfsp::Vector{Float64}`: Prescribed surface lift coefficient for each surface.
    - `t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc`: Work arrays for computation.

    **Outputs:**
    - `CLsurf::Vector{Float64}`: Lift coefficients for each surface.
    - `CL::Float64`: Sum of lift coefficients of all surfaces.
    - `CD::Float64`: Sum of induced drag coefficients of all surfaces.
    - `spanef::Float64`: Span efficiency of combined surfaces (``= (CL^2 / (π*AR))/CD``).

See [theory above](@ref trefftz) or Sections 2.14.7 and 3.8.1 of the [TASOPT Technical Desc](@ref dreladocs).
"""
function _trefftz_analysis(nsurf, trefftz_config::TrefftzPlaneConfig,
	wing, htail,
	Sref, bref,
	po, gammat, gammas, fLo,
	specifies_CL,CLsurfsp,t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)

      ifrst = zeros(Int, nsurf)
      ilast = zeros(Int, nsurf)

      CLsurf= zeros(Float64, nsurf)

      # Calculate total number of points
      # outboard point + inboard + points within fuselage + dummy between surfaces
      # Wing contribution
      isum = trefftz_config.wing_panels.n_outer_panels +
             trefftz_config.wing_panels.n_inner_panels +
             trefftz_config.wing_panels.n_image_panels + 1
      # Tail contribution (handle T-tail case where root_span == 0)
      tail_image_panels = (htail.layout.root_span == 0.0) ? 0 : trefftz_config.tail_panels.n_image_panels
      isum += trefftz_config.tail_panels.n_outer_panels +
              trefftz_config.tail_panels.n_inner_panels +
              tail_image_panels + 1

      if(isum > idim)
	      println("TREFFTZ: Passed array overflow. Increase idim to ",isum)
        exit()
      end

      if(isum > jdim)
	      println("TREFFTZ: Local array overflow. Increase jdim to ", isum)
        exit()
      end

      i::Int64 = 0

      @inbounds for isurf = 1:nsurf
          # Select the right surface and panel config
          surface = isurf == 1 ? wing : htail
          panels = isurf == 1 ? trefftz_config.wing_panels : trefftz_config.tail_panels
          root_contraction = isurf == 1 ? trefftz_config.wing_root_contraction : trefftz_config.tail_root_contraction

          # Handle T-tail special case (root_span = 0)
          if isurf == 2 && surface.layout.root_span == 0.0
              panels = SurfaceDiscretization(panels.n_outer_panels, panels.n_inner_panels, 0)
          end

          i = i + 1
          ifrst[isurf] = i

          # Generate points for this surface
          i = generate_trefftz_points_single_surface!(
              i, t, y, yp, z, zp, yc, ycp, zc, zcp, gc,
              surface,
              po[isurf], gammat[isurf], gammas[isurf],
              panels,
              trefftz_config,
              root_contraction
          )

          ilast[isurf] = i

          #---- dummy control point between surfaces
          yc[i] = 0.0
          zc[i] = 0.0
          gc[i] = 0.0
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
