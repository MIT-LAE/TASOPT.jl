"""
   aircraft_drag!(ac, imission, ip, computes_wing_direct)

Calculates aircraft `CD` components for operating point, ipoint.
If `computes_wing_direct` is `true`, computes wing `cdf`,`cdp` from airfoil database # `iairf`,
otherwise uses default values in para array. 
Called by [`_mission_iteration!()`](@ref TASOPT._mission_iteration!), [`size_aircraft!`](@ref TASOPT.size_aircraft!), and [`takeoff!`](@ref TASOPT.takeoff!).
Formerly, `cdsum!()`.

The total drag is computed by

```math
C_{D} = C_{D, i} + C_{D,fuse} + C_{D,wing} + C_{D,over} + C_{D,htail} + C_{D,vtail} + C_{D,strut} + C_{D,nace} + \\Delta C_{D,BLI,f} + \\Delta C_{D,BLI,w}
```

where:
- ``C_{D,i}`` (`CDi`) is the total induced drag including the wing and tail,
- ``C_{D,fuse}`` (`CDfuse`) is the fuselage profile drage computed by solving a boundary layer integral equation,
- ``C_{D,wing}`` (`CDwing`) is the wing profile drag (viscous + pressure) computed using airfoil data obtained from CFD,
- ``C_{D,over}`` (`CDover`) is the fuselage added CD due to lift carryover,
- ``C_{D,htail}`` (`CDhtail`) is the horizontal tail profile drag computed from assumed drag coefficients,
- ``C_{D,vtail}`` (`CDvtail`) is the vertical tail profile drag computed from assumed drag coefficients,
- ``C_{D,strut}`` (`CDstrut`) is the profile drag of the strut and, for a truss-braced wing, of the jury struts,
- ``C_{D,nace}`` (`CDnace`) is the nacelle profile drag,
- ``\\Delta C_{D,BLI,f}`` (`dCDBLIf`) is related to the boundary layer ingestion on the fuselage,
- and ``\\Delta C_{D,BLI,w}`` (`dCDBLIw`) is related to the boundary layer ingestion on the wing.


!!! details "🔃 Inputs and Outputs"
**Inputs:**
      - `ac::aircraft`: aircraft data storage object
      - `imission::Int64`: mission index
      - `ip::Int64`: operating point index
      - `computes_wing_direct::Bool`: Flag if drag should be computed with `wing_profiledrag_direct` (true) or if para values should be used (false).

      **Outputs:**
      - No explicit outputs. Computed drag values are saved to `para` of `aircraft` model.

See Section 2.14 of the [TASOPT Technical Desc](@ref dreladocs).
See also [`fuselage_drag!`](@ref), [`wing_profiledrag_direct`](@ref), [`wing_profiledrag_scaled`](@ref), [`cfturb`](@ref), and [`induced_drag!`](@ref).
"""
function aircraft_drag!(ac, imission::Int, ip::Int, computes_wing_direct::Bool; Ldebug=false)
      #Unpack data storage
      parg = ac.parg
      para = view(ac.para, :, ip, imission)
      pare = view(ac.pare, :, ip, imission)
      wing = ac.wing
      htail = ac.htail
      vtail = ac.vtail

      fSnace   = parg[igfSnace ]
      
      gammat   = wing.outboard.λ*para[iarclt]
      gammas   = wing.inboard.λ*para[iarcls]

      fCDwcen = 0.0
      fCDhcen = parg[igfCDhcen]
      fCDvcen = 0.0

      rVnace = parg[igrVnace]

      CL   = para[iaCL]
      Mach = para[iaMach]

      fduo = para[iafduo]
      fdus = para[iafdus]
      fdut = para[iafdut]

      DAfsurf = para[iaDAfsurf]
      DAfwake = para[iaDAfwake]
      PAfinf = para[iaPAfinf]

      fexcdw = para[iafexcdw]

#---- tail lift
      CLhtail = para[iaCLh]*htail.layout.S/wing.layout.S

#---- parameters for Re-scaling of cd's
      Reunit = para[iaReunit]
      Rerefw = para[iaRerefw]
      Rereft = para[iaRereft]
      aRexp  = para[iaaRexp]

#---- root shock unsweep constant
      rkSunsw = 0.5
      rkSunsh = 0.
      rkSunsv = 0.
#
#---- Re referenced to root chord
      Reco = Reunit*wing.layout.root_chord

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Strut lift is not carried by the wing
      CL_strut = strut_lift_share(wing, ac.options.trefftz_config) * (CL - CLhtail)

      if (computes_wing_direct) 
#----- integrated across span for CDwing
      wing_drag_components = wing_profiledrag_direct(wing,gammat,gammas,
                                    Mach,CL,CLhtail + CL_strut,Reco,
                                    aRexp, rkSunsw,fexcdw,
                                    fduo,fdus,fdut; imission, ip)
      
      # Will get rid of these when I start culling the para arrays.                                     
      clpo = wing_drag_components.clpo
      clps = wing_drag_components.clps
      clpt = wing_drag_components.clpt
      cdfw = wing_drag_components.CDfwing
      cdpw = wing_drag_components.CDpwing
      CDwing = wing_drag_components.CDwing
      CDover = wing_drag_components.CDover

      # Flag a wing outside the airfoil database at cruise
      if ipcruise1 <= ip <= ipcruisen && wing_drag_components.penalty_fraction > 0.0
            record_infeasibility!(ac, :airfoil_database_exceeded)
      end
		
       #if(Ldebug) write(*,*) '...exited wing_profiledrag_direct'

#----- store CD values
      para[iaCDwing] = CDwing
      para[iaCDover] = CDover
      para[iacdfw] = cdfw
      para[iacdpw] = cdpw

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
#----- just use stored values
      cdfw = para[iacdfw] * fexcdw
      cdpw = para[iacdpw] * fexcdw

	CDwing,CDover = wing_profiledrag_scaled(wing.layout.S,
	wing.layout.span,wing.layout.break_span,wing.layout.root_span,wing.outboard.λ,wing.inboard.λ,wing.layout.sweep,wing.layout.root_chord, 
	cdfw,cdpw,Reco,Rerefw,aRexp,rkSunsw,
       fCDwcen)
 
      para[iaCDwing] = CDwing
      para[iaCDover] = CDover

      clpo, clps, clpt = wing_section_cls(wing,gammat,gammas,
                              CL,CLhtail,
	                        fduo,fdus,fdut)

      end

#---- store local cl values
      para[iaclpo] = clpo
      para[iaclps] = clps
      para[iaclpt] = clpt

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- tail profile Cd's
      cdft = para[iacdft] * para[iafexcdt]
      cdpt = para[iacdpt] * para[iafexcdt]

#---- horizontal tail profile CD
      Recoh = Reunit*htail.layout.root_chord
	CDhtail,CDhover = wing_profiledrag_scaled(wing.layout.S,
	htail.layout.span,htail.layout.root_span,htail.layout.root_span,htail.outboard.λ,1.0,htail.layout.sweep,htail.layout.root_chord, 
	cdft,cdpt,Recoh,Rereft,aRexp,rkSunsh,
	fCDhcen)

#---- vertical tail profile CD
      Recov = Reunit*vtail.layout.root_chord
      CDvtail1,CDvover1 = wing_profiledrag_scaled(wing.layout.S,
	vtail.layout.span,vtail.layout.root_span,vtail.layout.root_span,vtail.outboard.λ,1.0,vtail.layout.sweep,vtail.layout.root_chord, 
	cdft,cdpt,Recov,Rereft,aRexp,rkSunsv,
      fCDvcen)
	
      CDvtail = CDvtail1*vtail.ntails

      para[iaCDhtail] = CDhtail
      para[iaCDvtail] = CDvtail

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- fuselage profile CD, using fuselage CDA from BLI calculation
      #note/reminder that the CDA is assumed constant across the entire flight
      CDfuse = PAfinf/wing.layout.S
      para[iaCDfuse] = CDfuse

#      if(Ldebug) write(*,*) 'nacelle CD...'
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- set new nacelle skin friction coefficient Cfnace
#-     (use fuselage excrescence factor)
      lnace = parg[iglnace]
      Reunit = para[iaReunit]
      if (Reunit == 0.0) 
       Cfnace = 0.
      else
       Renace = Reunit*lnace
       Cfnace = cfturb(Renace) * para[iafexcdf]
      end
      para[iaCfnace] = Cfnace

#---- set Mnac on outside of nacelle, using vortex sheet model of nacelle,
#-     (Mnac+M2)/2 = Mach*rVnace
      rVnLE = max( 2.0*rVnace - pare[ieM2] / max(Mach,0.001) , 0.0 )
      rVnsurf3 = 0.25*(rVnLE+rVnace)*(rVnLE^2+rVnace^2)
      CDnace = fSnace * Cfnace * rVnsurf3
      para[iaCDnace] = CDnace

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- strut and jury strut profile CD
      strut_exceeded = false
      if ac.options.has_strut
            strut = wing.strut
            rVstrut = strut.local_velocity_ratio
            # Main strut section cl from the lift of both struts on their effective area S cosΓ (area seen from above)
            if strut.S > 0.0 && rVstrut > 0.0
                  cl_strut = CL_strut * wing.layout.S / (strut.S * cos(strut.dihedral) * rVstrut^2)
            else
                  cl_strut = 0.0
            end
            mean_strut_chord = 0.5 * (strut.root_chord + strut.chord)
            main_strut = strut_member_drag(strut.S/wing.layout.S, mean_strut_chord, cl_strut,
                                           strut.thickness_to_chord, Mach, strut.cos_lambda,
                                           Reunit, rVstrut, aRexp, fexcdw, wing.airsection)
            CDstrut = main_strut.CD
            strut_exceeded = !isempty(main_strut.exceeded)
            if ac.options.opt_wing_config == WingConfig.TrussBraced
                  # Jury struts (one per side) are vertical, unswept and carry no lift
                  jury = wing.jury
                  jury_strut = strut_member_drag(2*jury.chord*jury.length/wing.layout.S, jury.chord, 0.0,
                                                 jury.thickness_to_chord, Mach, 1.0,
                                                 Reunit, rVstrut, aRexp, fexcdw, wing.airsection)
                  CDstrut += jury_strut.CD
                  # Jury cl is not flagged, it carries no lift
                  strut_exceeded |= :thickness in jury_strut.exceeded || :Mach in jury_strut.exceeded
            end
      else
            CDstrut = 0.0
      end
      if strut_exceeded && ipcruise1 <= ip <= ipcruisen
            record_infeasibility!(ac, :strut_airfoil_database_exceeded)
      end
      para[iaCDstrut] = CDstrut


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- induced CD
      induced_drag!(para, ac, ac.options.trefftz_config)
      CDi = para[iaCDi]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- apparent fuselage deltaCD from ingestion
      fBLIf = parg[igfBLIf]
      dCDBLIf = -fBLIf*DAfwake/wing.layout.S
      para[iadCDBLIf] = dCDBLIf

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- apparent wing deltaCD from ingestion
      fBLIw = parg[igfBLIw]
      CDWwake = CDwing * 0.15  # assume 15% of the wing dissipation is in wake
      dCDBLIw = -fBLIw*CDWwake
      para[iadCDBLIw] = dCDBLIw

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- total CD
      CD = CDi + CDfuse + CDwing + CDover + CDhtail + CDvtail + CDstrut + CDnace + dCDBLIf + dCDBLIw
      para[iaCD] = CD      

      return
end # aircraft_drag!

"""
      cfturb(Re)

Returns the total ``C_f`` for turbulent flat plate (one side) as a function of
``\\mathrm{Re}_l``

```math
C_{f, turb} = \\frac{0.523}{(\\log(0.06\\mathrm{Re}))^2} 
```

"""
function cfturb(Re)

#   cfturb = 0.427/(log10(re) - 0.407)**2.64   ! original Hoerner
#   cfturb = 0.310/(log10(re) - 0.407)**2.47   ! modified (weaker dependence on Re)

      cfturb = 0.523/(log(0.06*Re))^2           # White

      return cfturb
end # cfturb

"""
    strut_section_cd(cl, τ, Mach, cosΛ, Re_chord, aRexp, fexcd, air)

Section drag coefficient of a strut or jury strut member from the airfoil database `air` (see [`airfun`](@ref)),
using simple sweep theory: `cl/cosΛ^2` and `Mach*cosΛ` normal to the member. Below the lowest database `cl`, the
lookup is made at that `cl`. Drag is scaled by `(Re_chord/air.Re)^aRexp` and `fexcd`; the pressure part also by `cosΛ`.

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `cl::Float64`: Streamwise section lift coefficient of the member.
      - `τ::Float64`: Thickness to chord ratio.
      - `Mach::Float64`: Local Mach number at the member.
      - `cosΛ::Float64`: Cosine of the member sweep.
      - `Re_chord::Float64`: Reynolds number based on the member chord.
      - `aRexp::Float64`: Re-scaling exponent.
      - `fexcd::Float64`: Excrescence cd scaling factor.
      - `air::TASOPT.aerodynamics.airfoil`: Airfoil database.

      **Outputs (`NamedTuple`):**
      - `cd::Float64`, `cdf::Float64`, `cdp::Float64`: Total, friction and pressure section drag coefficients.
      - `exceeded::Vector{Symbol}`: Database bounds exceeded (`:cl`, `:thickness`, `:Mach`), see [`airfoil_database_exceedances`](@ref).
"""
function strut_section_cd(cl, τ, Mach, cosΛ, Re_chord, aRexp, fexcd, air::airfoil)
      cl_perp = cl / cosΛ^2
      M_perp = Mach * cosΛ

      cl_lookup = max(cl_perp, air.cl[1])
      cdf1, cdp1, _, _ = airfun(cl_lookup, τ, M_perp, air)

      Refac = (Re_chord / air.Re)^aRexp
      cdf = cdf1 * Refac * fexcd
      cdp = cdp1 * Refac * fexcd * cosΛ

      exceeded = airfoil_database_exceedances(cl_perp, τ, M_perp, air)
      return (cd = cdf + cdp, cdf = cdf, cdp = cdp, exceeded = exceeded)
end

"""
    strut_member_drag(area_ratio, chord, cl, τ, Mach, cosΛ, Reunit, rV, aRexp, fexcd, air)

Profile drag coefficient, referenced to the wing area, of a strut or jury strut member in a local flow of `rV`
times the freestream. Uses [`strut_section_cd`](@ref) scaled by `area_ratio` (member area / wing area) and `rV^3`.
Zero if `area_ratio` or `chord` is zero.

Returns a `NamedTuple` with `CD` and `exceeded` (database bounds exceeded).
"""
function strut_member_drag(area_ratio, chord, cl, τ, Mach, cosΛ, Reunit, rV, aRexp, fexcd, air::airfoil)
      (area_ratio == 0.0 || chord == 0.0) && return (CD = 0.0, exceeded = Symbol[])
      section = strut_section_cd(cl, τ, Mach*rV, cosΛ, Reunit*rV*chord, aRexp, fexcd, air)
      return (CD = area_ratio * section.cd * rV^3, exceeded = section.exceeded)
end
