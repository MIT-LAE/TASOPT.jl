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
- ``C_{D,htail}`` (`CDhtail`) is the horizontal tail profile drag computed in a similar manner with `CDwing`,
- ``C_{D,vtail}`` (`CDvtail`) is the vertical tail profile drag computed in a similar manner with `CDwing`,
- ``C_{D,strut}`` (`CDstrut`) is the strut profile drag, 
- ``C_{D,nace}`` (`CDnace`) is the nacelle profile drag,
- ``\\Delta C_{D,BLI,f}`` (`dCDBLIf`) is related to the boundary layer ingestion on the fuselage,
- and ``\\Delta C_{D,BLI,w}`` (`dCDBLIw`) is related to the boundary layer ingestion on the wing.


!!! details "🔃 Inputs and Outputs"
**Inputs:**
      - `ac::aircraft`: aircraft data storage object
      - `imission::Int64`: mission index
      - `computes_wing_direct::Bool`: Flag if drag should be computed with `wing_profiledrag_direct` (true) or if para values should be used (false).

      **Outputs:**
      - No explicit outputs. Computed drag values are saved to `para` of `aircraft` model.

See Section 2.14 of the [TASOPT Technical Desc](@ref dreladocs).
See also [`fuselage_drag!`](@ref), [`wing_profiledrag_direct`](@ref), [`wing_profiledrag_scaled`](@ref), [`cfturb`](@ref), and [`induced_drag!`](@ref).
"""
function aircraft_drag!(ac, imission, ip, computes_wing_direct; Ldebug=false)
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
      Rerefs = para[iaRerefs]
      aRexp  = para[iaaRexp]

#---- root shock unsweep constant
      rkSunsw = 0.5
      rkSunsh = 0.
      rkSunsv = 0.
#
#---- Re referenced to root chord
      Reco = Reunit*wing.layout.root_chord

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (computes_wing_direct) 
#----- integrated across span for CDwing
      clpo,clps,clpt,
	cdfw,cdpw,CDwing,CDover = wing_profiledrag_direct(wing,gammat,gammas,
                                    Mach,CL,CLhtail,Reco,
                                    aRexp, rkSunsw,fexcdw,
                                    fduo,fdus,fdut)
	
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
#---- strut profile CD
      if wing.has_strut
            cdfs = para[iacdfs]
            cdps = para[iacdps]
            rVstrut = wing.strut.local_velocity_ratio
            CDstrut = (wing.strut.S/wing.layout.S)*(cdfs + cdps*wing.strut.cos_lambda^3) * rVstrut^3
            para[iaCDstrut] = CDstrut
      else
            CDstrut = 0.0
            para[iaCDstrut] = 0.0
      end


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- induced CD
      induced_drag!(para, wing, htail)
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


