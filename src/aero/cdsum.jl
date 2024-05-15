"""
    cdsum!(pari,parg,para,pare,icdfun)

Calculates aircraft `CD` components for operating point, ipoint.
If `icdfun=1`, computes wing `cdf`,`cdp` from airfoil database # `iairf`,
otherwise uses default values in para array. Called by `mission!`, `wsize`, `takeoff!`, and `odperf!`.

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
- ``C_{D,strut}`` (`CDstrut`) is the struct profile drag, 
- ``C_{D,nace}`` (`CDnace`) is the nacelle profile drag,
- ``\\Delta C_{D,BLI,f}`` (`dCDBLIf`) is related to the boundary layer ingestion on the fuselage,
- and ``\\Delta C_{D,BLI,w}`` (`dCDBLIw`) is related to the boundary layer ingestion on the wing.


!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `pari::AbstractVector{Int64}`: Vector of `aircraft` model integer/flag parameters.
      - `parg::AbstractArray{Float64}`: Vector of `aircraft` model geometry parameters.
      - `para::AbstractArray{Float64}`: Vector of `aircraft` model aerodynamic parameters.
      - `pare::AbstractArray{Float64}`: Vector of `aircraft` model engine parameters.
      - `icdfun::Integer`: Flag if drag should be computed (=1) or if para values should be used (=0).

      **Outputs:**
      - No explicit outputs. Computed drag values are saved to `para` of `aircraft` model.

See Section 2.14 of the [TASOPT Technical Desc](@ref dreladocs).
See also [`trefftz1`](@ref), [`fusebl!`](@ref), [`surfcd2`](@ref), [`surfcd`](@ref), [`cfturb`](@ref), and `cditrp`.

!!! compat "Future Changes"
      In an upcoming revision, an `aircraft` struct and auxiliary indices will be passed in lieu of pre-sliced `par` arrays.

"""
function cdsum!(pari,parg,para,pare, wing, icdfun)

      Ldebug = false
#      Ldebug = true

      AR       = wing.layout.AR
      sweep    = wing.layout.sweep
      hboxo    = wing.layout.root_chord_thickness
      hboxs    = wing.layout.spanbreak_chord_thickness
      hboxt    = hboxs
      fSnace   = parg[igfSnace ]
      bo       = wing.layout.box_halfspan
      bs       = parg[igbs     ]
      boh      = parg[igboh    ]
      bov      = parg[igbov    ]

      lambdat  = wing.layout.λt
      lambdas  = wing.layout.λs
      gammat   = wing.layout.λt*para[iarclt]
      gammas   = wing.layout.λs*para[iarcls]

      lambdah  = parg[iglambdah]
      lambdav  = parg[iglambdav]
      sweeph   = parg[igsweeph ]
      sweepv   = parg[igsweepv ]
      cosLs    = parg[igcosLs  ]
      Sstrut   = parg[igSstrut ]

      co   = parg[igco ]
      coh  = parg[igcoh]
      cov  = parg[igcov]

      b    = parg[igb  ]
      bh   = parg[igbh ]
      bv   = parg[igbv ]

      S    = parg[igS  ]
      Sh   = parg[igSh ]
      Sv   = parg[igSv ]

      fLo  = parg[igfLo]
      fLt  = parg[igfLt]

      fCDwcen = 0.0
      fCDhcen = parg[igfCDhcen]
      fCDvcen = 0.0

      Astrut   = parg[igAstrut ]
      cstrut   = parg[igcstrut ]
      nvtail   = parg[ignvtail ]

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
      CLhtail = para[iaCLh]*Sh/S

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
      Reco = Reunit*co

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (icdfun==1) 
#----- integrated across span for CDwing

      # if(Ldebug) write(*,*) 'calling SURFCD2...'
      clpo,clps,clpt,
	cdfw,cdpw,CDwing,CDover = surfcd2(S, b,bs,bo,
	lambdat,lambdas,gammat,gammas,
	hboxo,hboxs,hboxt,
	Mach,sweep,co, 
	CL,CLhtail, fLo,fLt,
	Reco,aRexp, rkSunsw,fexcdw,
	fduo,fdus,fdut)
	
       #if(Ldebug) write(*,*) '...exited SURFCD2'

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

	CDwing,CDover = surfcd(S,
	b,bs,bo,lambdat,lambdas,sweep,co, 
	cdfw,cdpw,Reco,Rerefw,aRexp,rkSunsw,
       fCDwcen)
 
      para[iaCDwing] = CDwing
      para[iaCDover] = CDover

      clpo, clps, clpt = wingcl(b,bs,bo,
	lambdat,lambdas,gammat,gammas,
	sweep,AR,CL,CLhtail,fLo,fLt,
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
      Recoh = Reunit*coh
	CDhtail,CDhover = surfcd(S,
	bh,boh,boh,lambdah,1.0,sweeph,coh, 
	cdft,cdpt,Recoh,Rereft,aRexp,rkSunsh,
	fCDhcen)

#---- vertical tail profile CD
      Recov = Reunit*cov
      CDvtail1,CDvover1 = surfcd(S,
	bv,bov,bov,lambdav,1.0,sweepv,cov, 
	cdft,cdpt,Recov,Rereft,aRexp,rkSunsv,
      fCDvcen)
	
      CDvtail = CDvtail1*nvtail

      para[iaCDhtail] = CDhtail
      para[iaCDvtail] = CDvtail

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- fuselage profile CD, using fuselage CDA from BLI calculation
      CDfuse = PAfinf/S
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
      cdfs = para[iacdfs]
      cdps = para[iacdps]
      rVstrut = wing.strut_local_velocity_ratio
      CDstrut = (Sstrut/S)*(cdfs + cdps*cosLs^3) * rVstrut^3
      para[iaCDstrut] = CDstrut

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- induced CD
#      if(Ldebug) write(*,*) '...calling CDITRP...'

      cditrp(pari,parg,para, wing)
      CDi = para[iaCDi]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- apparent fuselage deltaCD from ingestion
      fBLIf = parg[igfBLIf]
      dCDBLIf = -fBLIf*DAfwake/S
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
      CD_components = [CDi  CDfuse  CDwing  CDover CDhtail  CDvtail  CDstrut 	CDnace dCDBLIf dCDBLIw]
      # println(CD_components)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#      if(Ldebug) write(*,*) '...exiting CDSUM...'
      # println("Total CD = ", CD)
      # println("PARA = ", para)
      return
end # cdsum


"""
      cditrp(pari,parg,para, wing)

Computes the induced drag via the Trefftz plane. Calls [`trefftz1`](@ref).

!!! details "🔃 Inputs and Outputs"
      **Inputs:**
      - `pari::AbstractVector{Int64}`: Vector of `aircraft` model integer/flag parameters.
      - `parg::AbstractArray{Float64}`: Vector of `aircraft` model geometry parameters.
      - `para::AbstractArray{Float64}`: Vector of `aircraft` model aerodynamic parameters.

      **Outputs:**
      - No explicit outputs. Computed induced drag value and span efficiency are saved to `para` of `aircraft` model.

!!! compat "Future Changes"
      In an upcoming revision, an `aircraft` struct and auxiliary indices will be passed in lieu of pre-sliced `par` arrays.

"""
function cditrp(pari,parg,para, wing)

      CL = para[iaCL]

      if (CL == 0.0) 
	    para[iaCDi]  = 0.
	    para[iaspaneff] = 1.0
       return
      end

      CLhtail = para[iaCLh]*parg[igSh]/parg[igS]
      # println("CLhtail: $(para[iaCLh]) $(parg[igSh]) $(parg[igS])")
      bref = parg[igb]
      Sref = parg[igS]

      Mach = para[iaMach]

      Lspec = true

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
      fLo = parg[igfLo]
#      fLo = 0.0

#---- span, wing-break span, wing-root span
      b[1]  = parg[igb]
      bs[1] = parg[igbs]
      bo[1] = wing.layout.box_halfspan

#---- span of wing-root streamline in Trefftz Plane
      bop[1] = wing.layout.box_halfspan * 0.2

      zcent[1]  = wing.layout.z_wing
      gammas[1] = wing.layout.λs*para[iarcls]
      gammat[1] = wing.layout.λt*para[iarclt]
      po[1]     = 1.0
      CLsurfsp[1] = CL - CLhtail

#---- velocity-change fractions at wing spanwise locations due to fuselage flow
      fduo = para[iafduo]
      fdus = para[iafdus]
      fdut = para[iafdut]

#---- horizontal tail wake parameters
      b[2]   = parg[igbh]
      bs[2]  = parg[igboh]
      bo[2]  = parg[igboh]
      bop[2] = parg[igboh]

      zcent[2] = parg[igzhtail]
      gammas[2] = 1.0
      gammat[2] = parg[iglambdah]
      po[2]     = 1.0
      CLsurfsp[2] = CLhtail


#---- number of surfaces  (wing, horizontal tail)
      nsurf = 2

#---- number of spanwise intervals
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

#     npout[1] = 40  # outer panel
#     npinn[1] = 12  # inner panel
#     npimg[1] = 6   # image inside fuselage
# 
#     npout[2] = 20  # outer panel
#     npinn[2] = 0   # inner panel
#     npimg[2] = 4   # image inside fuselage  (or inner panel if T-tail)

#      npout[1] = 160  # outer panel
#      npinn[1] = 48   # inner panel
#      npimg[1] = 24   # image inside fuselage
# 
#      npout[2] = 80  # outer panel
#      npinn[2] = 0   # inner panel
#      npimg[2] = 16  # image inside fuselage  (or inner panel if T-tail)

      ktip = 16
      #CLsurf = zeros(Float64, nsurf)
      # println("$nsurf, $npout, $npinn, $npimg, 
	# $Sref, $bref,
	# $b,$bs,$bo,$bop, $zcent,
	# $po,$gammat,gammas, $fLo, $ktip,
      # $Lspec,$CLsurfsp")


      CLsurf, CLtp, CDtp, sefftp = trefftz1(nsurf, npout, npinn, npimg, 
	Sref, bref,
	b,bs,bo,bop, zcent,
	po,gammat,gammas, fLo, ktip,
      Lspec,CLsurfsp, t, y, yp, z, zp, gw, yc, ycp, zc, zcp, gc, vc, wc, vnc)
      
      # println("$CLsurf, $CLtp, $CDtp, $sefftp")

      para[iaCDi] = CDtp
      para[iaspaneff] = sefftp

      return
end # cditrp

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


