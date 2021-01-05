"""
Calculates aircraft CD components for operating point ipoint
If icdfun=1, computes wing cdf,cdp from airfoil database # iairf,
otherwise uses default values in para array.
"""
function cdsum!(pari,parg,para,pare, icdfun)
#      implicit real (a-h,l-z)
#      include 'index.inc'
#      integer pari(iitotal)
#      real parg[igtotal], 
#	para[iatotal],
#	pare(ietotal)
#      integer icdfun
#
#      include 'airf.inc'
#      include 'constants.inc'
#
#      logical Ldebug
#
#      common /com_lu/ ilu


      Ldebug = false
#      Ldebug = true

      AR       = parg[igAR     ]
      sweep    = parg[igsweep  ]
      hboxo    = parg[ighboxo  ]
      hboxs    = parg[ighboxs  ]
      hboxt    = hboxs
      fSnace   = parg[igfSnace ]
      bo       = parg[igbo     ]
      bs       = parg[igbs     ]
      zs       = parg[igzs     ]
      boh      = parg[igboh    ]
      bov      = parg[igbov    ]
      hstrut   = parg[ighstrut ]


      Rfuse    = parg[igRfuse  ]
      dRfuse   = parg[igdRfuse ]
      wfb      = parg[igwfb    ]

      xnose    = parg[igxnose  ]
      xend     = parg[igxend   ]
      xblend1  = parg[igxblend1]
      xblend2  = parg[igxblend2]
      xhtail   = parg[igxhtail ]
      xvtail   = parg[igxvtail ]
      xwing    = parg[igxwing  ]

      lambdac  = parg[iglambdac]

      lambdat  = parg[iglambdat]
      lambdas  = parg[iglambdas]
      gammat   = parg[iglambdat]*para[iarclt]
      gammas   = parg[iglambdas]*para[iarcls]

      ARh      = parg[igARh    ]
      ARv      = parg[igARv    ]
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
      if(icdfun==1) 
#----- integrated across span for CDwing

      # if(Ldebug) write(*,*) 'calling SURFCD2...'
      clpo,clps,clpt,
	cdfw,cdpw,CDwing,CDover = surfcd2(S, b,bs,bo,
	lambdat,lambdas,gammat,gammas,
	hboxo,hboxs,hboxt,
	Mach,sweep,co, 
	CL,CLhtail, fLo,fLt,
	Reco,aRexp, rkSunsw,fexcdw,
	AMa,Acl,Aτ,ARe,A,
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

#      write(*,*) ip
#      write(*,*) 'cdf cdp', cdf, cdp

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -


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
#      write(*,*) 'CDhtail', CDhtail
#
#---- vertical tail profile CD
      Recov = Reunit*cov
      CDvtail1,CDvover1 = surfcd(S,
	bv,bov,bov,lambdav,1.0,sweepv,cov, 
	cdft,cdpt,Recov,Rereft,aRexp,rkSunsv,
      fCDvcen)
	
      CDvtail = CDvtail1*nvtail
#      write(*,*) 'CDvtail', CDvtail
#
      para[iaCDhtail] = CDhtail
      para[iaCDvtail] = CDvtail

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- fuselage profile CD from wetted area
#      ltot = xend - xnose
#      Rel = Reunit*ltot
#      Cfwet = cfturb(Rel) * para[iafexcdf]
#      call bodycd(S,
#     &  Rfuse,dRfuse,wfb,xnose,xblend1,xblend2,xend, Cfwet,
#     &  CDfuse)

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
      if(Reunit == 0.0) 
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
      rVstrut = parg[igrVstrut]
      CDstrut = (Sstrut/S)*(cdfs + cdps*cosLs^3) * rVstrut^3
      para[iaCDstrut] = CDstrut

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- induced CD
#      if(Ldebug) write(*,*) '...calling CDITRP...'

      cditrp(pari,parg,para)
      CDi = para[iaCDi]

#      spaneff = 0.8825
#      CDi = CL^2 / (pi*AR*spaneff)
#      para[iaCDi] = CDi
#      para[iaspaneff] = spaneff
#      write(*,*)
#      write(*,*) para[iaCDi], para[iaspaneff]
#      write(*,*) CDi, spaneff



#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- apparent fuselage deltaCD from ingestion
      fBLIf = parg[igfBLIf]
      dCDBLIf = -fBLIf*DAfwake/S

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- apparent wing deltaCD from ingestion
      fBLIw = parg[igfBLIw]
      CDWwake = CDwing * 0.15  # assume 15% of the wing dissipation is in wake
      dCDBLIw = -fBLIw*CDWwake

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#---- total CD
      CD = CDi + CDfuse + CDwing + CDover + CDhtail + CDvtail + CDstrut + CDnace + dCDBLIf + dCDBLIw
      para[iaCD] = CD      
      CD_components = [CDi  CDfuse  CDwing  CDover CDhtail  CDvtail  CDstrut 	CDnace dCDBLIf dCDBLIw]

#      tau = hboxo
#      rsb3 = 0.25*((0.5*clpo + tau+1.0)    + (tau+1.0)   )
#     &           *((0.5*clpo + tau+1.0)^2 + (tau+1.0)^2)
#      Cdiss = 0.0009 * 1.5
#      aksb = 0.25
#      CDsb = 2.0*Cdiss*(rsb3-1.0)*(co^2/S)*2.0*aksb
#      CD = CD + CDsb
#      para[iaCD] = CD
#c      write(*,'(1x,3f12.7)') CDsb, CD, CDsb/CD

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#      if(Ldebug) write(*,*) '...exiting CDSUM...'
      # println("Total CD = ", CD)
      # println("PARA = ", para)
      return
      end # cdsum


"""
cditrip calcualtes the induced drag from the treftz plane
"""
function cditrp(pari,parg,para)

#      include("index.inc") #include the array indices for pari, parg and para
      
#      integer pari(iitotal)
#      real parg[igtotal], 
#	para[iatotal]
#
#      include 'trp.inc'
#
#      include 'time.inc'
#
#      logical Lspec
#      common /com_lu/ ilu
#

     ifclose = pari[iifclose]

      CL = para[iaCL]

      if(CL == 0.0) 
	    para[iaCDi]  = 0.
	    para[iaspaneff] = 1.0
       return
      end

      CLhtail = para[iaCLh]*parg[igSh]/parg[igS]

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
      bo[1] = parg[igbo]

#---- span of wing-root streamline in Trefftz Plane
      bop[1] = parg[igbo] * 0.2

      zcent[1]  = parg[igzwing]
      gammas[1] = parg[iglambdas]*para[iarcls]
      gammat[1] = parg[iglambdat]*para[iarclt]
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
      if(bo[2] == 0.0) 
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
      CLsurf, CLtp, CDtp, sefftp = trefftz1(nsurf, npout, npinn, npimg, 
	Sref, bref,
	b,bs,bo,bop, zcent,
	po,gammat,gammas, fLo, ktip,
      Lspec,CLsurfsp)
      
      # println("$CLsurf, $CLtp, $CDtp, $sefftp")

      para[iaCDi] = CDtp
      para[iaspaneff] = sefftp

      return
end # cditrp

function cfturb(Re)
# -------------------------------------------------------------
#     Returns total Cf for turbulent flat plate, versus Re_l
# -------------------------------------------------------------
#cc   cfturb = 0.427/(log10(re) - 0.407)**2.64   ! original Hoerner
#cc   cfturb = 0.310/(log10(re) - 0.407)**2.47   ! modified (weaker dependence on Re)
#
      cfturb = 0.523/(log(0.06*Re))^2           # White
#
      return cfturb
end # cfturb


