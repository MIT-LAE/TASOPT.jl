# Inputs for testing runs
nmisx = 1
pari = zeros(Int64, iitotal)
parg = zeros(Float64, igtotal)
parm = zeros(Float64, (imtotal, nmisx))
para = zeros(Float64, (iatotal, iptotal, nmisx))
pare = zeros(Float64, (ietotal, iptotal, nmisx))

ft_to_m = 0.3048
in_to_m = 0.0254
nmi_to_m = 1852.0
deg_to_rad = π/180.0
lbf_to_N = 4.448222
kts_to_mps = 0.51444

pari[iifuel]   = 1 # 1 = H2 120 MJ/kg
pari[iifwcen ] = 1 # Fuel in center box
pari[iiwplan ] = 1 # wing cantilever with engine 
pari[iiengloc] = 2 # engines on fuselage
pari[iiengwgt] = 2 # advanced tech for eng weight
pari[iiBLIc  ] = 0 # core in clean flow
# pari[iifclose] = 0 # fuse tapers to a point
pari[iifclose] = 1 # fuse end is flat
# pari[iiHTsize] = 1 # set Sh via Vh 
pari[iiHTsize] = 2 # set Sh via igCLhCGfwd at max-forward CG during landing

# pari[iixwmove] = 0 # ixwmove   fix wing position 
# pari[iixwmove] = 1 # ixwmove   move wing to get CLh=CLhspec in cruise 
pari[iixwmove] = 2 # ixwmove   move wing to get min static margin = SMmin

pari[iiVTsize] = 1

 
pax = 220
seat_pitch = 30.0 * in_to_m
seat_width = 19.0 * in_to_m
aisle_halfwidth = 10.0 * in_to_m # per CFR § 25.815 

# parm[imwOpt   , :]  .=
parm[imRange  , :]  .= 3000.0 * nmi_to_m       # m
parm[imWpay   , :]  .= pax * 215.0 * lbf_to_N  # [N]
parm[imaltTO  , :]  .= 0.0
parm[imT0TO   , :]  .= 288.0 # dK
parm[imgamVCB , :]  .=  3.0 * π/180.0
parm[imgamVDE1, :]  .= -3.0 * π/180.0
parm[imgamVDEn, :]  .= -3.0 * π/180.0
parm[imthCB   , :]  .= 40.0 * π/180.0


# Geometry parameters
    # parg[igRange   ] = parm[imRange, 1]
    # Take-off and initial climb parameters
        parg[igcdefan  ] = 0.5       # cdefan    CDA_fan /A_fan  of dead engine in engine-out climb           
        parg[igCDgear  ] = 0.015     # CDgear    CDA_gear/S      during climb           
        parg[igCDspoil ] = 0.10      # CDspoiler CDA_spoiler/S   during braking            
        parg[igmuroll  ] = 0.025     # muroll    rolling-resistance coefficient         
        parg[igmubrake ] = 0.35      # mubrake   braking-resistance coefficient      
        parg[ighobst   ] = 35.0   * ft_to_m
        parg[iglBFmax  ] = 8000.0 * ft_to_m

        # Noise stuff
        parg[iggtocmin ] = 0.015 #  gtocmin  specified min top-of-climb gradient
        parg[igdBSLmax ] = 90.0  # dBSLmax  max dBA for sideline   (not used)
        parg[igdBCBmax ] = 75.0  # dBSLmax  max dBA for cutback    (not used)


    # Sizing load parameters
        parg[igNlift   ] = 3.0           # Nlift  max vertical load factor for wing bending loads
        parg[igNland   ] = 6.0           # Nland  max vertical load factor for fuse bending loads
        parg[igVne     ] = 280.0 * kts_to_mps # Vne    never-exceed IAS, for tail loads

        cabinPressureAlt = 8000.0 # Altitude for cabin pressure [ft]
        _, p_cabin, _, _, _ = atmos(cabinPressureAlt*ft_to_m/1000.0)

        parg[igpcabin  ] = p_cabin

    # Wing 
        parg[igsweep   ] =  27.567         # λ wing
        parg[igAR      ] =  10.4411        # Aspect ratio
        parg[igbmax    ] = 117.5  * ft_to_m # Max span for span constraint ICAO Code D/ FAA Group IV
        # parg[igbmax    ] = 171.0  * ft_to_m # Max span for span constraint ICAO Code E/ FAA Group V
        
        parg[iglambdas ] = 0.8784
        parg[iglambdat ] = 0.1503

        parg[igfLo     ] = -0.3   # fLo   fuselage lift carryover loss factor
        parg[igfLt     ] = -0.05  # fLt   tip lift rolloff factor

        parg[igzs      ] = 154.0 * in_to_m

        parg[igbo      ] = 2 * (120.0 * in_to_m) # 2 × wing centerbox halfspan
        parg[igetas    ] = 0.29 # ηs panel break eta location  (strut-attach if iwplan=2)
        parg[igrVstrut ] = 1.0   # Strut local/free sream velocity

        # Structural box
            parg[igwbox    ] = 0.50   # wbox    box width/c
            parg[ighboxo   ] = 0.13 # hboxo   box height/c  (airfoil t/c) at root
            parg[ighboxs   ] = 0.10 # hboxs   box height/c  (airfoil t/c) at break and tip
            parg[igrh      ] = 0.75   # rh      web-height/hbox ratio
            parg[igXaxis   ] = 0.40   # Xaxis   spar box axis x/c location
            parg[ighstrut  ] = 0.15   # hstrut  strut t/c    (used only if iwplan=2)
        # Weight fractions of flight surfaces and secondary wing components, 
        # as fractions of dry total wing weight
            parg[igfflap   ] = 0.20
            parg[igfslat   ] = 0.10
            parg[igfaile   ] = 0.04
            parg[igflete   ] = 0.10
            parg[igfribs   ] = 0.15
            parg[igfspoi   ] = 0.02
            parg[igfwatt   ] = 0.03

    # Tails
        # Horizontal Tail
            parg[igVh      ] = 1.45 # HT volume coeff
            parg[igCLhCGfwd] = -0.305
            parg[igCLhspec ] = -0.02
            parg[igSMmin   ] = 0.05
            parg[igdepsda  ] = 0.60
            parg[igARh     ] = 6.0
            parg[iglambdah ] = 0.25
            parg[igsweeph  ] = 25.0
            parg[igboh     ] = 2 * (2.5 * ft_to_m) # 2 × half span
            parg[igfCDhcen ] = 0.1
            parg[igCLhmax  ] = 2.0
            
            parg[igwboxh   ] = 0.50
            parg[ighboxh   ] = 0.14
            parg[igrhh     ] = 0.75

            parg[igfhadd   ] = 0.30

            parg[igsweeph ] = 8.0#parg[igsweep] #Copy wing values

            parg[igCLhNrat ] = -0.5 #CLh/ CLmax
    
        # Vertical Tail
            parg[ignvtail  ] = 2.0

            parg[igVv      ] = 0.10
            parg[igCLveout ] =  0.453
            parg[igARv     ] = 2.0
            parg[iglambdav ] = 0.30
            parg[igsweepv  ] = 25.0
            parg[igbov     ] = 0.0
            parg[igCLvmax  ] = 2.6
            
            parg[igwboxv   ] = 0.50
            parg[ighboxv   ] = 0.14
            parg[igrhv     ] = 0.75

            parg[igfvadd   ] = 0.40


        parg[igdCLnda  ] = 3.8

    # Cabin and Fuselage and PowerTrain stuff

        parg[igRfuse   ] = 120.0 * in_to_m 
        parg[igdRfuse  ] = 15.0 * in_to_m
        parg[igwfb     ] =  0.0 * in_to_m
        parg[ignfweb   ] =  1.0
        parg[ighfloor  ] =  5.0 * in_to_m

        parg[iganose   ] = 1.65
        parg[igbtail   ] = 2.0

        seats_per_row = 2*Int(parg[igRfuse] ÷ (seat_width + aisle_halfwidth/3))
        rows = ceil(pax ÷ seats_per_row)
        lcabin = rows * seat_pitch 
        println("Seats per row = $seats_per_row, rows = $rows, lcabin = $(lcabin/ft_to_m) ft")

        parg[igxnose   ] =   0.0 * ft_to_m
        parg[igxblend1 ] =  20.0 * ft_to_m
        parg[igxshell1 ] =  17.0 * ft_to_m
        parg[igxshell2 ] = parg[igxshell1] + lcabin + 30.0*ft_to_m # 10ft for cockpit + 2 ends* 10 ft for galley + lavatory 

        ltank = 20.0 * ft_to_m
        parg[igxftank ]  = parg[igxshell2] + ltank/2
        parg[igxblend2 ] = parg[igxshell2] + ltank 

        ltshaft = 9.0 * ft_to_m # length of T406 ~ 6.5 ft + 2.5 ft margin
        lgen    = 5.0 * ft_to_m 
        parg[igxtshaft]  = parg[igxblend2] + ltshaft/2
        parg[igxgen   ]  = parg[igxblend2] + ltshaft + lgen/2
        parg[igxcat   ]  = parg[igxgen   ]

        parg[igxconend ] = parg[igxgen] + lgen/2 + 5.0*ft_to_m
        parg[igxend    ] = parg[igxconend] + 10.0*ft_to_m # 10 ft margin/ other things not accounted for

        parg[igxwbox   ] =  65.0 * ft_to_m  # x location of wing box
        parg[igxhbox   ] = parg[igxconend ] #- 3*ft_to_m
        parg[igxvbox   ] = parg[igxconend ] #- 5*ft_to_m

        parg[igxinv   ]  =  60.0 * ft_to_m
        parg[igxmot   ]  =  57.0 * ft_to_m
        parg[igxfan   ]  =  55.0 * ft_to_m
        
        parg[igzwing   ] = -5.5 * ft_to_m
        parg[igzhtail  ] =  0.0 * ft_to_m
        # parg[igzhtail  ] =  13.0 * ft_to_m

        parg[igneng    ] =  parpt[ipt_nfan] # Represents ducted fans + motors for TE config

        parg[igxeng    ] = 52.0 * ft_to_m
        parg[igyeng    ] = 50.0 * ft_to_m
        # parg[igneng    ] =  2.0

        parg[iglambdac ] =  0.3 # Tail cone taper ratio
        
        parg[igfstring ] = 0.35
        parg[igfframe  ] = 0.25
        parg[igffadd   ] = 0.20

        parg[igWfix    ] = 3000.0 * 4.45  # cockpit, pilots etc converted to [N]
        parg[igxfix    ] =    7.0 * ft_to_m

        # Insulation weight per length or area
            parg[igWpwindow] = 145.0 * 3.0 #[N/m]
            parg[igWppinsul] =  22.0       #[N/m²]
            parg[igWppfloor] =  60.0       #[N/m²]
        
        # Bending moment inertia due to HT and VT
        parg[igrMh     ] = 0.4
        parg[igrMv     ] = 0.7

        parg[igCMVf1   ] =  2390.0 * 0.0283  # CMVf1  fuselage moment volume derivative  d(Mfuse/q)/dCL
        parg[igCLMf0   ] =  0.185            # CLMf1  CL where Mfuse = 0

        para[iafduo, :, :] .= 0.018    # fduo   fuselage velocity overspeed at wing root
        para[iafdus, :, :] .= 0.014    # fdus   fuselage velocity overspeed at wing break
        para[iafdut, :, :] .= 0.0045   # fdut   fuselage velocity overspeed at wing tip

    # Landing gear weight fractions and locations
        parg[igxhpesys ] =  62.0 * ft_to_m   #  xhpesys   hyd/pneu/ele system location
        parg[igxlgnose ] =  14.0 * ft_to_m   #  xlgnose   nose LG location
        parg[igdxlgmain] =   1.0 * ft_to_m   # dxlgmain   main LG offset behind wing lift centroid

        parg[igfhpesys ] =  0.010     # fhpesys    Whpesys/WMTO
        parg[igflgnose ] =  0.011     # flgnose    Wlgnose/WMTO
        parg[igflgmain ] =  0.044     # flgmain    Wlgmain/WMTO

        parg[igxapu    ] = parg[igxconend] * ft_to_m # xapu      APU location
        parg[igfapu    ] = 0.035          # fapu   Wapu/Wpay     APU weight fraction

        parg[igfreserve] = 0.20  # freserve Wfreserve/Wburn
        parg[igfpadd   ] = 0.35  # fpadd    Wpadd/Wpay    other payload-proportional fraction
        parg[igfseat   ] = 0.10  # fseat    Wseat/Wpay    seat weight fraction
        parg[igfeadd   ] = 0.10  # feadd    Weadd/Wbare   engine accessories, fuel system fraction 
        parg[igfpylon  ] = 0.10  # fpylon   Wpylon/We+a+n engine pylon weight fraction   

    # Allowable stresses
        parg[igsigfac  ] = 1.0   #  sigfac   convenient multiplier on all the stress values below                      
                                                
        parg[igsigskin ] = 15000.0 / 0.000145   # sigskin   fuselage pressurization skin stress                      
        parg[igsigbend ] = 30000.0 / 0.000145   # sigbend   fuselage bending skin+stringer stress                      
                                                    
        parg[igsigcap  ] = 30000.0 / 0.000145   # sigcap    wing,tail bending caps                      
        parg[igtauweb  ] = 20000.0 / 0.000145   # tauweb    wing,tail shear webs                      
        parg[igsigstrut] = 30000.0 / 0.000145   # sigstrut  strut       
        
        parg[igrEshell ] = 1.0                  # rEshell   Ebend/Eskin  ratio
        parg[igEcap    ] =  10.0e6 / 0.000145   # Ecap     wing sparcap
        parg[igEstrut  ] =  10.0e6 / 0.000145   # Estrut   strut

        parg[igrhoskin ] =  2700.0  #  rhoskin     fuselage skin
        parg[igrhobend ] =  2700.0  #  rhobend     fuselage bending stringers 
        parg[igrhocap  ] =  2700.0  #  rhocap  	wing, tail bending caps	 
        parg[igrhoweb  ] =  2700.0  #  rhoweb  	wing, tail shear webs	 
        parg[igrhostrut] =  2700.0  #  rhostrut	strut   

        parg[igrhofuel]  = ρmix(0.1, 1.5)#817.0 #kg/m³
        # if pari[iifuel] == 0
        #     parg[igrhofuel]  = 817.0 #kg/m³
        #     pare[iehfuel] = 43.0 #MJ/kg
        # elseif pari[iifuel] == 1
        #     # parg[igrhofuel] = ρmix(0.1, 1.5)
        #     pare[iehfuel]= 120.0
        # end
        
    # Nacelle Drag stuff
        parg[igrSnace  ] = 16.0   # rSnace   nacelle+pylon wetted area/fan area  Snace/Afan
        parg[igrVnace  ] =  1.02  # rVnace   nacelle local/freesteam velocity ratio

    parg[igrWfmax  ] = 0.90

# Aerodynamic parameters

para[iaalt, ipcruise1, :] .=  33000.0 * ft_to_m # Cruise altitude [m] Max fuel max payload

# Takeoff and initial climb and descent
    para[iaalt,    ipstatic:ipcutback, :] .= parm[imaltTO]
    para[iaclpmax, ipstatic:ipcutback, :] .= 2.25 # clpmax   wing max cl_perp  = CLmax/cos(sweep)^2
    
    para[iaalt,    ipclimb1, :] .= parm[imaltTO]
    para[iaclpmax, ipclimb1, :] .= 2.25 # clpmax   wing max cl_perp  = CLmax/cos(sweep)^2
    
    para[iaalt,    ipdescentn, :] .= parm[imaltTO]
    para[iaclpmax, ipdescentn, :] .= 2.25 # clpmax   wing max cl_perp  = CLmax/cos(sweep)^2
    
# Cruise 
    para[iaCL  , ipclimb1+1:ipdescentn-1, :] .= 0.57067
    para[iaMach, ipclimbn:ipdescent1  , :] .= 0.80

# Wing span load parameters
    # Takeoff
        para[iarcls, 1:ipclimb1, :] .= 1.1     #  rcls    break/root cl ratio = cls/clo
        para[iarclt, 1:ipclimb1, :] .= 0.6     #  rclt    tip  /root cl ratio = clt/clo
        para[iacmpo, 1:ipclimb1, :] .= -0.20   #  cmpo    root  cm
        para[iacmps, 1:ipclimb1, :] .= -0.20   #  cmps    break cm
        para[iacmpt, 1:ipclimb1, :] .= -0.02   #  cmpt    tip   cm

    # Clean climb cruise descent and for wing structure sizing
        para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= 1.238   #  rcls    
        para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= 0.90    #  rclt    
        para[iacmpo, ipclimb1+1 : ipdescentn-1, :] .= -0.06   #  cmpo    
        para[iacmps, ipclimb1+1 : ipdescentn-1, :] .= -0.06   #  cmps    
        para[iacmpt, ipclimb1+1 : ipdescentn-1, :] .= -0.06   #  cmpt   
   
    # Landing, forward CG tail sizing case
        para[iarcls, ipdescentn, :] .= 1.1     #  rcls  
        para[iarclt, ipdescentn, :] .= 0.5     #  rclt  
        para[iacmpo, ipdescentn, :] .= -0.35   #  cmpo      
        para[iacmps, ipdescentn, :] .= -0.35   #  cmps      
        para[iacmpt, ipdescentn, :] .= -0.02   #  cmpt     

# Wing and tail cd's
    para[iacdfw  , 1 : iptotal, :] .=  0.0085  #  cdfw    wing profile cd for low speed (takeoff, initial climb)
    para[iacdpw  , 1 : iptotal, :] .=  0.0035  #  cdpw    
    para[iaRerefw, 1 : iptotal, :] .=  20.0e6  #  Rerefw
                                        
    para[iacdft  , 1 : iptotal, :] .=  0.0060  #  cdft    tail profile cd
    para[iacdpt  , 1 : iptotal, :] .=  0.0030  #  cdpt    
    para[iaRereft, 1 : iptotal, :] .=  10.0e6  #  Rereft  
                                       
    para[iacdfs  , 1 : iptotal, :] .=  0.0085  #  cdfs    strut profile cd (not used if there's no strut)
    para[iacdps  , 1 : iptotal, :] .=  0.0035  #  cdps    
    para[iaRerefs, 1 : iptotal, :] .=  1.0e6   #  Rerefs  
                                       
    para[iaaRexp , 1 : iptotal, :] .=  -0.15   #  aRexp   exponent for Re-scaling:  CD = cd * (Re/Re_ref)^aRexp
                                       
    para[iafexcdw, 1 : iptotal, :] .=  1.02     #  fexcdw   # wing excrescence drag factor
    para[iafexcdt, 1 : iptotal, :] .=  1.02     #  fexcdt   # tail excrescence drag factor
    para[iafexcdf, 1 : iptotal, :] .=  1.03     #  fexcdf   # fuse excrescence drag factor

# Engine parameters

pare[ieTt4, 1:iptotal, :] .= 3200.0 # [R]
pare[ieTt4, ipstatic:iptakeoff, :] .= 3300.0 #[R]

parg[igfTt4CL1] = 0.2
parg[igfTt4CLn] = 0.2