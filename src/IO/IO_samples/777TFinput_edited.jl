# Inputs for testing runs
nmisx = 1
pari = zeros(Int64, iitotal)
parg = zeros(Float64, igtotal)
parm = zeros(Float64, (imtotal, nmisx))
para = zeros(Float64, (iatotal, iptotal, nmisx))
pare = zeros(Float64, (ietotal, iptotal, nmisx))
parad = zeros(Float64, (iatotal, iptotal, nmisx))
pared = zeros(Float64, (ietotal, iptotal, nmisx))

ft_to_m = 0.3048
in_to_m = 0.0254
nmi_to_m = 1852.0
deg_to_rad = π/180.0
lbf_to_N = 4.448222
kts_to_mps = 0.51444
hp_to_W    = 745.7

# parg[igrhofuel]  = 817.0 #kg/m³
# parg[igLHVfuel]  = 43.2 # MJ/kg

pari[iiopt]    = 0 # 0 run sizing loop only; 1 run optimization procedure
pari[iifuel]   = 24 # 1 = H2 120 MJ/kg; 2= JetA; 3 = specied by composition (JANAF)
pari[iifwing]   = 1 # 0 = all fuel stored in tanks; 1 = all fuel stored in wings
pari[iifwcen]   = 1 # Fuel in center box 
pari[iiengtype] = 1  # 0 = Turboelectric engine; 1 = Turbofan engine
# pari[iicompfuse] = 0 # 0 = Aluminum fuselage; 1 = Composite fuselage
# pari[iicompwing] = 0 # 0 = Aluminum wingbox; 1 = Composite wingbox
# pari[iibleedeng] = 0 # 0 = customer bleed active; 1 = no bleed architecture
# pari[iiactype]  = 777 # 777 = 777 size aircraft; 787 = 787 size aircraft

# Boundary Layer Ingestion (BLI)
        # parg[igfBLIf] = 0.5
        parg[igfBLIf] = 0.0


pari[iiwplan ] = 1 # wing cantilever with engine 
pari[iiengloc] = 1 # engines on fuselage
pari[iiengwgt] = 2 # advanced tech for eng weight
pari[iiBLIc  ] = 1 # core in clean flow
pari[iifclose] = 0 # 0 = fuse tapers to a point; 1 = fuse end is flat

pari[iiHTsize] = 2 # 1 = set Sh via Vh; 2 = set Sh via igCLhCGfwd at max-forward CG during landing

# pari[iixwmove] = 0 # ixwmove   fix wing position 
# pari[iixwmove] = 1 # ixwmove   move wing to get CLh=CLhspec in cruise 
pari[iixwmove] = 2 # 0 = fix wing position; 1 = move wing to get Clh=CLhspec in cruise; 2 = move wing to get min static margin (SMmin)

pari[iiVTsize] = 1


pax = 370
paxmax = 670
seat_pitch = 30.0 * in_to_m  
seat_width = 19.0 * in_to_m
aisle_halfwidth = 10.0 * in_to_m # per CFR § 25.815 

# parm[imwOpt   , :]  .=
parm[imRange  , :]  .= 7870.0 * nmi_to_m       # m
parm[imWpay   , :]  .= pax * 230.0 * lbf_to_N  # [N]
parg[igWpaymax   ]   = paxmax * 230.0 * lbf_to_N 
parm[imaltTO  , :]  .= 0.0
parm[imT0TO   , :]  .= 288.0 # dK
parm[imgamVCB , :]  .=  3.0 * π/180.0
parm[imgamVDE1, :]  .= -3.0 * π/180.0
parm[imgamVDEn, :]  .= -3.0 * π/180.0
parm[imthCB   , :]  .= 40.0 * π/180.0


# Fuel Parameters
if pari[iifuel] == 1
    parg[igrhofuel]  = ρmix(0.1, 1.5)
    parg[igLHVfuel]  = 120 # MJ/kg
 elseif pari[iifuel] == 2
    parg[igrhofuel]  = 817.0 #kg/m³
     parg[igLHVfuel]  = 43.0 # MJ/kg
 elseif pari[iifuel] == 3

    """
    // Setup Thermodynamic package
    // +-----------------+----------+----+----+-------------------+
    // |       Fuel      |   Fuel   |  C |  H |  O | hfuel [J/mol] | |
    // +-----------------+----------+----+----+----------------------+
    // |     Hydrogen    |    H2    |  0 |  2 |    |         0       |
    // +-----------------+----------+----+----+----------------------+
    // | Vaporized Jet-A | Jet-A(g) | 12 | 23 |    |   -249657.0     |
    // +-----------------+----------+----+----+----------------------+
    // |     Methane     |    CH4   |  1 |  4 |    |    -74600.0     |
    // +-----------------+----------+----+----+----------------------+
    // |    Napthalene   |   C10H8  | 10 |  8 |    |   -150580.0     |
    // +-----------------+----------+----+----+----------------------+
    // |    Ethanol      |  C2H5OH  | 2  |  6 | 1  |   -234950.0     |
    // +-----------------+----------+----+----+----------------------+
    // |    Butanol      |  C4H10O  |  4 | 10 | 1  |   -281400.0     |
    // +-----------------+----------+----+----+----------------------+
    """
    
    # Octane
    # n_carbon   =  8.0
    # n_hydrogen =  18.0
    # n_oxygen   =  0.0
    # n_nitrogen =  0.0
    # hfuel            =  -208759.0
    # parg[igrhofuel]   =  814.0 #kg/m³

    # Butanol
    # n_carbon   =  4.0
    # n_hydrogen =  10.0
    # n_oxygen   =  1.0
    # n_nitrogen =  0.0
    # hfuel            =  -281400.0
    # parg[igrhofuel]   =  810.0 #kg/m³

    # Pentanol
    n_carbon   =  5.0
    n_hydrogen =  12.0
    n_oxygen   =  1.0
    n_nitrogen =  0.0
    hfuel            =  -298000.0
    parg[igrhofuel]   =  814.0 #kg/m³

    # JP
    # n_carbon   =  12.5
    # n_hydrogen =  n_carbon * 1.906349
    # n_oxygen   =  0.0
    # n_nitrogen =  0.0
    # parg[igrhofuel]   =  817.0 #kg/m³
    # hfuel            =  -249657.0

    # H2
    # n_carbon   =  0.0
    # n_hydrogen =  2.0
    # n_oxygen   =  0.0
    # n_nitrogen =  0.0
    # parg[igrhofuel]   =  ρmix(0.1, 1.5) #kg/m³
    # hfuel            = 0.0

    hf_CO2 = -393510 #J/mol
    hf_H2O = -241826 #J/mol (vapor)
    MW_C = 12.011e-3 #kg/mol
    MW_H = 1.0008e-3 #kg/mol
    MW_O = 15.99940e-3 #kg/mol
    MW_N = 14.0067e-3 #kg/mol

    #Calcualte the molecular weight of the fuel and mass fractions of C and H
    MW = n_carbon*MW_C + n_hydrogen*MW_H + n_oxygen*MW_O + n_nitrogen*MW_N #kg/mol
    YC = n_carbon*MW_C/MW
    YH = n_hydrogen*MW_H/MW
    YO = n_oxygen*MW_O/MW
    YN = n_nitrogen*MW_N/MW

    #Calculate the LHV - this assumes fuel is of the form CnHm + m/4 O2 --> nCO2 + m/2 H2O
    #Assume fuel is CnHmOz --> nCO2 + m/2 H2O 
    LHV = (abs(n_carbon*hf_CO2 + n_hydrogen/2*hf_H2O - hfuel)/MW)/(10^6);

    parg[ighfuel] = hfuel
    parg[ig_YC] = YC
    parg[ig_YH] = YH
    parg[ig_YO] = YO
    parg[ig_YN] = YN
    parg[ig_MW] = MW
    parg[igTfuel] = 250.0 #Temperature of fuel in K
    parg[igLHVfuel] = LHV


end


# Geometry parameters
    # parg[igRange   ] = parm[imRange, 1]
    # Take-off and initial climb parameters
        parg[igcdefan  ] = 0.5       # cdefan    CDA_fan /A_fan  of dead engine in engine-out climb           
        parg[igCDgear  ] = 0.015     # CDgear    CDA_gear/S      during climb           
        parg[igCDspoil ] = 0.10      # CDspoiler CDA_spoiler/S   during braking            
        parg[igmuroll  ] = 0.025     # muroll    rolling-resistance coefficient         
        parg[igmubrake ] = 0.35      # mubrake   braking-resistance coefficient      
        parg[ighobst   ] = 35.0   * ft_to_m
        parg[iglBFmax  ] = 8800.0 * ft_to_m

        # Noise stuff
        parg[iggtocmin ] = 0.015 #  gtocmin  specified min top-of-climb gradient
        parg[igdBSLmax ] = 90.0  # dBSLmax  max dBA for sideline   (not used)
        parg[igdBCBmax ] = 75.0  # dBSLmax  max dBA for cutback    (not used)


    # Sizing load parameters
        parg[igNlift   ] = 3.0           # Nlift  max vertical load factor for wing bending loads
        parg[igNland   ] = 6.0           # Nland  max vertical load factor for fuse bending loads
        parg[igVne     ] = 280.0 * kts_to_mps # Vne    never-exceed IAS, for tail loads

        # if pari[iicompfuse] == 0
        #     cabinPressureAlt = 8000.0 # Altitude for cabin pressure [ft]
        # elseif pari[iicompfuse] == 1
        cabinPressureAlt = 6000.0 # Altitude for cabin pressure [ft]
        # end
        _, p_cabin, _, _, _ = atmos(cabinPressureAlt*ft_to_m/1000.0)

        parg[igpcabin  ] = p_cabin
        
        # Wing 
        parg[igsweep   ] =  32.583         # λ wing
        # parg[igsweep   ] =  30.0
        parg[igAR      ] =  8.455        # Aspect ratio
        # parg[igAR      ] =  9.1        # Aspect ratio

        parg[igbmax    ] = 213.0  * ft_to_m # Max span for span constraint ICAO Code D/ FAA Group IV
        # parg[igbmax    ] = 117.5*1.10  * ft_to_m # Max span for span constraint ICAO Code D/ FAA Group IV + 10% foldable wing tips based on the 777x dimensions 64.82 m ⇾ 71.75 m
        # parg[igbmax    ] = 171.0  * ft_to_m # Max span for span constraint ICAO Code E/ FAA Group V
        
        parg[igzwing   ] = -7.0 * ft_to_m
        
        parg[iglambdas ] = 0.78
        parg[iglambdat ] = 0.175

        parg[igfLo     ] = -0.3   # fLo   fuselage lift carryover loss factor
        parg[igfLt     ] = -0.05  # fLt   tip lift rolloff factor

        parg[igzs      ] = 154.0 * in_to_m #for strut- not used

        # parg[igbo      ] = 2 * (120.0 * in_to_m) # 2 × wing centerbox halfspan
        parg[igetas    ] = 0.32 # ηs panel break eta location  (strut-attach if iwplan=2)
        parg[igrVstrut ] = 1.0   # Strut local/free sream velocity

        # Structural box
            parg[igwbox    ] = 0.50   # wbox    box width/c
            parg[ighboxo   ] = 0.14466 # hboxo   box height/c  (airfoil t/c) at root
            parg[ighboxs   ] = 0.1395 # hboxs   box height/c  (airfoil t/c) at break and tip
            parg[igrh      ] = 0.75   # rh      web-height/hbox ratio
            parg[igXaxis   ] = 0.40   # Xaxis   spar box axis x/c location
            parg[ighstrut  ] = 0.15   # hstrut  strut t/c    (used only if iwplan=2)
            # Weight fractions of flight surfaces and secondary wing components, 
            # as fractions of dry total wing weight
            parg[igfflap   ] = 0.20
            parg[igfslat   ] = 0.10
            parg[igfaile   ] = 0.04
            parg[igflete   ] = 0.10
            parg[igfribs   ] = 0.15 #+ 0.01 #Adding this for folding wingtips
            parg[igfspoi   ] = 0.02
            parg[igfwatt   ] = 0.03

    # Tails
        # Horizontal Tail
            parg[igVh      ] = 0.85 # HT volume coeff
            parg[igCLhCGfwd] = -1.0
            parg[igCLhspec ] = -0.02
            parg[igSMmin   ] = 0.15
            parg[igdepsda  ] = 0.60
            parg[igARh     ] = 4.8
            parg[iglambdah ] = 0.32
            parg[igsweeph  ] = 33.0
            parg[igboh     ] = 2 * (5. * ft_to_m) # 2 × half span
            parg[igfCDhcen ] = 1.0
            parg[igCLhmax  ] = 2.0
            
            parg[igwboxh   ] = 0.50
            parg[ighboxh   ] = 0.14
            parg[igrhh     ] = 0.75

            parg[igfhadd   ] = 0.30

            # parg[igsweeph ] = 8.0#parg[igsweep] #Copy wing values

            parg[igCLhNrat ] = -0.5 #CLh/ CLmax
    
        # Vertical Tail
            parg[ignvtail  ] = 1.0

            parg[igVv      ] = 0.06
            parg[igCLveout ] = 1.0
            parg[igARv     ] = 2.35
            parg[iglambdav ] = 0.25
            parg[igsweepv  ] = 28.0
            parg[igbov     ] = 0.0
            parg[igCLvmax  ] = 2.6
            
            parg[igwboxv   ] = 0.50
            parg[ighboxv   ] = 0.14
            parg[igrhv     ] = 0.75

            parg[igfvadd   ] = 0.40

# What is up with this
        # parg[igdCLnda  ] = 1.0# 3.8
        parg[igdCLnda  ] = 3.8

    # Cabin and Fuselage and PowerTrain stuff

        parg[igRfuse   ] = 122.0 * in_to_m 
        parg[igbo      ] = 2 * (120 * in_to_m) # 2 × wing centerbox halfspan
        parg[igdRfuse  ] =  0.0 * in_to_m
        parg[igwfb     ] =  0.0 * in_to_m
        parg[ignfweb   ] =  1.0
        parg[ighfloor  ] =  8.0 * in_to_m

        parg[iganose   ] = 1.65
        parg[igbtail   ] = 2.0

        seats_per_row = 2*Int(parg[igRfuse] ÷ (seat_width + aisle_halfwidth/3))
        rows = Int(ceil(pax / seats_per_row))
        lcabin = rows * seat_pitch 
        println("Seats per row = $seats_per_row, rows = $rows, lcabin = $(lcabin/ft_to_m) ft")

        parg[igxnose   ] =   0.0 * ft_to_m
        parg[igxblend1 ] =  40.0 * ft_to_m
        parg[igxshell1 ] =  38.0 * ft_to_m
        # parg[igxshell2 ] = parg[igxshell1] + lcabin + 20.0*ft_to_m + 2*seat_pitch # 2 ends* 10 ft (for galley (6ft) + lavatory (4ft length) ) + space for emergency_exit
        parg[igxshell2 ] = 204.0 * ft_to_m

        if pari[iifwing] == 1
            ltank = 0.0
        end
        # ltank = 10.7
        parg[iglftankin] = ltank
        # ltank = 58.0 * ft_to_m 
        # parg[igxftank ]  = parg[igxshell2] + ltank/2 + 1.0*ft_to_m #(buffer)
        # parg[igxblend2 ] = parg[igxftank]  + ltank/2
        # parg[igxblend2 ] = 171.0 * ft_to_m

        # parg[iglftank] = 0.0
        # parg[igxftankfront ]  = 0.0
        # parg[igxftankaft ]  = 0.0

         
        if pari[iiengtype] == 0
            ltshaft = 9.0 * ft_to_m # length of T46 ~ 6.5 ft + 2.5 ft margin
            lgen    = 5.0 * ft_to_m  
            parg[igxtshaft]  = parg[igxblend2] + ltshaft/2
            parg[igxgen   ]  = parg[igxblend2] + ltshaft + lgen/2
            parg[igxcat   ]  = parg[igxgen   ]
            parg[igxeng    ] = parg[igxtshaft]

        elseif pari[iiengtype] == 1
            ltshaft = 0.0
            lgen = 0.0
            parg[igxtshaft]  = 0.0
            parg[igxgen   ]  = 0.0
            parg[igxcat   ]  = 0.0  
        end


        # parg[igxconend ] = parg[igxblend2] +  5.0*ft_to_m
        # parg[igxend    ] = parg[igxconend] + 5.0*ft_to_m # 5 ft margin/ other things not accounted for

        parg[igxconend ] = 235.0 * ft_to_m
        parg[igxend    ] = 242.0 * ft_to_m

        parg[igxwbox   ] =  114.0 * ft_to_m  # x location of wing box\
        parg[igxeng    ] =  102.0 * ft_to_m 
        
        parg[igxhbox   ] = parg[igxconend ] - 10*ft_to_m
        parg[igxvbox   ] = parg[igxconend ] - 23*ft_to_m

        # parg[igxinv   ]  =  60.0 * ft_to_m
        # parg[igxmot   ]  =  parg[igxwbox] # 57.0 * ft_to_m
        # parg[igxfan   ]  =  parg[igxwbox] # 55.0 * ft_to_m
        
        parg[igzhtail  ] =  9.0 * ft_to_m
        # parg[igzhtail  ] =  13.0 * ft_to_m
        # parg[igzhtail  ] =  7.0 #5.0 * ft_to_m

        # parg[igneng    ] =  parpt[ipt_nfan] # Represents ducted fans + motors for TE config
        parg[igneng    ] = 2.0
    
        parg[igyeng    ] = 31.5 * ft_to_m
        # parg[igneng    ] =  2.0

        parg[iglambdac ] =  0.3 # Tail cone taper ratio
        
        parg[igfstring ] = 0.34
        parg[igfframe  ] = 0.24
        parg[igffadd   ] = 0.20

        parg[igWfix    ] = 3000.0 * 4.45  # cockpit, pilots etc converted to [N]
        parg[igxfix    ] =    10.0 * ft_to_m

        # Insulation weight per length or area
            parg[igWpwindow] = 145.0 * 3.0 #[N/m]
            parg[igWppinsul] =  40.0       #[N/m²]
            parg[igWppfloor] =  60.0       #[N/m²]
        
        # Bending moment inertia due to HT and VT
        parg[igrMh     ] = 0.4
        parg[igrMv     ] = 0.7

        # parg[igCMVf1   ] =  7470.0 * 0.0283  # CMVf1  fuselage moment volume derivative  d(Mfuse/q)/dCL (0.0283 is conversion from ft³ to m³)
        # parg[igCMVf1   ] =  127.0  #m³ based on cab vol of 800 m³ CMVf1 ≈ 2𝒱/(∂Cl/∂α), where 𝒱 = fuselage volume
        # parg[igCMVf1   ] =  2390.0 *0.0283 #60.0 
        
        #What is this too (igCMVf1 value originally at 83)
        parg[igCMVf1   ] =  7470.0 * 0.0283 #60.0  
        # parg[igCMVf1   ] =  83.0 #60.0  
        parg[igCLMf0   ] =  0.185            # CLMf1  CL where Mfuse = 0

        #What, why did these values slightly change?? (originals now in #)
        para[iafduo, :, :] .= 0.019# 0.019    # fduo   fuselage velocity overspeed at wing root
        para[iafdus, :, :] .= 0.014#0.011    # fdus   fuselage velocity overspeed at wing break
        para[iafdut, :, :] .= 0.004#0.004   # fdut   fuselage velocity overspeed at wing tip

    # Landing gear weight fractions and locations
        parg[igxhpesys ] = 125.0 * ft_to_m   #  xhpesys   hyd/pneu/ele system location
        parg[igxlgnose ] =  28.0 * ft_to_m   #  xlgnose   nose LG location
        parg[igdxlgmain] =   3.0 * ft_to_m   # dxlgmain   main LG offset behind wing lift centroid

        parg[igfhpesys ] =  0.010     # fhpesys    Whpesys/WMTO
        parg[igflgnose ] =  0.010     # flgnose    Wlgnose/WMTO
        parg[igflgmain ] =  0.040     # flgmain    Wlgmain/WMTO

        parg[igxapu    ] = parg[igxconend]-3.0 * ft_to_m # xapu      APU location
        parg[igfapu    ] = 0.035          # fapu   Wapu/Wpay     APU weight fraction

        parg[igfreserve] = 0.05  # freserve Wfreserve/Wburn
        parg[igfpadd   ] = 0.35  # fpadd    Wpadd/Wpay    other payload-proportional fraction
        parg[igfseat   ] = 0.10  # fseat    Wseat/Wpay    seat weight fraction
        parg[igfeadd   ] = 0.10  # feadd    Weadd/Wbare   engine accessories, fuel system fraction 

        # if pari[iibleedeng] == 1
        #     parg[igfeadd   ] = 0.10 * 0.90
        #     parg[igfhpesys ] =  0.010 * 0.95
        # end


        parg[igfpylon  ] = 0.05  # fpylon   Wpylon/We+a+n engine pylon weight fraction   

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
        
        # Composite Properties
            # # Lamina Properties
            # parg[igE11] = 163000e6 # Young's modulus along fiber orientation
            # parg[igE22] = 7800e6# Young's moduls perpendicular to fiber orientation
            # # parg[igVf] = 0.6 # fraction of volume in lamina comprised of carbon fiber
            # parg[igTTfx] = 2297e6/(1.5*1.1) # ultimate strength along fiber orientation divide by factor of safety and factor accounting for stress concentrations
            # parg[igTTfy] = 69.96e6/(1.5*1.1) # ultimate strength perpendicular to fiber orientation divide by factor of safety and factor accounting for stress concentrations
            # parg[igCCfx] = 1600e6/(1.5*1.1) 
            # parg[igCCfy] = 160.6e6/(1.5*1.1) 
            # parg[igtaucomp] = 63.4e6/(1.5*1.1) # max in-plane shear stress 
            # parg[igGf] = 4000e6 # shear modulus of fiber
            # parg[igrhof] = 1600.0 # fiber density
            # parg[ignuf] = 0.26 # poisson's ratio for fiber

            # # Matrix Properties
            # # parg[igEm] = 3300e6 # Young's modulus for matrix
            # # parg[igGm] = 2200e6 # shear modulus of matrix
            # # parg[igrhom] = 1290.0 # matrix density
            # # parg[ignum] = 0.37 # poisson's ratio for matrix

            # parg[igcompfac] = 0.5 #factor used to characterize imperfections in composite manufacturing when calculating ultimate stress

            # parg[igWpplsp] = 0.2 * 9.81 #Lightning strile protection weight in N/m^2

    # Nacelle Drag stuff
        parg[igrSnace  ] = 12.0   # rSnace   nacelle+pylon wetted area/fan area  Snace/Afan
        parg[igrVnace  ] =  1.02  # rVnace   nacelle local/freesteam velocity ratio

    parg[igrWfmax  ] = 0.98

    

# Aerodynamic parameters

para[iaalt, ipcruise1, :] .=  32000.0 * ft_to_m # Cruise altitude [m] Max fuel max payload

# Takeoff and initial climb and descent
    para[iaalt,    ipstatic:ipcutback, :] .= parm[imaltTO]
    para[iaclpmax, ipstatic:ipcutback, :] .= 2.8 # clpmax   wing max cl_perp  = CLmax/cos(sweep)^2
    
    para[iaalt,    ipclimb1, :] .= parm[imaltTO]
    para[iaclpmax, ipclimb1, :] .= 2.8 # clpmax   wing max cl_perp  = CLmax/cos(sweep)^2
    
    para[iaalt,    ipdescentn, :] .= parm[imaltTO]
    para[iaclpmax, ipdescentn, :] .= 2.8 # clpmax   wing max cl_perp  = CLmax/cos(sweep)^2
    
# Cruise 
    para[iaCL  , ipclimb1+1:ipdescentn-1, :] .= 0.4702
    para[iaMach, ipclimbn:ipdescent1  , :] .= 0.84

# Wing span load parameters
    # Takeoff
        para[iarcls, 1:ipclimb1, :] .= 1.1     #  rcls    break/root cl ratio = cls/clo
        para[iarclt, 1:ipclimb1, :] .= 0.6     #  rclt    tip  /root cl ratio = clt/clo
        para[iacmpo, 1:ipclimb1, :] .= -0.30   #  cmpo    root  cm
        para[iacmps, 1:ipclimb1, :] .= -0.30   #  cmps    break cm
        para[iacmpt, 1:ipclimb1, :] .= -0.05   #  cmpt    tip   cm

    # Clean climb cruise descent and for wing structure sizing
        para[iarcls, ipclimb1+1 : ipdescentn-1, :] .= 1.132   #  rcls    
        para[iarclt, ipclimb1+1 : ipdescentn-1, :] .= 1.0266    #  rclt    
        para[iacmpo, ipclimb1+1 : ipdescentn-1, :] .= -0.10   #  cmpo    
        para[iacmps, ipclimb1+1 : ipdescentn-1, :] .= -0.10   #  cmps    
        para[iacmpt, ipclimb1+1 : ipdescentn-1, :] .= -0.10   #  cmpt   
   
    # Landing, forward CG tail sizing case
        para[iarcls, ipdescentn, :] .= 1.0     #  rcls  
        para[iarclt, ipdescentn, :] .= 0.5     #  rclt  
        para[iacmpo, ipdescentn, :] .= -0.40   #  cmpo      
        para[iacmps, ipdescentn, :] .= -0.40   #  cmps      
        para[iacmpt, ipdescentn, :] .= -0.05   #  cmpt     

# Wing and tail cd's
    para[iacdfw  , 1 : iptotal, :] .=  0.0085  #  cdfw    wing profile cd for low speed (takeoff, initial climb)
    para[iacdpw  , 1 : iptotal, :] .=  0.0035  #  cdpw    
    para[iaRerefw, 1 : iptotal, :] .=  20.0e6  #  Rerefw
                                        
    para[iacdft  , 1 : iptotal, :] .=  0.0060  #  cdft    tail profile cd
    para[iacdpt  , 1 : iptotal, :] .=  0.0035  #  cdpt    
    para[iaRereft, 1 : iptotal, :] .=  10.0e6  #  Rereft  
                                       
    para[iacdfs  , 1 : iptotal, :] .=  0.0085  #  cdfs    strut profile cd (not used if there's no strut)
    para[iacdps  , 1 : iptotal, :] .=  0.0035  #  cdps    
    para[iaRerefs, 1 : iptotal, :] .=  1.0e6   #  Rerefs  
                                       
    para[iaaRexp , 1 : iptotal, :] .=  -0.15   #  aRexp   exponent for Re-scaling:  CD = cd * (Re/Re_ref)^aRexp
                                       
    para[iafexcdw, 1 : iptotal, :] .=  1.02     #  fexcdw   # wing excrescence drag factor
    para[iafexcdt, 1 : iptotal, :] .=  1.02     #  fexcdt   # tail excrescence drag factor
    para[iafexcdf, 1 : iptotal, :] .=  1.030    #  fexcdf   # fuse excrescence drag factor

# Engine parameters

# pare[ieTt4, 1:iptotal, :] .= 3000.0 # [R]
pare[ieTt4, 1:iptotal, :] .= 2630.0 # [R]
# pare[ieTt4, ipstatic:iptakeoff, :] .= 3200.0 #[R]
pare[ieTt4, ipstatic:iptakeoff, :] .= 2900.0 #[R]


parg[igfTt4CL1] = 0.2
parg[igfTt4CLn] = 0.2

parg[igfanPCT] =   100.0

mofftpax = 0.008
mofftmMTO = 0.0085 * 0.001
Pofftpax  = 200.0
PofftmMTO = 0.5

# mofftpax = 0.008 / 2.0
# mofftmMTO = 0.0085 * 0.001 / 2.0
# Pofftpax  = 200.0 / 2.0
# PofftmMTO = 0.5 / 2.0


parg[igmofWpay] = mofftpax / (230.0 * lbf_to_N)  #This later gets multiplied by Wpay, resulting in mofftpax*pax
parg[igmofWMTO] = mofftmMTO  / 9.81

parg[igPofWpay] = Pofftpax / (230.0 * lbf_to_N)
parg[igPofWMTO] = PofftmMTO / 9.81

parpt = zeros(Union{Int64, Float64}, ipt_total)

# 777-300ER
parpt[ipt_pifan]   = 1.57
parpt[ipt_piLPC]   = 1.72
parpt[ipt_piHPC]   = 25.0
parpt[ipt_Tt41 ]   = 2630.0





