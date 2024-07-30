
#---- indices for integer variables and flags
      iifuel   =  1  # index specifying the fuel (see gasfun.f)    0 = JetA 1 = LH2      
      iifwcen  =  2  # 0 = no fuel in center box  1 = fuel in center box
      iiwplan  =  3  # 0 = bare wing 1 = cantilever with engine 2 = strut braced
      iiengloc =  4  # 1 = engines on wing 2 = engines on fuselage "Options.engine_location"  
      iiengwgt =  5  # 1 = basic tech eng. weight 2 = advanced tech eng. weight 
      iiBLIc   =  6  # 0 = core in clean flow 1 = core ingests KE defect   
      iifclose =  7  # 0 = fuse tapers to point 1 = tapers to edge 
      iixwmove = 8  # 0 = wing centroid fixed 1 = move wing to balance
      iifwing  = 9  # 0 = no fuel in wing (for LH2) 1 = fuel in wing
      iicalcCMVf = 10
      iiengtype  = 11  # 0 = Turboelectric engine; 1 = Turbofan engine; "prop_sys_arch"
      iiopt      = 12  # 0 = run sizing loop only; 1 = run optimization procedure "optimize"
      iiaircraftclass = 13 # 737 = 737 size aircraft; 777 = 777 size aircraft
      iinftanks = 14 # number of fuel tanks in fuselage
      iiengmodel = 15 #engine model 
      iidoubledeck = 16 #double decker flag
      iitotal    = 16

#---- indices for missions to be examined 
#-      (different Ranges Payloads ...)
      imwOpt     =  1 
      imRange    =  2
      imWpay     =  3 
      imWperpax  =  4
      imaltTO    =  5
      imT0TO     =  6
      imWTO      =  7
      imWfuel    =  8
      imPFEI     =  9
      imV1       = 10
      imV2       = 11
      imtTO      = 12
      imFTO      = 13
      iml1       = 14
      imlTO      = 15
      imlBF      = 16
      imlCB      = 17
      imgamVTO   = 18
      imgamVBF   = 19
      imgamVCB   = 20
      imgamVDE1  = 21
      imgamVDEn  = 22
      imthCB     = 23
      imxCB      = 24
      imzCB      = 25
      imxFO      = 26
      imzFO      = 27
      imdBSL     = 28
      imdBCB     = 29
      imdBFO     = 30
      imfexcdw   = 31
      imfexcdt   = 32
      imfexcdf   = 33
      imDeltaTatm= 34
      imLHVfuel  = 35
      imtotal    = 35

#---- indices for geometry (airframe) variables and other sizing variables
      const igFOpt     =   1
      const igPFEI     =   2
      const igRange    =   3
      const igWMTO     =   4
      const igWpay     =   5
### Good Riddance ###
      const igWfuel    =   6
      const igWfmax    =   7
      const igrWfmax   =   8
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igWebare   =  9
      const igWnace    =  10
      const igWeng     =  11
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igdxWfuel  =  12
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igxCGfwd   =  13
      const igxCGaft   =  14
      const igfreserve =  15
### Good Riddance ###
### Good Riddance ###
      const igfeadd    =  16
      const igfpylon   =  17
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igflgnose  =  18
      const igflgmain  =  19
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igNlift    =  20
      const igNland    =  21
      const igVne      =  22
      const igneng     =  23
      const igGearf    =  24
      const igfTt4CL1  =  25
      const igfTt4CLn  =  26
      const igHTRf     =  27
      const igHTRlc    =  28
      const igHTRhc    =  29
      const igrSnace   =  30
      const igrVnace   =  31
### Good Riddance ###
      const igfSnace   =  32
      const igpcabin   =  33
      const igdeltap   =  34
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igxeng     =  35     
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igxlgnose  = 36
      const igdxlgmain = 37
      const igyeng     = 38
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igcabVol   = 39
### Good Riddance ###
### Good Riddance ###
      const igrpayfwd  = 40
      const igrpayaft  = 41
      const igxNP      = 42
      const igCMVf1    = 43
      const igCLMf0    = 44
### Good Riddance ###
      const igdCLnda   = 45
      const igdCLhdCL  = 46
      const igdCLndCL  = 47
      const igCLhspec  = 48
### Good Riddance ###
      const igCLveout  = 49
### Good Riddance ###
### Good Riddance ###
      const igfCDhcen  = 50
### Good Riddance ###
      const igrMh      = 51
      const igrMv      = 52
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igsigfac   = 53
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igrhofuel  = 54
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igfLo      = 55
      const igfLt      = 56
      const igfLn      = 57
      const igcma      = 58
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igfBLIw    = 59
      const igfBLIf    = 60
      const igdfan     = 61
      const igdlcomp   = 62
      const igdhcomp   = 63
      const iglnace    = 64
      const igA5       = 65
      const igA7       = 66
      const igTmetal   = 67
      const igcdefan   = 68
      const igCDgear   = 69
      const igCDspoil  = 70
      const igmuroll   = 71
      const igmubrake  = 72
      const ighobst    = 73
      const iglBFmax   = 74
### Good Riddance ###
      const iggtocmin  = 75
      const igdBSLmax  = 76
      const igdBCBmax  = 77
      const igmofWpay  = 78
      const igmofWMTO  = 79
      const igPofWpay  = 80
      const igPofWMTO  = 81

      # Adding indices to store TE weights and locations
      const igWtshaft  = 82
      const igWgen     = 83
      const igWinv     = 84
      const igWmot     = 85
      const igWfan     = 86
      const igWftank   = 87

      const igxtshaft  = 88
      const igxgen     = 89
      const igxinv     = 90
      const igxmot     = 91
      const igxfan     = 92
      const igxftank   = 93
      
      const igxcables  = 94
      const igWcables  = 95

      const igxcat     = 96
      const igWcat     = 97

      const igWtesys   = 98
      const igxWtesys  = 99

      const iglftank   = 100
      const igWinsftank= 101
      const igxWftank  = 102
      
      const igRftank   = 103
      const igWc3des   = 104

      const igdaftfan  = 105
      const iglnaceaft = 106
### Good Riddance ###
      const igneout    = 107
      const igyeout    = 108
      const igyeinn    = 109
      const iglftankin = 110 #lftank input
      const igLHVfuel  = 111
      const igWfburn   = 112

      const igWaftfan  = 113
      const igWfanGB   = 114
      const igWaftfanGB   = 115
      const igWrect    = 116
      const igWtms     = 117

      const igfanPCT   = 118

      const iggamVtarget = 119

      const igWpaymax    = 120

      const ig_YC = 121
      const ig_YH = 122
      const ig_YO = 123
      const ig_YN = 124
      const ig_MW = 125
      const ighfuel      = 126
      const igTfuel      = 127

      const igxftankaft  = 128
      const igdxcabin    = 129
      const igdxeng2wbox = 130

      const igxWfuel = 131

      const igseatpitch = 132
      const igseatwidth = 133
      const igaislehalfwidth = 134
      const igfloordist = 135

      const igWfvent = 136

      const igtotal    = 136

# indices for turbo-electric systems - really just the electrical machines
      const ite_ratSM   = 1
      const ite_ratAg   = 2
      const ite_ratM    = 3
      const ite_ratSp   = 4
      const ite_ratW    = 5
      const ite_ratShft = 6
      const ite_z       = 7
      const ite_kServ   = 8
      const ite_ratSd   = 9
      const ite_kpf     = 10
      const ite_wSd     = 11
      const ite_Nshrt   = 12
      const ite_hRS     = 13
      const ite_Tarm    = 14
      const ite_kst     = 15
      const ite_p        = 16
      const ite_Br       = 17
      const ite_BSat     = 18
      const ite_mu0      = 19
      const ite_rhoMag   = 20
      const ite_rhoCu    = 21
      const ite_rhoFe    = 22
      const ite_rhoSteel = 23
      const ite_tauMax   = 24
      const ite_sigCu    = 25
      const ite_thetaCu  = 26
      const ite_psi      = 27
      const ite_pb0      = 28
      const ite_Bb0      = 29
      const ite_fb0      = 30
      const ite_epsb     = 31   
      const ite_epsf     = 32 
      const ite_lambda   = 33
      const ite_kw       = 34
      const ite_ks       = 35
      const ite_BAg      = 36
      const ite_rRoti    = 37
      const ite_hM       = 38
      const ite_lRot     = 39
      const ite_hAg      = 40
      const ite_areaArm  = 41
      const ite_Bsbi     = 42
      const ite_Bt       = 43
      const ite_mSBI     = 44
      const ite_mTeeth   = 45
      const ite_lArm     = 46
      const ite_NSz      = 47
      const ite_Wpmsm    = 48
      const ite_k1       = 49  # k1, k2, k3 are inverter constants for efficiency calcs
      const ite_k2       = 50
      const ite_k3       = 51
      const ite_Pinvdes  = 52  # Stores design inverter power
      const ite_SPdes    = 53  # Specific power at design
      const ite_effdes   = 54  

      const ite_total    = 54


# indices for powertrain
      const ipt_nfan    = 1
      const ipt_ngen    = 2
      const ipt_nTshaft = 3

      const ipt_pifan   = 4
      const ipt_piLPC   = 5
      const ipt_piHPC   = 6

      const ipt_ARgen       = 7
      const ipt_sigAgGen    = 8
      const ipt_ratSplitGen = 9

      const ipt_ARmot       = 10
      const ipt_sigAgMot    = 11
      const ipt_ratSplitMot = 12
      
      const ipt_Tt41     = 13

      const ipt_cpsi    = 14
      const ipt_wcat    = 15
      const ipt_lcat    = 16
      const ipt_deNOx   = 17

      const ipt_Wfan     = 18
      const ipt_Wmot     = 19
      const ipt_Winv     = 20
      const ipt_Wcables  = 21
      const ipt_Wgen     = 22
      const ipt_Wtshaft  = 23
      const ipt_Wnacelle = 24
      const ipt_Wcatalyst= 25
      const ipt_Wpttotal = 26
      
      const ipt_NdesGen  = 27
      const ipt_NdesFan  = 28
      const ipt_NdesMot  = 29
      
      const ipt_fb0      = 30
      const ipt_epsb     = 31   
      const ipt_epsf     = 32 
      const ipt_lambda   = 33
      const ipt_kw       = 34
      const ipt_ks       = 35
      const ipt_BAg      = 36
      const ipt_rRoti    = 37
      const ipt_hM       = 38
      const ipt_lRot     = 39
      const ipt_hAg      = 40
      const ipt_areaArm  = 41
      const ipt_Bsbi     = 42
      const ipt_Bt       = 43
      const ipt_mSBI     = 44
      const ipt_mTeeth   = 45
      const ipt_lArm     = 46
      const ipt_NSz      = 47
      const ipt_Wpmsm    = 48
      const ipt_k1       = 49  # k1, k2, k3 are inverter constants for efficiency calcs
      const ipt_k2       = 50
      const ipt_k3       = 51
      const ipt_Pinvdes  = 52  # Stores design inverter power

      const ipt_time_NPSS  = 53
      const ipt_calls_NPSS = 54
      const ipt_FanGR      = 55

      const ipt_Ptshaft   = 56
      const ipt_Fnsplit   = 57 # Thrust split between fuse mounted fans and wing podded fans

      const ipt_Rcable    = 58
      const ipt_Vcable    = 59
      const ipt_lcable    = 60
      const ipt_sigcon    = 61 # Conductance of conductor
      const ipt_alphacon  = 62 # Temp coeff
      const ipt_rholcable = 63 # linear density kg/m
      const ipt_rhocon    = 64
      const ipt_Jmax      = 65
      const ipt_rhoins    = 66
      const ipt_Emax      = 67
      const ipt_kpf       = 68
      
      const ipt_Pelec_mot = 69
      const ipt_Pelec_gen = 70


      const ipt_total     = 70

#---- indices for aero variables at each mission point
      iaalt     =  1
      iatime    =  2
      iaRange   =  3
      iafracW   =  4
      iaWbuoy   =  5
      iaMach    =  6
      iaReunit  =  7
      iagamV    =  9
      iaCL      = 10
      iaCD      = 11
      iaCDi     = 12
      iaDAfsurf = 13
      iaDAfwake = 14
      iaKAfTE   = 15
      iaPAfinf  = 16
      iaCDfuse  = 17
      iaCDwing  = 18
      iaCDhtail = 19
      iaCDvtail = 20
      iaCDnace  = 21
      iaCDstrut = 22
      iaCDover  = 23
      iaCfnace  = 24
      iarcls    = 25
      iarclt    = 26
      iaclpo    = 27
      iaclps    = 28
      iaclpt    = 29
      iacmpo    = 30
      iacmps    = 31
      iacmpt    = 32
      iafduo    = 33
      iafdus    = 34
      iafdut    = 35
      iaclpmax  = 36
      iaCLh     = 37
      iaCMw0    = 38
      iaCMw1    = 39
      iaCMh0    = 40
      iaCMh1    = 41
      iaCMwing  = 42
      iaCMfuse  = 43
      iaCMtail  = 44
      iaxCG     = 45
      iaxCP     = 46
      iaxNP     = 47
      iafexcdw  = 48
      iafexcdt  = 49
      iafexcdf  = 50
      iacdfw    = 51
      iacdpw    = 52
      iaspaneff = 53
      iaRerefw  = 54
      iaaRexp   = 55
      iacdft    = 56
      iacdpt    = 57
      iaRereft  = 58
      iacdfs    = 59
      iacdps    = 60
      iaRerefs  = 61
      iadCDBLIf = 62
      iadCDBLIw = 63
      iaROC     = 64
      iatotal   = 64 

#---- indices for engine variables at each mission point
      iehfuel =   1
      ieTfuel =   2
      ieff    =   3
      iepid   =   4
      iepib   =   5
      iepifn  =   6
      iepitn  =   7
      ieBPR   =   8
      ieepolf =   9
      ieepollc=  10
      ieepolhc=  11
      ieepolht=  12
      ieepollt=  13
      ieetab  =  14
      iepifK  =  15
      ieepfK  =  16
      ieNf    =  17
      ieN1    =  18
      ieN2    =  19
      ieNbf   =  20
      ieNblc  =  21
      ieNbhc  =  22
      iembf   =  23
      iemblc  =  24
      iembhc  =  25
      iepif   =  26
      iepilc  =  27
      iepihc  =  28
      ieNbfD  =  29
      ieNblcD =  30
      ieNbhcD =  31
      ieNbhtD =  32
      ieNbltD =  33
      iembfD  =  34
      iemblcD =  35
      iembhcD =  36
      iembhtD =  37
      iembltD =  38
      iepifD  =  39
      iepilcD =  40
      iepihcD =  41
      iepihtD =  42
      iepiltD =  43
      ieM2    =  44
      ieM25   =  45
      ieM0    =  46
      iep0    =  47
      iea0    =  48
      ierho0  =  49
      iemu0   =  50
      ieT0    =  51
      ieu0    =  52
      ieTt0   =  53
      ieht0   =  54
      iept0   =  55
      iecpt0  =  56
      ieRt0   =  57
      ieTt18  =  58
      ieht18  =  59
      iept18  =  60
      iecpt18 =  61
      ieRt18  =  62
      ieTt19  =  63
      ieht19  =  64
      iept19  =  65
      iecpt19 =  66
      ieRt19  =  67
      ieTt2   =  68
      ieht2   =  69
      iept2   =  70
      iecpt2  =  71
      ieRt2   =  72
      ieTt21  =  73
      ieht21  =  74
      iept21  =  75
      iecpt21 =  76
      ieRt21  =  77
      ieTt25  =  78
      ieht25  =  79 
      iept25  =  80 
      iecpt25 =  81 
      ieRt25  =  82 
      ieTt3   =  83
      ieht3   =  84
      iept3   =  85
      iecpt3  =  86
      ieRt3   =  87
      ieTt4   =  88
      ieht4   =  89
      iept4   =  90
      iecpt4  =  91
      ieRt4   =  92
      ieTt41  =  93
      ieht41  =  94
      iept41  =  95
      iecpt41 =  96
      ieRt41  =  97
      ieTt45  =  98
      ieht45  =  99
      iept45  = 100
      iecpt45 = 101
      ieRt45  = 102
      ieTt49  = 103
      ieht49  = 104
      iept49  = 105
      iecpt49 = 106
      ieRt49  = 107
      ieTt5   = 108
      ieht5   = 109
      iept5   = 110
      iecpt5  = 111
      ieRt5   = 112
      ieTt7   = 113
      ieht7   = 114
      iept7   = 115
      iecpt7  = 116
      ieRt7   = 117
      ieTt9   = 118
      iept9   = 119
      iep2    = 120
      ieT2    = 121
      ieR2    = 122
      iecp2   = 123
      ieu2    = 124
      ieA2    = 125
      iep25   = 126
      ieT25   = 127
      ieR25   = 128
      iecp25  = 129
      ieu25   = 130
      ieA25   = 131
      iep5    = 132
      ieT5    = 133
      ieR5    = 134
      iecp5   = 135
      ieu5    = 136
      ieA5    = 137
      iep6    = 138
      ieT6    = 139
      ieR6    = 140
      iecp6   = 141
      ieu6    = 142
      ieA6    = 143
      iep7    = 144
      ieT7    = 145
      ieR7    = 146
      iecp7   = 147
      ieu7    = 148
      ieA7    = 149
      iep8    = 150
      ieT8    = 151
      ieR8    = 152
      iecp8   = 153
      ieu8    = 154
      ieA8    = 155
      ieu9    = 156
      ieA9    = 157
      ieepf   = 158
      ieeplc  = 159
      ieephc  = 160
      ieepht  = 161
      ieeplt  = 162
      ieetaf  = 163
      ieetalc = 164
      ieetahc = 165
      ieetaht = 166
      ieetalt = 167
      iemcore = 168
      iemofft = 169
      iePofft = 170
      iePhiinl= 171
      ieKinl  = 172
      ieepsl  = 173
      ieepsh  = 174
      ieFe    = 175 
      ieFsp   = 176
      ieTSFC  = 177
      ieA5fac = 178
      ieA7fac = 179
      iedTstrk= 180
      ieStA   = 181
      ieMtexit= 182
      ieM4a   = 183
      ieruc   = 184
      ieefilm = 185
      ietfilm = 186
      iefc    = 187
      ieepsc1 = 188
      ieepsc2 = 189
      ieepsc3 = 190
      ieepsc4 = 191
      ieTmet1 = 192
      ieTmet2 = 193
      ieTmet3 = 194
      ieTmet4 = 195

      iedeNOx = 196
      iemdotf = 197

      iePLH2  = 198
      ieyg    = 199

      ieemot     = 200
      ieeinv     = 201
      ieecable   = 202
      ieegen     = 203
      ieethermal = 204

      iePinmot     = 205
      iePininv     = 206
      iePincable   = 207
      iePingen     = 208
      iePinthermal = 209

      ieEINOx1     = 210
      ieEINOx2     = 211
      ieFAR        = 212
      ieOPR        = 213
      ieWc3        = 214
      
      # Heat rejection:
      ieHrejmot   = 215
      ieHrejinv   = 216
      ieHrejcab   = 217
      ieHrejgen   = 218
      ieHrejtot   = 219

      ieHexcess   = 220
#Guesses for NPSS offdes:
      iegsFnsplit = 221
      iegsMotShP  = 222
      iegsmdotf   = 223
      iegsWin     = 224
      iegsRlineF  = 225
      iegsBPR     = 226
      iegsRlinelc = 227
      iegsRlinehc = 228
      iegsPRhtrb  = 229
      iegsPRltrb  = 230
      iegsNmechH  = 231
      iegsGBtrq   = 232
      iegsNmechL  = 233
      iegsNmechF  = 234
      iegsPodWin  = 235

      iegsPodRlineF   = 236
      iegsPodGBtrq    = 237
      iegsPodMotNmech = 238
      iegsPodFanNmech = 239

#Heat exchanger variables
      ieDi = 240
      ieTft = 241
      iefrecirc = 242
      ierecircT = 243
      iehlat = 244
      iePreCorder = 245
      iePreCepsilon = 246
      iePreCMp = 247
      iePreCDeltah = 248
      iePreCDeltap = 249
      ieInterCorder = 250
      ieInterCepsilon = 251 
      ieInterCMp = 252
      ieInterCDeltah = 253
      ieInterCDeltap = 254
      ieRegenorder = 255
      ieRegenepsilon = 256 
      ieRegenMp = 257
      ieRegenDeltah = 258
      ieRegenDeltap = 259
      ieTurbCorder = 260
      ieTurbCepsilon = 261 
      ieTurbCMp = 262
      ieTurbCDeltah = 263
      ieTurbCDeltap = 264

      ietotal = 264

#---- max number of blade rows allowed by ieepsc* and ieTmet* indices above
# integer ncrowx
ncrowx = ieTmet1 - ieepsc1

#---- indices for indexing mission points
      ipstatic   =  1
      iprotate   =  2
      iptakeoff  =  3
      ipcutback  =  4
      ipclimb1   =  5
      ipclimb2   =  6
      ipclimb3   =  7
      ipclimb4   =  8
      ipclimb5   =  9
      ipcruise1  = 10
      ipcruise2  = 11
      ipdescent1 = 12
      ipdescent2 = 13
      ipdescent3 = 14
      ipdescent4 = 15
      ipdescent5 = 16
      iptest     = 17
      iptotal    = 17 
      
      ipclimbn   = ipclimb5      # last climb   point
      ipcruisen  = ipcruise2     # last cruise  point
      ipdescentn = ipdescent5    # last descent point

#---- indices for indexing optimization parameters (to be sequenced over)
      isRange  =  1
      isMach   =  2
      isNmax   =  3
      issigfac =  4
      isCL     =  5
      isAR     =  6
      issweep  =  7
      isetas   =  8
      isT4CR   =  9
      isT4TO   = 10
      isTmetal = 11
      isOPR    = 12
      isFPR    = 13
      islBFmax = 14
      isbmax   = 15
      isalt    = 16
      istotal  = 16  

#---- indices for indexing optimization variables      
      ioCL    =  1
      ioAR    =  2
      iosweep =  3
      iohboxo =  4
      iohboxs =  5
      iolams  =  6
      iolamt  =  7
      iorcls  =  8
      iorclt  =  9
      ioFPR   = 10
      ioBPR   = 11
      ioalt   = 12
      ioT4CR  = 13
      ioT4TO  = 14
      ioOPR   = 15
      iototal = 15

#=
      character*8 cpars(istotal)
      character*8 cparo(iototal)
      character*2 cplab(iptotal)

#---- keywords for ij sequence parameters 
#-     (must not have any leading blanks must have at least one trailing blank)
      data cpars /
      'Range   '
      'Mach    '
      'Nmax    '
      'sigfac  '
      'CL      '
      'AR      '
      'sweep   '
      'etas    '
      'Tt4CR   '
      'Tt4TO   '
      'Tmetal  '
      'OPR     '
      'FPR     '
      'lBFmax  '
      'bmax    '
      'alt     ' /

#---- keywords for optimization variables
#-     (must not have any leading blanks must have at least one trailing blank)
      data cparo /
      'CL      '
      'AR      '
      'sweep   '
      'hboxo   '
      'hboxs   '
      'lambdas '
      'lambdat '
      'rcls    '
      'rclt    '
      'FPR     '
      'BPR     '
      'alt     '
      'Tt4CR   '
      'Tt4TO   '
      'OPR     ' /

#---- printout labels for mission points
      data cplab /
      'ST'
      'RO'
      'TO'
      'CB'
      'B1'
      'B2'
      'B3'
      'B4'
      'B5'
      'C1'
      'C2'
      'D1'
      'D2'
      'D3'
      'D4'
      'D5'
      'TE'  /
=#
