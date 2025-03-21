
#---- indices for integer variables and flags
      const iifuel = 1  # index specifying the fuel (see gasfun.f)    0 = JetA 1 = LH2      
      const iifwcen = 2  # 0 = no fuel in center box  1 = fuel in center box
      const iiwplan = 3  # 0 = bare wing 1 = cantilever with engine 2 = strut braced
      const iiengloc = 4  # 1 = engines on wing 2 = engines on fuselage "Options.engine_location"  
      const iiengwgt = 5  # 1 = basic tech eng. weight 2 = advanced tech eng. weight 
      const iiBLIc = 6  # 0 = core in clean flow 1 = core ingests KE defect   
      const iifclose = 7  # 0 = fuse tapers to point 1 = tapers to edge 
      const iixwmove = 8  # 0 = wing centroid fixed 1 = move wing to balance
      const iifwing = 9  # 0 = no fuel in wing (for LH2) 1 = fuel in wing
      const iicalcCMVf = 10
      const iiengtype = 11  # 0 = Turboelectric engine; 1 = Turbofan engine; "prop_sys_arch"
      const iiopt = 12  # 0 = run sizing loop only; 1 = run optimization procedure "optimize"
      const iiaircraftclass = 13 # 737 = 737 size aircraft; 777 = 777 size aircraft
      const iinftanks = 14 # number of fuel tanks in fuselage
      const iitotal = 14

#---- indices for missions to be examined 
#-      (different Ranges Payloads ...)
      const imwOpt = 1 
      const imRange = 2
      const imWpay = 3 
      const imWperpax = 4
      const imaltTO = 5
      const imT0TO = 6
      const imWTO = 7
      const imWfuel = 8
      const imPFEI = 9
      const imV1 = 10
      const imV2 = 11
      const imtTO = 12
      const imFTO = 13
      const iml1 = 14
      const imlTO = 15
      const imlBF = 16
      const imlCB = 17
      const imgamVTO = 18
      const imgamVBF = 19
      const imgamVCB = 20
      const imgamVDE1 = 21
      const imgamVDEn = 22
      const imthCB = 23
      const imxCB = 24
      const imzCB = 25
      const imxFO = 26
      const imzFO = 27
      const imdBSL = 28
      const imdBCB = 29
      const imdBFO = 30
      const imfexcdw = 31
      const imfexcdt = 32
      const imfexcdf = 33
      const imDeltaTatm = 34
      const imWfvent = 35
      const imtotal = 35

#---- indices for geometry (airframe) variables and other sizing variables
### Good Riddance ###
### Good Riddance ###
      const igRange    =   1
      const igWMTO     =   2
      const igWpay     =   3
      const igWfuel    =   4
      const igWfmax    =   5
      const igrWfmax   =   6
      const igWebare   =  7
      const igWnace    =  8
      const igWeng     =  9
      const igdxWfuel  =  10
      const igxCGfwd   =  11
      const igxCGaft   =  12
      const igfreserve =  13
      const igfeadd    =  14
      const igfpylon   =  15
### Good Riddance ###
### Good Riddance ###
      const igNlift    =  16
      const igNland    =  17
      const igVne      =  18
      const igneng     =  19
      const igGearf    =  20
      const igfTt4CL1  =  21
      const igfTt4CLn  =  22
      const igHTRf     =  23
      const igHTRlc    =  24
      const igHTRhc    =  25
      const igrSnace   =  26
      const igrVnace   =  27
      const igfSnace   =  28
      const igpcabin   =  29
      const igdeltap   =  30
      const igxeng     =  31     
### Good Riddance ###
### Good Riddance ###
      const igyeng     = 32
      const igcabVol   = 33
      const igrpayfwd  = 34
      const igrpayaft  = 35
      const igxNP      = 36
      const igCMVf1    = 37
      const igCLMf0    = 38
      const igdCLnda   = 39
      const igdCLhdCL  = 40
      const igdCLndCL  = 41
      const igCLhspec  = 42
      const igCLveout  = 43
      const igfCDhcen  = 44
      const igrMh      = 45
      const igrMv      = 46
      const igsigfac   = 47
      const igrhofuel  = 48
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igfBLIw    = 49
      const igfBLIf    = 50
      const igdfan     = 51
      const igdlcomp   = 52
      const igdhcomp   = 53
      const iglnace    = 54
      const igA5       = 55
      const igA7       = 56
      const igTmetal   = 57
      const igcdefan   = 58
      const igCDgear   = 59
      const igCDspoil  = 60
      const igmuroll   = 61
      const igmubrake  = 62
      const ighobst    = 63
      const iglBFmax   = 64
      const iggtocmin  = 65
      const igdBSLmax  = 66
      const igdBCBmax  = 67
      const igmofWpay  = 68
      const igmofWMTO  = 69
      const igPofWpay  = 70
      const igPofWMTO  = 71

      # Adding indices to store TE weights and locations
      const igWtshaft  = 72
      const igWgen     = 73
      const igWinv     = 74
      const igWmot     = 75
      const igWfan     = 76
      const igWftank   = 77

      const igxtshaft  = 78
      const igxgen     = 79
      const igxinv     = 80
      const igxmot     = 81
      const igxfan     = 82
      const igxftank   = 83
      
      const igxcables  = 84
      const igWcables  = 85

      const igxcat     = 86
      const igWcat     = 87

      const igWtesys   = 88
      const igxWtesys  = 89

      const iglftank   = 90
      const igWinsftank= 91
      const igxWftank  = 92
      
      const igRftank   = 93
      const igWc3des   = 94

      const igdaftfan  = 95
      const iglnaceaft = 96

      const igneout    = 97
      const igyeout    = 98
      const igyeinn    = 99
      const iglftankin = 100 #lftank input
      const igLHVfuel  = 101
      const igWfburn   = 102

      const igWaftfan  = 103
      const igWfanGB   = 104
      const igWaftfanGB   = 105
      const igWrect    = 106
      const igWtms     = 107

      const igfanPCT   = 108

      const iggamVtarget = 109

      const igWpaymax    = 110

      const ig_YC = 111
      const ig_YH = 112
      const ig_YO = 113
      const ig_YN = 114
      const ig_MW = 115
      const ighfuel      = 116
      const igTfuel      = 117

      const igxftankaft  = 118
      const igdxeng2wbox = 119

      const igxWfuel = 120

### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###

      const igWfvent = 121
      const igdHXPreC = 122
      const igdHXInterC = 123
      const igdHXRegen = 124
      const igHXaddmassfrac = 125
      const igWHXs = 126
      const igtotal    = 126

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
      const iaalt     =  1
      const iatime    =  2
      const iaRange   =  3
      const iafracW   =  4
      const iaWbuoy   =  5
      const iaMach    =  6
      const iaReunit  =  7
      const iagamV    =  9
      const iaCL      = 10
      const iaCD      = 11
      const iaCDi     = 12
      const iaDAfsurf = 13
      const iaDAfwake = 14
      const iaKAfTE   = 15
      const iaPAfinf  = 16
      const iaCDfuse  = 17
      const iaCDwing  = 18
      const iaCDhtail = 19
      const iaCDvtail = 20
      const iaCDnace  = 21
      const iaCDstrut = 22
      const iaCDover  = 23
      const iaCfnace  = 24
      const iarcls    = 25
      const iarclt    = 26
      const iaclpo    = 27
      const iaclps    = 28
      const iaclpt    = 29
      const iacmpo    = 30
      const iacmps    = 31
      const iacmpt    = 32
      const iafduo    = 33
      const iafdus    = 34
      const iafdut    = 35
      const iaclpmax  = 36
      const iaCLh     = 37
      const iaCMw0    = 38
      const iaCMw1    = 39
      const iaCMh0    = 40
      const iaCMh1    = 41
      const iaCMwing  = 42
      const iaCMfuse  = 43
      const iaCMtail  = 44
      const iaxCG     = 45
      const iaxCP     = 46
      const iaxNP     = 47
      const iafexcdw  = 48
      const iafexcdt  = 49
      const iafexcdf  = 50
      const iacdfw    = 51
      const iacdpw    = 52
      const iaspaneff = 53
      const iaRerefw  = 54
      const iaaRexp   = 55
      const iacdft    = 56
      const iacdpt    = 57
      const iaRereft  = 58
      const iacdfs    = 59
      const iacdps    = 60
      const iaRerefs  = 61
      const iadCDBLIf = 62
      const iadCDBLIw = 63
      const iaROC     = 64
      const iatotal   = 64 

#---- indices for engine variables at each mission point
      const iehfuel = 1
      const ieTfuel = 2
      const ieff = 3
      const iepid = 4
      const iepib = 5
      const iepifn = 6
      const iepitn = 7
      const ieBPR = 8
      const ieepolf = 9
      const ieepollc = 10
      const ieepolhc = 11
      const ieepolht = 12
      const ieepollt = 13
      const ieetab = 14
      const iepifK = 15
      const ieepfK = 16
      const ieNf = 17
      const ieN1 = 18
      const ieN2 = 19
      const ieNbf = 20
      const ieNblc = 21
      const ieNbhc = 22
      const iembf = 23
      const iemblc = 24
      const iembhc = 25
      const iepif = 26
      const iepilc = 27
      const iepihc = 28
      const ieNbfD = 29
      const ieNblcD = 30
      const ieNbhcD = 31
      const ieNbhtD = 32
      const ieNbltD = 33
      const iembfD = 34
      const iemblcD = 35
      const iembhcD = 36
      const iembhtD = 37
      const iembltD = 38
      const iepifD = 39
      const iepilcD = 40
      const iepihcD = 41
      const iepihtD = 42
      const iepiltD = 43
      const ieM2 = 44
      const ieM25 = 45
      const ieM0 = 46
      const iep0 = 47
      const iea0 = 48
      const ierho0 = 49
      const iemu0 = 50
      const ieT0 = 51
      const ieu0 = 52
      const ieTt0 = 53
      const ieht0 = 54
      const iept0 = 55
      const iecpt0 = 56
      const ieRt0 = 57
      const ieTt18 = 58
      const ieht18 = 59
      const iept18 = 60
      const iecpt18 = 61
      const ieRt18 = 62
      const ieTt19 = 63
      const ieht19 = 64
      const iept19 = 65
      const iecpt19 = 66
      const ieRt19 = 67
      const ieTt2 = 68
      const ieht2 = 69
      const iept2 = 70
      const iecpt2 = 71
      const ieRt2 = 72
      const ieTt21 = 73
      const ieht21 = 74
      const iept21 = 75
      const iecpt21 = 76
      const ieRt21 = 77
      const ieTt25 = 78
      const ieht25 = 79 
      const iept25 = 80 
      const iecpt25 = 81 
      const ieRt25 = 82 
      const ieTt3 = 83
      const ieht3 = 84
      const iept3 = 85
      const iecpt3 = 86
      const ieRt3 = 87
      const ieTt4 = 88
      const ieht4 = 89
      const iept4 = 90
      const iecpt4 = 91
      const ieRt4 = 92
      const ieTt41 = 93
      const ieht41 = 94
      const iept41 = 95
      const iecpt41 = 96
      const ieRt41 = 97
      const ieTt45 = 98
      const ieht45 = 99
      const iept45 = 100
      const iecpt45 = 101
      const ieRt45 = 102
      const ieTt49 = 103
      const ieht49 = 104
      const iept49 = 105
      const iecpt49 = 106
      const ieRt49 = 107
      const ieTt5 = 108
      const ieht5 = 109
      const iept5 = 110
      const iecpt5 = 111
      const ieRt5 = 112
      const ieTt7 = 113
      const ieht7 = 114
      const iept7 = 115
      const iecpt7 = 116
      const ieRt7 = 117
      const ieTt9 = 118
      const iept9 = 119
      const iep2 = 120
      const ieT2 = 121
      const ieR2 = 122
      const iecp2 = 123
      const ieu2 = 124
      const ieA2 = 125
      const iep25 = 126
      const ieT25 = 127
      const ieR25 = 128
      const iecp25 = 129
      const ieu25 = 130
      const ieA25 = 131
      const iep5 = 132
      const ieT5 = 133
      const ieR5 = 134
      const iecp5 = 135
      const ieu5 = 136
      const ieA5 = 137
      const iep6 = 138
      const ieT6 = 139
      const ieR6 = 140
      const iecp6 = 141
      const ieu6 = 142
      const ieA6 = 143
      const iep7 = 144
      const ieT7 = 145
      const ieR7 = 146
      const iecp7 = 147
      const ieu7 = 148
      const ieA7 = 149
      const iep8 = 150
      const ieT8 = 151
      const ieR8 = 152
      const iecp8 = 153
      const ieu8 = 154
      const ieA8 = 155
      const ieu9 = 156
      const ieA9 = 157
      const ieepf = 158
      const ieeplc = 159
      const ieephc = 160
      const ieepht = 161
      const ieeplt = 162
      const ieetaf = 163
      const ieetalc = 164
      const ieetahc = 165
      const ieetaht = 166
      const ieetalt = 167
      const iemcore = 168
      const iemofft = 169
      const iePofft = 170
      const iePhiinl = 171
      const ieKinl = 172
      const ieepsl = 173
      const ieepsh = 174
      const ieFe = 175 
      const ieFsp = 176
      const ieTSFC = 177
      const ieA5fac = 178
      const ieA7fac = 179
      const iedTstrk = 180
      const ieStA = 181
      const ieMtexit = 182
      const ieM4a = 183
      const ieruc = 184
      const ieefilm = 185
      const ietfilm = 186
      const iefc = 187
      const ieepsc1 = 188
      const ieepsc2 = 189
      const ieepsc3 = 190
      const ieepsc4 = 191
      const ieTmet1 = 192
      const ieTmet2 = 193
      const ieTmet3 = 194
      const ieTmet4 = 195

      const iedeNOx = 196
      const iemdotf = 197

      const iePLH2 = 198
      const ieyg = 199

      const ieemot = 200
      const ieeinv = 201
      const ieecable = 202
      const ieegen = 203
      const ieethermal = 204

      const iePinmot = 205
      const iePininv = 206
      const iePincable = 207
      const iePingen = 208
      const iePinthermal = 209

      const ieEINOx1 = 210
      const ieEINOx2 = 211
      const ieFAR = 212
      const ieOPR = 213
      const ieWc3 = 214
      
      # Heat rejection:
      const ieHrejmot = 215
      const ieHrejinv = 216
      const ieHrejcab = 217
      const ieHrejgen = 218
      const ieHrejtot = 219

      const ieHexcess = 220
#Guesses for NPSS offdes:
      const iegsFnsplit = 221
      const iegsMotShP = 222
      const iegsmdotf = 223
      const iegsWin = 224
      const iegsRlineF = 225
      const iegsBPR = 226
      const iegsRlinelc = 227
      const iegsRlinehc = 228
      const iegsPRhtrb = 229
      const iegsPRltrb = 230
      const iegsNmechH = 231
      const iegsGBtrq = 232
      const iegsNmechL = 233
      const iegsNmechF = 234
      const iegsPodWin = 235

      const iegsPodRlineF = 236
      const iegsPodGBtrq = 237
      const iegsPodMotNmech = 238
      const iegsPodFanNmech = 239

#Heat exchanger variables
      const ieDi = 240
      const ieTft = 241
      const iefrecirc = 242
      const ierecircT = 243
      const iehvap = 244
      const iePreCorder = 245
      const iePreCepsilon = 246
      const iePreCMp = 247
      const iePreCDeltah = 248
      const iePreCDeltap = 249
      const ieInterCorder = 250
      const ieInterCepsilon = 251 
      const ieInterCMp = 252
      const ieInterCDeltah = 253
      const ieInterCDeltap = 254
      const ieRegenorder = 255
      const ieRegenepsilon = 256 
      const ieRegenMp = 257
      const ieRegenDeltah = 258
      const ieRegenDeltap = 259
      const ieTurbCorder = 260
      const ieTurbCepsilon = 261 
      const ieTurbCMp = 262
      const ieTurbCDeltah = 263
      const ieTurbCDeltap = 264
      const ieRadiatorepsilon = 265
      const ieRadiatorMp = 266
      const ieRadiatorDeltah = 267
      const ieRadiatorDeltap = 268
      const iehvapcombustor = 269
      const ieHXminTwall = 270
      const iemfuel = 271

      const ietotal = 271

#---- max number of blade rows allowed by ieepsc* and ieTmet* indices above
# integer ncrowx
ncrowx = ieTmet1 - ieepsc1

#---- indices for indexing mission points
      const ipstatic = 1
      const iprotate = 2
      const iptakeoff = 3
      const ipcutback = 4
      const ipclimb1 = 5
      const ipclimb2 = 6
      const ipclimb3 = 7
      const ipclimb4 = 8
      const ipclimb5 = 9
      const ipcruise1 = 10
      const ipcruise2 = 11
      const ipdescent1 = 12
      const ipdescent2 = 13
      const ipdescent3 = 14
      const ipdescent4 = 15
      const ipdescent5 = 16
      const iptest = 17
      const iptotal = 17 
      
      const ipclimbn   = ipclimb5      # last climb   point
      const ipcruisen  = ipcruise2     # last cruise  point
      const ipdescentn = ipdescent5    # last descent point

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
      iototal = 15;