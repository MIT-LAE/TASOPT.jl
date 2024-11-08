
#---- indices for integer variables and flags
      iifuel   =  1  # index specifying the fuel (see gasfun.f)    0 = JetA 1 = LH2      
      iifwcen  =  2  # 0 = no fuel in center box  1 = fuel in center box
      iiwplan  =  3  # 0 = bare wing 1 = cantilever with engine 2 = strut braced
      iiengloc =  4  # 1 = engines on wing 2 = engines on fuselage "Options.engine_location"  
      iiengwgt =  5  # 1 = basic tech eng. weight 2 = advanced tech eng. weight 
      iiBLIc   =  6  # 0 = core in clean flow 1 = core ingests KE defect   
      iifclose =  7  # 0 = fuse tapers to point 1 = tapers to edge 
      iiHTsize =  8  # 1 = set Sh via Vh 2 = set Sh via CLh at forward CG trim
      iiVTsize =  9  # 1 = set Sv via Vv 2 = set Sv via CLv at engine-out trim
      iixwmove = 10  # 0 = wing centroid fixed 1 = move wing to balance
      iifwing  = 11  # 0 = no fuel in wing (for LH2) 1 = fuel in wing
      iicalcCMVf = 12
      iiengtype  = 13  # 0 = Turboelectric engine; 1 = Turbofan engine; "prop_sys_arch"
      iiopt      = 14  # 0 = run sizing loop only; 1 = run optimization procedure "optimize"
      iiaircraftclass = 15 # 737 = 737 size aircraft; 777 = 777 size aircraft
      iinftanks = 16 # number of fuel tanks in fuselage
      iiengmodel = 17 #engine model 
      iidoubledeck = 18 #double decker flag
      iitotal    = 18

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
      const igWweb     =  9
      const igWcap     =  10
      const igWwing    =  11
      const igWebare   =  12
      const igWnace    =  13
      const igWeng     =  14
      const igWhtail   =  15
      const igWvtail   =  16
      const igWstrut   =  17
### Good Riddance ###
      const igdxWfuel  =  18
      const igdxWwing  =  19
      const igdxWstrut =  20
      const igdxWhtail =  21
      const igdxWvtail =  22
      const igWinn     =  23
      const igWout     =  24
      const igdyWinn   =  25
      const igdyWout   =  26
      const igxCGfwd   =  27
      const igxCGaft   =  28
      const igfreserve =  29
### Good Riddance ###
### Good Riddance ###
      const igfeadd    =  30
      const igfpylon   =  31
### Good Riddance ###
      const igfflap    =  32
      const igfslat    =  33
      const igfaile    =  34
      const igflete    =  35
      const igfribs    =  36
      const igfspoi    =  37
      const igfwatt    =  38
      const igfhadd    =  39
      const igfvadd    =  40
### Good Riddance ###
### Good Riddance ###
      const igflgnose  =  41
      const igflgmain  =  42
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igNlift    =  43
      const igNland    =  44
      const igVne      =  45
      const igneng     =  46
      const igGearf    =  47
      const igfTt4CL1  =  48
      const igfTt4CLn  =  49
      const igHTRf     =  50
      const igHTRlc    =  51
      const igHTRhc    =  52
      const igrSnace   =  53
      const igrVnace   =  54
      const igrVstrut  =  55
      const igfSnace   =  56
      const igpcabin   =  57
      const igdeltap   =  58
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
      const igxhtail   =  59 
      const igxvtail   =  60
      const igxeng     =  61     
      const igxwing    =  62
      const igxwbox    =  63
      const igxhbox    =  64
      const igxvbox    =  65
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igxlgnose  = 66
      const igdxlgmain = 67
      const igyeng     = 68
      const igzwing    = 69
      const igzhtail   = 70
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
### Good Riddance ###
      const igcabVol   = 71
      const igcosLs    = 72
      const igSstrut   = 73
      const igrpayfwd  = 74
      const igrpayaft  = 75
      const igxNP      = 76
      const igCMVf1    = 77
      const igCLMf0    = 78
      const igdepsda   = 79
      const igdCLnda   = 80
      const igdCLhdCL  = 81
      const igdCLndCL  = 82
      const igCLhspec  = 83
      const igCLhCGfwd = 84
      const igCLveout  = 85
      const igCLhmax   = 86
      const igCLvmax   = 87
      const igfCDhcen  = 88
      const igSMmin    = 89
      const igrMh      = 90
      const igrMv      = 91
      const igXaxis    = 92
      const igwbox     = 93
      const ighboxo    = 94
      const ighboxs    = 95
      const igrh       = 96
      const igwboxh    = 97
      const ighboxh    = 98
      const igrhh      = 99
      const igwboxv    = 100
      const ighboxv    = 101
      const igrhv      = 102
      const igsigfac   = 103
### Good Riddance ###
### Good Riddance ###
      const igsigcap   = 104
      const igtauweb   = 105
      const igsigstrut = 106
### Good Riddance ###
      const igEcap     = 107
      const igEstrut   = 108
### Good Riddance ###
### Good Riddance ###
      const igrhocap   = 109
      const igrhoweb   = 110
      const igrhostrut = 111
      const igrhofuel  = 112
      const igrcls     = 113
      const igrclt     = 114
      const igCLhNrat  = 115
      const igSomax    = 116
      const igMomax    = 117
      const igSsmax    = 118
      const igMsmax    = 119
      const igtbcapo   = 120
      const igtbwebo   = 121
      const igtbcaps   = 122
      const igtbwebs   = 123
      const igtbcaph   = 124
      const igtbwebh   = 125
      const igtbcapv   = 126
      const igtbwebv   = 127
      const igEIco     = 128
      const igEIno     = 129
      const igGJo      = 130
      const igEIcs     = 131
      const igEIns     = 132
      const igGJs      = 133
      const igEIch     = 134
      const igEInh     = 135
      const igGJh      = 136
      const igEIcv     = 137
      const igEInv     = 138
      const igGJv      = 139
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
      const igfLo      = 140
      const igfLt      = 141
      const igfLn      = 142
      const igcma      = 143
      const igAR       = 144
      const igS        = 145
      const igb        = 146
      const igbo       = 147
      const igbs       = 148
      const igetas     = 149
      const iglambdat  = 150
      const iglambdas  = 151
      const igco       = 152
      const igsweep    = 153
      const igVh       = 154
      const igARh      = 155
      const igSh       = 156
      const igbh       = 157
      const igboh      = 158
      const iglambdah  = 159
      const igcoh      = 160
      const igsweeph   = 161
      const igVv       = 162
      const igARv      = 163
      const igSv       = 164
      const igbv       = 165
      const igbov      = 166
      const iglambdav  = 167
      const igcov      = 168
      const igsweepv   = 169
      const ignvtail   = 170
      const igzs       = 171
      const ighstrut   = 172
      const igAstrut   = 173
      const igcstrut   = 174
      const igfBLIw    = 175
      const igfBLIf    = 176
      const igdfan     = 177
      const igdlcomp   = 178
      const igdhcomp   = 179
      const iglnace    = 180
      const igA5       = 181
      const igA7       = 182
      const igTmetal   = 183
      const igcdefan   = 184
      const igCDgear   = 185
      const igCDspoil  = 186
      const igmuroll   = 187
      const igmubrake  = 188
      const ighobst    = 189
      const iglBFmax   = 190
      const igbmax     = 191
      const iggtocmin  = 192
      const igdBSLmax  = 193
      const igdBCBmax  = 194
      const igmofWpay  = 195
      const igmofWMTO  = 196
      const igPofWpay  = 197
      const igPofWMTO  = 198

      # Adding indices to store TE weights and locations
      const igWtshaft  = 199
      const igWgen     = 200
      const igWinv     = 201
      const igWmot     = 202
      const igWfan     = 203
      const igWftank   = 204

      const igxtshaft  = 205
      const igxgen     = 206
      const igxinv     = 207
      const igxmot     = 208
      const igxfan     = 209
      const igxftank   = 210
      
      const igxcables  = 211
      const igWcables  = 212

      const igxcat     = 213
      const igWcat     = 214

      const igWtesys   = 215
      const igxWtesys  = 216

      const iglftank   = 217
      const igWinsftank= 218
      const igxWftank  = 219
      
      const igRftank   = 220
      const igWc3des   = 221

      const igdaftfan  = 222
      const iglnaceaft = 223
### Good Riddance ###
      const igneout    = 224
      const igyeout    = 225
      const igyeinn    = 226
      const iglftankin = 227 #lftank input
      const igLHVfuel  = 228
      const igWfburn   = 229

      const igWaftfan  = 230
      const igWfanGB   = 231
      const igWaftfanGB   = 232
      const igWrect    = 233
      const igWtms     = 234

      const igfanPCT   = 235

      const iggamVtarget = 236

      const igWpaymax    = 237

      const ig_YC = 238
      const ig_YH = 239
      const ig_YO = 240
      const ig_YN = 241
      const ig_MW = 242
      const ighfuel      = 243
      const igTfuel      = 244

      const igxftankaft  = 245
      const igdxeng2wbox = 246

      const igxWfuel = 247

      const igseatpitch = 248
      const igseatwidth = 249
      const igaislehalfwidth = 250
      const igfloordist = 251

      const igWfvent = 252
      const igdHXPreC = 253
      const igdHXInterC = 254
      const igdHXRegen = 255
      const igHXaddmassfrac = 256
      const igWHXs = 257

      const igtotal    = 257

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
      iehvap = 244
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
      ieRadiatorepsilon = 265
      ieRadiatorMp = 266
      ieRadiatorDeltah = 267
      ieRadiatorDeltap = 268
      iehvapcombustor = 269

      ietotal = 269

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
