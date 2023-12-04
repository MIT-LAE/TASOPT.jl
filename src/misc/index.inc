
#---- indices for integer variables and flags
      iifuel   =  1  # index specifying the fuel (see gasfun.f)    0 = JetA 1 = LH2      
      iifwcen  =  2  # 0 = no fuel in center box  1 = fuel in center box
      iiwplan  =  3  # 0 = bare wing 1 = cantilever with engine 2 = strut braced
      iiengloc =  4  # 1 = engines on wing 2 = engines on fuselage   
      iiengwgt =  5  # 1 = basic tech eng. weight 2 = advanced tech eng. weight 
      iiBLIc   =  6  # 0 = core in clean flow 1 = core ingests KE defect   
      iifclose =  7  # 0 = fuse tapers to point 1 = tapers to edge 
      iiHTsize =  8  # 1 = set Sh via Vh 2 = set Sh via CLh at forward CG trim
      iiVTsize =  9  # 1 = set Sv via Vv 2 = set Sv via CLv at engine-out trim
      iixwmove = 10  # 0 = wing centroid fixed 1 = move wing to balance
      iifwing  = 11  # 0 = no fuel in wing (for LH2) 1 = fuel in wing
      iicalcCMVf = 12
      iiengtype  = 13  # 0 = Turboelectric engine; 1 = Turbofan engine
      iiopt      = 14  # 0 = run sizing loop only; 1 = run optimization procedure
      iiaircraftclass = 15 # 737 = 737 size aircraft; 777 = 777 size aircraft
      iiengmodel = 15
      iitotal    = 16

#---- indices for missions to be examined 
#-      (different Ranges Payloads ...)
      imwOpt     =  1 
      imRange    =  2
      imWpay     =  3 
      imaltTO    =  4
      imT0TO     =  5
      imWTO      =  6 
      imWfuel    =  7 
      imPFEI     =  8 
      imV1       =  9
      imV2       = 10
      imtTO      = 11
      imFTO      = 12
      iml1       = 13
      imlTO      = 14
      imlBF      = 15
      imlCB      = 16
      imgamVTO   = 17
      imgamVBF   = 18
      imgamVCB   = 19
      imgamVDE1  = 20
      imgamVDEn  = 21
      imthCB     = 22
      imxCB      = 23
      imzCB      = 24
      imxFO      = 25
      imzFO      = 26
      imdBSL     = 27
      imdBCB     = 28
      imdBFO     = 29
      imfexcdw   = 30
      imfexcdt   = 31
      imfexcdf   = 32
      imtotal    = 32 

#---- indices for geometry (airframe) variables and other sizing variables
      const igFOpt     =   1
      const igPFEI     =   2
      const igRange    =   3
      const igWMTO     =   4
      const igWpay     =   5
      const igWfix     =   6
      const igWfuel    =   7
      const igWfmax    =   8
      const igrWfmax   =   9
      const igWshell   =  10
      const igWwindow  =  11
      const igWinsul   =  12
      const igWfloor   =  13
      const igWcone    =  14
      const igWhbend   =  15
      const igWvbend   =  16
      const igWfuse    =  17
      const igWweb     =  18
      const igWcap     =  19
      const igWwing    =  20
      const igWebare   =  21
      const igWnace    =  22
      const igWeng     =  23
      const igWhtail   =  24
      const igWvtail   =  25
      const igWstrut   =  26
      const igxWfuse   =  27
      const igdxWfuel  =  28
      const igdxWwing  =  29
      const igdxWstrut =  30
      const igdxWhtail =  31
      const igdxWvtail =  32
      const igWinn     =  33
      const igWout     =  34
      const igdyWinn   =  35
      const igdyWout   =  36
      const igxCGfwd   =  37
      const igxCGaft   =  38
      const igfreserve =  39
      const igfpadd    =  40
      const igfseat    =  41
      const igfeadd    =  42
      const igfpylon   =  43
      const igfnace    =  44
      const igfflap    =  45
      const igfslat    =  46
      const igfaile    =  47
      const igflete    =  48
      const igfribs    =  49
      const igfspoi    =  50
      const igfwatt    =  51
      const igfhadd    =  52
      const igfvadd    =  53
      const igfapu     =  54
      const igfhpesys  =  55
      const igflgnose  =  56
      const igflgmain  =  57
      const igfstring  =  58
      const igfframe   =  59
      const igffadd    =  60
      const igWpwindow =  61
      const igWppinsul =  62
      const igWppfloor =  63
      const igNlift    =  64
      const igNland    =  65
      const igVne      =  66
      const igneng     =  67
      const igGearf    =  68
      const igfTt4CL1  =  69
      const igfTt4CLn  =  70
      const igHTRf     =  71
      const igHTRlc    =  72
      const igHTRhc    =  73
      const igrSnace   =  74
      const igrVnace   =  75
      const igrVstrut  =  76
      const igfSnace   =  77
      const igpcabin   =  78
      const igdeltap   =  79
      const iganose    =  80 
      const igbtail    =  81
      const igxnose    =  82
      const igxend     =  83
      const igxblend1  =  84
      const igxblend2  =  85
      const igxshell1  =  86
      const igxshell2  =  87
      const igxconend  =  88
      const igxhbend   =  89
      const igxvbend   =  90
      const igxhtail   =  91 
      const igxvtail   =  92
      const igxeng     =  93     
      const igxwing    =  94
      const igxwbox    =  95
      const igxhbox    =  96
      const igxvbox    =  97
      const igxfix     =  98
      const igxapu     =  99
      const igxhpesys  = 100
      const igxlgnose  = 101
      const igdxlgmain = 102
      const igyeng     = 103
      const igzwing    = 104
      const igzhtail   = 105
      const ignfweb    = 106
      const igwfb      = 107
      const igRfuse    = 108
      const igdRfuse   = 109
      const ighfloor   = 110
      const iglambdac  = 111
      const igcabVol   = 112
      const igcosLs    = 113
      const igSstrut   = 114
      const igrpayfwd  = 115
      const igrpayaft  = 116
      const igxNP      = 117
      const igCMVf1    = 118
      const igCLMf0    = 119
      const igdepsda   = 120
      const igdCLnda   = 121
      const igdCLhdCL  = 122
      const igdCLndCL  = 123
      const igCLhspec  = 124
      const igCLhCGfwd = 125
      const igCLveout  = 126
      const igCLhmax   = 127
      const igCLvmax   = 128
      const igfCDhcen  = 129
      const igSMmin    = 130
      const igrMh      = 131
      const igrMv      = 132
      const igXaxis    = 133
      const igwbox     = 134
      const ighboxo    = 135
      const ighboxs    = 136
      const igrh       = 137
      const igwboxh    = 138
      const ighboxh    = 139
      const igrhh      = 140
      const igwboxv    = 141
      const ighboxv    = 142
      const igrhv      = 143
      const igsigfac   = 144
      const igsigskin  = 145
      const igsigbend  = 146
      const igsigcap   = 147
      const igtauweb   = 148
      const igsigstrut = 149
      const igrEshell  = 150
      const igEcap     = 151
      const igEstrut   = 152
      const igrhoskin  = 153
      const igrhobend  = 154
      const igrhocap   = 155
      const igrhoweb   = 156
      const igrhostrut = 157
      const igrhofuel  = 158
      const igrcls     = 159
      const igrclt     = 160
      const igCLhNrat  = 161
      const igSomax    = 162
      const igMomax    = 163
      const igSsmax    = 164
      const igMsmax    = 165
      const igtbcapo   = 166
      const igtbwebo   = 167
      const igtbcaps   = 168
      const igtbwebs   = 169
      const igtbcaph   = 170
      const igtbwebh   = 171
      const igtbcapv   = 172
      const igtbwebv   = 173
      const igEIco     = 174
      const igEIno     = 175
      const igGJo      = 176
      const igEIcs     = 177
      const igEIns     = 178
      const igGJs      = 179
      const igEIch     = 180
      const igEInh     = 181
      const igGJh      = 182
      const igEIcv     = 183
      const igEInv     = 184
      const igGJv      = 185
      const igtskin    = 186
      const igtcone    = 187
      const igtfweb    = 188
      const igtfloor   = 189
      const igEIhshell = 190
      const igEIhbend  = 191
      const igEIvshell = 192
      const igEIvbend  = 193
      const igGJshell  = 194
      const igGJcone   = 195
      const igfLo      = 196
      const igfLt      = 197
      const igfLn      = 198
      const igcma      = 199
      const igAR       = 200
      const igS        = 201
      const igb        = 202
      const igbo       = 203
      const igbs       = 204
      const igetas     = 205
      const iglambdat  = 206
      const iglambdas  = 207
      const igco       = 208
      const igsweep    = 209
      const igVh       = 210
      const igARh      = 211
      const igSh       = 212
      const igbh       = 213
      const igboh      = 214
      const iglambdah  = 215
      const igcoh      = 216
      const igsweeph   = 217
      const igVv       = 218
      const igARv      = 219
      const igSv       = 220
      const igbv       = 221
      const igbov      = 222
      const iglambdav  = 223
      const igcov      = 224
      const igsweepv   = 225
      const ignvtail   = 226
      const igzs       = 227
      const ighstrut   = 228
      const igAstrut   = 229
      const igcstrut   = 230
      const igfBLIw    = 231
      const igfBLIf    = 232
      const igdfan     = 233
      const igdlcomp   = 234
      const igdhcomp   = 235
      const iglnace    = 236
      const igA5       = 237
      const igA7       = 238
      const igTmetal   = 239
      const igcdefan   = 240
      const igCDgear   = 241
      const igCDspoil  = 242
      const igmuroll   = 243
      const igmubrake  = 244
      const ighobst    = 245
      const iglBFmax   = 246
      const igbmax     = 247
      const iggtocmin  = 248
      const igdBSLmax  = 249
      const igdBCBmax  = 250
      const igmofWpay  = 251
      const igmofWMTO  = 252
      const igPofWpay  = 253
      const igPofWMTO  = 254

      # Adding indices to store TE weights and locations
      const igWtshaft  = 255
      const igWgen     = 256
      const igWinv     = 257
      const igWmot     = 258
      const igWfan     = 259
      const igWftank   = 260

      const igxtshaft  = 261
      const igxgen     = 262
      const igxinv     = 263
      const igxmot     = 264
      const igxfan     = 265
      const igxftank   = 266
      
      const igxcables  = 267
      const igWcables  = 268

      const igxcat     = 269
      const igWcat     = 270

      const igWtesys   = 271
      const igxWtesys  = 272

      const iglftank   = 273
      const igWinsftank= 274
      const igxWftank  = 275
      
      const igRftank   = 276
      const igWc3des   = 277

      const igdaftfan  = 278
      const iglnaceaft = 279
      const igfuseVol  = 280
      const igneout    = 281
      const igyeout    = 282
      const igyeinn    = 283
      const iglftankin = 284 #lftank input
      const igLHVfuel  = 285
      const igWfburn   = 286

      const igWaftfan  = 287
      const igWfanGB   = 288
      const igWaftfanGB   = 289
      const igWrect    = 290
      const igWtms     = 291

      const igfanPCT   = 292

      const iggamVtarget = 293

      const igWpaymax    = 294

      const ig_YC = 295
      const ig_YH = 296
      const ig_YO = 297
      const ig_YN = 298
      const ig_MW = 299
      const ighfuel      = 300
      const igTfuel      = 301

      const igtotal    = 301

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
      ieFe    = 175 # now stores total thrust & NOT thrust per engine
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
      iegsPRhtrb  = 239
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
      iePreCMh = 247
      iePreCDeltah = 248
      iePreCDeltap = 249
      ieInterCorder = 250
      ieInterCepsilon = 251 
      ieInterCMh = 252
      ieInterCDeltah = 253
      ieInterCDeltap = 254
      ieRegenorder = 255
      ieRegenepsilon = 256 
      ieRegenMh = 257
      ieRegenDeltah = 258
      ieRegenDeltap = 259
      ieTurbCorder = 260
      ieTurbCepsilon = 261 
      ieTurbCMh = 262
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
