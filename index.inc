
#---- indices for integer variables and flags
      iifuel   =  1  # index specifying the fuel (see gasfun.f)          
      iifwcen  =  2  # 0 = no fuel in center box  1 = fuel in center box
      iiwplan  =  3  # 0 = bare wing 1 = cantilever with engine 2 = strut braced
      iiengloc =  4  # 1 = engines on wing 2 = engines on fuselage   
      iiengwgt =  5  # 1 = basic tech eng. weight 2 = advanced tech eng. weight 
      iiBLIc   =  6  # 0 = core in clean flow 1 = core ingests KE defect   
      iifclose =  7  # 0 = fuse tapers to point 1 = tapers to edge 
      iiHTsize =  8  # 1 = set Sh via Vh 2 = set Sh via CLh at forward CG trim
      iiVTsize =  9  # 1 = set Sv via Vv 2 = set Sv via CLv at engine-out trim
      iixwmove = 10  # 0 = wing centroid fixed 1 = move wing to balance
      iitotal  = 10 

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
      igFOpt     =   1
      igPFEI     =   2
      igRange    =   3
      igWMTO     =   4
      igWpay     =   5
      igWfix     =   6
      igWfuel    =   7
      igWfmax    =   8
      igrWfmax   =   9
      igWshell   =  10
      igWwindow  =  11
      igWinsul   =  12
      igWfloor   =  13
      igWcone    =  14
      igWhbend   =  15
      igWvbend   =  16
      igWfuse    =  17
      igWweb     =  18
      igWcap     =  19
      igWwing    =  20
      igWebare   =  21
      igWnace    =  22
      igWeng     =  23
      igWhtail   =  24
      igWvtail   =  25
      igWstrut   =  26
      igxWfuse   =  27
      igdxWfuel  =  28
      igdxWwing  =  29
      igdxWstrut =  30
      igdxWhtail =  31
      igdxWvtail =  32
      igWinn     =  33
      igWout     =  34
      igdyWinn   =  35
      igdyWout   =  36
      igxCGfwd   =  37
      igxCGaft   =  38
      igfreserve =  39
      igfpadd    =  40
      igfseat    =  41
      igfeadd    =  42
      igfpylon   =  43
      igfnace    =  44
      igfflap    =  45
      igfslat    =  46
      igfaile    =  47
      igflete    =  48
      igfribs    =  49
      igfspoi    =  50
      igfwatt    =  51
      igfhadd    =  52
      igfvadd    =  53
      igfapu     =  54
      igfhpesys  =  55
      igflgnose  =  56
      igflgmain  =  57
      igfstring  =  58
      igfframe   =  59
      igffadd    =  60
      igWpwindow =  61
      igWppinsul =  62
      igWppfloor =  63
      igNlift    =  64
      igNland    =  65
      igVne      =  66
      igneng     =  67
      igGearf    =  68
      igfTt4CL1  =  69
      igfTt4CLn  =  70
      igHTRf     =  71
      igHTRlc    =  72
      igHTRhc    =  73
      igrSnace   =  74
      igrVnace   =  75
      igrVstrut  =  76
      igfSnace   =  77
      igpcabin   =  78
      igdeltap   =  79
      iganose    =  80 
      igbtail    =  81
      igxnose    =  82
      igxend     =  83
      igxblend1  =  84
      igxblend2  =  85
      igxshell1  =  86
      igxshell2  =  87
      igxconend  =  88
      igxhbend   =  89
      igxvbend   =  90
      igxhtail   =  91 
      igxvtail   =  92
      igxeng     =  93     
      igxwing    =  94
      igxwbox    =  95
      igxhbox    =  96
      igxvbox    =  97
      igxfix     =  98
      igxapu     =  99
      igxhpesys  = 100
      igxlgnose  = 101
      igdxlgmain = 102
      igyeng     = 103
      igzwing    = 104
      igzhtail   = 105
      ignfweb    = 106
      igwfb      = 107
      igRfuse    = 108
      igdRfuse   = 109
      ighfloor   = 110
      iglambdac  = 111
      igcabVol   = 112
      igcosLs    = 113
      igSstrut   = 114
      igrpayfwd  = 115
      igrpayaft  = 116
      igxNP      = 117
      igCMVf1    = 118
      igCLMf0    = 119
      igdepsda   = 120
      igdCLnda   = 121
      igdCLhdCL  = 122
      igdCLndCL  = 123
      igCLhspec  = 124
      igCLhCGfwd = 125
      igCLveout  = 126
      igCLhmax   = 127
      igCLvmax   = 128
      igfCDhcen  = 129
      igSMmin    = 130
      igrMh      = 131
      igrMv      = 132
      igXaxis    = 133
      igwbox     = 134
      ighboxo    = 135
      ighboxs    = 136
      igrh       = 137
      igwboxh    = 138
      ighboxh    = 139
      igrhh      = 140
      igwboxv    = 141
      ighboxv    = 142
      igrhv      = 143
      igsigfac   = 144
      igsigskin  = 145
      igsigbend  = 146
      igsigcap   = 147
      igtauweb   = 148
      igsigstrut = 149
      igrEshell  = 150
      igEcap     = 151
      igEstrut   = 152
      igrhoskin  = 153
      igrhobend  = 154
      igrhocap   = 155
      igrhoweb   = 156
      igrhostrut = 157
      igrhofuel  = 158
      igrcls     = 159
      igrclt     = 160
      igCLhNrat  = 161
      igSomax    = 162
      igMomax    = 163
      igSsmax    = 164
      igMsmax    = 165
      igtbcapo   = 166
      igtbwebo   = 167
      igtbcaps   = 168
      igtbwebs   = 169
      igtbcaph   = 170
      igtbwebh   = 171
      igtbcapv   = 172
      igtbwebv   = 173
      igEIco     = 174
      igEIno     = 175
      igGJo      = 176
      igEIcs     = 177
      igEIns     = 178
      igGJs      = 179
      igEIch     = 180
      igEInh     = 181
      igGJh      = 182
      igEIcv     = 183
      igEInv     = 184
      igGJv      = 185
      igtskin    = 186
      igtcone    = 187
      igtfweb    = 188
      igtfloor   = 189
      igEIhshell = 190
      igEIhbend  = 191
      igEIvshell = 192
      igEIvbend  = 193
      igGJshell  = 194
      igGJcone   = 195
      igfLo      = 196
      igfLt      = 197
      igfLn      = 198
      igcma      = 199
      igAR       = 200
      igS        = 201
      igb        = 202
      igbo       = 203
      igbs       = 204
      igetas     = 205
      iglambdat  = 206
      iglambdas  = 207
      igco       = 208
      igsweep    = 209
      igVh       = 210
      igARh      = 211
      igSh       = 212
      igbh       = 213
      igboh      = 214
      iglambdah  = 215
      igcoh      = 216
      igsweeph   = 217
      igVv       = 218
      igARv      = 219
      igSv       = 220
      igbv       = 221
      igbov      = 222
      iglambdav  = 223
      igcov      = 224
      igsweepv   = 225
      ignvtail   = 226
      igzs       = 227
      ighstrut   = 228
      igAstrut   = 229
      igcstrut   = 230
      igfBLIw    = 231
      igfBLIf    = 232
      igdfan     = 233
      igdlcomp   = 234
      igdhcomp   = 235
      iglnace    = 236
      igA5       = 237
      igA7       = 238
      igTmetal   = 239
      igcdefan   = 240
      igCDgear   = 241
      igCDspoil  = 242
      igmuroll   = 243
      igmubrake  = 244
      ighobst    = 245
      iglBFmax   = 246
      igbmax     = 247
      iggtocmin  = 248
      igdBSLmax  = 249
      igdBCBmax  = 250
      igmofWpay  = 251
      igmofWMTO  = 252
      igPofWpay  = 253
      igPofWMTO  = 254
      igtotal    = 254 

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
      iatotal   = 61 

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
      ietotal = 195 

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
