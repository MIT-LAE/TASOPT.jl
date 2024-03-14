
@testset "structural components" begin

# Fuselage weight
  gee, Nland, Wfix, Wpay, Wpadd, Wseat, Wapu, Wengtail,
  ifwing, nftanks, xblend1, xblend2,
  Waftfuel, Wftank, ltank, xftankaft, tank_placement,
  fstring, fframe, ffadd, deltap,
  Wpwindow, Wppinsul, Wppfloor,
  Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
  bv, lambdav, nvtail,
  Rfuse, dRfuse, wfb, nfweb, lambdac,
  xnose, xshell1, xshell2, xconend,
  xhtail, xvtail,
  xwing, xwbox, cbox,
  xfix, xapu, xeng,
  hfloor,
  sigskin, sigbend, rhoskin, rhobend,
  Eskin, Ebend, Gskin = [9.810,
    6.0, 13350.00, 172215.0, 60275.25, 17221.50, 6027.5250, 0.00,
    1,
    0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    "",
    0.34999999999999998,
    0.25000000000000000,
    0.20000000000000001,
    56108.434242330870,
    435.00000000000000,
    22.000000000000000,
    60.000000000000000,
    11637.065832885963,
    7847.1938058748428,
    0.40000000000000002,
    0.69999999999999996,
    1077432.1225689780,
    912092.01487599779,
    7.4348760062795316,
    0.29999999999999999,
    1.0000000000000000,
    1.9558000000000000,
    0.38100000000000001,
    0.0000000000000000,
    1.0000000000000000,
    0.29999999999999999,
    0.0000000000000000,
    5.1816000000000004,
    31.089600000000001,
    35.661600000000000,
    36.104115424580897,
    34.950334202732712,
    18.941262714160111,
    16.361970523376034,
    2.9420828049786860,
    2.1335999999999999,
    36.576000000000001,
    15.849600000000001,
    0.12700000000000000,
    103448275.86206897,
    206896551.72413793,
    2700.0000000000000,
    2700.0000000000000,
    68965517241.379318,
    68965517241.379318,
    26525198938.992046]

  Waftfuel = 0.0
  xfuel = 0.0
  out = TASOPT.fusew(Nland, Wfix, Wpay, Wpadd, Wseat, Wapu, Wengtail, 
  ifwing, nftanks, xblend1, xblend2,
  Waftfuel, Wftank, ltank, xftankaft,tank_placement,
    fstring, fframe, ffadd, deltap,
    Wpwindow, Wppinsul, Wppfloor,
    Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
    bv, lambdav, nvtail,
    Rfuse, dRfuse, wfb, nfweb, lambdac,
    xnose, xshell1, xshell2, xconend,
    xhtail, xvtail,
    xwing, xwbox, cbox,
    xfix, xapu, xeng, xfuel,
    hfloor,
    sigskin, sigbend, rhoskin, rhobend,
    Eskin, Ebend, Gskin)

  fort_out = [1.0607897983477732E-003,
    9.9550382519854471E-004,
    0.0000000000000000,
    7.0928467266495859E-004,
    30.064013766112691,
    30.406344539593768,
    2941018305.5863328,
    17628813266.493961,
    2321230662.5304923,
    13387246823.029541,
    1573471398.5658336,
    1476632598.2324898,
    20976.978379674187,
    4611.2450539189012,
    11269.980000000001,
    4779.8514160257109,
    11808.103388796697,
    8991.3441092588018,
    5384.5872422420944,
    164696.36458991640,
    2962392.0525438893,
    414.54271063740043]

@test all(isapprox.(out, fort_out))
# end Fuselage weight

# surfdx:
  b = 35.486921629195265
  bs = 10.113772664320649
  bo = 3.6067999999999998
  lambdat = 0.25000000000000000
  lambdas = 0.69999999999999996
  sweep = 26.000000000000000
  dxwing = 2.5792921907840762
  macco = 0.68475212288241949

dx_out, macco_out = TASOPT.structures.surfdx(b,bs,bo,lambdat,lambdas,sweep)
@test dx_out ≈ dxwing
@test macco_out ≈ macco 

# end surfdx


#surfw:
  gee = 9.8100000000000005
  po = 110091.58394892939
  b = 35.486921629195265
  bs = 10.113772664320649
  bo = 3.6067999999999998
  co = 5.8841656099573720
  zs = 3.9116000000000000
  lambdat = 0.25000000000000000
  lambdas = 0.69999999999999996
  gammat = 0.22500000000000001
  gammas = 0.86659999999999993
  Nlift = 3.0000000000000000
  iwplan = 1
  Weng1 = 27312.133624204889
  Winn = 48907.685542960695
  Wout = 70390.482266430423
  dyWinn = 70296.379288596174
  dyWout = 315316.83291120041
  sweep = 26.000000000000000
  wbox = 0.50000000000000000
  hboxo = 0.12680000000000000
  hboxs = 0.12659999999999999
  rh = 0.75000000000000000
  fLt = -5.0000000000000003E-002
  tauweb = 137931034.48275861
  sigcap = 206896551.72413793
  sigstrut = 206896551.72413793
  Ecap = 68965517241.379318
  Eweb = 68965517241.379318
  Gcap = 26525198938.992046
  Gweb = 26525198938.992046
  rhoweb = 2700.0000000000000
  rhocap = 2700.0000000000000
  rhostrut = 2700.0000000000000
  rhofuel = 817.00000000000000
  fort_Ss = 549317.34029970923
  fort_Ms = 2919107.4513926785
  fort_tbwebs = 1.5302067790179924E-003
  fort_tbcaps = 6.4712234567801050E-003
  fort_EIcs = 256525134.96391767
  fort_EIns = 1981689904.7351539
  fort_GJs = 220465338.41465056
  fort_So = 654948.64283556491
  fort_Mo = 4759686.3416206865
  fort_tbwebo = 8.9257473924635666E-004
  fort_tbcapo = 3.4237821474538827E-003
  fort_EIco = 598779822.13510609
  fort_EIno = 4420836188.7699680
  fort_GJo = 529475637.81714487
  fort_Astrut = 0.0000000000000000
  fort_lstrutp = 0.0000000000000000
  fort_cosLs = 1.0000000000000000
  fort_Wscen = 5943.1718341806891
  fort_Wsinn = 9074.6000905588244
  fort_Wsout = 17146.749196179870
  fort_dxWsinn = 6361.5817898457608
  fort_dxWsout = 64671.501022901110
  fort_dyWsinn = 13043.175582241660
  fort_dyWsout = 76809.512837201284
  fort_Wfcen = 27272.842046818590
  fort_Wfinn = 31647.789098991288
  fort_Wfout = 39316.171146500550
  fort_dxWfinn = 22186.101515424907
  fort_dxWfout = 148286.75531592543
  fort_dyWfinn = 45488.249166743939
  fort_dyWfout = 176118.27862155536
  fort_Wweb = 2855.5435405652265
  fort_Wcap = 61473.498701273536
  fort_Wstrut = 0.0000000000000000
  fort_dxWweb = 6133.4918937186821
  fort_dxWcap = 135932.67373177505
  fort_dxWstrut = 0.0000000000000000

neout = 0
neinn = 1
dyeout = 0 
dyeinn = 0.5*(bs - bo)

Ss,Ms,tbwebs,tbcaps,EIcs,EIns,GJs,
So,Mo,tbwebo,tbcapo,EIco,EIno,GJo,
Astrut,lsp,cosLs,
Wscen,Wsinn,Wsout,dxWsinn,dxWsout,dyWsinn,dyWsout,
Wfcen,Wfinn,Wfout,dxWfinn,dxWfout,dyWfinn,dyWfout,
Wweb,  Wcap,  Wstrut,
dxWweb,dxWcap,dxWstrut = TASOPT.surfw(po,b,bs,bo,co,zs,
lambdat,lambdas,gammat,gammas,
Nlift,iwplan,Weng1,neout, dyeout, neinn, dyeinn,
Winn,Wout,dyWinn,dyWout,
sweep,wbox,hboxo,hboxs,rh, fLt,
tauweb,sigcap,sigstrut,Ecap,Eweb,Gcap,Gweb,
rhoweb,rhocap,rhostrut,rhofuel)


@test fort_Ss ≈ Ss 
@test fort_Ms ≈ Ms 
@test fort_tbwebs ≈ tbwebs 
@test fort_tbcaps ≈ tbcaps 
@test fort_EIcs ≈ EIcs 
@test fort_EIns ≈ EIns 
@test fort_GJs ≈ GJs 
@test fort_So ≈ So 
@test fort_Mo ≈ Mo 
@test fort_tbwebo ≈ tbwebo 
@test fort_tbcapo ≈ tbcapo 
@test fort_EIco ≈ EIco 
@test fort_EIno ≈ EIno 
@test fort_GJo ≈ GJo 
@test fort_Astrut ≈ Astrut 
@test fort_lstrutp ≈ lsp 
@test fort_cosLs ≈ cosLs 
@test fort_Wscen ≈ Wscen 
@test fort_Wsinn ≈ Wsinn 
@test fort_Wsout ≈ Wsout 
@test fort_dxWsinn ≈ dxWsinn 
@test fort_dxWsout ≈ dxWsout 
@test fort_dyWsinn ≈ dyWsinn 
@test fort_dyWsout ≈ dyWsout 
@test fort_Wfcen ≈ Wfcen 
@test fort_Wfinn ≈ Wfinn 
@test fort_Wfout ≈ Wfout 
@test fort_dxWfinn ≈ dxWfinn 
@test fort_dxWfout ≈ dxWfout 
@test fort_dyWfinn ≈ dyWfinn 
@test fort_dyWfout ≈ dyWfout 
@test fort_Wweb ≈ Wweb 
@test fort_Wcap ≈ Wcap 
@test fort_Wstrut ≈ Wstrut 
@test fort_dxWweb ≈ dxWweb 
@test fort_dxWcap ≈ dxWcap 
@test fort_dxWstrut ≈ dxWstrut 
#end surfw

#tailpo:
  Sh = 42.443587259163984
  ARh = 6.0000000000000000
  lambdah = 0.25000000000000000
  qne = 12692.519555311032
  CLhmax = 2.0000000000000000
  fort_bh = 15.958117794871169
  fort_coh = 4.2554980786323124
  fort_poh = 108025.98516125829

bh, coh, poh = TASOPT.structures.tailpo(Sh, ARh, lambdah, qne, CLhmax)

@test fort_bh ≈ bh
@test fort_coh ≈ coh
@test fort_poh ≈ poh
#end tailpo:

#TODO write tests for H2 tanks 
end