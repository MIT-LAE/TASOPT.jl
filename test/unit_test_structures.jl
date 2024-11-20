@testset "structural components" begin

  # Fuselage
  include(joinpath(TASOPT.__TASOPTroot__, "../test/default_fuselage.jl"))
  fuselage = ac_test.fuselage

  Nland, Wpaymax, Wengtail, 
  nftanks,
  Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
  Δp,
  Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
  bv, λv, nvtail,
  xhtail, xvtail,
  xwing, xwbox, cbox,
  xeng = [6.0, 219964.5779, 0.0, 0, 0.0, 0.0, 0.0, 0.0, "", 
          56016.16406368774, 14401.51302958861, 9645.666840197786, 
          0.4, 0.7, 1.2379954317075564e6, 1.0469874186878382e6, 7.958908725623409, 
          0.3, 1.0, 36.21248821998368, 35.0505846520266, 18.465757410492177, 
          15.688376236863938, 3.119653923461824, 14.164376236863939 ]
  
  cabVol = TASOPT.fusew!(fuselage, Nland, Wpaymax, Wengtail, 
                          nftanks,
                          Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
                          Δp,
                          Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
                          bv, λv, nvtail,
                          xhtail, xvtail,
                          xwing, xwbox, cbox,
                          xeng)
    
    out = [fuselage.skin.thickness, fuselage.cone.thickness, 
            fuselage.layout.thickness_webs, fuselage.floor.thickness, 
            fuselage.bendingmaterial_h.weight.x, fuselage.bendingmaterial_v.weight.x,
            fuselage.shell.EIh, fuselage.bendingmaterial_h.EIh,
            fuselage.shell.EIv, fuselage.bendingmaterial_v.EIh,
            fuselage.shell.GJ , fuselage.cone.GJ,
            fuselage.shell.weight.W, fuselage.cone.weight.W, 
            fuselage.window.W, fuselage.insulation.W, 
            fuselage.floor.weight.W,
            fuselage.bendingmaterial_h.weight.W, fuselage.bendingmaterial_v.weight.W,
            fuselage.weight,fuselage.moment,
            cabVol]

  ## Old FORTRAN outputs just temporarily retained. Better representation of 
  # material properties and sizing consideration (τ instead of σ and proper G vals) means 
  # the new calculations are "better" and more accurate.
  # fort_out = [0.0010592914302080562, 0.0012235633080289908, 0.0, 
  #             0.0009061565981319232, 31.85596336970722, 31.098559122927245, 
  #             2.361319618160982e9, 2.4686136573868637e10, 2.3173490053921146e9, 
  #             2.4686136573868637e10, 1.570840175144195e9, 1.8144415656291723e9, 
  #             20947.348347295865, 5667.632920626953, 11269.980000000001, 
  #             5018.7283576257105, 13415.7630120004, 14314.44885191115, 
  #             7381.15959985673, 198042.54737081684, 3.620266695235501e6, 414.54271063740043]
  
  # for i in eachindex(out)
  #   @test out[i] ≈ fort_out[i] rtol=1e-3
  # end

  test_out = [0.0010592914302080562, 0.0012235812333327315, 0.0,
    0.0009061565981319232, 31.855960860329677, 31.098559122927245,
    2.361319618160982e9, 2.4686144041160843e10, 2.3173490053921146e9,
    2.4686144041160843e10, 1.570840175144195e9, 1.6421579148805373e9,
    20947.348347295865, 5667.715951918366, 11269.980000000001,
    5018.7283576257105, 13415.7630120004, 14314.452366401405,
    7381.155204479751, 198042.62952122153, 3.620269431679788e6,
    414.54271063740043]
    
@test all(isapprox.(out, test_out))

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

# #Test fuselage layout updating for tank

parg = zeros(igtotal)
fuselage.layout.x_start_cylinder = 6.096
fuselage.layout.x_end_cylinder = 29.5656
fuselage.layout.x_pressure_shell_aft = 31.0896
parg[igdxeng2wbox] = 1.5239999999999991
fuselage.APU.r = [36.576, 0.0, 0.0]
fuselage.layout.x_end = 37.7952
fuselage.layout.x_cone_end = 35.6616
parg[igxhbox ] = 34.8996
parg[igxvbox ] = 33.528
parg[igxwbox] = 16.04432532088372
parg[igxeng] = parg[igxwbox] - parg[igdxeng2wbox]
fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 

pari = zeros(iitotal)
pari[iinftanks] = 1

# parg_orig = deepcopy(parg)
# deleteat!(parg_orig, parg_orig .== 0)

#Update fuel tank length and check changes
parg[iglftank] = 5.0
TASOPT.update_fuse!(fuselage, pari, parg)

# parg_check = [43.40480000000001, 6.096, 35.175200000000004, 36.699200000000005, 41.27120000000001, 14.52032532088372, 16.04432532088372, 40.50920000000001, 39.137600000000006, 42.18560000000001, 21.660776608000003, 5.0, 23.4696, 1.5239999999999991, 18.716634144] 
update_fuse_out = [fuselage.layout.x_end_cylinder, 
fuselage.layout.x_pressure_shell_aft, 
fuselage.layout.x_cone_end, 
fuselage.APU.x,
fuselage.layout.x_end,
fuselage.HPE_sys.x,
parg[igxhbox],
parg[igxvbox],
parg[igxeng]]

update_fuse_out_test = [35.175200000000004, 36.699200000000005, 24.326234144000004, 42.18560000000001, 43.40480000000001, 12.767380728136962, 40.50920000000001, 39.137600000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))


#Return to original points?
pari[iinftanks] = 0.0
parg[iglftank] = 0.0
TASOPT.update_fuse!(fuselage, pari, parg)

update_fuse_out = [fuselage.layout.x_end_cylinder, 
fuselage.layout.x_pressure_shell_aft, 
fuselage.layout.x_cone_end, 
fuselage.APU.x,
fuselage.layout.x_end,
fuselage.HPE_sys.x,
parg[igxhbox],
parg[igxvbox],
parg[igxeng]]

update_fuse_out_test = [29.5656, 31.0896, 18.716634144, 36.57600000000001, 37.79520000000001, 9.82323826413696, 34.89960000000001, 33.528000000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))

# #Test cabin resizing
# parg = zeros(igtotal)
# fuselage.layout.x_start_cylinder = 6.096
# fuselage.layout.x_end_cylinder = 29.5656
# fuselage.layout.x_pressure_shell_aft = 31.0896
# parg[igdxcabin] = 23.4696
# parg[igdxeng2wbox] = 1.5239999999999991
# fuselage.APU.r = [36.576, 0.0, 0.0]
# fuselage.layout.x_end = 37.7952
# fuselage.layout.x_cone_end = 35.6616
# parg[igxhbox ] = 34.8996
# parg[igxvbox ] = 33.528
# parg[igxwbox] = 16.04432532088372
# parg[igxeng] = parg[igxwbox] - parg[igdxeng2wbox]
# fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 

# parg[igseatpitch] = 0.762
# parg[igseatwidth] = 0.4826
# parg[igaislehalfwidth] = 0.254
# parg[igWpaymax] = 219964.5779
# fuselage.layout.radius

# pari = zeros(iitotal)
# pari[iidoubledeck] = 0

# parm = zeros(imtotal)
# parm[imWperpax,1] = 956.36773

# fuse_tank = TASOPT.fuselage_tank()

# TASOPT.update_fuse_for_pax!(pari, parg, fuselage, fuse_tank)

# parg_check = [47.091600000000014, 6.096, 38.86200000000001, 40.38600000000001, 44.95800000000001, 18.460895740194808, 19.98489574019481, 44.19600000000001, 42.82440000000001, 45.87240000000001, 23.595756720000004, 1.9558, 219964.5779, 32.76600000000001, 1.5239999999999991, 0.762, 0.4826, 0.254]
# parg_nz = deepcopy(parg)
# deleteat!(parg_nz, parg_nz .== 0)
# for (i,item) in enumerate(parg_nz) #For every nonzero element in parg
#   @test parg_nz[i] ≈ parg_check[i]
# end

#Test landing gear sizing
ac = load_default_model()

#Test simple sizing based on mass fractions
WMTO = 80e4
ac.parg[igWMTO] = WMTO
ac.landing_gear.main_gear.overall_mass_fraction = 0.04
ac.landing_gear.nose_gear.overall_mass_fraction = 0.01

TASOPT.landing_gear_size!(ac)
@test ac.landing_gear.main_gear.weight.W ≈ 0.04 * 80e4
@test ac.landing_gear.nose_gear.weight.W ≈ 0.01 * 80e4

#Test models based on historical data
ac.landing_gear.model = "historical_correlations"
ac.landing_gear.tailstrike_angle = 10*pi/180
ac.landing_gear.wing_dihedral_angle = 6*pi/180
ac.landing_gear.engine_ground_clearance = 7*in_to_m
ac.landing_gear.nose_gear.number_struts = 1
ac.landing_gear.nose_gear.wheels_per_strut = 2
ac.landing_gear.main_gear.number_struts = 2
ac.landing_gear.main_gear.wheels_per_strut = 2

ac.parg[igxCGaft] = 30
ac.pare[ieu0,iprotate,1] = 70
ac.parg[igdfan] = 1.5
ac.parg[igb] = 40
TASOPT.landing_gear_size!(ac)

@test ac.landing_gear.main_gear.weight.W ≈ 25733.572498747944
@test ac.landing_gear.nose_gear.weight.W ≈ 4247.610705426578
end