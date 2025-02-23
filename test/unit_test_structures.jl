@testset "structural components" begin

  # Fuselage
  include(joinpath(TASOPT.__TASOPTroot__, "../test/default_structures.jl"))
  fuselage = ac_test.fuselage

  Nland, Wpaymax, Wengtail, 
                          nftanks,
                          Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
                          Δp,
                          Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
                          bv, λv, nvtail,
                          xhtail, xvtail,
                          xwing, xwbox, cbox,
                          xeng = [6.0, 219964.5779, 0.0, 0, 0.0, 0.0, 0.0, 0.0, "", 56016.1756316446, 14400.815502104786, 
        9645.22120844135, 0.4, 0.7, 1.2379754637499652e6, 1.0469704911343623e6, 7.958844386073947, 
        0.3, 1.0, 36.212475021199154, 35.05057234350367, 18.465770971525686, 15.688410884477133, 
        3.119635043953234, 14.164410884477133 ]
  
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

  test_out = [0.0010592916489634418, 0.0012236016395195024, 0.0, 0.0009061565981319232, 31.85585727200865, 31.098482101397437, 2.361320105799576e9, 2.4685708558261826e10, 2.317349483950291e9, 2.4685708558261826e10, 1.5708404995400493e9, 1.6421853018493836e9, 20947.352673154768, 5667.810474837753, 11269.980000000001, 5018.7283576257105, 13415.796495703155, 14314.079534858021, 7380.701144174144, 198041.9349618536, 3.6202561572463284e6, 414.54271063740043]

@test all(isapprox.(out, test_out))

# end Fuselage weight

# calculate_centroid_offset:
  b = 35.486921629195265
  bs = 10.113772664320649
  bo = 3.6067999999999998
  lambdat = 0.25000000000000000
  lambdas = 0.69999999999999996
  sweep = 26.000000000000000
  dxwing = 2.5792921907840762
  macco = 0.68475212288241949

dx_out, macco_out = TASOPT.structures.calculate_centroid_offset(b,bs,bo,lambdat,lambdas,sweep)
@test dx_out ≈ dxwing
@test macco_out ≈ macco 

# end calculate_centroid_offset

#size_wing_section
sigfac = 1.0
test_tbwebs, test_tbcaps, test_Abcaps, test_Abwebs = 0.001541944450552547, 0.006480816113833214, 0.006480816113833214, 0.00029281525115992866
tbwebs, tbcaps, Abcaps, Abwebs = TASOPT.structures.size_wing_section!(wing.outboard, wing.sweep, sigfac)
@test test_tbwebs ≈ tbwebs
@test test_tbcaps ≈ tbcaps
@test test_Abcaps ≈ Abcaps
@test test_Abwebs ≈ Abwebs
#end size_wing_section!

#size_cap
sigcap = 2.0684848484848484e8
moment = 3.374158220198118e6
toc = 0.14
width_to_chord = 0.5
h_rms = 0.1287568768389997
cs = 6.12218798433346
cosL = 0.9063077870366499

test_tbcap = 0.0016511998748901524
test_Abcap = 0.0016511998748901524
tbcap, Abcap = TASOPT.structures.size_cap(sigcap, moment, wing.outboard.cross_section.thickness_to_chord,
  wing.outboard.cross_section.width_to_chord, h_rms, cs, cosL)

@test test_tbcap ≈ tbcap
@test test_Abcap ≈ Abcap
#end size_cap

# size_web
tauweb = 1.378913257881327e8
shear_load = 1.0413949053835782e6
cp = 5.548586639725954
web_height = 0.10500000000000001

test_tbweb = 0.0011681388369913896
test_Abweb = 0.00024530915576819183

tbweb, Abweb = TASOPT.structures.size_web(tauweb, shear_load, cp, web_height)

@test test_tbweb ≈ tbweb
@test test_Abweb ≈ Abweb
#end size_web

#calc_sparbox_internal_area
h_avgs = 0.12833333333333333
tbcaps = 0.0018311280070005542
tbwebs = 0.0011681388369913902
Abfuels = TASOPT.structures.calc_sparbox_internal_area(wing.inboard.cross_section.width_to_chord, h_avgs, tbcaps, tbwebs)
@test Abfuels ≈ 0.06204427240513358
#end calc_sparbox_internal_area

#calc_wing_weights:
  po = 114119.45308868506
  gammat = 0.225
  gammas = 0.8665999999999999
  Nlift = 3.0000000000000000
  Weng1 = 30083.769249783934
  fLt = -0.05
  sigfac = 1.0
  rhofuel = 817.00000000000000

  fort_Wwing = 128938.6987239125
  fort_Wsinn = 11628.586033197067
  fort_Wsout = 20426.086030894225
  fort_dyWsinn = 18223.86914912334
  fort_dyWsout = 96817.0136533246
  fort_Wfcen =  30442.713635875174
  fort_Wfinn = 38583.76809013594
  fort_Wfout = 46690.92699159976
  fort_dxWfinn = 29491.719272607046
  fort_dxWfout = 188722.8051493076
  fort_dyWfinn =  60466.98531940439
  fort_dyWfout = 221308.97271238975
  fort_lstrutp = 0.0


Wwing,Wsinn,Wsout,
        dyWsinn,dyWsout,
        Wfcen,Wfinn,Wfout,
        dxWfinn,dxWfout,
        dyWfinn,dyWfout,lstrutp = TASOPT.calc_wing_weights!(wing, po, gammat, gammas,
                                            Nlift, Weng1, 0, 0.0, 0, 0.0,
                                            sigfac, rhofuel)


@test fort_Wwing ≈ Wwing 
@test fort_Wsinn ≈ Wsinn 
@test fort_Wsout ≈ Wsout 
@test fort_dyWsinn ≈ dyWsinn 
@test fort_dyWsout ≈ dyWsout 
@test fort_Wfcen ≈ Wfcen 
@test fort_Wfinn ≈ Wfinn 
@test fort_Wfout ≈ Wfout 
@test fort_dxWfinn ≈ dxWfinn 
@test fort_dxWfout ≈ dxWfout 
@test fort_dyWfinn ≈ dyWfinn 
@test fort_dyWfout ≈ dyWfout 
@test fort_lstrutp ≈ lstrutp 
#end calc_wing_weights

#calc_wing_weights for Htail
poh = 115893.98734144184
λhs = 1.0
fLt = -0.05
tauwebh = 1.378913257881327e8
σcaph = 2.0684848484848484e8
surft_f_out = [14366.067634789782, 14032.558269851817, 0.0011577052661293624, 0.0023921269535798137, 1.8915676188667163e8, 1.258557904500963e9, 1.7895336288389182e8]

TASOPT.calc_wing_weights!(htail, poh, htail.outboard.λ, λhs,
0.0, 0.0, 0, 0.0, 0, 0.0,
sigfac, rhofuel)

surft_out = [htail.weight, htail.dxW, htail.outboard.webs.thickness, htail.outboard.caps.thickness, htail.outboard.EI[1], htail.outboard.EI[4], htail.outboard.GJ]

@test all(isapprox.(surft_out, surft_f_out))
#end surft


#tailpo:
  Sh = 42.443587259163984
  ARh = 6.0000000000000000
  lambdah = 0.25000000000000000
  qne = 12692.519555311032
  CLhmax = 2.0000000000000000
  fort_bh = 15.958117794871169
  fort_coh = 4.2554980786323124
  fort_poh = 108025.98516125829

poh,htail_span = TASOPT.aerodynamics.tailpo!(htail, Sh, qne)

@test fort_bh ≈ htail_span
@test fort_coh ≈ htail.layout.root_chord
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
htail.layout.box_x = 34.8996
vtail.layout.box_x = 33.528
wing.layout.box_x = 16.04432532088372
parg[igxeng] = wing.layout.box_x - parg[igdxeng2wbox]
fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 

fuselage_fueltank_count = 1

#Update fuel tank length and check changes
parg[iglftank] = 5.0
TASOPT.update_fuse!(fuselage, wing, htail, vtail, parg, fuselage_fueltank_count)

update_fuse_out = [fuselage.layout.x_end_cylinder, 
fuselage.layout.x_pressure_shell_aft, 
fuselage.layout.x_cone_end, 
fuselage.APU.x,
fuselage.layout.x_end,
fuselage.HPE_sys.x,
htail.layout.box_x,
vtail.layout.box_x,
parg[igxeng]]

update_fuse_out_test = [35.175200000000004, 36.699200000000005, 24.326234144000004, 42.18560000000001, 43.40480000000001, 12.767380728136962, 40.50920000000001, 39.137600000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))


#Return to original points?
fuselage_fueltank_count = 0
parg[iglftank] = 0.0
TASOPT.update_fuse!(fuselage, wing, htail, vtail, parg, fuselage_fueltank_count)

update_fuse_out = [fuselage.layout.x_end_cylinder, 
fuselage.layout.x_pressure_shell_aft, 
fuselage.layout.x_cone_end, 
fuselage.APU.x,
fuselage.layout.x_end,
fuselage.HPE_sys.x,
htail.layout.box_x,
vtail.layout.box_x,
parg[igxeng]]

update_fuse_out_test = [29.5656, 31.0896, 18.716634144, 36.57600000000001, 37.79520000000001, 9.82323826413696, 34.89960000000001, 33.528000000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))

#Test cabin resizing
parg = zeros(igtotal)
fuselage.layout.x_start_cylinder = 6.096
fuselage.layout.x_end_cylinder = 29.5656
fuselage.layout.x_pressure_shell_aft = 31.0896
fuselage.layout.l_cabin_cylinder = 23.4696
parg[igdxeng2wbox] = 1.5239999999999991
fuselage.APU.r = [36.576, 0.0, 0.0]
fuselage.layout.x_end = 37.7952
fuselage.layout.x_cone_end = 35.6616
htail.layout.box_x = 34.8996
vtail.layout.box_x = 33.528
wing.layout.box_x = 16.04432532088372
parg[igxeng] = wing.layout.box_x - parg[igdxeng2wbox]
fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 

fuselage.cabin.seat_pitch = 0.762
fuselage.cabin.seat_width = 0.4826
fuselage.cabin.aisle_halfwidth = 0.254
parg[igWpaymax] = 219964.5779
fuselage.layout.cross_section.radius = 2.5 #Change radius to 2.5 m

fuse_tank = TASOPT.fuselage_tank()
has_wing_fuel = false

TASOPT.update_fuse_for_pax!(has_wing_fuel, parg, fuselage, fuse_tank, wing, htail, vtail)

parg_check = [14.584924835954398, 219964.5779, 1.5239999999999991]

parg_nz = deepcopy(parg)
deleteat!(parg_nz, parg_nz .== 0)
for (i,item) in enumerate(parg_nz) #For every nonzero element in parg
  @test parg_nz[i] ≈ parg_check[i]
end

@test fuselage.layout.x_pressure_shell_aft ≈ 31.24200000000001
@test fuselage.layout.x_cone_end ≈ 18.86903414400001
@test fuselage.layout.x_end ≈ 37.94760000000001

#Test minimum radius calculation
Rmin = TASOPT.structures.find_minimum_radius_for_seats_per_row(5, ac_test)
@test Rmin ≈ 1.7113052179793784
end