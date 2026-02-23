@testset "structural components" begin

  # Fuselage
  include(joinpath(TASOPT.__TASOPTroot__, "../test/structures_params.jl"))
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

# --------------------------------------------------------------
# Structural sizing tests using the Wing and WingSection types
# --------------------------------------------------------------
# ------------------------------
# Wing
# ------------------------------
wing = ac_test.wing
wing.inboard.webs.weight = Weight(W = 3441.16315609358662186423) 
wing.outboard.webs.weight = Weight(W = 3441.16315609358662186423) 
wing.inboard.caps.weight = Weight(W = 74843.53209370737022254616) 
wing.outboard.caps.weight = Weight(W = 74843.53209370737022254616) 
wing.inboard.caps.material = TASOPT.materials.StructuralAlloy("TASOPT-Al",
        max_avg_stress = 1.1,
        safety_factor = 1.5)
wing.outboard.caps.material = TASOPT.materials.StructuralAlloy("TASOPT-Al",
        max_avg_stress = 1.1,
        safety_factor = 1.5)
wing.inboard.caps.material = TASOPT.materials.StructuralAlloy("TASOPT-Al",
        max_avg_stress = 1.1,
        safety_factor = 1.5)
wing.inboard.webs.material = TASOPT.materials.StructuralAlloy("TASOPT-Al",
        max_avg_stress = 1.1,
        safety_factor = 1.5)
wing.outboard.webs.material = TASOPT.materials.StructuralAlloy("TASOPT-Al",
        max_avg_stress = 1.1,
        safety_factor = 1.5)
wing.weight = 128386.90028286514279898256 
wing.strut.weight = 0.00000000000000000000 
wing.dxW = 299510.30434786697151139379 
wing.strut.dxW = 0.00000000000000000000 
wing.inboard.weight = 49780.59713745497720083222 
wing.outboard.weight = 70752.93408993091725278646 
wing.inboard.dyW = 78014.22168667706137057394 
wing.outboard.dyW = 335359.78341273276600986719 
wing.weight_frac_flap = 0.20000000000000001110 
wing.weight_frac_slat = 0.10000000000000000555 
wing.weight_frac_ailerons = 0.04000000000000000083 
wing.weight_frac_leading_trailing_edge = 0.10000000000000000555 
wing.weight_frac_ribs = 0.14999999999999999445 
wing.weight_frac_spoilers = 0.02000000000000000042 
wing.weight_frac_attachments = 0.02999999999999999889 
wing.strut.local_velocity_ratio = 0.00000000000000000000 
wing.layout.x = 18.47304210589511797025 
wing.layout.box_x = 15.69813728625496729308 
wing.layout.z = -1.67640000000000011227 
wing.strut.cos_lambda = 1.00000000000000000000 
wing.strut.S = 0.00000000000000000000 
wing.layout.spar_box_x_c = 0.40000000000000002220 
wing.inboard.cross_section.width_to_chord = 0.50000000000000000000 
wing.inboard.cross_section.web_to_box_height = 0.75000000000000000000 
wing.inboard.cross_section.thickness_to_chord = 0.12679999999999999605 
wing.outboard.cross_section.thickness_to_chord = 0.12659999999999999032 
wing.layout.max_span = 35.81400000000000005684 
wing.strut.thickness_to_chord = 0.00000000000000000000 
wing.strut.z = 0.00000000000000000000 
wing.outboard.moment = 3476534.48820107197389006615 
wing.outboard.max_shear_load = 621301.47081772319506853819 
wing.outboard.GJ = 252453751.39864748716354370117 
wing.outboard.EI[4] = 2502969300.33896875381469726562 
wing.outboard.EI[1] = 323733378.60267031192779541016 
wing.outboard.caps.thickness = 0.00648081611396861268 
wing.inboard.moment = 6074518.19739722646772861481 
wing.inboard.max_shear_load = 759371.76371910830494016409 
wing.inboard.GJ = 637210633.19837903976440429688 
wing.inboard.EI[4] = 5973671886.48499298095703125000 
wing.inboard.EI[1] = 809365031.46897768974304199219 
wing.inboard.caps.thickness = 0.00369099037929488533 
wing.inboard.webs.thickness = 0.00092200085359745961 
wing.outboard.webs.thickness = 0.00154194445053246480 
wing.layout.S = 139.59886994635940027365 
wing.layout.root_span = 3.60679999999999978400 
wing.layout.ηs = 0.28499999999999997558 
wing.inboard.λ = 0.69999999999999995559 
wing.outboard.λ = 0.25000000000000000000 
wing.layout.root_chord = 6.23487346145193299662 
wing.layout.span= 37.54928210310058034338 
wing.layout.sweep = 26.00000000000000000000 
wing.layout.AR = 10.09999999999999964473 
wing.fuse_lift_carryover = -0.29999999999999998890 
wing.tip_lift_loss = -0.05000000000000000278 
wing.inboard.co = 6.23487346145193299662 
wing.outboard.co = 4.36441142301635309764 
wing.mean_aero_chord = 4.25957824083664071679 
# ------------------------------
# Htail
# ------------------------------
htail = ac_test.htail
htail.weight = 14366.06764700812527735252 
htail.dxW = 14021.14191882828163215891 
htail.weight_fraction_added = 0.29999999999999998890 
htail.layout.box_x = 34.89959999999999951115 
htail.layout.z = 0.00000000000000000000 
htail.downwash_factor = 0.59999999999999997780 
htail.CL_max_fwd_CG = -0.11279688449577944531 
htail.CL_max = 2.00000000000000000000 
htail.SM_min = 0.05000000000000000278 
htail.layout.x = 36.21118914551904310883 
htail.outboard.cross_section.thickness_to_chord = 0.14000000000000001332 
htail.CL_CLmax = -0.50000000000000000000 
htail.opt_sizing = TASOPT.TailSizing.FixedVh
htail.volume = 1.44999999999999884537 
htail.outboard.GJ = 178813858.18902274966239929199 
htail.outboard.EI[4] = 1257539450.85079479217529296875 
htail.outboard.EI[1] = 189008140.11557617783546447754 
htail.layout.sweep = 26.00000000000000000000 
htail.layout.root_chord = 4.55405723428494280114 
htail.outboard.λ = 0.25000000000000000000 
htail.layout.root_span = 1.52400000000000002132 
htail.layout.span = 17.07771462856853617041 
htail.layout.AR = 6.00000000000000000000 
htail.layout.S = 48.60805615580396477071 
htail.outboard.cross_section.width_to_chord = 0.50000000000000000000 
htail.outboard.cross_section.web_to_box_height = 0.75000000000000000000 
htail.layout.ηs = htail.layout.root_span/htail.layout.span 
htail.strut.cos_lambda = 1.00000000000000000000 
htail.inboard.moment = 1738109.50381608842872083187 
htail.outboard.moment = 1738109.50381608842872083187 
htail.inboard.max_shear_load = 561216.24125756067223846912 
htail.outboard.max_shear_load = 561216.24125756067223846912 
htail.outboard.webs.thickness = 0.00115679564788597486 
htail.inboard.webs.weight.W = 1019.53692057405589821428 
htail.inboard.caps.weight.W = 10031.28435002535297826398 
htail.inboard.webs.thickness = 0.00115679564788597486 
htail.inboard.caps.thickness = 0.00239017403624703501 
htail.outboard.webs.thickness = 0.00115679564788597486 
htail.outboard.caps.thickness = 0.00239017403624703501 
htail.inboard.GJ = 178813858.18902274966239929199 
htail.outboard.co = htail.layout.root_chord*htail.inboard.λ 
htail.inboard.co = htail.layout.root_chord 
# ------------------------------
# Vtail
# ------------------------------
vtail = ac_test.vtail
vtail.weight = 9622.75029359092150116339 
vtail.dxW = 12001.92569704716697742697 
vtail.weight_fraction_added = 0.40000000000000002220 
vtail.layout.box_x = 33.52799999999999869260 
vtail.CL_max = 2.60000000000000008882 
vtail.layout.x = 35.04938902182657045614 
vtail.outboard.cross_section.thickness_to_chord = 0.14000000000000001332 
vtail.ntails = 1.00000000000000000000 
vtail.volume = 0.10000000000000000555 
vtail.outboard.GJ = 500649964.25236284732818603516 
vtail.outboard.EI[4] = 3484110272.10832309722900390625 
vtail.outboard.EI[1] = 495174831.98658269643783569336 
vtail.layout.sweep = 25.00000000000000000000 
vtail.layout.root_chord = 6.11742991142690240025 
vtail.outboard.λ = 0.29999999999999998890 
vtail.layout.span = 7.95265888485497285387 
vtail.layout.AR = 2.00000000000000000000 
vtail.layout.S = 31.62239166943137291810 
vtail.opt_sizing = TASOPT.TailSizing.FixedVv
vtail.dxW = 12001.92569704716697742697 
vtail.outboard.cross_section.width_to_chord = 0.50000000000000000000 
vtail.outboard.cross_section.web_to_box_height = 0.75000000000000000000 
vtail.layout.ηs = vtail.layout.root_span/vtail.layout.span 
vtail.strut.cos_lambda = 1.00000000000000000000 
vtail.inboard.moment = 3366297.29689422156661748886 
vtail.outboard.moment = 3366297.29689422156661748886 
vtail.inboard.max_shear_load = 1039776.82288908748887479305 
vtail.outboard.max_shear_load = 1039776.82288908748887479305 
vtail.outboard.webs.thickness = 0.00116813883699138938 
vtail.inboard.webs.weight.W = 1624.03782910894824453862 
vtail.inboard.caps.weight.W = 12122.74830915780330542475 
vtail.inboard.webs.thickness = 0.00116813883699138938 
vtail.inboard.caps.thickness = 0.00183112800700055423 
vtail.outboard.webs.thickness = 0.00116813883699138938 
vtail.outboard.caps.thickness = 0.00183112800700055423 
vtail.inboard.GJ = 500649964.25236284732818603516 
vtail.outboard.co = vtail.layout.root_chord*vtail.inboard.λ 
vtail.inboard.co = vtail.layout.root_chord 
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

#wing_weights:
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
        dyWfinn,dyWfout,lstrutp = TASOPT.wing_weights!(wing, po, gammat, gammas,
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
#end wing_weights

#wing_weights for Htail
poh = 115893.98734144184
λhs = 1.0
fLt = -0.05
tauwebh = 1.378913257881327e8
σcaph = 2.0684848484848484e8
surft_f_out = [14366.067634789782, 14032.558269851817, 0.0011577052661293624, 0.0023921269535798137, 1.8915676188667163e8, 1.258557904500963e9, 1.7895336288389182e8]

TASOPT.wing_weights!(htail, poh, htail.outboard.λ, λhs,
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

poh,htail_span = TASOPT.aerodynamics.tail_loading!(htail, Sh, qne)

@test fort_bh ≈ htail_span
@test fort_coh ≈ htail.layout.root_chord
@test fort_poh ≈ poh
#end tailpo:

# #Test fuselage layout updating for tank
ac = deepcopy(ac_test)

ac.fuselage.layout.x_start_cylinder = 6.096
ac.fuselage.layout.x_end_cylinder = 29.5656
ac.fuselage.layout.x_pressure_shell_aft = 31.0896
ac.parg[igdxeng2wbox] = 1.5239999999999991
ac.fuselage.APU.r = [36.576, 0.0, 0.0]
ac.fuselage.layout.x_end = 37.7952
ac.fuselage.layout.x_cone_end = 35.6616
ac.htail.layout.box_x = 34.8996
ac.vtail.layout.box_x = 33.528
ac.wing.layout.box_x = 16.04432532088372
ac.parg[igxeng] = ac.wing.layout.box_x - ac.parg[igdxeng2wbox]
ac.fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 
ac.fuselage.cabin.front_seat_offset = 10.0 * ft_to_m
ac.fuselage.cabin.rear_seat_offset = 0.0

ac.fuse_tank.tank_count = 1

#Update fuel tank length and check changes
ac.parg[iglftank] = 5.0
TASOPT.update_fuse!(ac, 1)

update_fuse_out = [ac.fuselage.layout.x_end_cylinder, 
ac.fuselage.layout.x_pressure_shell_aft, 
ac.fuselage.layout.x_cone_end, 
ac.fuselage.APU.x,
ac.fuselage.layout.x_end,
ac.fuselage.HPE_sys.x,
ac.htail.layout.box_x,
ac.vtail.layout.box_x,
ac.parg[igxeng]]

update_fuse_out_test = [35.175200000000004, 36.699200000000005, 24.326234144000004, 42.18560000000001, 43.40480000000001, 12.767380728136962, 40.50920000000001, 39.137600000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))


#Return to original points?
ac.fuse_tank.tank_count = 0
ac.parg[iglftank] = 0.0
TASOPT.update_fuse!(ac, 1)

update_fuse_out = [ac.fuselage.layout.x_end_cylinder, 
ac.fuselage.layout.x_pressure_shell_aft, 
ac.fuselage.layout.x_cone_end, 
ac.fuselage.APU.x,
ac.fuselage.layout.x_end,
ac.fuselage.HPE_sys.x,
ac.htail.layout.box_x,
ac.vtail.layout.box_x,
ac.parg[igxeng]]

update_fuse_out_test = [29.5656, 31.0896, 18.716634144, 36.57600000000001, 37.79520000000001, 9.82323826413696, 34.89960000000001, 33.528000000000006, 14.52032532088372]
@test all(isapprox.(update_fuse_out, update_fuse_out_test))

#Test cabin resizing
ac = TASOPT.load_default_model()

ac.parg = zeros(igtotal)
ac.fuselage.cabin.exit_limit = 189
ac.fuselage.layout.x_start_cylinder = 6.096
ac.fuselage.layout.x_end_cylinder = 29.5656
ac.fuselage.layout.x_pressure_shell_aft = 31.0896
ac.fuselage.layout.l_cabin_cylinder = 23.4696
ac.parg[igdxeng2wbox] = 1.5239999999999991
ac.fuselage.APU.r = [36.576, 0.0, 0.0]
ac.fuselage.layout.x_end = 37.7952
ac.fuselage.layout.x_cone_end = 35.6616
ac.htail.layout.box_x = 34.8996
ac.vtail.layout.box_x = 33.528
ac.wing.layout.box_x = 16.04432532088372
ac.parg[igxeng] = wing.layout.box_x - ac.parg[igdxeng2wbox]
ac.fuselage.layout.x_cone_end = fuselage.layout.x_cone_end * 0.52484 
ac.fuselage.cabin.front_seat_offset = 10.0 * ft_to_m
ac.fuselage.cabin.rear_seat_offset = 0.0

ac.fuselage.cabin.seat_pitch = 0.762
ac.fuselage.cabin.seat_width = 0.4826
ac.fuselage.cabin.aisle_halfwidth = 0.254
ac.parg[igWpaymax] = 219964.5779
ac.fuselage.layout.cross_section.radius = 2.5 #Change radius to 2.5 m

ac.options.is_doubledecker = false

TASOPT.update_fuse_for_pax!(ac)

parg_check = [14.584924835954398, 219964.5779, 1.5239999999999991]

parg_nz = deepcopy(ac.parg)
deleteat!(parg_nz, parg_nz .== 0)
for (i,item) in enumerate(parg_nz) #For every nonzero element in parg
  @test parg_nz[i] ≈ parg_check[i]
end

@test ac.fuselage.layout.x_pressure_shell_aft ≈ 31.24200000000001
@test ac.fuselage.layout.x_cone_end ≈ 18.86903414400001
@test ac.fuselage.layout.x_end ≈ 37.94760000000001

#Test minimum radius calculation
Rmin = TASOPT.structures.find_minimum_radius_for_seats_per_row(5, ac_test)
@test Rmin ≈ 1.7113052179793784

#Test landing gear sizing
ac = load_default_model()

#Test simple sizing based on mass fractions
WMTO = 80e4
ac.parg[igWMTO] = WMTO
ac.landing_gear.main_gear.overall_mass_fraction = 0.04
ac.landing_gear.nose_gear.overall_mass_fraction = 0.01

TASOPT.size_landing_gear!(ac)
@test ac.landing_gear.main_gear.weight.W ≈ 0.04 * 80e4
@test ac.landing_gear.nose_gear.weight.W ≈ 0.01 * 80e4

#Test models based on historical data
ac.landing_gear.model = "historical_correlations"
ac.landing_gear.tailstrike_angle = 10*pi/180
ac.landing_gear.wing_dihedral_angle = 6*pi/180
ac.landing_gear.engine_ground_clearance = 20*in_to_m
ac.landing_gear.nose_gear.number_struts = 1
ac.landing_gear.nose_gear.wheels_per_strut = 2
ac.landing_gear.main_gear.number_struts = 2
ac.landing_gear.main_gear.wheels_per_strut = 2
ac.landing_gear.main_gear.y_offset_halfspan_fraction = 0.2

ac.parg[igxCGaft] = 30
ac.pare[ieu0,iprotate,1] = 70
ac.parg[igdfan] = 1.5
ac.wing.layout.span = 40

ac.landing_gear.main_gear.weight = TASOPT.Weight(y = ac.landing_gear.main_gear.y_offset_halfspan_fraction * ac.wing.layout.span / 2)

TASOPT.size_landing_gear!(ac)

@test ac.landing_gear.main_gear.weight.W ≈ 31787.399917196755
@test ac.landing_gear.nose_gear.weight.W ≈ 4854.384534273811
end