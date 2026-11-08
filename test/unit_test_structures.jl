using Random

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

  test_out = [0.0010592916489634418, 0.0012236016395195024, 0.0, 0.0009061565981319232, 18.465770971525686, 22.676674681482936, 2.361320105799576e9, 2.4685708558261826e10, 2.317349483950291e9, 2.4685708558261826e10, 1.5708404995400493e9, 1.6421853018493836e9, 20947.352673154768, 5667.810474837753, 11269.980000000001, 5018.7283576257105, 13415.796495703155, 14314.082074792976, 7380.9251468929015, 198042.1615045073, 3.6202612837849665e6, 414.54271063740043]

# @test all(isapprox.(out, test_out))
for i in 1:length(out)
  @test out[i] ≈ test_out[i]
end

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
wing.strut.area = 0.00000000000000000000 
wing.strut.axial_force = 0.00000000000000000000 
wing.strut.length = 0.00000000000000000000 
wing.strut.chord = 0.00000000000000000000 
wing.strut.box_width_to_chord = 0.50000000000000000000 
wing.strut.web_to_box_height = 0.75000000000000000000 
wing.strut.end_fixity_K = 1.00000000000000000000 
wing.strut.dihedral = 0.00000000000000000000 
wing.strut.axial_force_compression = 0.00000000000000000000 
wing.strut.buckling_load = 0.00000000000000000000 
wing.strut.EI_min = 0.00000000000000000000 
wing.strut.CL = 0.00000000000000000000 
wing.strut.lift_fraction = 0.00000000000000000000 
wing.strut.weight_frac_added = 0.59999999999999997780 
wing.strut.sizing_criterion = :none
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
sd = TASOPT.structures.size_wing_section!(wing.outboard, wing.sweep, sigfac)
@test sd isa TASOPT.structures.WingSectionDimensions
@test sd.tbweb ≈ 0.001541944450552547
@test sd.tbcap ≈ 0.006480816113833214
@test sd.Abweb ≈ 0.00029281525115992866
@test sd.Abcap ≈ 0.006480816113833214
@test_throws MethodError begin; a, b, c, d = sd; end
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

# The strut branch is selected by the has_strut keyword (cantilever by default)
wing_cantilever = deepcopy(wing)
wing_braced = deepcopy(wing)
wing_braced.strut.z = 154 * in_to_m
wing_braced.strut.thickness_to_chord = 0.1 # the dumped strut has no box height
TASOPT.wing_weights!(wing_cantilever, po, gammat, gammas, Nlift, Weng1, 0, 0.0, 0, 0.0, sigfac, rhofuel)
TASOPT.wing_weights!(wing_braced, po, gammat, gammas, Nlift, Weng1, 0, 0.0, 0, 0.0, sigfac, rhofuel; has_strut = true)
@test wing_cantilever.strut.cos_lambda == 1.0
@test wing_braced.strut.cos_lambda < 1.0
@test wing_braced.strut.weight > 0
@test wing_braced.inboard.moment == wing_braced.outboard.moment

# Engine at the span break ηs: root shear and moment should match
# TASOPT 2.16 surfw.f (Drela eqs. 177-178 with R = N*We), i.e. moment arm (b/2)(ηs-ηo)
wing_noeng = deepcopy(wing)
wing_eng = deepcopy(wing)
b, etas, etao = wing.layout.span, wing.layout.ηs, wing.layout.ηo
Winn, dyWinn = wing.inboard.weight, wing.inboard.dyW
dyeinn = 0.5*b*(etas - etao)
TASOPT.wing_weights!(wing_noeng, po, gammat, gammas, Nlift, Weng1, 0, 0.0, 0, 0.0, sigfac, rhofuel)
TASOPT.wing_weights!(wing_eng, po, gammat, gammas, Nlift, Weng1, 0, 0.0, 1, dyeinn, sigfac, rhofuel)
Ss, Ms = wing_eng.outboard.max_shear_load, wing_eng.outboard.moment
fort_So = Ss - Nlift*Weng1 + 0.25*po*b*(1.0+gammas)*(etas-etao) - Nlift*Winn
fort_Mo = Ms + (Ss - Nlift*Weng1)*0.5*b*(etas-etao) +
          (1.0/24.0)*po*b^2*(1.0+2.0*gammas)*(etas-etao)^2 - Nlift*dyWinn
@test fort_Mo > Ms # Mo = max(Mo, Ms) clamp is inactive
@test wing_eng.inboard.max_shear_load ≈ fort_So
@test wing_eng.inboard.moment ≈ fort_Mo
@test wing_noeng.inboard.moment - wing_eng.inboard.moment ≈ Nlift*Weng1*dyeinn

@testset "wing_weights! panel inertia update (design_wing_fuel)" begin
    args = (po, gammat, gammas, Nlift, Weng1, 0, 0.0, 1, dyeinn, sigfac, rhofuel)
    fwadd = TASOPT.wing_additional_weight(wing)

    # Default (no design_wing_fuel): panel fields keep their input values
    w_ref = deepcopy(wing)
    ref = TASOPT.wing_weights!(w_ref, args...)
    @test w_ref.inboard.weight == wing.inboard.weight
    @test w_ref.outboard.weight == wing.outboard.weight
    @test w_ref.inboard.dyW == wing.inboard.dyW
    @test w_ref.outboard.dyW == wing.outboard.dyW
    Wsinn_r, Wsout_r, dyWsinn_r, dyWsout_r, Wfcen_r, Wfinn_r, Wfout_r, _, _, dyWfinn_r, dyWfout_r = ref[2:12]
    design_wing_fuel = 0.6 * 2.0 * (Wfcen_r + Wfinn_r + Wfout_r) # N

    panels = Dict{Bool, NTuple{4, Float64}}()
    for has_centerbox_fuel in (true, false)
        w_new = deepcopy(wing)
        out = TASOPT.wing_weights!(w_new, args...; design_wing_fuel = design_wing_fuel,
                                   has_wing_fuel = true, has_centerbox_fuel = has_centerbox_fuel)
        @test out == ref # return values do not depend on the new keywords

        Wfmax = 2.0 * ((has_centerbox_fuel ? Wfcen_r : 0.0) + Wfinn_r + Wfout_r)
        rfmax = design_wing_fuel / Wfmax
        @test w_new.inboard.weight  == Wsinn_r * (1.0 + fwadd) + rfmax * Wfinn_r
        @test w_new.outboard.weight == Wsout_r * (1.0 + fwadd) + rfmax * Wfout_r
        @test w_new.inboard.dyW     == dyWsinn_r * (1.0 + fwadd) + rfmax * dyWfinn_r
        @test w_new.outboard.dyW    == dyWsout_r * (1.0 + fwadd) + rfmax * dyWfout_r
        panels[has_centerbox_fuel] = (w_new.inboard.weight, w_new.outboard.weight,
                                      w_new.inboard.dyW, w_new.outboard.dyW)
    end
    # Without centerbox fuel the same fuel mass fills a smaller tank, so rfmax is larger
    @test panels[true] != panels[false]
    @test panels[false][1] > panels[true][1]

    # No wing fuel: rfmax = 0, panels carry only structure weight
    w_dry = deepcopy(wing)
    out_dry = TASOPT.wing_weights!(w_dry, args...; design_wing_fuel = design_wing_fuel, has_wing_fuel = false)
    @test out_dry == ref
    @test w_dry.inboard.weight  == Wsinn_r * (1.0 + fwadd)
    @test w_dry.outboard.weight == Wsout_r * (1.0 + fwadd)
    @test w_dry.inboard.dyW     == dyWsinn_r * (1.0 + fwadd)
    @test w_dry.outboard.dyW    == dyWsout_r * (1.0 + fwadd)
end

@testset "wing_case_loads matches the section loads of wing_weights!" begin
    dyeout = 0.25 * b * (1.0 - etas)
    configs = (
        "no engines"                      => (We = Weng1,      neout = 0, dyeout = 0.0,    neinn = 0, dyeinn = 0.0),
        "inboard engine at ηs"            => (We = Weng1,      neout = 0, dyeout = 0.0,    neinn = 1, dyeinn = dyeinn),
        "outboard and heavy inboard engine" => (We = 10*Weng1, neout = 1, dyeout = dyeout, neinn = 1, dyeinn = dyeinn),
    )
    for (label, c) in configs
        @testset "$label" begin
            w_sized = deepcopy(wing)
            w_loads = deepcopy(wing)
            TASOPT.wing_weights!(w_sized, po, gammat, gammas, Nlift, c.We, c.neout, c.dyeout,
                                 c.neinn, c.dyeinn, sigfac, rhofuel)
            loads = TASOPT.wing_case_loads(TASOPT.structures.Cantilever(), w_loads, po, gammat, gammas,
                                           Nlift, c.We, c.neout, c.dyeout, c.neinn, c.dyeinn)
            @test loads isa TASOPT.structures.WingCaseLoads
            @test loads.S_root  == w_sized.inboard.max_shear_load
            @test loads.M_root  == w_sized.inboard.moment
            @test loads.S_break == w_sized.outboard.max_shear_load
            @test loads.M_break == w_sized.outboard.moment
        end
    end

    # The heavy inboard engine drives the inboard loads onto the outboard ones (max(So, outboard) clamp)
    w_heavy = deepcopy(wing)
    heavy = TASOPT.wing_case_loads(TASOPT.structures.Cantilever(), w_heavy, po, gammat, gammas,
                                   Nlift, 10*Weng1, 1, dyeout, 1, dyeinn)
    @test heavy.S_root == heavy.S_break
    @test heavy.M_root == heavy.M_break
end

@testset "wing_weights! with a section loads envelope" begin
    args = (po, gammat, gammas, Nlift, Weng1, 0, 0.0, 1, dyeinn, sigfac, rhofuel)
    w_ref = deepcopy(wing)
    ref = TASOPT.wing_weights!(w_ref, args...)
    section_loads(w) = (S_root = w.inboard.max_shear_load, M_root = w.inboard.moment,
                        S_break = w.outboard.max_shear_load, M_break = w.outboard.moment)
    computed = section_loads(w_ref)
    thicknesses(w) = (w.inboard.caps.thickness, w.outboard.caps.thickness,
                      w.inboard.webs.thickness, w.outboard.webs.thickness)

    # No envelope leaves the section loads and the returns unchanged
    w_none = deepcopy(wing)
    @test TASOPT.wing_weights!(w_none, args...; section_loads_envelope = nothing) == ref
    @test section_loads(w_none) == computed
    @test thicknesses(w_none) == thicknesses(w_ref)

    # An envelope equal to the computed loads changes nothing
    w_same = deepcopy(wing)
    @test TASOPT.wing_weights!(w_same, args...; section_loads_envelope = computed) == ref
    @test section_loads(w_same) == computed
    @test thicknesses(w_same) == thicknesses(w_ref)

    # Larger envelope loads replace the section loads and thicken the caps and webs
    larger = map(x -> 1.5 * x, computed)
    w_larger = deepcopy(wing)
    TASOPT.wing_weights!(w_larger, args...; section_loads_envelope = larger)
    @test section_loads(w_larger) == larger
    @test all(thicknesses(w_larger) .> thicknesses(w_ref))
end
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
ac.wing.η_engine = ac.wing.layout.ηs #Engine at the span break

ac.landing_gear.main_gear.weight = TASOPT.Weight(y = ac.landing_gear.main_gear.y_offset_halfspan_fraction * ac.wing.layout.span / 2)

TASOPT.size_landing_gear!(ac)

@test ac.landing_gear.main_gear.weight.W ≈ 31787.399917196755
@test ac.landing_gear.nose_gear.weight.W ≈ 4854.384534273811

@testset "landing gear length follows the engine spanwise location" begin
    gear_ac = deepcopy(ac)
    lg = gear_ac.landing_gear
    half_span = gear_ac.wing.layout.span / 2
    ηs = gear_ac.wing.layout.ηs
    η_outboard = 0.6

    nose_length(η_engine) = (gear_ac.wing.η_engine = η_engine; TASOPT.size_landing_gear!(gear_ac); lg.nose_gear.length)

    l_inboard = nose_length(ηs)
    l_outboard = nose_length(η_outboard)

    # The engine clearance, not the tailstrike, sets the nose gear length
    x_lg = gear_ac.parg[igxCGaft] + lg.main_gear.distance_CG_to_landing_gear
    l_tailstrike = (gear_ac.fuselage.layout.x_end - x_lg) * tan(lg.tailstrike_angle) - 2 * gear_ac.fuselage.layout.radius
    l_clearance = lg.engine_ground_clearance + gear_ac.parg[igdfan] - ηs * half_span * tan(lg.wing_dihedral_angle)
    @test l_clearance > l_tailstrike
    @test l_inboard ≈ l_clearance

    # Moving the engine outboard shortens the gear by the dihedral rise of the wing
    @test l_inboard - l_outboard ≈ (η_outboard - ηs) * half_span * tan(lg.wing_dihedral_angle)
    @test lg.main_gear.length - l_outboard ≈ lg.main_gear.weight.y * tan(lg.wing_dihedral_angle)

    # Without dihedral the engine location does not matter
    lg.wing_dihedral_angle = 0.0
    @test nose_length(ηs) == nose_length(η_outboard)
end

@testset "wing_engine_arms" begin
    local wing = deepcopy(ac.wing) # local: the tests below use the global wing
    wing.layout.span = 40.0
    wing.layout.root_span = 4.0 # ηo = 0.1
    wing.layout.ηs = 0.3

    # Engine inboard of the span break
    wing.η_engine = 0.2
    arms = TASOPT.structures.wing_engine_arms(wing)
    @test (arms.neinn, arms.neout, arms.dyeout) == (1, 0, 0.0)
    @test arms.dyeinn ≈ 2.0

    # Engine at the span break: the arms are the inboard engine arms of the span break
    wing.η_engine = wing.layout.ηs
    arms = TASOPT.structures.wing_engine_arms(wing)
    @test (arms.neinn, arms.neout, arms.dyeout) == (1, 0, 0.0)
    @test arms.dyeinn == 0.5 * wing.layout.span * (wing.layout.ηs - wing.layout.ηo)
    @test arms.dyeinn ≈ 4.0

    # Engine outboard of the span break
    wing.η_engine = 0.8
    arms = TASOPT.structures.wing_engine_arms(wing)
    @test (arms.neinn, arms.dyeinn, arms.neout) == (0, 0.0, 1)
    @test arms.dyeout ≈ 10.0

    # An engine inboard of the wing root station has no inboard arm
    wing.η_engine = 0.05
    @test TASOPT.structures.wing_engine_arms(wing).dyeinn == 0.0
end

@testset "size_wing_section! property tests" begin
    # Re-use wing.outboard from the setup above which already has max_shear_load
    # and moment set. size_wing_section! only side-effects EI and GJ, so calling it
    # multiple times with different sigfac is safe (loads don't change).
    rng = MersenneTwister(42)

    @testset "Web shear: t_web·c_perp²·sigfac is invariant" begin
        # t_web = V / (2·h_web·τmax·c_perp²), so
        # t_web · c_perp² · sigfac = V / (2·h_web·τmax_material) = const.
        for _ in 1:10
            sf1 = 0.5 + 1.5 * rand(rng)
            sf2 = 0.5 + 1.5 * rand(rng)
            sw1 = 5.0 + 30.0 * rand(rng)
            sw2 = 5.0 + 30.0 * rand(rng)
            s1 = TASOPT.structures.size_wing_section!(wing.outboard, sw1, sf1)
            s2 = TASOPT.structures.size_wing_section!(wing.outboard, sw2, sf2)
            @test s1.tbweb * cosd(sw1)^2 * sf1 ≈ s2.tbweb * cosd(sw2)^2 * sf2 rtol = 1e-12
        end
    end

    @testset "Cap bending: t_cap is monotone increasing with sweep angle" begin
        # con (in calc_cap_thickness) ∝ M / (σmax · cosΛ⁴), so as sweep increases cosΛ shrinks
        # and con grows. So the cap must be thicker to resist the same bending moment
        # on the shorter effective structural chord. A more swept wing always needs
        # a heavier spar cap for equal loads.
        for _ in 1:10
            sw_lo = 5.0  + 15.0 * rand(rng)
            sw_hi = sw_lo + 5.0 + 10.0 * rand(rng)
            sf    = 0.5  +  1.5 * rand(rng)
            s_lo = TASOPT.structures.size_wing_section!(wing.outboard, sw_lo, sf)
            s_hi = TASOPT.structures.size_wing_section!(wing.outboard, sw_hi, sf)
            @test s_hi.tbcap > s_lo.tbcap
            @test s_hi.Abcap > s_lo.Abcap
        end
    end
end

@testset "strut_force_compatibility hand calculations" begin
    strut_force_compatibility = TASOPT.structures.strut_force_compatibility
    L, EI0 = 6.0, 2.0e8   # bay length [m], constant bending stiffness [N m^2]
    EI_const(x) = EI0
    f_ww = L^3 / (3EI0)   # flexibility of the clamped bay at x = L for a unit vertical load

    @testset "tip load: the rigid strut takes the whole load" begin
        Q = 3.0e5
        M0_tip(x) = Q * (L - x)
        @test strut_force_compatibility(M0_tip, EI_const, L, 0.0) ≈ Q rtol = 1e-12
        # M0 is proportional to (L - x), so the result does not depend on the stiffness distribution
        @test strut_force_compatibility(M0_tip, x -> EI0 * (1 + x / L), L, 0.0) ≈ Q rtol = 1e-12
    end

    @testset "uniform load: R = 3 w L / 8" begin
        w = 5.0e4
        M0_uniform(x) = w * (L - x)^2 / 2
        @test strut_force_compatibility(M0_uniform, EI_const, L, 0.0; npts = 201) ≈ 3w * L / 8 rtol = 1e-4
        @test strut_force_compatibility(M0_uniform, EI_const, L, 0.0) ≈ 3w * L / 8 rtol = 2e-3
    end

    @testset "strut flexibility lowers the strut reaction" begin
        Q = 3.0e5
        M0_tip(x) = Q * (L - x)
        R(f_strut) = strut_force_compatibility(M0_tip, EI_const, L, f_strut)
        # R = Q f_ww / (f_ww + f_strut)
        @test R(f_ww) ≈ Q / 2 rtol = 2e-3
        @test R(1e6 * f_ww) < 1e-5 * Q
        @test R(1e-6 * f_ww) ≈ Q rtol = 1e-5
        @test R(0.1 * f_ww) > R(f_ww) > R(10 * f_ww)
        # a downward load gives a downward reaction
        @test strut_force_compatibility(x -> -Q * (L - x), EI_const, L, f_ww) ≈ -R(f_ww) rtol = 1e-12
    end
end

@testset "size_wing_section! axial force in the caps" begin
    st = TASOPT.structures
    section = wing.outboard
    sweep, sigfac = wing.sweep, 1.0
    dims(d) = (d.tbweb, d.tbcap, d.Abweb, d.Abcap)

    @testset "axial force is the equivalent moment M + |P| h_rms cp / 2" begin
        d0 = st.size_wing_section!(section, sweep, sigfac)
        h_avg, h_rms = st.get_average_sparbox_heights(section.cross_section)
        cp = section.co * cosd(sweep)
        P = 0.2 * section.moment / (h_rms * cp / 2)   # the axial force adds 20 % to the moment

        for P_axial in (P, -P)
            d_axial = st.size_wing_section!(section, sweep, sigfac; P_axial)
            @test d_axial.tbcap > d0.tbcap
            @test d_axial.tbweb == d0.tbweb

            equivalent = deepcopy(section)
            equivalent.moment = section.moment + abs(P_axial) * h_rms * cp / 2
            d_equivalent = st.size_wing_section!(equivalent, sweep, sigfac)
            @test d_axial.tbcap ≈ d_equivalent.tbcap rtol = 1e-12
            @test d_axial.Abcap ≈ d_equivalent.Abcap rtol = 1e-12
            @test d_axial.tbweb == d_equivalent.tbweb
        end
    end
end

@testset "fusew! property-based tests" begin
    # Helper: run fusew! on a fresh deep copy of fuselage, overriding specific load params.
    # All geometry/material fields come from the already-run `fuselage` base state all the way on the top.
    function call_fusew(; Lhmax_=Lhmax, Lvmax_=Lvmax,
                          Whtail_=Whtail, Wvtail_=Wvtail,
                          xwbox_=xwbox)
        fuse_ = deepcopy(fuselage)
        TASOPT.fusew!(fuse_, Nland, Wpaymax, Wengtail,
                      nftanks, Waftfuel, Wftank_single, ltank, xftank_fuse, tank_placement,
                      Δp, Whtail_, Wvtail_, rMh, rMv, Lhmax_, Lvmax_,
                      bv, λv, nvtail, xhtail, xvtail, xwing, xwbox_, cbox, xeng)
        return fuse_
    end

    rng = MersenneTwister(42)

    @testset "Larger vertical tail load requires more vertical bending material" begin
        # The vertical bending material is sized by the vertical tail torsion Qv ∝ Lvmax.
        for _ in 1:10
            Lv_lo = Lvmax * (0.3 + rand(rng))          # ∈ [0.3, 1.3] × baseline
            Lv_hi = Lv_lo * (1.2 + rand(rng))          # ∈ [1.2, 2.2] × Lv_lo
            f_lo = call_fusew(Lvmax_=Lv_lo)
            f_hi = call_fusew(Lvmax_=Lv_hi)
            @test f_hi.bendingmaterial_v.weight.W ≥ f_lo.bendingmaterial_v.weight.W
        end
    end

    @testset "Larger horizontal tail load requires more horizontal bending material" begin
        for _ in 1:10
            Lh_lo = Lhmax + 1e6 * rand(rng)
            Lh_hi = Lh_lo * (1.5 + rand(rng))
            f_lo = call_fusew(Lhmax_=Lh_lo)
            f_hi = call_fusew(Lhmax_=Lh_hi)
            @test f_hi.bendingmaterial_h.weight.W ≥ f_lo.bendingmaterial_h.weight.W
        end
    end

    @testset "Bending material never reduces fuselage stiffness" begin
        # Physics: the sizing routine can only add material to the shell; it never removes it.
        # Therefore EI_total ≥ EI_shell for both bending axes across all load combinations.
        for _ in 1:10
            Lh = 1e5 + 2e6 * rand(rng)
            Lv = 1e5 + 2e6 * rand(rng)
            f = call_fusew(Lhmax_=Lh, Lvmax_=Lv)
            @test f.shell.EIh + f.bendingmaterial_h.EIh ≥ f.shell.EIh
            @test f.shell.EIv + f.bendingmaterial_h.EIv ≥ f.shell.EIv
        end
    end

    @testset "Wing-box symmetry preserves horizontal bending weight and stiffness" begin
        for _ in 1:10
            δ = 0.5 + 3.0 * rand(rng)
            xwbox_fwd = xwing - δ  # dxwing = +δ
            xwbox_aft = xwing + δ  # dxwing = -δ  (mirror image)
            f_fwd = call_fusew(xwbox_=xwbox_fwd)
            f_aft = call_fusew(xwbox_=xwbox_aft)
            @test f_fwd.bendingmaterial_h.weight.W ≈ f_aft.bendingmaterial_h.weight.W rtol=1e-12
            @test f_fwd.bendingmaterial_h.EIh       ≈ f_aft.bendingmaterial_h.EIh       rtol=1e-12
        end
    end
end

@testset "fusew! strut-attach reinforcement is the wing-style fraction of the strut weight" begin
    # Strut weight 40 kN, wing-style fraction 0.03: the reinforcement is 0.03 * 40e3 = 1200 N at the wing-box x.
    # The reinforcement only adds to the fuselage weight and moment; it does not size the other components.
    function call_fusew_with_attach(; f_strut_weight, strut_weight)
        fuse_ = deepcopy(fuselage)
        TASOPT.fusew!(fuse_, Nland, Wpaymax, Wengtail,
                      nftanks, Waftfuel, Wftank_single, ltank, xftank_fuse, tank_placement,
                      Δp, Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
                      bv, λv, nvtail, xhtail, xvtail, xwing, xwbox, cbox, xeng;
                      strut_weight = strut_weight,
                      strut_weight_frac_attachments = f_strut_weight)
        return fuse_
    end
    attach(; kwargs...) = call_fusew_with_attach(; f_strut_weight = 0.03, strut_weight = 40e3, kwargs...)
    base = attach(strut_weight = 0.0)

    @testset "no strut weight (cantilever): no reinforcement, whatever the fraction" begin
        @test base.strut_attach.W == 0
        @test base.strut_attach_criterion === :none
        # the weight and moment are those of a fuselage without the keywords
        plain = deepcopy(fuselage)
        TASOPT.fusew!(plain, Nland, Wpaymax, Wengtail,
                      nftanks, Waftfuel, Wftank_single, ltank, xftank_fuse, tank_placement,
                      Δp, Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
                      bv, λv, nvtail, xhtail, xvtail, xwing, xwbox, cbox, xeng)
        @test plain.strut_attach.W == 0
        @test plain.strut_attach_criterion === :none
        @test base.weight == plain.weight
        @test base.moment == plain.moment
    end

    @testset "the reinforcement is the wing-style fraction times the strut weight, at the wing-box x" begin
        fuse_ = attach()
        @test fuse_.strut_attach.W ≈ 0.03 * 40e3
        @test fuse_.strut_attach.x == xwbox
        @test fuse_.strut_attach_criterion === :strut_weight

        # proportional to both the fraction and the strut weight
        @test attach(strut_weight = 20e3).strut_attach.W ≈ 0.03 * 20e3
        @test attach(f_strut_weight = 0.06).strut_attach.W ≈ 0.06 * 40e3
        @test attach(strut_weight = 20e3).strut_attach_criterion === :strut_weight
    end

    @testset "a zero wing-style fraction gives no reinforcement" begin
        fuse_ = attach(f_strut_weight = 0.0)
        @test fuse_.strut_attach.W == 0
        @test fuse_.strut_attach_criterion === :none
        @test fuse_.weight == base.weight
        @test fuse_.moment == base.moment
    end

    @testset "the reinforcement adds its weight and its moment about the wing-box x to the fuselage" begin
        for fuse_ in (attach(), attach(strut_weight = 10e3))
            W_attach = fuse_.strut_attach.W
            @test W_attach > 0
            @test fuse_.weight ≈ base.weight + W_attach rtol = 1e-12
            @test fuse_.moment ≈ base.moment + W_attach * xwbox rtol = 1e-12
            # the other components do not depend on the reinforcement
            @test fuse_.shell.weight.W == base.shell.weight.W
            @test fuse_.bendingmaterial_h.weight.W == base.bendingmaterial_h.weight.W
        end
    end
end

@testset "Phase 0 robustness: engine warning, tail guard, cap guard, minimum gauge" begin
    st = TASOPT.structures

    @testset "check_engine_convergence_failure warns instead of throwing" begin
        pare_g = zeros(ietotal, 3)
        @test_logs TASOPT.engine.check_engine_convergence_failure(pare_g) # converged: no log
        pare_g[ieConvFail, 2] = 1.0
        @test_logs (:warn, r"failed to converge") TASOPT.engine.check_engine_convergence_failure(pare_g)
    end

    @testset "tail_loading! rejects negative area" begin
        htail_g = deepcopy(ac_test.htail)
        po_g, b_g = TASOPT.aerodynamics.tail_loading!(htail_g, 42.443587259163984, 12692.519555311032)
        @test isfinite(po_g) && isfinite(b_g) && b_g > 0
        @test_throws "negative" TASOPT.aerodynamics.tail_loading!(htail_g, -1.0, 12692.519555311032)
    end

    # Cap/web inputs (same values as the size_cap / size_web tests above)
    σ_g, hb_g, wb_g, c_g, cosL_g = 2.0684848484848484e8, 0.14, 0.5, 6.12218798433346, 0.9063077870366499
    hrms_g = 0.12
    # moment that gives a chosen con = moment*6h̄/w̄/(c^3 σmax cosΛ^4)
    moment_for(con) = con * wb_g * c_g^3 * σ_g * cosL_g^4 / (6hb_g)
    τ_g, shear_g, cperp_g, hweb_g = 1.378913257881327e8, 1.0413949053835782e6, 5.548586639725954, 0.105

    @testset "calc_cap_thickness feasibility guard" begin
        con = 0.5 * hrms_g^3
        @test st.calc_cap_thickness(σ_g, moment_for(con), hb_g, wb_g, hrms_g, c_g, cosL_g) ≈
              0.5 * (hrms_g - cbrt(hrms_g^3 - con))

        # moment_for round-trips con with ~1e-16 relative error, so con == h_rms^3 exactly
        # cannot be reproduced through the moment; test just above it instead.
        for con in (2e-3, hrms_g^3 * (1 + 1e-12))
            m = moment_for(con)
            @test_throws r"infeasible" st.calc_cap_thickness(σ_g, m, hb_g, wb_g, hrms_g, c_g, cosL_g)
            @test_throws r"infeasible" st.size_cap(σ_g, m, hb_g, wb_g, hrms_g, c_g, cosL_g)
        end
    end

    @testset "size_cap / size_web minimum gauge" begin
        m = moment_for(0.5 * hrms_g^3)
        cap0 = st.size_cap(σ_g, m, hb_g, wb_g, hrms_g, c_g, cosL_g)
        web0 = st.size_web(τ_g, shear_g, cperp_g, hweb_g)

        # t_min = 0 reproduces the no-keyword call exactly
        @test st.size_cap(σ_g, m, hb_g, wb_g, hrms_g, c_g, cosL_g; t_min = 0.0) == cap0
        @test st.size_web(τ_g, shear_g, cperp_g, hweb_g; t_min = 0.0) == web0

        # t_min above the required thickness governs, and Ab follows the returned t
        t_hi = 2 * max(cap0[1], web0[1])
        t, Ab = st.size_cap(σ_g, m, hb_g, wb_g, hrms_g, c_g, cosL_g; t_min = t_hi)
        @test t == t_hi
        @test Ab ≈ 2 * t_hi * wb_g
        t, Ab = st.size_web(τ_g, shear_g, cperp_g, hweb_g; t_min = t_hi)
        @test t == t_hi
        @test Ab ≈ 2 * t_hi * hweb_g

        # t_min below the required thickness leaves the result unchanged
        @test st.size_cap(σ_g, m, hb_g, wb_g, hrms_g, c_g, cosL_g; t_min = 0.5 * cap0[1]) == cap0
        @test st.size_web(τ_g, shear_g, cperp_g, hweb_g; t_min = 0.5 * web0[1]) == web0
    end

    @testset "size_wing_section! defaults unchanged, t_cap_min passes through" begin
        # uses wing.outboard from the setup above (loads already set)
        sd0 = st.size_wing_section!(wing.outboard, wing.sweep, 1.0)
        EI0, GJ0 = copy(wing.outboard.EI), wing.outboard.GJ
        sd1 = st.size_wing_section!(wing.outboard, wing.sweep, 1.0; t_cap_min = 0.0, t_web_min = 0.0)
        @test (sd1.tbweb, sd1.tbcap, sd1.Abweb, sd1.Abcap) == (sd0.tbweb, sd0.tbcap, sd0.Abweb, sd0.Abcap)
        @test wing.outboard.EI == EI0 && wing.outboard.GJ == GJ0

        sd2 = st.size_wing_section!(wing.outboard, wing.sweep, 1.0; t_cap_min = 2 * sd0.tbcap)
        @test sd2.tbcap == 2 * sd0.tbcap
        @test sd2.tbweb == sd0.tbweb
    end

    @testset "size_wing_section! absolute minimum gauge of the section" begin
        section = deepcopy(wing.outboard)
        cp = section.co * cosd(wing.sweep)
        cross = section.cross_section
        web_height = cross.web_to_box_height * cross.thickness_to_chord

        # the floor fields default to zero, so sizing is unchanged
        @test section.min_cap_thickness == 0.0
        @test section.min_web_thickness == 0.0
        sd0 = st.size_wing_section!(section, wing.sweep, 1.0)
        EI0, GJ0 = copy(section.EI), section.GJ

        @testset "floor above the required thickness governs" begin
            cap_floor = 2 * sd0.tbcap * cp
            web_floor = 2 * sd0.tbweb * cp
            section.min_cap_thickness = cap_floor
            section.min_web_thickness = web_floor
            sd = st.size_wing_section!(section, wing.sweep, 1.0)
            @test sd.tbcap == cap_floor / cp
            @test sd.tbweb == web_floor / cp
            @test sd.Abcap ≈ 2 * sd.tbcap * cross.width_to_chord
            @test sd.Abweb ≈ 2 * sd.tbweb * web_height
        end

        @testset "each floor acts on its own member" begin
            section.min_cap_thickness = 2 * sd0.tbcap * cp
            section.min_web_thickness = 0.0
            sd = st.size_wing_section!(section, wing.sweep, 1.0)
            @test sd.tbcap == section.min_cap_thickness / cp
            @test sd.tbweb == sd0.tbweb

            section.min_cap_thickness = 0.0
            section.min_web_thickness = 2 * sd0.tbweb * cp
            sd = st.size_wing_section!(section, wing.sweep, 1.0)
            @test sd.tbcap == sd0.tbcap
            @test sd.tbweb == section.min_web_thickness / cp
        end

        @testset "floor below the required thickness changes nothing" begin
            section.min_cap_thickness = 0.5 * sd0.tbcap * cp
            section.min_web_thickness = 0.5 * sd0.tbweb * cp
            sd = st.size_wing_section!(section, wing.sweep, 1.0)
            @test (sd.tbweb, sd.tbcap, sd.Abweb, sd.Abcap) == (sd0.tbweb, sd0.tbcap, sd0.Abweb, sd0.Abcap)
            @test section.EI == EI0 && section.GJ == GJ0
        end
    end
end

end
