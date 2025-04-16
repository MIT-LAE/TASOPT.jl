ac_test = load_default_model()
# ------------------------------
# Fuselage
# ------------------------------
fuse = ac_test.fuselage
Weight = TASOPT.structures.Weight
fuse.n_decks = 1.00000000000000000000 
fuse.shell.weight = Weight(W = 20947.92842800837024697103 ) 
fuse.shell.weight.r = [ 16.97078113922183817408 ,0.0,0.0] 
fuse.window.W = 11269.98000000000138243195 
fuse.window.r = [ 18.13560000000000016485 ,0.0,0.0] 
fuse.window_W_per_length = 435.00000000000000000000 
fuse.insulation.W = 5018.72835762571048690006 
fuse.insulation.r = [ 18.13560000000000016485 ,0.0,0.0] 
fuse.insulation_W_per_area = 22.00000000000000000000 
fuse.floor.weight.W = 13415.79649570315450546332 
fuse.floor_W_per_area = 60.00000000000000000000 
fuse.cone.weight = Weight(W = 5654.60591629168902727542 ) 
fuse.cone.weight.r = [ 33.37559999999999860165 ,0.0,0.0] 
fuse.bendingmaterial_h.weight = Weight(W = 14269.45818342337406647857 ) 
fuse.bendingmaterial_v.weight = Weight(W = 7354.26902517856979102362 ) 
fuse.weight = 197958.25268773091374896467 
fuse.moment = 3618523.01813895627856254578 
fuse.volume = 422.77743761859198912134 
fuse.weight_frac_stringers = 0.34999999999999997780 
fuse.weight_frac_frame = 0.25000000000000000000 
fuse.weight_frac_skin_addl = 0.20000000000000001110 
fuse.layout.nose_radius = 1.64999999999999991118 
fuse.layout.tail_radius = 2.00000000000000000000 
fuse.layout.l_cabin_cylinder = 23.46959999999999979536 
fuse.layout.x_nose = 0.00000000000000000000 
fuse.layout.x_end = 37.79520000000000123919 
fuse.layout.x_start_cylinder = 6.09600000000000008527 
fuse.layout.x_end_cylinder = 29.56559999999999988063 
fuse.layout.x_pressure_shell_fwd = 5.18160000000000042775 
fuse.layout.x_pressure_shell_aft = 31.08960000000000079012 
fuse.layout.x_cone_end = 35.66159999999999996589 
fuse.bendingmaterial_h.weight.r = [ 31.84857592130468972869 ,0.0,0.0] 
fuse.bendingmaterial_v.weight.r = [ 31.09103979825993491204 ,0.0,0.0] 
fuse.layout.cross_section.radius = 1.95579999999999998295 
fuse.layout.cross_section.bubble_lower_downward_shift = 0.38100000000000000533 
fuse.layout.floor_depth = 0.12700000000000000178 
fuse.layout.taper_tailcone = 0.29999999999999998890 
fuse.ratio_young_mod_fuse_bending = 1.00000000000000000000 
fuse.skin.thickness = 0.00105932076444728292 
fuse.cone.thickness = 0.00122075095854528443 
fuse.layout.thickness_webs = 0.00000000000000000000 
fuse.floor.thickness = 0.00090615659813192320 
fuse.shell.EIh = 2361385008.58385801315307617188 
fuse.bendingmaterial_h.EIh = 24633362601.06134414672851562500 
fuse.bendingmaterial_v.EIh = 24633362601.06134414672851562500 
fuse.shell.EIv = 2317413178.16662979125976562500 
fuse.bendingmaterial_h.EIv = 16793530425.73552703857421875000 
fuse.bendingmaterial_v.EIv = 16793530425.73552703857421875000 
fuse.shell.GJ = 1570883675.35591268539428710938 
fuse.cone.GJ = 1638359427.27964949607849121094 
fuse.APU.W = 7698.76022650000140856719 
fuse.APU.r = [36.57600000000000051159,0.0,0.0] 
fuse.seat.W = 21996.45779000000038649887 
fuse.fixed.W = 13344.66600000000107684173 
fuse.fixed.r = [2.13359999999999994102,0.0,0.0] 
fuse.HPE_sys.r = [18.89760000000000061959,0.0,0.0] 
fuse.HPE_sys.W = 0.01000000000000000021 
fuse.added_payload.W = 76987.60226500000862870365 
fuse.cabin.design_pax = 180.00000000000000000000 
fuse.cabin.exit_limit = 189.00000000000000000000 
fuse.cabin.seat_pitch = 0.76200000000000001066 
fuse.cabin.seat_width = 0.48259999999999997344 
fuse.cabin.seat_height = 1.14300000000000001599 
fuse.cabin.aisle_halfwidth = 0.25400000000000000355 
fuse.cabin.floor_distance = 0.00000000000000000000 
fuse.cabin.cabin_width_main = 0.00000000000000000000 
fuse.cabin.cabin_width_top = 0.00000000000000000000 
fuse.cabin.seats_abreast_main = 0 
fuse.cabin.seats_abreast_top = 0 
fuse.cabin.floor_angle_main = 0.00000000000000000000 
fuse.cabin.floor_angle_top = 0.00000000000000000000 
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
htail.opt_sizing = "fixed_Vh"
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
vtail.opt_sizing = "fixed_Vv"
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
