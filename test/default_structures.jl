ac_test = load_default_model()
# ------------------------------
# Fuselage
# ------------------------------
fuse = ac_test.fuselage
Weight = TASOPT.structures.Weight
fuse.n_decks = 1.00000000000000000000 
fuse.shell.weight = Weight(W = 20947.35267315476812655106 ) 
fuse.shell.weight.r = [ 16.97078113922184527951 ,0.0,0.0] 
fuse.window.W = 11269.98000000000138243195 
fuse.window.r = [ 18.13560000000000016485 ,0.0,0.0] 
fuse.window_W_per_length = 435.00000000000000000000 
fuse.insulation.W = 5018.72835762571048690006 
fuse.insulation.r = [ 18.13560000000000016485 ,0.0,0.0] 
fuse.insulation_W_per_area = 22.00000000000000000000 
fuse.floor.weight.W = 13415.79649570315450546332 
fuse.floor_W_per_area = 60.00000000000000000000 
fuse.cone.weight = Weight(W = 5667.81047483775273576612 ) 
fuse.cone.weight.r = [ 33.37559999999999860165 ,0.0,0.0] 
fuse.bendingmaterial_h.weight = Weight(W = 14314.08207479297561803833 ) 
fuse.bendingmaterial_v.weight = Weight(W = 7380.92514863496762700379 ) 
fuse.weight = 198042.16150624936562962830 
fuse.moment = 3620261.28382447082549333572 
fuse.volume = 422.77743761859198912134 
fuse.weight_frac_stringers = 0.34999999999999997780 
fuse.weight_frac_frame = 0.25000000000000000000 
fuse.weight_frac_skin_addl = 0.20000000000000001110 
fuse.layout.nose_radius = 1.64999999999999991118 
fuse.layout.tail_radius = 2.00000000000000000000 
fuse.layout.x_nose = 0.00000000000000000000 
fuse.layout.x_end = 37.79520000000000123919 
fuse.layout.x_start_cylinder = 6.09600000000000008527 
fuse.layout.x_end_cylinder = 29.56559999999999988063 
fuse.layout.x_pressure_shell_fwd = 5.18160000000000042775 
fuse.layout.x_pressure_shell_aft = 31.08960000000000079012 
fuse.layout.x_cone_end = 35.66159999999999996589 
fuse.bendingmaterial_h.weight.r = [ 31.85585727200864880615 ,0.0,0.0] 
fuse.bendingmaterial_v.weight.r = [ 31.09848210139743684977 ,0.0,0.0] 
fuse.layout.cross_section.radius = 1.95579999999999998295 
fuse.layout.cross_section.bubble_lower_downward_shift = 0.38100000000000000533 
fuse.layout.floor_depth = 0.12700000000000000178 
fuse.layout.taper_tailcone = 0.29999999999999998890 
fuse.ratio_young_mod_fuse_bending = 1.00000000000000000000 
fuse.skin.thickness = 0.00105929164896344176 
fuse.cone.thickness = 0.00122360163951950240 
fuse.layout.thickness_webs = 0.00000000000000000000 
fuse.floor.thickness = 0.00090615659813192320 
fuse.shell.EIh = 2361320105.79957580566406250000 
fuse.bendingmaterial_h.EIh = 24685708558.26182556152343750000 
fuse.bendingmaterial_v.EIh = 24685708558.26182556152343750000 
fuse.shell.EIv = 2317349483.95029115676879882812 
fuse.bendingmaterial_h.EIv = 16837510581.35239791870117187500 
fuse.bendingmaterial_v.EIv = 16837510581.35239791870117187500 
fuse.shell.GJ = 1570840499.54004931449890136719 
fuse.cone.GJ = 1642185301.84938359260559082031 
fuse.APU.W = 7698.76022650000140856719 
fuse.APU.r = [36.57600000000000051159,0.0,0.0] 
fuse.seat.W = 21996.45779000000038649887 
fuse.fixed.W = 13344.66600000000107684173 
fuse.fixed.r = [2.13359999999999994102,0.0,0.0] 
fuse.HPE_sys.r = [18.89760000000000061959,0.0,0.0] 
fuse.HPE_sys.W = 0.01000000000000000021 
fuse.added_payload.W = 76987.60226500000862870365 
# ------------------------------
# Wing
# ------------------------------
wing = ac_test.wing
wing.inboard.webs.weight = Weight(W = 3591.78987989722190832254) 
wing.outboard.webs.weight = Weight(W = 3591.78987989722190832254) 
wing.inboard.caps.weight = Weight(W = 75100.79461642308160662651) 
wing.outboard.caps.weight = Weight(W = 75100.79461642308160662651) 
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
wing.weight = 129055.83864577353233471513 
wing.strut.weight = 0.00000000000000000000 
wing.dxW = 300492.41100064787315204740 
wing.strut.dxW = 0.00000000000000000000 
wing.inboard.weight = 49973.26554599140945356339 
wing.outboard.weight = 70883.07125228032236918807 
wing.inboard.dyW = 78397.50028844068583566695 
wing.outboard.dyW = 336207.94655290641821920872 
wing.weight_frac_flap = 0.20000000000000001110 
wing.weight_frac_slat = 0.10000000000000000555 
wing.weight_frac_ailerons = 0.04000000000000000083 
wing.weight_frac_leading_trailing_edge = 0.10000000000000000555 
wing.weight_frac_ribs = 0.14999999999999999445 
wing.weight_frac_spoilers = 0.02000000000000000042 
wing.weight_frac_attachments = 0.02999999999999999889 
wing.strut.local_velocity_ratio = 1.00000000000000000000 
wing.layout.x = 18.46577097152568569527 
wing.layout.box_x = 15.68841088447713261189 
wing.layout.z = -1.67640000000000011227 
wing.strut.cos_lambda = 1.00000000000000000000 
wing.strut.S = 0.00000000000000000000 
wing.layout.spar_box_x_c = 0.40000000000000002220 
wing.layout.box_width = 0.50000000000000000000 
wing.inboard.layout.chord_thickness = 0.12679999999999999605 
wing.outboard.layout.chord_thickness = 0.12659999999999999032 
wing.layout.hweb_to_hbox = 0.75000000000000000000 
wing.layout.b_max = 35.81400000000000005684 
wing.strut.thickness_to_chord = 0.14999999999999999445 
wing.strut.z = 3.91159999999999996589 
wing.outboard.moment = 3483196.13483797945082187653 
wing.outboard.max_shear_load = 622048.31389977096114307642 
wing.outboard.web_cap.GJ = 279675706.51725339889526367188 
wing.outboard.web_cap.EI_normal = 2509473607.51622676849365234375 
wing.outboard.web_cap.EI_bending = 324582402.63841980695724487305 
wing.outboard.caps.thickness = 0.00647934795725599311 
wing.inboard.moment = 6112466.36944914143532514572 
wing.inboard.max_shear_load = 850343.35977413994260132313 
wing.inboard.web_cap.GJ = 740693528.54392838478088378906 
wing.inboard.web_cap.EI_normal = 6101999745.09585475921630859375 
wing.inboard.web_cap.EI_bending = 816023778.80346381664276123047 
wing.inboard.caps.thickness = 0.00370725840496802239 
wing.inboard.webs.thickness = 0.00103100060251482169 
wing.outboard.webs.thickness = 0.00154162299420744823 
wing.layout.S = 139.79117197415885698319 
wing.outboard.layout.b = 3.60679999999999978400 
wing.ηs = 0.28499999999999997558 
wing.inboard.layout.λ = 0.69999999999999995559 
wing.outboard.layout.λ = 0.25000000000000000000 
wing.layout.chord = 6.23927008790646819847 
wing.inboard.layout.b = 10.70891372784236317273 
wing.layout.b= 37.57513588716619068464 
wing.layout.sweep = 26.00000000000000000000 
wing.layout.AR = 10.09999999999999964473 
# ------------------------------
# Htail
# ------------------------------
htail = ac_test.htail
htail.weight = 14400.81550210478599183261 
htail.dxW = 14069.61121976982394699007 
htail.weight_fraction_added = 0.29999999999999998890 
htail.layout.box_x = 34.89959999999999951115 
htail.layout.z = 0.00000000000000000000 
htail.downwash_factor = 0.59999999999999997780 
htail.CL_max_fwd_CG = -0.11129374433574321102 
htail.CL_max = 2.00000000000000000000 
htail.SM_min = 0.05000000000000000278 
htail.layout.x = 36.21247502119915395724 
htail.layout.box_width = 0.50000000000000000000 
htail.outboard.chord_thickness = 0.14000000000000001332 
htail.layout.hweb_to_hbox = 0.75000000000000000000 
htail.thickness_web = 0.00115688496649875012 
htail.move_wingbox = 2.00000000000000000000 
htail.CL_CLmax = -0.50000000000000000000 
htail.size = 1.00000000000000000000 
htail.volume = 1.44999999999999862332 
htail.thickness_cap = 0.00239055785595602283 
htail.GJ = 198224681.32888934016227722168 
htail.EI_normal = 1261677460.09846663475036621094 
htail.EI_bending = 189632296.67963379621505737305 
htail.layout.sweep = 26.00000000000000000000 
htail.layout.chord = 4.55763820582986856067 
htail.outboard.λ = 0.25000000000000000000 
htail.outboard.b = 1.52400000000000002132 
htail.layout.b = 17.09114327186200554820 
htail.layout.AR = 6.00000000000000000000 
htail.layout.S = 48.68452972321898641894 
# ------------------------------
# Vtail
# ------------------------------
vtail = ac_test.vtail
vtail.weight = 9645.22120844134951767046 
vtail.dxW = 12039.30921083126122539397 
vtail.weight_fraction_added = 0.40000000000000002220 
vtail.layout.box_x = 33.52799999999999869260 
vtail.CL_max = 2.60000000000000008882 
vtail.layout.x = 35.05057234350366712761 
vtail.layout.box_width = 0.50000000000000000000 
vtail.outboard.chord_thickness = 0.14000000000000001332 
vtail.layout.hweb_to_hbox = 0.75000000000000000000 
vtail.ntails = 1.00000000000000000000 
vtail.volume = 0.10000000000000000555 
vtail.thickness_web = 0.00116813883699139025 
vtail.thickness_cap = 0.00183112800700054035 
vtail.GJ = 554905788.63482689857482910156 
vtail.EI_normal = 3494962554.28956699371337890625 
vtail.EI_bending = 496717199.07202911376953125000 
vtail.layout.sweep = 25.00000000000000000000 
vtail.layout.chord = 6.12218798928765028933 
vtail.outboard.λ = 0.29999999999999998890 
vtail.layout.b = 7.95884438607394706366 
vtail.layout.AR = 2.00000000000000000000 
vtail.layout.S = 31.67160198087038835979 
vtail.size = 1.00000000000000000000 
vtail.dxW = 12039.30921083126122539397 