ac_test = load_default_model()
# ------------------------------
# Fuselage
# ------------------------------
fuse_test = ac_test.fuselage
Weight = TASOPT.structures.Weight
fuse_test.n_decks = 1.00000000000000000000 
fuse_test.shell.weight = Weight(W = 20947.34828228715196019039 ) 
fuse_test.shell.weight.r = [ 16.97078113922183462137 ,0.0,0.0] 
fuse_test.window.W = 11269.98000000000138243195 
fuse_test.window.r = [ 18.13560000000000016485 ,0.0,0.0] 
fuse_test.window_W_per_length = 435.00000000000000000000 
fuse_test.insulation.W = 5018.72835762571048690006 
fuse_test.insulation.r = [ 18.13560000000000016485 ,0.0,0.0] 
fuse_test.insulation_W_per_area = 22.00000000000000000000 
fuse_test.floor.weight.W = 13415.76301200039961258881 
fuse_test.floor_W_per_area = 60.00000000000000000000 
fuse_test.cone.weight = Weight(W = 5667.72197228375807753764 ) 
fuse_test.cone.weight.r = [ 33.3756 ,0.0,0.0] 
fuse_test.bendingmaterial_h.weight = Weight(W = 14314.45400479989075392950 ) 
fuse_test.bendingmaterial_v.weight = Weight(W = 7381.15991952439981105272 ) 
fuse_test.weight = 198042.64183002131176181138 
fuse_test.moment = 3620269.85519108083099126816 
fuse_test.volume = 422.77743761859198912134 
fuse_test.weight_frac_stringers = 0.34999999999999997780 
fuse_test.weight_frac_frame = 0.25000000000000000000 
fuse_test.weight_frac_skin_addl = 0.20000000000000001110 
fuse_test.layout.nose_radius = 1.64999999999999991118 
fuse_test.layout.tail_radius = 2.00000000000000000000 
fuse_test.layout.x_nose = 0.00000000000000000000 
fuse_test.layout.x_end = 37.79520000000000123919 
fuse_test.layout.x_start_cylinder = 6.09600000000000008527 
fuse_test.layout.x_end_cylinder = 29.56559999999999988063 
fuse_test.layout.x_pressure_shell_fwd = 5.18160000000000042775 
fuse_test.layout.x_pressure_shell_aft = 31.08960000000000079012 
fuse_test.layout.x_cone_end = 35.66159999999999996589 
fuse_test.bendingmaterial_h.weight.r = [ 31.85596402453565900714 ,0.0,0.0] 
fuse_test.bendingmaterial_v.weight.r = [ 31.09856247291204311978 ,0.0,0.0] 
fuse_test.layout.cross_section.radius = 1.95579999999999998295 
fuse_test.layout.cross_section.bubble_lower_downward_shift = 0.38100000000000000533 
fuse_test.layout.floor_depth = 0.12700000000000000178 
fuse_test.layout.taper_tailcone = 0.29999999999999998890 
fuse_test.ratio_young_mod_fuse_bending = 1.00000000000000000000 
fuse_test.skin.thickness = 0.00105929142692061514 
fuse_test.layout.cross_section.skin_thickness = fuse_test.skin.thickness
fuse_test.cone.thickness = 0.00122358253304608554 
fuse_test.layout.thickness_webs = 0.00000000000000000000 
fuse_test.floor.thickness = 0.00090615659813192320 
fuse_test.shell.EIh = 2361319610.83278226852416992188 
fuse_test.bendingmaterial_h.EIh = 24686146267.78833389282226562500 
fuse_test.bendingmaterial_v.EIh = 24686146267.78833389282226562500 
fuse_test.shell.EIv = 2317348998.20037508010864257812 
fuse_test.bendingmaterial_h.EIv = 16837922621.63796615600585937500 
fuse_test.bendingmaterial_v.EIv = 16837922621.63796615600585937500 
fuse_test.shell.GJ = 1570840170.26919627189636230469 
fuse_test.cone.GJ = 1642159659.21472048759460449219 
fuse_test.APU.W = 7698.76022650000140856719 
fuse_test.APU.r = [36.57600000000000051159,0.0,0.0] 
fuse_test.seat.W = 21996.45779000000038649887 
fuse_test.fixed.W = 13344.66600000000107684173 
fuse_test.fixed.r = [2.13359999999999994102,0.0,0.0] 
fuse_test.HPE_sys.r = [18.89760000000000061959,0.0,0.0] 
fuse_test.HPE_sys.W = 0.01000000000000000021 
fuse_test.added_payload.W = 76987.60226500000862870365 