using TOML
export save_aircraft_model

"""
    save_aircraft_model(ac::TASOPT.aircraft=TASOPT.read_aircraft_model(), 
    datafile=joinpath(TASOPT.__TASOPTroot__, "IO/default_output.toml"),
    save_output::Bool=false)

Converts an aircraft model into a dictionary and writes 
it to a TOML file. Values to be written are explicitly set
following the default_input.toml. All values are written in SI units.

This save operation makes add'l* assumptions about parameter repetition. Namely:
The same value is applied for all flight segments/points for:
    - parm[] parameters
    - excrescence_drag_factors, wing overspeeds, wing/stabilizer Re_refs
The same value is applied for all missions and flight segments for:
    - parg[], pare[], and pari[] parameters
    - fuel temperature

Said value is the first entry in the corresponding array axis, 
except for some aero parameters where other points are more relevant (e.g., "Cruise" "Takeoff").

*and modifiable

!!! note "Deviating from default"
    Extending `read_input.jl` and `save_model.jl` is recommended for models deviating appreciably 
    from the default functionality. Thorough knowledge of the model is required.
"""
function save_aircraft_model(ac::TASOPT.aircraft=TASOPT.read_aircraft_model(), 
    datafile=joinpath(TASOPT.__TASOPTroot__, "IO/IO_samples/default_output.toml"),
    save_output::Bool=false)

    #get parameter arrays from aircraft struct
    ac_i = ac.pari      #options parameters
    ac_g = ac.parg      #geometry   "
    ac_m = ac.parm      #mission    "
    ac_a = ac.para      #aero       "
    ac_e = ac.pare      #engine     "

    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    fuselage = ac.fuselage  #fuselage     "
    landing_gear = ac.landing_gear

    #dictionaries to map some selections (e.g., ints) to outputs
    propsysarch = Dict(0 => "te", 1 => "tf")
    engloc = Dict(1 => "wing", 2 => "fuselage")
    fueltype = Dict(1 => "LH2", 24 => "JET-A")
    taperfuse = Dict(0=>"point", 1=>"edge")
    engweightmodel = Dict(0 => "md", 1 =>"basic", 2=>"advanced")

    #Save everything in a dict() of dicts()
    d_out = Dict()
    #operations on ac
    #TODO: reverse operations of read_input.jl
    # General description data
    d_desc = Dict()
    d_desc["name"] = ac.name
    d_desc["description"] = ac.description
    d_desc["sized"] = ac.sized
    d_out["Aircraft Description"] = d_desc

    #Options------------------------
    d_opt = Dict()
        d_opt["optimize"] = ac_i[iiopt]
        d_opt["prop_sys_arch"] = propsysarch[ac_i[iiengtype]]
        d_opt["engine_location"] = engloc[ac_i[iiengloc]]
    d_out["Options"] = d_opt
    #--end options----------------

    #Fuel------------------------

    d_fuel = Dict()
        d_fuel["fuel_type"] = fueltype[ac_i[iifuel]]
        d_fuel["fuel_in_wing"] = ac_i[iifwing]
        d_fuel["fuel_in_wingcen"] = ac_i[iifwcen]
        d_fuel["fuel_usability_factor"] = ac_g[igrWfmax]

        d_fuel["fuel_temp"] = ac_e[ieTfuel,1,1]             

        d_fuel["fuel_density"] = ac_g[igrhofuel]
    d_out["Fuel"] = d_fuel
    #--end fuel----------------

    #Mission------------------------
    d_miss = Dict()

        d_miss["N_missions"] = size(ac_m,2)
        
        d_miss["range"] = ac_m[imRange,:]
        d_miss["weight_per_pax"] = ac_m[imWperpax, :]
        d_miss["pax"] = ac_m[imWpay,:] ./ ac_m[imWperpax, :]
        d_miss["max_pax"] = ac_g[igWpaymax] ./ ac_m[imWperpax, :]
        d_miss["fuel_reserves"] = ac_g[igfreserve]
        d_miss["Vne"] = ac_g[igVne]
        d_miss["Nlift"] = ac_g[igNlift]

    # mission: Takeoff
    d_to = Dict()
        d_to["takeoff_alt"] = ac_m[imaltTO, :]
        d_to["braking_resistance_coeff"] = ac_g[igmubrake]
        d_to["rolling_resistance_coeff"] = ac_g[igmuroll]
        d_to["takeoff_obstacle_height"] = ac_g[ighobst]
        d_to["CD_dead_engine"] = ac_g[igcdefan]
        d_to["CD_landing_gear"] = ac_g[igCDgear]
        d_to["CD_spoilers"] = ac_g[igCDspoil]
        d_to["max_balanced_field_length"] = ac_g[iglBFmax]
        d_to["Nland"] = ac_g[igNland]

        d_to["takeoff_T"] = ac_m[imT0TO, :]
        d_to["CL_max_perp"] = ac_a[iaclpmax, iptakeoff, 1]
    d_miss["Takeoff"] = d_to

    # mission: Climb
    d_climb = Dict()
        d_climb["minimum_top-of-climb_gradient"] = ac_g[iggtocmin]
    d_miss["Climb"] = d_climb

    # mission: Cruise
    d_crz = Dict()
        d_crz["cruise_alt"] = ac_a[iaalt, ipcruise1, :]
        d_crz["cruise_mach"] = ac_a[iaMach, ipclimbn, :]
        d_crz["cruise_CL"] = ac_a[iaCL, ipclimb1+1, :]
    d_miss["Cruise"] = d_crz

    # mission: Descent
    d_desc = Dict()
        d_desc["descent_angle_top-of-descent"] = ac_m[imgamVDE1,:]
        d_desc["descent_angle_bottom-of-descent"] = ac_m[imgamVDEn, :]
    d_miss["Descent"] = d_desc

    d_out["Mission"] = d_miss #add to output dictionary
    #--end Mission----------------

    #Fuselage------------------------
    d_fuse = Dict()
        d_fuse["cabin_pressure"] = ac_g[igpcabin]

    #aero 
    d_fuse_aero = Dict()
        d_fuse_aero["excrescence_drag_factor"] = ac_a[iafexcdf,1,:]
        d_fuse_aero["BLI_frac"] = ac_g[igfBLIf]

        d_fuse_aero["wingroot_fuse_overspeed"] = ac_a[iafduo,1,:]
        d_fuse_aero["wingbreak_fuse_overspeed"] = ac_a[iafdus,1,:]
        d_fuse_aero["wingtip_fuse_overspeed"] = ac_a[iafdut,1,:]

        d_fuse_aero["fuse_moment_volume_deriv"] = ac_g[igCMVf1]
        d_fuse_aero["CL_zero_fuse_moment"] = ac_g[igCLMf0]
    d_fuse["Aero"] = d_fuse_aero

    #weight
    d_fuse_weights = Dict()
        d_fuse_weights["frame"] = fuselage.weight_frac_frame
        d_fuse_weights["stringer"] = fuselage.weight_frac_stringers
        d_fuse_weights["additional"] = fuselage.weight_frac_skin_addl
        d_fuse_weights["fixed_weight"] = fuselage.fixed.W

        d_fuse_weights["window_per_length"] = fuselage.window_W_per_length
        d_fuse_weights["window_insul_per_area"] = fuselage.insulation_W_per_area
        d_fuse_weights["floor_weight_per_area"] =fuselage.floor_W_per_area

        d_fuse_weights["HPE_sys_weight_fraction"] = fuselage.HPE_sys.W

        
        d_fuse_weights["add_payload_weight_fraction"] = fuselage.added_payload.W/ac_g[igWpaymax] 
    d_fuse["Weights"] = d_fuse_weights


    #geometry
    d_fuse_geom = Dict()
        d_fuse_geom["radius"] = fuselage.layout.cross_section.radius
        d_fuse_geom["dRadius"] = fuselage.layout.cross_section.bubble_lower_downward_shift
        d_fuse_geom["y_offset"] = fuselage.layout.bubble_center_y_offset
        d_fuse_geom["floor_depth"] = fuselage.layout.floor_depth
        d_fuse_geom["Nwebs"] = fuselage.layout.n_webs

        d_fuse_geom["a_nose"] = fuselage.layout.nose_radius
        d_fuse_geom["b_tail"] = fuselage.layout.tail_radius

        d_fuse_geom["taper_fuse_to"] = taperfuse[ac_i[iifclose]]

        d_fuse_geom["tailcone_taper"] = fuselage.layout.taper_tailcone
        d_fuse_geom["HT_load_fuse_bend_relief"] = ac_g[igrMh]
        d_fuse_geom["VT_load_fuse_bend_relief"] = ac_g[igrMv]

        d_fuse_geom["calculate_cabin_length"] = false #Use final fuselage parameters by default
        if ac_i[iidoubledeck] == 1
            d_fuse_geom["double_decker"] = true
        else
            d_fuse_geom["double_decker"] = false
        end
      
        d_fuse_geom["seat_pitch"] = fuselage.cabin.seat_pitch
        d_fuse_geom["seat_width"] = fuselage.cabin.seat_width
        d_fuse_geom["aisle_halfwidth"] = fuselage.cabin.aisle_halfwidth  

        d_fuse_geom["x_nose_tip"] = fuselage.layout.x_nose
        d_fuse_geom["x_pressure_shell_fwd"] = fuselage.layout.x_pressure_shell_fwd
        d_fuse_geom["x_start_cylinder"] = fuselage.layout.x_start_cylinder
        d_fuse_geom["x_end_cylinder"] = fuselage.layout.x_end_cylinder
        d_fuse_geom["x_pressure_shell_aft"] = fuselage.layout.x_pressure_shell_aft
        d_fuse_geom["x_cone_end"] = fuselage.layout.x_cone_end
        d_fuse_geom["x_end"] = fuselage.layout.x_end

        d_fuse_geom["x_APU"] = fuselage.APU.x
        d_fuse_geom["x_HPE_sys"] = fuselage.HPE_sys.x

        d_fuse_geom["x_fixed_weight"] = fuselage.fixed.x

        d_fuse_geom["x_engines"] = ac_g[igxeng]
        d_fuse_geom["y_critical_engines"] = ac_g[igyeng]
    d_fuse["Geometry"] = d_fuse_geom

    d_out["Fuselage"] = d_fuse
    #--end Fuselage----------------

    #Wing ------------------------
    d_wing = Dict()
        d_wing["planform"] = wing.planform
        d_wing["strut_braced_wing"] = wing.has_strut

        d_wing["sweep"] = wing.layout.sweep
        d_wing["AR"] = wing.layout.AR
        d_wing["maxSpan"] = wing.layout.max_span

        d_wing["inner_panel_taper_ratio"] = wing.inboard.λ
        d_wing["outer_panel_taper_ratio"] = wing.outboard.λ
        d_wing["panel_break_location"] = wing.layout.ηs

        d_wing["center_box_halfspan"] = wing.layout.root_span/2
        d_wing["box_width_to_chord"] = wing.inboard.cross_section.width_to_chord
        d_wing["root_thickness_to_chord"] = wing.inboard.cross_section.thickness_to_chord
        d_wing["spanbreak_thickness_to_chord"] = wing.outboard.cross_section.thickness_to_chord
        d_wing["hweb_to_hbox"] = wing.inboard.cross_section.web_to_box_height
        d_wing["spar_box_x_c"] = wing.layout.spar_box_x_c

        d_wing["x_wing_box"] = wing.layout.box_x
        d_wing["z_wing"] = wing.layout.z

    # Strut details (only used if strut_braced_wing is True)
        if d_wing["strut_braced_wing"]
            d_wing["z_strut"] = wing.strut.z
            d_wing["strut_toc"] = wing.strut.thickness_to_chord
            d_wing["strut_local_velocity_ratio"] = wing.strut.local_velocity_ratio
        end

    #Aero , for multiple segments
    d_wing_aero = Dict()
        d_wing_aero["fuselage_lift_carryover_loss_factor"] = wing.fuse_lift_carryover
        d_wing_aero["wing_tip_lift_rolloff_factor"] = wing.tip_lift_loss

        d_wing_aero["lowspeed_cdf"] = ac_a[iacdfw, 1,:]
        d_wing_aero["lowspeed_cdp"] = ac_a[iacdpw, 1,:]
        d_wing_aero["Re_ref"] = ac_a[iaRerefw, 1,:]

        d_wing_aero["strut_lowspeed_cdf"] = ac_a[iacdfs, 1,:]
        d_wing_aero["strut_lowspeed_cdp"] = ac_a[iacdps, 1,:]
        d_wing_aero["strut_Re_ref"] = ac_a[iaRerefs, 1,:]

        d_wing_aero["Reynolds_scaling"] = ac_a[iaaRexp, 1,:]
        d_wing_aero["excrescence_drag_factor"] = ac_a[iafexcdw, 1,:]
        d_wing_aero["BLI_frac"] = ac_g[igfBLIw]
      
    d_wing_aero_to = Dict()
        d_wing_aero_to["cls_clo"] = ac_a[iarcls, 1,:]
        d_wing_aero_to["clt_clo"] = ac_a[iarclt, 1,:]
        d_wing_aero_to["cm_o"] = ac_a[iacmpo, 1,:]
        d_wing_aero_to["cm_s"] = ac_a[iacmps, 1,:]
        d_wing_aero_to["cm_t"] = ac_a[iacmpt, 1,:]
    d_wing_aero["Takeoff"] = d_wing_aero_to
    
    d_wing_aero_clm = Dict()
        d_wing_aero_clm["cls_clo"] = ac_a[iarcls, ipclimb1+1,:]
        d_wing_aero_clm["clt_clo"] = ac_a[iarclt, ipclimb1+1,:]
        d_wing_aero_clm["cm_o"] = ac_a[iacmpo, ipclimb1+1,:]
        d_wing_aero_clm["cm_s"] = ac_a[iacmps, ipclimb1+1,:]
        d_wing_aero_clm["cm_t"] = ac_a[iacmpt, ipclimb1+1,:]
    d_wing_aero["Climb"] = d_wing_aero_clm

    d_wing_aero_la = Dict()
        d_wing_aero_la["cls_clo"] = ac_a[iarcls, ipdescentn,:]
        d_wing_aero_la["clt_clo"] = ac_a[iarclt, ipdescentn,:]
        d_wing_aero_la["cm_o"] = ac_a[iacmpo, ipdescentn,:]
        d_wing_aero_la["cm_s"] = ac_a[iacmps, ipdescentn,:]
        d_wing_aero_la["cm_t"] = ac_a[iacmpt, ipdescentn,:]
    d_wing_aero["Landing"] = d_wing_aero_la

    d_wing["Aero"] = d_wing_aero

    #Weight
    d_wing_weight = Dict()
        d_wing_weight["flap"] = wing.weight_frac_flap
        d_wing_weight["slat"] = wing.weight_frac_slat
        d_wing_weight["aileron"] = wing.weight_frac_ailerons
        d_wing_weight["leading_trailing_edge"] = wing.weight_frac_leading_trailing_edge
        d_wing_weight["ribs"] = wing.weight_frac_ribs
        d_wing_weight["spoilers"] = wing.weight_frac_spoilers
        d_wing_weight["attachments"] = wing.weight_frac_attachments
    d_wing["Weightfracs"] = d_wing_weight

    d_out["Wing"] = d_wing
    #--end Wing----------------


    #Stabilizers------------------------
    d_stab = Dict()

        d_stab["lowspeed_cdf"] = ac_a[iacdft,1,:]
        d_stab["lowspeed_cdp"] = ac_a[iacdpt,1,:]
        d_stab["Re_ref"] = ac_a[iaRereft,1,:]

        d_stab["excrescence_drag_factor"] = ac_a[iafexcdt,1,:]


    #Horz tail
    d_stab_htail  = Dict()
        d_stab_htail["AR_Htail"] = htail.layout.AR
        d_stab_htail["taper"] = htail.outboard.λ
        d_stab_htail["sweep"] = htail.layout.sweep
        d_stab_htail["center_box_halfspan"] = htail.layout.root_span/2

        d_stab_htail["x_Htail"] = htail.layout.box_x
        d_stab_htail["z_Htail"] = htail.layout.z

        d_stab_htail["max_tail_download"] = htail.CL_CLmax
    
        if htail.size == 1
            d_stab_htail["HTsize"] = "vh"
            d_stab_htail["Vh"] = htail.volume
        elseif htail.size == 2
            d_stab_htail["HTsize"] = "maxforwardcg"
            d_stab_htail["CLh_at_max_forward_CG"] = htail.CL_max_fwd_CG
        end

        d_stab_htail["move_wingbox"] = htail.move_wingbox

        d_stab_htail["SM_min"] = htail.SM_min

        d_stab_htail["CLh_spec"] = ac_g[igCLhspec]

        d_stab_htail["downwash_factor"] = htail.downwash_factor
        d_stab_htail["nacelle_lift_curve_slope"] = ac_g[igdCLnda]

        d_stab_htail["CD_Htail_from_center"] = ac_g[igfCDhcen]
        d_stab_htail["CLh_max"] = htail.CL_max

        d_stab_htail["added_weight_fraction"] = htail.weight_fraction_added

        d_stab_htail["box_width_to_chord"] = htail.outboard.cross_section.width_to_chord
        d_stab_htail["box_height_chord"] = htail.outboard.cross_section.thickness_to_chord
        d_stab_htail["web_height_hbox"] = htail.outboard.cross_section.web_to_box_height
    d_stab["Htail"] = d_stab_htail

    #Vertical tail
    d_stab_vtail = Dict()
        d_stab_vtail["AR_Vtail"] = vtail.layout.AR
        d_stab_vtail["taper"] = vtail.outboard.λ
        d_stab_vtail["sweep"] = vtail.layout.sweep
        d_stab_vtail["center_box_halfspan"] = vtail.layout.root_span
        d_stab_vtail["x_Vtail"] = vtail.layout.box_x
        d_stab_vtail["number_Vtails"] = vtail.ntails

        if vtail.size == 1
            d_stab_vtail["VTsize"] = "vv"
            d_stab_vtail["Vv"] = vtail.volume
        elseif avtail.size == 2
            d_stab_vtail["VTsize"] = "oei"
            d_stab_vtail["CLv_at_engine_out"] = ac_g[igCLveout]
        end

        d_stab_vtail["CLv_max"] = vtail.CL_max
        d_stab_vtail["added_weight_fraction"] = vtail.weight_fraction_added
        d_stab_vtail["box_width_to_chord"] = vtail.outboard.cross_section.width_to_chord
        d_stab_vtail["box_height_chord"] = vtail.outboard.cross_section.thickness_to_chord
        d_stab_vtail["web_height_hbox"] = vtail.outboard.cross_section.web_to_box_height

    d_stab["Vtail"] = d_stab_vtail

    d_out["Stabilizers"] = d_stab
    #--end Stabilizers----------------

    #Landing gear------------------------
    d_lg  = Dict()
        d_lg["x_nose_landing_gear"] = landing_gear.nose_gear.weight.r[1]
        d_lg["x_main_landing_gear_offset"] = landing_gear.main_gear.distance_CG_to_landing_gear
        d_lg["LG_nose_weight_fraction"] = landing_gear.nose_gear.overall_mass_fraction
        d_lg["LG_main_weight_fraction"] = landing_gear.main_gear.overall_mass_fraction
    
    d_out["Landing_gear"] = d_lg
    #--end Landing gear----------------

    #Structures------------------------
    d_struct = Dict()
        d_struct["stress_factor"] = ac_g[igsigfac]

        d_struct["sigma_fuse_skin"] = fuselage.skin.σ
        d_struct["sigma_fuse_bending"] = fuselage.bendingmaterial_h.σ

        d_struct["sigma_caps"] = wing.inboard.caps.σ
        d_struct["tau_webs"] = wing.inboard.webs.σ

        d_struct["sigma_struts"] = wing.strut.material.σmax
        
        d_struct["fuse_shell_modulus_ratio"] = fuselage.ratio_young_mod_fuse_bending
        
        d_struct["E_wing_spar_cap"] = wing.inboard.caps.material.E
        d_struct["E_struts"] = wing.strut.material.E
        
        # Structural material densities
        d_struct["skin_density"] = fuselage.skin.ρ
        d_struct["fuse_stringer_density"] = fuselage.bendingmaterial_h.ρ
        d_struct["wing_tail_cap_density"] = wing.inboard.caps.ρ
        d_struct["wing_tail_web_density"] = wing.inboard.webs.ρ
        d_struct["strut_density"] = wing.strut.material.ρ
    d_out["Structures"] = d_struct
    #--end Structures----------------

    #Propulsion ------------------------
    d_prop = Dict()
        d_prop["number_of_engines"] = ac_g[igneng]
        d_prop["T_max_metal"] = ac_g[igTmetal]
        d_prop["Tt4_frac_bottom_of_climb"] = ac_g[igfTt4CL1]
        d_prop["Tt4_frac_top_of_climb"] = ac_g[igfTt4CLn]
        
        d_prop["Tt4_cruise"] = ac_e[ieTt4,ipcruise1,:]
        d_prop["Tt4_takeoff"] =  ac_e[ieTt4,ipstatic,:]

        d_prop["core_in_clean_flow"] = !Bool(ac_i[iiBLIc])
            #expression negates 0 to 1 and vice versa (), see read_input.jl
        
    #Turbomachinery
    d_prop_turb = Dict()
        d_prop_turb["BPR"] = ac_e[ieBPR, 1, 1]
        d_prop_turb["Fan_PR"] = ac_e[iepif, 1, 1]
        d_prop_turb["LPC_PR"] = ac_e[iepilc, 1, 1]
        d_prop_turb["OPR"] = d_prop_turb["LPC_PR"]*ac_e[iepihc, 1, 1]
        
        d_prop_turb["diffuser_PR"] = ac_e[iepid, 1, 1]
        d_prop_turb["burner_PR"] = ac_e[iepib, 1, 1]
        d_prop_turb["fan_nozzle_PR"] = ac_e[iepifn, 1, 1]
        d_prop_turb["core_nozzle_PR"] = ac_e[iepitn, 1, 1]
        
        d_prop_turb["fan_eta_poly"] = ac_e[ieepolf, 1, 1]
        d_prop_turb["LPC_eta_poly"] = ac_e[ieepollc, 1, 1]
        d_prop_turb["HPC_eta_poly"] = ac_e[ieepolhc, 1, 1]
        d_prop_turb["HPT_eta_poly"] = ac_e[ieepolht, 1, 1]
        d_prop_turb["LPT_eta_poly"] = ac_e[ieepollt, 1, 1]
        
        d_prop_turb["FPR0"] = ac_e[iepifK, 1, 1]
        d_prop_turb["Kf_polyeff"] = ac_e[ieepfK, 1, 1]
        
        d_prop_turb["M2"] = ac_e[ieM2, 1, 1]
        d_prop_turb["M25"] = ac_e[ieM25, 1, 1]
        
        d_prop_turb["low_spool_loss"] = ac_e[ieepsl, 1, 1]
        d_prop_turb["high_spool_loss"] = ac_e[ieepsh, 1, 1]

        d_prop_turb["gear_ratio"] = ac_g[igGearf]
        d_prop_turb["HTR_fan"] = ac_g[igHTRf]
        d_prop_turb["HTR_LPC"] = ac_g[igHTRlc]
        d_prop_turb["HTR_HPC"] = ac_g[igHTRhc]
    d_prop["Turbomachinery"] = d_prop_turb

    #Combustor
    d_prop_comb = Dict()
        d_prop_comb["combustion_efficiency"] = ac_e[ieetab,1,1]
    d_prop["Combustor"] = d_prop_comb

    #Cooling
    d_prop_cool = Dict()
        d_prop_cool["hot_streak_T_allowance"] = ac_e[iedTstrk,1,1]
        d_prop_cool["M_turbine_blade_exit"] = ac_e[ieMtexit,1,1]
        d_prop_cool["St"] = ac_e[ieStA,1,1]

        d_prop_cool["e_film_cooling"] = ac_e[ieefilm,1,1]
        d_prop_cool["t_film_cooling"] = ac_e[ietfilm,1,1]
        
        d_prop_cool["M41"] = ac_e[ieM4a,1,1]
        d_prop_cool["cooling_air_V_ratio"] = ac_e[ieruc,1,1]
    d_prop["Cooling"] = d_prop_cool

    #Offtakes
    d_prop_offt = Dict()
        #TODO: check if htese expressions are correct (Wpax)
        #Also, check matrix slicing thing here too
        d_prop_offt["LPC_mass_offtake_per_pax"] = ac_g[igmofWpay]*ac_m[imWperpax, 1]
        d_prop_offt["LPC_mass_offtake_per_max_mass"] = ac_g[igmofWMTO]*gee

        d_prop_offt["Low_spool_power_offtake_per_pax"] = ac_g[igPofWpay]*ac_m[imWperpax, 1]
        d_prop_offt["Low_spool_power_offtake_per_max_mass"] = ac_g[igPofWMTO]*gee

        d_prop_offt["Tt_offtake_air"] = ac_e[ieTt9,1,1]
        d_prop_offt["Pt_offtake_air"] = ac_e[iept9,1,1]
    d_prop["Offtakes"] = d_prop_offt

    #Nozzles
    d_prop_nozz = Dict()
        #core nozzle
        d_prop_cnoz = Dict()
            d_prop_cnoz["static"] = ac_e[ieA5fac, ipstatic,1]
            d_prop_cnoz["rotation"] = ac_e[ieA5fac, iprotate,1]
            d_prop_cnoz["cutback"] = ac_e[ieA5fac, ipcutback,1]
            d_prop_cnoz["climbstart"] = ac_e[ieA5fac,ipclimb1,1]
            d_prop_cnoz["climbend"] = ac_e[ieA5fac,ipclimbn,1]
            d_prop_cnoz["descentstart"] = ac_e[ieA5fac,ipdescent1,1]
            d_prop_cnoz["descentend"] = ac_e[ieA5fac,ipdescentn,1]
        d_prop_nozz["core_nozzle_area"] = d_prop_cnoz
        #fan nozzle
        d_prop_fnoz = Dict()
            d_prop_fnoz["static"] = ac_e[ieA7fac, ipstatic,1]
            d_prop_fnoz["rotation"] = ac_e[ieA7fac, iprotate,1]
            d_prop_fnoz["cutback"] = ac_e[ieA7fac, ipcutback,1]
            d_prop_fnoz["climbstart"] = ac_e[ieA7fac,ipclimb1,1]
            d_prop_fnoz["climbend"] = ac_e[ieA7fac,ipclimbn,1]
            d_prop_fnoz["descentstart"] = ac_e[ieA7fac,ipdescent1,1]
            d_prop_fnoz["descentend"] = ac_e[ieA7fac,ipdescentn,1]
        d_prop_nozz["fan_nozzle_area"] = d_prop_fnoz
    d_prop["Nozzles"] = d_prop_nozz
    
    #Nacelles
    d_prop_nace = Dict()
        d_prop_nace["nacelle_pylon_wetted_area_ratio"] = ac_g[igrSnace]
        d_prop_nace["nacelle_local_velocity_ratio"] = ac_g[igrVnace]        
    d_prop["Nacelles"] = d_prop_nace

    #Weight
    d_prop_weight = Dict()
        d_prop_weight["engine_access_weight_fraction"] = ac_g[igfeadd]
        d_prop_weight["pylon_weight_fraction"] = ac_g[igfpylon]
        d_prop_weight["weight_model"] = engweightmodel[ac_i[iiengwgt]]
    d_prop["Weight"] = d_prop_weight

    d_out["Propulsion"] = d_prop
    #--end Propulsion----------------

    #Output to TOML file
    #pre-proc: replace repetitive vectors with single value
    d_out = make_dict_singletons(d_out)

    #write any preamble to the toml
    #TODO: copy paste template stuff? remember to mention units

    open(datafile, "w") do io

        #uncomment for debug:
        # print_nested_dict(d_out)

      #write model dictionary out
        TOML.print(io,d_out)

        if save_output
            println("Functions for saving sized params not yet implemented.")
            if ac.sized[1]
                #TODO: sized aircraft output

            else
                @warn ac.name * " is not sized. Outputs will not be saved."
            end #if
        end #if
    end #open()
end #function

"""
    print_nested_dict(dict, indent = "    ")

Prints dictionary contents to console, including any nested `dict`s.
Useful for debugging.

"""
function print_nested_dict(dict, indent = "    ")
    for (key, value) in dict
        if isa(value, Dict)
            println("$indent$key (Type: Dict)")
            print_nested_dict(value, "$indent  ")
        else
            value_type = typeof(value)
            println("$indent$key (Type: $value_type): $value")
        end
    end
end

"""
    make_dict_singletons(dict::Dict{K, V})


replaces Dict values that are vectors with identical elements 
with a single value ("singleton") for simplicity in processing/readability in output

returns a new Dict. treats nested Dicts recursively

"""
function make_dict_singletons(dict::Dict{K, V}) where {K, V}
    new_dict = Dict{K, V}()
    #iterate over all dict entries
    for (key, value) in dict 
        #if element is a dict, call this fxn recursively
        if isa(value, Dict)  
            new_dict[key] = make_dict_singletons(value)
        #if element is a vector and has identical elements, save singleton
        #AND not ac.sized (exception needed bc ac.sized needs to be a vector for mutability)
        elseif isa(value, Vector) && all(x -> x == value[1], value) && key != "sized"
            new_dict[key] = value[1]
        else
            new_dict[key] = value
        end
    end
    return new_dict
end






# Store label names
iglabels = ["igFOpt     ", "igPFEI     ","igRange    ","igWMTO     ","igWpay     ","igWfix     ","igWfuel    ","igWfmax    ","igrWfmax   ","igWshell   ","igWwindow  ","igWinsul   ","igWfloor   ","igWcone    ","igWhbend   ","igWvbend   ","igWfuse    ","igWweb     ","igWcap     ","igWwing    ","igWebare   ","igWnace    ","igWeng     ","igWhtail   ","igWvtail   ","igWstrut   ","igxWfuse   ","igdxWfuel  ","igdxWwing  ","igdxWstrut ","igdxWhtail ","igdxWvtail ","igWinn     ","igWout     ","igdyWinn   ","igdyWout   ","igxCGfwd   ","igxCGaft   ","igfreserve ","igfpadd    ","igfseat    ","igfeadd    ","igfpylon   ","igfnace    ","igfflap    ","igfslat    ","igfaile    ","igflete    ","igfribs    ","igfspoi    ","igfwatt    ","igfhadd    ","igfvadd    ","igfapu     ","igfhpesys  ","igflgnose  ","igflgmain  ","igfstring  ","igfframe   ","igffadd    ","igWpwindow ","igWppinsul ","igWppfloor ","igNlift    ","igNland    ","igVne      ","igneng     ","igGearf    ","igfTt4CL1  ","igfTt4CLn  ","igHTRf     ","igHTRlc    ","igHTRhc    ","igrSnace   ","igrVnace   ","igrVstrut  ","igfSnace   ","igpcabin   ","igdeltap   ","iganose    ","igbtail    ","igxnose    ","igxend     ","igxblend1  ","igxblend2  ","igxshell1  ","igxshell2  ","igxconend  ","igxhbend   ","igxvbend   ","igxhtail   ","igxvtail   ","igxeng     ","igxwing    ","igxwbox    ","igxhbox    ","igxvbox    ","igxfix     ","igxapu     ","igxhpesys  ","igxlgnose  ","igdxlgmain ","igyeng     ","igzwing    ","igzhtail   ","ignfweb    ","igwfb      ","igRfuse    ","igdRfuse   ","ighfloor   ","iglambdac  ","igcabVol   ","igcosLs    ","igSstrut   ","igrpayfwd  ","igrpayaft  ","igxNP      ","igCMVf1    ","igCLMf0    ","igdepsda   ","igdCLnda   ","igdCLhdCL  ","igdCLndCL  ","igCLhspec  ","igCLhCGfwd ","igCLveout  ","igCLhmax   ","igCLvmax   ","igfCDhcen  ","igSMmin    ","igrMh      ","igrMv      ","igXaxis    ","igwbox     ","ighboxo    ","ighboxs    ","igrh       ","igwboxh    ","ighboxh    ","igrhh      ","igwboxv    ","ighboxv    ","igrhv      ","igsigfac   ","igsigskin  ","igsigbend  ","igsigcap   ","igtauweb   ","igsigstrut ","igrEshell  ","igEcap     ","igEstrut   ","igrhoskin  ","igrhobend  ","igrhocap   ","igrhoweb   ","igrhostrut ","igrhofuel  ","igrcls     ","igrclt     ","igCLhNrat  ","igSomax    ","igMomax    ","igSsmax    ","igMsmax    ","igtbcapo   ","igtbwebo   ","igtbcaps   ","igtbwebs   ","igtbcaph   ","igtbwebh   ","igtbcapv   ","igtbwebv   ","igEIco     ","igEIno     ","igGJo      ","igEIcs     ","igEIns     ","igGJs      ","igEIch     ","igEInh     ","igGJh      ","igEIcv     ","igEInv     ","igGJv      ","igtskin    ","igtcone    ","igtfweb    ","igtfloor   ","igEIhshell ","igEIhbend  ","igEIvshell ","igEIvbend  ","igGJshell  ","igGJcone   ","igfLo      ","igfLt      ","igfLn      ","igcma      ","igAR       ","igS        ","igb        ","igbo       ","igbs       ","igetas     ","iglambdat  ","iglambdas  ","igco       ","igsweep    ","igVh       ","igARh      ","igSh       ","igbh       ","igboh      ","iglambdah  ","igcoh      ","igsweeph   ","igVv       ","igARv      ","igSv       ","igbv       ","igbov      ","iglambdav  ","igcov      ","igsweepv   ","ignvtail   ","igzs       ","ighstrut   ","igAstrut   ","igcstrut   ","igfBLIw    ","igfBLIf    ","igdfan     ","igdlcomp   ","igdhcomp   ","iglnace    ","igA5       ","igA7       ","igTmetal   ","igcdefan   ","igCDgear   ","igCDspoil  ","igmuroll   ","igmubrake  ","ighobst    ","iglBFmax   ","igbmax     ","iggtocmin  ","igdBSLmax  ","igdBCBmax  ","igmofWpay  ","igmofWMTO  ","igPofWpay  ","igPofWMTO  ","igWtshaft  ","igWgen     ","igWinv     ","igWmot     ","igWfan     ","igWftank   ","igxtshaft  ","igxgen     ","igxinv     ","igxmot     ","igxfan     ","igxftank   ","igxcables  ","igWcables  ","igxcat     ","igWcat     ","igWtesys   ","igxWtesys  ","iglftank   ","igWinsftank","igxWftank  ","igRftank   ","igWc3des   ", "igdaftfan", "lnaceaft", "igfuseVol", "igneout", "igyeout", "igyeinn", "iglftankin", "igLHVfuel", "igWfburn", "igWaftfan", "igWfanGB", "igWaftfanGB", "igWrect", "igWtms"] 
function savemodel(fname, pari, parg, parm, para, pare, parpt, parmot, pargen)
    open(fname, "w") do io

        @printf(io, "# ------------------------------\n")
        @printf(io, "# Flags  - stored in pari array:\n")
        @printf(io, "# ------------------------------\n")
        for (i,val) in enumerate(pari)
            @printf(io, "pari[%d] = %d \n", i, val )
        end

        @printf(io, "# --------------------------------\n")
        @printf(io, "# Geometry - stored in parg array:\n")
        @printf(io, "# --------------------------------\n")
        for (i,val) in enumerate(parg)
            @printf(io, "parg[%d] = %20.20f # %s\n", i, val, i<291 ? iglabels[i] : "" )
        end

        @printf(io, "# --------------------------------\n")
        @printf(io, "# Mission  - stored in parm array:\n")
        @printf(io, "# --------------------------------\n")
        for (i,val) in enumerate(parm)
            @printf(io, "parm[%d] = %20.20f \n", i, val)
        end

        @printf(io, "# --------------------------------\n")
        @printf(io, "# Aero     - stored in para array:\n")
        @printf(io, "# --------------------------------\n")
        l = size(para)[1]
        m = size(para)[2]
        for i = 1:l
            @printf(io, "para[%d, :] .= [", i)
            for j = 1:m
                @printf(io, "%20.20f, ", para[i, j])
            end
            @printf(io, "]\n")
        end
        @printf(io, "# --------------------------------\n")
        @printf(io, "# Engine   - stored in pare array:\n")
        @printf(io, "# --------------------------------\n")
        l = size(pare)[1]
        m = size(pare)[2]
        for i = 1:l
            @printf(io, "pare[%d, :] .= [", i)
            for j = 1:m
                @printf(io, "%20.20f, ", pare[i, j])
            end
            @printf(io, "]\n")
        end

        @printf(io, "# ---------------------------------\n")
        @printf(io, "# Powertrain-stored in parpt array:\n")
        @printf(io, "# ---------------------------------\n")
        for (i,val) in enumerate(parpt)
            @printf(io, "parpt[%d] = %20.20f \n", i, val)
        end

        @printf(io, "# ---------------------------------\n")
        @printf(io, "# Motor   - stored in parmot array:\n")
        @printf(io, "# ---------------------------------\n")
        for (i,val) in enumerate(parmot)
            @printf(io, "parmot[%d] = %20.20f \n", i, val)
        end
        @printf(io, "# ---------------------------------\n")
        @printf(io, "# Generator-stored in pargen array:\n")
        @printf(io, "# ---------------------------------\n")
        for (i,val) in enumerate(pargen)
            @printf(io, "pargen[%d] = %20.20f \n", i, val)
        end
    end

end

function reset_regression_test(ac)
    wing = ac.wing
    htail = ac.htail
    vtail = ac.vtail
    open(joinpath(__TASOPTroot__,"../test/default_sized.jl"), "w") do io
        @printf(io, "pari = zeros(Int64, iitotal)\n")
        @printf(io, "parg = zeros(Float64, igtotal)\n")
        @printf(io, "parm = zeros(Float64, imtotal)\n")
        @printf(io, "para = zeros(Float64, (iatotal, iptotal))\n")
        @printf(io, "pare = zeros(Float64, (ietotal, iptotal))\n \n")
        @printf(io, "# ------------------------------\n")
        @printf(io, "# Flags  - stored in pari array:\n")
        @printf(io, "# ------------------------------\n")
        for (i,val) in enumerate(ac.pari)
            @printf(io, "pari[%d] = %d \n", i, val )
        end

        @printf(io, "# --------------------------------\n")
        @printf(io, "# Geometry - stored in parg array:\n")
        @printf(io, "# --------------------------------\n")
        for (i,val) in enumerate(ac.parg)
            @printf(io, "parg[%d] = %20.20f\n", i, val )
        end

        @printf(io, "# --------------------------------\n")
        @printf(io, "# Mission  - stored in parm array:\n")
        @printf(io, "# --------------------------------\n")
        l = size(ac.parm)[1]
        for i = 1:l
            @printf(io, "parm[%d] = %20.20f \n", i, ac.parm[i,1])
        end

        @printf(io, "# --------------------------------\n")
        @printf(io, "# Aero     - stored in para array:\n")
        @printf(io, "# --------------------------------\n")
        l = size(ac.para)[1]
        m = size(ac.para)[2]
        for i = 1:l
            @printf(io, "para[%d, :] .= [", i)
            for j = 1:m
                @printf(io, "%20.20f, ", ac.para[i, j,1])
            end
            @printf(io, "]\n")
        end
        @printf(io, "# --------------------------------\n")
        @printf(io, "# Engine   - stored in pare array:\n")
        @printf(io, "# --------------------------------\n")
        l = size(ac.pare)[1]
        m = size(ac.pare)[2]
        for i = 1:l
            @printf(io, "pare[%d, :] .= [", i)
            for j = 1:m
                @printf(io, "%20.20f, ", ac.pare[i, j,1])
            end
            @printf(io, "]\n")
        end
    end
    open(joinpath(__TASOPTroot__,"../test/default_structures.jl"), "w") do io
        @printf(io,"ac_test = load_default_model()\n")
        @printf(io, "# ------------------------------\n")
        @printf(io, "# Fuselage\n")
        @printf(io, "# ------------------------------\n")
        @printf(io, "fuse = ac_test.fuselage\n")
        @printf(io, "Weight = TASOPT.structures.Weight\n")
        @printf(io, "fuse.n_decks = %20.20f \n", ac.fuselage.n_decks)
        @printf(io, "fuse.shell.weight = Weight(W = %20.20f ) \n", ac.fuselage.shell.weight.W)
        @printf(io, "fuse.shell.weight.r = [ %20.20f ,0.0,0.0] \n", ac.fuselage.shell.weight.x)
        @printf(io,"fuse.window.W = %20.20f \n", ac.fuselage.window.W)
        @printf(io,"fuse.window.r = [ %20.20f ,0.0,0.0] \n", ac.fuselage.window.x)
        @printf(io,"fuse.window_W_per_length = %20.20f \n", ac.fuselage.window_W_per_length)
        @printf(io,"fuse.insulation.W = %20.20f \n", ac.fuselage.insulation.W)
        @printf(io,"fuse.insulation.r = [ %20.20f ,0.0,0.0] \n", ac.fuselage.insulation.x)
        @printf(io,"fuse.insulation_W_per_area = %20.20f \n", ac.fuselage.insulation_W_per_area)
        @printf(io,"fuse.floor.weight.W = %20.20f \n", ac.fuselage.floor.weight.W)
        @printf(io,"fuse.floor_W_per_area = %20.20f \n", ac.fuselage.floor_W_per_area)
        @printf(io,"fuse.cone.weight = Weight(W = %20.20f ) \n", ac.fuselage.cone.weight.W)
        @printf(io,"fuse.cone.weight.r = [ %20.20f ,0.0,0.0] \n", ac.fuselage.cone.weight.x)
        @printf(io,"fuse.bendingmaterial_h.weight = Weight(W = %20.20f ) \n", ac.fuselage.bendingmaterial_h.weight.W)
        @printf(io,"fuse.bendingmaterial_v.weight = Weight(W = %20.20f ) \n", ac.fuselage.bendingmaterial_v.weight.W)
        @printf(io,"fuse.weight = %20.20f \n", ac.fuselage.weight)
        @printf(io,"fuse.moment = %20.20f \n", ac.fuselage.moment)
        @printf(io,"fuse.volume = %20.20f \n", ac.fuselage.volume)
        @printf(io,"fuse.weight_frac_stringers = %20.20f \n", ac.fuselage.weight_frac_stringers)
        @printf(io,"fuse.weight_frac_frame = %20.20f \n", ac.fuselage.weight_frac_frame)
        @printf(io,"fuse.weight_frac_skin_addl = %20.20f \n", ac.fuselage.weight_frac_skin_addl)
        @printf(io,"fuse.layout.nose_radius = %20.20f \n", ac.fuselage.layout.nose_radius)
        @printf(io,"fuse.layout.tail_radius = %20.20f \n", ac.fuselage.layout.tail_radius)
        @printf(io,"fuse.layout.l_cabin_cylinder = %20.20f \n", ac.fuselage.layout.l_cabin_cylinder)
        @printf(io,"fuse.layout.x_nose = %20.20f \n", ac.fuselage.layout.x_nose)
        @printf(io,"fuse.layout.x_end = %20.20f \n", ac.fuselage.layout.x_end)
        @printf(io,"fuse.layout.x_start_cylinder = %20.20f \n", ac.fuselage.layout.x_start_cylinder)
        @printf(io,"fuse.layout.x_end_cylinder = %20.20f \n", ac.fuselage.layout.x_end_cylinder)
        @printf(io,"fuse.layout.x_pressure_shell_fwd = %20.20f \n", ac.fuselage.layout.x_pressure_shell_fwd)
        @printf(io,"fuse.layout.x_pressure_shell_aft = %20.20f \n", ac.fuselage.layout.x_pressure_shell_aft)
        @printf(io,"fuse.layout.x_cone_end = %20.20f \n", ac.fuselage.layout.x_cone_end)
        @printf(io,"fuse.bendingmaterial_h.weight.r = [ %20.20f ,0.0,0.0] \n", ac.fuselage.bendingmaterial_h.weight.x)
        @printf(io,"fuse.bendingmaterial_v.weight.r = [ %20.20f ,0.0,0.0] \n", ac.fuselage.bendingmaterial_v.weight.x)
        @printf(io,"fuse.layout.cross_section.radius = %20.20f \n", ac.fuselage.layout.cross_section.radius)
        @printf(io,"fuse.layout.cross_section.bubble_lower_downward_shift = %20.20f \n", ac.fuselage.layout.cross_section.bubble_lower_downward_shift)
        @printf(io,"fuse.layout.floor_depth = %20.20f \n", ac.fuselage.layout.floor_depth)
        @printf(io,"fuse.layout.taper_tailcone = %20.20f \n", ac.fuselage.layout.taper_tailcone)
        @printf(io,"fuse.ratio_young_mod_fuse_bending = %20.20f \n", ac.fuselage.ratio_young_mod_fuse_bending)
        @printf(io,"fuse.skin.thickness = %20.20f \n", ac.fuselage.skin.thickness)
        @printf(io,"fuse.cone.thickness = %20.20f \n", ac.fuselage.cone.thickness)
        @printf(io,"fuse.layout.thickness_webs = %20.20f \n", ac.fuselage.layout.thickness_webs)
        @printf(io,"fuse.floor.thickness = %20.20f \n", ac.fuselage.floor.thickness)
        @printf(io,"fuse.shell.EIh = %20.20f \n", ac.fuselage.shell.EIh)
        @printf(io,"fuse.bendingmaterial_h.EIh = %20.20f \n", ac.fuselage.bendingmaterial_h.EIh)
        @printf(io,"fuse.bendingmaterial_v.EIh = %20.20f \n", ac.fuselage.bendingmaterial_v.EIh)
        @printf(io,"fuse.shell.EIv = %20.20f \n", ac.fuselage.shell.EIv)
        @printf(io,"fuse.bendingmaterial_h.EIv = %20.20f \n", ac.fuselage.bendingmaterial_h.EIv)
        @printf(io,"fuse.bendingmaterial_v.EIv = %20.20f \n", ac.fuselage.bendingmaterial_v.EIv)
        @printf(io,"fuse.shell.GJ = %20.20f \n", ac.fuselage.shell.GJ)
        @printf(io,"fuse.cone.GJ = %20.20f \n", ac.fuselage.cone.GJ)
        @printf(io,"fuse.APU.W = %20.20f \n", ac.fuselage.APU.W)
        @printf(io,"fuse.APU.r = [%20.20f,0.0,0.0] \n", ac.fuselage.APU.x)
        @printf(io,"fuse.seat.W = %20.20f \n", ac.fuselage.seat.W)
        @printf(io,"fuse.fixed.W = %20.20f \n", ac.fuselage.fixed.W)
        @printf(io,"fuse.fixed.r = [%20.20f,0.0,0.0] \n", ac.fuselage.fixed.x)
        @printf(io,"fuse.HPE_sys.r = [%20.20f,0.0,0.0] \n", ac.fuselage.HPE_sys.x)
        @printf(io,"fuse.HPE_sys.W = %20.20f \n", ac.fuselage.HPE_sys.W)
        @printf(io,"fuse.added_payload.W = %20.20f \n", ac.fuselage.added_payload.W)

        @printf(io, "fuse.cabin.design_pax = %20.20f \n", ac.fuselage.cabin.design_pax)
        @printf(io, "fuse.cabin.exit_limit = %20.20f \n", ac.fuselage.cabin.exit_limit)
        @printf(io, "fuse.cabin.seat_pitch = %20.20f \n", ac.fuselage.cabin.seat_pitch)
        @printf(io, "fuse.cabin.seat_width = %20.20f \n", ac.fuselage.cabin.seat_width)
        @printf(io, "fuse.cabin.seat_height = %20.20f \n", ac.fuselage.cabin.seat_height)
        @printf(io, "fuse.cabin.aisle_halfwidth = %20.20f \n", ac.fuselage.cabin.aisle_halfwidth)
        @printf(io, "fuse.cabin.floor_distance = %20.20f \n", ac.fuselage.cabin.floor_distance)
        @printf(io, "fuse.cabin.cabin_width_main = %20.20f \n", ac.fuselage.cabin.cabin_width_main)
        @printf(io, "fuse.cabin.cabin_width_top = %20.20f \n", ac.fuselage.cabin.cabin_width_top)
        @printf(io, "fuse.cabin.seats_abreast_main = %d \n", ac.fuselage.cabin.seats_abreast_main)
        @printf(io, "fuse.cabin.seats_abreast_top = %d \n", ac.fuselage.cabin.seats_abreast_top)
        @printf(io, "fuse.cabin.floor_angle_main = %20.20f \n", ac.fuselage.cabin.floor_angle_main)
        @printf(io, "fuse.cabin.floor_angle_top = %20.20f \n", ac.fuselage.cabin.floor_angle_top)

        @printf(io, "# ------------------------------\n")
        @printf(io, "# Wing\n")
        @printf(io, "# ------------------------------\n")
        @printf(io, "wing = ac_test.wing\n")
        @printf(io, "wing.inboard.webs.weight = Weight(W = %20.20f) \n", wing.inboard.webs.weight.W)
        @printf(io, "wing.outboard.webs.weight = Weight(W = %20.20f) \n", wing.outboard.webs.weight.W)
        @printf(io, "wing.inboard.caps.weight = Weight(W = %20.20f) \n", wing.inboard.caps.weight.W)
        @printf(io, "wing.outboard.caps.weight = Weight(W = %20.20f) \n", wing.outboard.caps.weight.W)
        @printf(io, "wing.inboard.caps.material = TASOPT.materials.StructuralAlloy(\"TASOPT-Al\",
        max_avg_stress = 1.1,
        safety_factor = 1.5)\n")
        @printf(io, "wing.outboard.caps.material = TASOPT.materials.StructuralAlloy(\"TASOPT-Al\",
        max_avg_stress = 1.1,
        safety_factor = 1.5)\n")
        @printf(io, "wing.inboard.caps.material = TASOPT.materials.StructuralAlloy(\"TASOPT-Al\",
        max_avg_stress = 1.1,
        safety_factor = 1.5)\n")
        @printf(io, "wing.inboard.webs.material = TASOPT.materials.StructuralAlloy(\"TASOPT-Al\",
        max_avg_stress = 1.1,
        safety_factor = 1.5)\n")
        @printf(io, "wing.outboard.webs.material = TASOPT.materials.StructuralAlloy(\"TASOPT-Al\",
        max_avg_stress = 1.1,
        safety_factor = 1.5)\n")


        @printf(io, "wing.weight = %20.20f \n", wing.weight)
        @printf(io, "wing.strut.weight = %20.20f \n", wing.strut.weight)
        @printf(io, "wing.dxW = %20.20f \n", wing.dxW)
        @printf(io, "wing.strut.dxW = %20.20f \n", wing.strut.dxW)
        @printf(io, "wing.inboard.weight = %20.20f \n", wing.inboard.weight)
        @printf(io, "wing.outboard.weight = %20.20f \n", wing.outboard.weight)
        @printf(io, "wing.inboard.dyW = %20.20f \n", wing.inboard.dyW)
        @printf(io, "wing.outboard.dyW = %20.20f \n", wing.outboard.dyW)
        @printf(io, "wing.weight_frac_flap = %20.20f \n", wing.weight_frac_flap)
        @printf(io, "wing.weight_frac_slat = %20.20f \n", wing.weight_frac_slat)
        @printf(io, "wing.weight_frac_ailerons = %20.20f \n", wing.weight_frac_ailerons)
        @printf(io, "wing.weight_frac_leading_trailing_edge = %20.20f \n", wing.weight_frac_leading_trailing_edge)
        @printf(io, "wing.weight_frac_ribs = %20.20f \n", wing.weight_frac_ribs)
        @printf(io, "wing.weight_frac_spoilers = %20.20f \n", wing.weight_frac_spoilers)
        @printf(io, "wing.weight_frac_attachments = %20.20f \n", wing.weight_frac_attachments)
        @printf(io, "wing.strut.local_velocity_ratio = %20.20f \n", wing.strut.local_velocity_ratio)
        @printf(io, "wing.layout.x = %20.20f \n", wing.layout.x)
        @printf(io, "wing.layout.box_x = %20.20f \n", wing.layout.box_x)
        @printf(io, "wing.layout.z = %20.20f \n", wing.layout.z)
        @printf(io, "wing.strut.cos_lambda = %20.20f \n", wing.strut.cos_lambda)
        @printf(io, "wing.strut.S = %20.20f \n", wing.strut.S)
        @printf(io, "wing.layout.spar_box_x_c = %20.20f \n", wing.layout.spar_box_x_c)
        @printf(io, "wing.inboard.cross_section.width_to_chord = %20.20f \n", wing.inboard.cross_section.width_to_chord)
        @printf(io, "wing.inboard.cross_section.web_to_box_height = %20.20f \n", wing.inboard.cross_section.web_to_box_height)
        @printf(io, "wing.inboard.cross_section.thickness_to_chord = %20.20f \n", wing.inboard.cross_section.thickness_to_chord)
        @printf(io, "wing.outboard.cross_section.thickness_to_chord = %20.20f \n", wing.outboard.cross_section.thickness_to_chord)
        @printf(io, "wing.layout.max_span = %20.20f \n", wing.layout.max_span)
        @printf(io, "wing.strut.thickness_to_chord = %20.20f \n", wing.strut.thickness_to_chord)
        @printf(io, "wing.strut.z = %20.20f \n", wing.strut.z)
        @printf(io, "wing.outboard.moment = %20.20f \n", wing.outboard.moment)
        @printf(io, "wing.outboard.max_shear_load = %20.20f \n", wing.outboard.max_shear_load)
        @printf(io, "wing.outboard.GJ = %20.20f \n", wing.outboard.GJ)
        @printf(io, "wing.outboard.EI[4] = %20.20f \n", wing.outboard.EI[4])
        @printf(io, "wing.outboard.EI[1] = %20.20f \n", wing.outboard.EI[1])
        @printf(io, "wing.outboard.caps.thickness = %20.20f \n", wing.outboard.caps.thickness)
        @printf(io, "wing.inboard.moment = %20.20f \n", wing.inboard.moment)
        @printf(io, "wing.inboard.max_shear_load = %20.20f \n", wing.inboard.max_shear_load)
        @printf(io, "wing.inboard.GJ = %20.20f \n", wing.inboard.GJ)
        @printf(io, "wing.inboard.EI[4] = %20.20f \n", wing.inboard.EI[4])
        @printf(io, "wing.inboard.EI[1] = %20.20f \n", wing.inboard.EI[1])
        @printf(io, "wing.inboard.caps.thickness = %20.20f \n", wing.inboard.caps.thickness)
        @printf(io, "wing.inboard.webs.thickness = %20.20f \n", wing.inboard.webs.thickness)
        @printf(io, "wing.outboard.webs.thickness = %20.20f \n", wing.outboard.webs.thickness)
        @printf(io, "wing.layout.S = %20.20f \n", wing.layout.S)
        @printf(io, "wing.layout.root_span = %20.20f \n", wing.layout.root_span)
        @printf(io, "wing.layout.ηs = %20.20f \n", wing.layout.ηs)
        @printf(io, "wing.inboard.λ = %20.20f \n", wing.inboard.λ)
        @printf(io, "wing.outboard.λ = %20.20f \n", wing.outboard.λ)
        @printf(io, "wing.layout.root_chord = %20.20f \n", wing.layout.root_chord)
        @printf(io, "wing.layout.span= %20.20f \n", wing.layout.span)
        @printf(io, "wing.layout.sweep = %20.20f \n", wing.layout.sweep)
        @printf(io, "wing.layout.AR = %20.20f \n", wing.layout.AR)
        @printf(io, "wing.fuse_lift_carryover = %20.20f \n", wing.fuse_lift_carryover)
        @printf(io, "wing.tip_lift_loss = %20.20f \n", wing.tip_lift_loss)
        @printf(io, "wing.inboard.co = %20.20f \n", wing.inboard.co)
        @printf(io, "wing.outboard.co = %20.20f \n", wing.outboard.co)
        @printf(io, "wing.mean_aero_chord = %20.20f \n", wing.mean_aero_chord)

        @printf(io, "# ------------------------------\n")
        @printf(io,"# Htail\n")
        @printf(io, "# ------------------------------\n")
        @printf(io,"htail = ac_test.htail\n")
        @printf(io, "htail.weight = %20.20f \n", htail.weight)
        @printf(io, "htail.dxW = %20.20f \n", htail.dxW)
        @printf(io, "htail.weight_fraction_added = %20.20f \n", htail.weight_fraction_added)
        @printf(io, "htail.layout.box_x = %20.20f \n", htail.layout.box_x)
        @printf(io, "htail.layout.z = %20.20f \n", htail.layout.z)
        @printf(io, "htail.downwash_factor = %20.20f \n", htail.downwash_factor)
        @printf(io, "htail.CL_max_fwd_CG = %20.20f \n", htail.CL_max_fwd_CG)
        @printf(io, "htail.CL_max = %20.20f \n", htail.CL_max)
        @printf(io, "htail.SM_min = %20.20f \n", htail.SM_min)
        @printf(io, "htail.layout.x = %20.20f \n", htail.layout.x)
        @printf(io, "htail.outboard.cross_section.thickness_to_chord = %20.20f \n", htail.outboard.cross_section.thickness_to_chord)
        @printf(io, "htail.move_wingbox = %20.20f \n", htail.move_wingbox)
        @printf(io, "htail.CL_CLmax = %20.20f \n", htail.CL_CLmax)
        @printf(io, "htail.size = %20.20f \n", htail.size)
        @printf(io, "htail.volume = %20.20f \n", htail.volume)
        @printf(io, "htail.outboard.GJ = %20.20f \n", htail.outboard.GJ)
        @printf(io, "htail.outboard.EI[4] = %20.20f \n", htail.outboard.EI[4])
        @printf(io, "htail.outboard.EI[1] = %20.20f \n", htail.outboard.EI[1])
        @printf(io, "htail.layout.sweep = %20.20f \n", htail.layout.sweep)
        @printf(io, "htail.layout.root_chord = %20.20f \n", htail.layout.root_chord)
        @printf(io, "htail.outboard.λ = %20.20f \n", htail.outboard.λ)
        @printf(io, "htail.layout.root_span = %20.20f \n", htail.layout.root_span)
        @printf(io, "htail.layout.span = %20.20f \n", htail.layout.span)
        @printf(io, "htail.layout.AR = %20.20f \n", htail.layout.AR)
        @printf(io, "htail.layout.S = %20.20f \n", htail.layout.S)
        @printf(io, "htail.outboard.cross_section.width_to_chord = %20.20f \n", htail.outboard.cross_section.width_to_chord)
        @printf(io, "htail.outboard.cross_section.web_to_box_height = %20.20f \n", htail.outboard.cross_section.web_to_box_height)
        @printf(io, "htail.layout.ηs = htail.layout.root_span/htail.layout.span \n")
        @printf(io, "htail.strut.cos_lambda = %20.20f \n", htail.strut.cos_lambda)
        @printf(io, "htail.inboard.moment = %20.20f \n", htail.inboard.moment)
        @printf(io, "htail.outboard.moment = %20.20f \n", htail.outboard.moment)
        @printf(io, "htail.inboard.max_shear_load = %20.20f \n", htail.inboard.max_shear_load)
        @printf(io, "htail.outboard.max_shear_load = %20.20f \n", htail.outboard.max_shear_load)
        @printf(io, "htail.outboard.webs.thickness = %20.20f \n", htail.outboard.webs.thickness)
        @printf(io, "htail.inboard.webs.weight.W = %20.20f \n", htail.inboard.webs.weight.W)
        @printf(io, "htail.inboard.caps.weight.W = %20.20f \n", htail.inboard.caps.weight.W)
        @printf(io, "htail.inboard.webs.thickness = %20.20f \n", htail.inboard.webs.thickness)
        @printf(io, "htail.inboard.caps.thickness = %20.20f \n", htail.inboard.caps.thickness)
        @printf(io, "htail.outboard.webs.thickness = %20.20f \n", htail.inboard.webs.thickness)
        @printf(io, "htail.outboard.caps.thickness = %20.20f \n", htail.inboard.caps.thickness)
        @printf(io, "htail.inboard.GJ = %20.20f \n", htail.inboard.GJ)
        @printf(io, "htail.outboard.co = htail.layout.root_chord*htail.inboard.λ \n")
        @printf(io, "htail.inboard.co = htail.layout.root_chord \n")
        
        @printf(io, "# ------------------------------\n")
        @printf(io,"# Vtail\n")
        @printf(io, "# ------------------------------\n")
        @printf(io,"vtail = ac_test.vtail\n")
        @printf(io, "vtail.weight = %20.20f \n", vtail.weight)
        @printf(io, "vtail.dxW = %20.20f \n", vtail.dxW)
        @printf(io, "vtail.weight_fraction_added = %20.20f \n", vtail.weight_fraction_added)
        @printf(io, "vtail.layout.box_x = %20.20f \n", vtail.layout.box_x)
        @printf(io, "vtail.CL_max = %20.20f \n", vtail.CL_max)
        @printf(io, "vtail.layout.x = %20.20f \n", vtail.layout.x)
        @printf(io, "vtail.outboard.cross_section.thickness_to_chord = %20.20f \n", vtail.outboard.cross_section.thickness_to_chord)
        @printf(io, "vtail.ntails = %20.20f \n", vtail.ntails)
        @printf(io, "vtail.volume = %20.20f \n", vtail.volume)
        @printf(io, "vtail.outboard.GJ = %20.20f \n", vtail.outboard.GJ)
        @printf(io, "vtail.outboard.EI[4] = %20.20f \n", vtail.outboard.EI[4])
        @printf(io, "vtail.outboard.EI[1] = %20.20f \n", vtail.outboard.EI[1])
        @printf(io, "vtail.layout.sweep = %20.20f \n", vtail.layout.sweep)
        @printf(io, "vtail.layout.root_chord = %20.20f \n", vtail.layout.root_chord)
        @printf(io, "vtail.outboard.λ = %20.20f \n", vtail.outboard.λ)
        @printf(io, "vtail.layout.span = %20.20f \n", vtail.layout.span)
        @printf(io, "vtail.layout.AR = %20.20f \n", vtail.layout.AR)
        @printf(io, "vtail.layout.S = %20.20f \n", vtail.layout.S)
        @printf(io, "vtail.size = %20.20f \n", vtail.size)
        @printf(io, "vtail.dxW = %20.20f \n", vtail.dxW)
        @printf(io, "vtail.outboard.cross_section.width_to_chord = %20.20f \n", vtail.outboard.cross_section.width_to_chord)
        @printf(io, "vtail.outboard.cross_section.web_to_box_height = %20.20f \n", vtail.outboard.cross_section.web_to_box_height)
        @printf(io, "vtail.layout.ηs = vtail.layout.root_span/vtail.layout.span \n")
        @printf(io, "vtail.strut.cos_lambda = %20.20f \n", vtail.strut.cos_lambda)
        @printf(io, "vtail.inboard.moment = %20.20f \n", vtail.inboard.moment)
        @printf(io, "vtail.outboard.moment = %20.20f \n", vtail.outboard.moment)
        @printf(io, "vtail.inboard.max_shear_load = %20.20f \n", vtail.inboard.max_shear_load)
        @printf(io, "vtail.outboard.max_shear_load = %20.20f \n", vtail.outboard.max_shear_load)
        @printf(io, "vtail.outboard.webs.thickness = %20.20f \n", vtail.outboard.webs.thickness)
        @printf(io, "vtail.inboard.webs.weight.W = %20.20f \n", vtail.inboard.webs.weight.W)
        @printf(io, "vtail.inboard.caps.weight.W = %20.20f \n", vtail.inboard.caps.weight.W)
        @printf(io, "vtail.inboard.webs.thickness = %20.20f \n", vtail.inboard.webs.thickness)
        @printf(io, "vtail.inboard.caps.thickness = %20.20f \n", vtail.inboard.caps.thickness)
        @printf(io, "vtail.outboard.webs.thickness = %20.20f \n", vtail.inboard.webs.thickness)
        @printf(io, "vtail.outboard.caps.thickness = %20.20f \n", vtail.inboard.caps.thickness)
        @printf(io, "vtail.inboard.GJ = %20.20f \n", vtail.inboard.GJ)
        @printf(io, "vtail.outboard.co = vtail.layout.root_chord*vtail.inboard.λ \n")
        @printf(io, "vtail.inboard.co = vtail.layout.root_chord \n")
    end
end



#NOTE:  the following functions are not currently used in the codebase, 
#       but are retained for potential future use with saving models in a less-manual, yet human-readable way
"""
    fix_dict_for_toml(dict::Dict)

Prepares a `dict` for output into .toml by:
    - replacing any multi-dimensional arrays with nested arrays (TOML library compatibility restriction),
    - and recursively applying this for any nested `dict`s.
    - Also replaces structs with dictionaries, then recursively applying this function.
"""
function fix_dict_for_toml(dict::Dict)
    #deep copy of dictionary
    dict = deepcopy(dict)
    for (key, value) in dict
        #if it's an array, convert to nested vectors
        if isa(value, Array) && ndims(value) >= 2
            dict[key] = array2nestedvectors(value)
        #if it's another dict, pass it through this fxn
        elseif isa(value, Dict)
            dict[key] = fix_dict_for_toml(value)
        #if it's a TOML-compatible value, assign directly
        #TODO: convert structs to dictionaries for output 
        elseif (value isa TOML.Internals.Printer.TOMLValue)
            dict[key] = value
        else #if not TOML-compatible
            @warn "TOML cannot output the following `aircraft` field and will skip it: "*String(key)
            # println("Field value: ", value)
            delete!(dict, key)
        end
    end
    return dict
end

function struct2dict(obj)
    type = typeof(obj)
    # println("T = $type")
    d = Dict()
    for key ∈ fieldnames(type)
        # println("field = $key")
        if fieldtype(type, key) <: Union{AbstractArray,Number,String}
            # println("pushing...")
            push!(d, key => getfield(obj, key))
        else
            # println("recursing...")
            push!(d, key => struct2dict(getfield(obj, key)))
        end
    end
    return d
end