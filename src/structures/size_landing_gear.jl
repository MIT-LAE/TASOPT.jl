function landing_gear_size!(ac)
    parg = ac.parg
    pare = ac.pared
    fuse = ac.fuselage
    landing_gear = ac.landing_gear
    nose_gear = landing_gear.nose_gear
    main_gear = landing_gear.main_gear
    
    WMTO = parg[igWMTO]

    if lowercase(landing_gear.model) == "mass_fractions"
        Wlgnose = WMTO * nose_gear.overall_mass_fraction
        Wlgmain = WMTO * main_gear.overall_mass_fraction
        
    elseif lowercase(landing_gear.model) == "historical_correlations"
        #Model from Raymer (2012). Aircraft Design: A Conceptual Approach
        load_factor = 4.5 #Assumed 3 gear load factor times 1.5 ultimate load factor
        Vstall = pare[ieu0,iprotate] #Stall speed, for correlation

        #Calculate landing gear length
        tailstrike_angle = landing_gear.tailstrike_angle
        x_lg = main_gear.weight.r[1] #Main gear longitudinal position
        l_tailstrike = (fuse.layout.x_end - x_lg) * tan(tailstrike_angle) - 2*fuse.layout.radius

        yeng = parg[igetas] * parg[igb ] / 2
        ground_clearance = landing_gear.engine_ground_clearance
        dihedral_angle = landing_gear.wing_dihedral_angle
        Deng = parg[igdfan]
        l_clearance = ground_clearance + Deng - yeng * tan(dihedral_angle)
        lgmain_length = max(l_clearance, l_tailstrike)

        lgmain_nshock_struts = main_gear.number_struts
        lgmain_nwheels = main_gear.number_struts * main_gear.wheels_per_strut

        lgnose_length = lgmain_length #Assume same length as main gear #TODO drop this assumption
        lgnose_nwheels = nose_gear.number_struts * nose_gear.wheels_per_strut

        Wlgmain = 0.0106 * (WMTO * lb_N)^0.888 * (load_factor)^0.25 * (lgmain_length / in_to_m)^0.4 * (lgmain_nwheels)^0.321 * 
            (lgmain_nshock_struts)^(-0.5) * (Vstall/ kts_to_mps)^0.1 / lb_N #Model from Raymer (2012)

        Wlgnose = 0.032 * (WMTO * lb_N)^0.646 * (load_factor)^0.2 * (lgnose_length / in_to_m)^0.5 * 
            (lgnose_nwheels)^0.45 / lb_N #Model from Raymer (2012)

        #Store lengths
        nose_gear.length = lgnose_length
        main_gear.length = lgmain_length
    end
    
    nose_gear.weight = Weight(W = Wlgnose, x = nose_gear.weight.r[1])
    main_gear.weight = Weight(W = Wlgmain, x = parg[igxCGaft] + main_gear.distance_CG_to_landing_gear)
end