"""
    size_landing_gear!(ac)
Function to calculate the landing gear mass and geometric properties.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: structure with aircraft parameters
    **Outputs:**
    No direct outputs; parameters in `ac` are modified.
"""
function size_landing_gear!(ac)
    parg, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear = unpack_ac_components(ac)
    nose_gear = landing_gear.nose_gear
    main_gear = landing_gear.main_gear

    WMTO = parg[igWMTO]

    #Two models of landing gear mass currently supported
    if lowercase(landing_gear.model) == "mass_fractions"
        #This model simply uses fixed fractions of the MTOW
        Wlgnose = WMTO * nose_gear.overall_mass_fraction
        Wlgmain = WMTO * main_gear.overall_mass_fraction

    elseif lowercase(landing_gear.model) == "historical_correlations"
        #Model from Raymer (2012). Aircraft Design: A Conceptual Approach, based on historical-data correlations
        load_factor = 4.5 #Assumed 3 gear load factor times 1.5 ultimate load factor
        Vstall = ac.pare[ieu0,iprotate, 1] #Stall speed, for correlation

        #Calculate landing gear length
        #First calculate required length to avoid tailstrike
        tailstrike_angle = landing_gear.tailstrike_angle
        x_lg = parg[igxCGaft] + main_gear.distance_CG_to_landing_gear #Main gear longitudinal position
        l_tailstrike = (fuse.layout.x_end - x_lg) * tan(tailstrike_angle) - 2*fuse.layout.radius

        #Next, find the nose length that gives desired engine ground clearance
        yeng = wing.layout.ηs * wing.layout.span / 2
        ground_clearance = landing_gear.engine_ground_clearance
        dihedral_angle = landing_gear.wing_dihedral_angle
        Deng = parg[igdfan]
        l_clearance = ground_clearance + Deng - yeng * tan(dihedral_angle)

        #Calculate nose gear length
        lgnose_length = max(l_clearance, l_tailstrike) #Choose the maximum of the tailstrike and clearance lengths

        lgmain_nshock_struts = main_gear.number_struts
        lgmain_nwheels = main_gear.number_struts * main_gear.wheels_per_strut

        lgmain_length = lgnose_length + main_gear.weight.y * tan(dihedral_angle) #Main gear length
        lgnose_nwheels = nose_gear.number_struts * nose_gear.wheels_per_strut

        Wlgmain = 0.0106 * (WMTO * lb_N)^0.888 * (load_factor)^0.25 * (lgmain_length / in_to_m)^0.4 * (lgmain_nwheels)^0.321 * 
            (lgmain_nshock_struts)^(-0.5) * (Vstall/ kts_to_mps)^0.1 / lb_N #Model from Raymer (2012)

        Wlgnose = 0.032 * (WMTO * lb_N)^0.646 * (load_factor)^0.2 * (lgnose_length / in_to_m)^0.5 * 
            (lgnose_nwheels)^0.45 / lb_N #Model from Raymer (2012)

        #Store lengths and mass fractions
        nose_gear.length = lgnose_length
        main_gear.length = lgmain_length

        nose_gear.overall_mass_fraction = Wlgnose / WMTO
        main_gear.overall_mass_fraction = Wlgmain / WMTO
    end

    #Store mass and location of landing gear
    nose_gear.weight = Weight(W = Wlgnose, x = nose_gear.weight.r[1])
    main_gear.weight = Weight(W = Wlgmain, x = parg[igxCGaft] + main_gear.distance_CG_to_landing_gear, 
        y = wing.span / 2 * main_gear.y_offset_halfspan_fraction)
end 