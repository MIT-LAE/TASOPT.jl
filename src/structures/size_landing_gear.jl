function landing_gear_size!(ac)
    parg = ac.parg
    landing_gear = ac.landing_gear
    nose_gear = landing_gear.nose_gear
    main_gear = landing_gear.main_gear
    
    WMTO = parg[igWMTO]
    
    nose_gear.weight = Weight(W = WMTO * parg[igflgnose], x = parg[igxlgnose])

    delxw = parg[igxwing] - parg[igxwbox]
    main_gear.weight = Weight(W = WMTO * parg[igflgmain], x = parg[igxwbox] + delxw + parg[igdxlgmain])
    # load_factor = 1.5*2.5
    # deg = 3.5 * (pi/180)
    # lgmain_length = max(1.1 * parg[igdfan]/2.0, (parg[igxend]/ 2.0) * tan(deg))
    # lgmain_nwheels = 6
    # lgmain_nshock_struts = 2
    # Vstall = pare[ieu0,iprotate]

    # lgnose_length = lgmain_length
    # lgnose_nwheels = 2.0

    # Wlgmain = 0.0106 *(parg[igWMTO] / lbf_to_N)^0.888 * (load_factor)^0.25 * (lgmain_length / in_to_m)^0.4 * (lgmain_nwheels)^0.321 * (lgmain_nshock_struts)^(-0.5) * 
    # (Vstall/ kts_to_mps)^0.1

    # parg[igflgmain] = (Wlgmain * 4.44822)/parg[igWMTO]

    # Wlgnose = 2.24*0.032 * (parg[igWMTO] / lbf_to_N)^0.646 * (load_factor)^0.2 * (lgnose_length / in_to_m)^0.5 * (lgnose_nwheels)^0.45

    # parg[igflgnose] = (Wlgnose * 4.44822)/parg[igWMTO]
end