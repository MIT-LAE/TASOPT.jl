using Printf
"""
    _size_aircraft!(ac; itermax=35,
    wrlx1=0.5, wrlx2=0.9, wrlx3=0.5, initwgt=false, initializes_engine=true, 
    iairf=1, Ldebug=false, printiter=true)

Main weight sizing function. Calls on various sub-functions to calculate weight of fuselage, wings, tails, etc.,
and iterates until the MTOW converges to within a specified tolerance. Formerly, `wsize()`.

!!! warning
    `_size_aircraft!()` Should not be called directly by users, instead use `size_aircraft!()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - Array of flags that control design choices - fuel types, where to store fuel, etc.
    - Geometric and structural parameters - dimensions primarily
    - Aerodynamic parameters - CL, CD, KE dissipation, etc.
    - Mission-specific parameters - alt, mach, P, T etc.
    - Engine-specific parameters 

    **Outputs:**
    - No explicit outputs. Computed quantities are saved to `par` arrays of `aircraft` model.
"""
function _size_aircraft!(ac; itermax=35,
    wrlx1=0.5, wrlx2=0.9, wrlx3=0.5, initwgt=false, initializes_engine=true, 
    iairf=1, Ldebug = false, printiter=true)

    # Unpack data storage arrays and components
    imission = 1 #Design mission
    parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, eng, landing_gear  = unpack_ac(ac, imission)

    # Initialize variables
    ichoke5 = zeros(iptotal)
    ichoke7 = zeros(iptotal)

    # Weight convergence settings
    tolerW = 1.0e-8
    errw = 1.0
    fsum = 0.0
    ifirst = true

    # Extract flags
    ifuel = options.ifuel

    # Calculate sea level temperature for takeoff conditions
    altTO = parm[imaltTO]
    T_std, _, _, _, _ = atmos(altTO / 1e3)
    ΔTatmos = parm[imT0TO] - T_std
    parm[imDeltaTatm] = ΔTatmos

    # Set atmospheric conditions for different flight phases
    set_ambient_conditions!(ac, ipcruise1)
    set_ambient_conditions!(ac, iprotate, 0.25)
    set_ambient_conditions!(ac, ipclimbn)

    # Calculate fuselage boundary layer development
    fuselage_drag!(fuse, parm, para, ipcruise1)

    # Extract and set constant values for all mission points
    broadcast_fuselage_drag!(para, ipcruise1)

    #Calculate fuel lower heating value for PFEI
    parg[igLHVfuel] = fuelLHV(ifuel)

    # Unpack and set design mission parameters
    Rangetot = parm[imRange]
    Wpay = parm[imWpay]

    Wpaymax = parg[igWpaymax] # Max payload in offdesign mission (used to size structures)
    # if Wpay or Wpaymax is unset
    if (Wpaymax == 0)
        println("Max payload weight was not set, setting Wpaymax = Wpay")
        Wpaymax = parg[igWpaymax] = max(Wpay, Wpaymax)
    end

    parg[igRange] = Rangetot
    parg[igWpay] = Wpay

    # Extract weight fractions and factors
    fwadd = wing_additional_weight(wing)
    flgnose = landing_gear.nose_gear.overall_mass_fraction
    flgmain = landing_gear.main_gear.overall_mass_fraction
    fLo =  wing.fuse_lift_carryover
    fLt =  wing.tip_lift_loss

    # Extract layout parameters
    xhbox = htail.layout.box_x
    xvbox = vtail.layout.box_x
    xeng = parg[igxeng]

    # Fuselage-bending inertial relief factors
    rMh = parg[igrMh]
    rMv = parg[igrMv]

    # Tail surface taper ratios and strut parameters
    λhs = λvs = 1.0
    tohstrut = 0.05


    # Load factors and dynamic pressure
    Nlift = parg[igNlift]
    Nland = parg[igNland]
    Vne = parg[igVne]
    qne = 0.5 * ρSL * Vne^2

    # Engine parameters
    neng = parg[igneng]
    yeng = parg[igyeng]

    # Fuel tank parameters (initialized here, updated in loop)
    nftanks = fuse_tank.tank_count
    xfuel = ltank = 0.0

    # Set up fuel storage parameters (wing vs fuselage)
    setup_fuel_storage!(options, fuse, fuse_tank, parg, pare)

    # -------------------------------------------------------
    ## Initial guess section [Section 3.2 of TASOPT docs]
    # -------------------------------------------------------
    if initwgt == 0
        # Initialize weights, geometry, and parameters for first sizing run
        initialize_sizing_loop!(ac)

        # Extract working variables set by initialization
        Sh = htail.layout.S
        Sv = vtail.layout.S
        bv = sqrt(Sv * vtail.layout.AR)
        cbox = 0.0

        # Variables needed for first-iteration WMTO calculation in the loop
        Whtail = htail.weight
        Wvtail = vtail.weight
        Wwing = wing.weight
        Wstrut = wing.strut.weight
        ffuel = parg[igWfuel] / parg[igWMTO]

    else #Second iteration onwards use previously calculated values

        # Extract layout parameters
        bv = vtail.layout.span
        cbox = wing.layout.root_chord * wing.inboard.cross_section.width_to_chord
        xwing, xhtail, xvtail = wing.layout.x, htail.layout.x, vtail.layout.x
        xhbox, xvbox = htail.layout.box_x, vtail.layout.box_x

        # Extract weights
        Whtail, Wvtail = htail.weight, vtail.weight
        Wwing, Wstrut = wing.weight, wing.strut.weight
        Weng = parg[igWeng]
        Winn, Wout = wing.inboard.weight, wing.outboard.weight
        Wtesys = parg[igWtesys]

        # Extract weight moments
        dyWinn, dyWout = wing.inboard.dyW, wing.outboard.dyW

        # Calculate weight fractions
        WMTO = parg[igWMTO]
        feng = Weng / WMTO
        ffuel = parg[igWfuel] / WMTO

        # Extract other parameters
        Sh, Sv = htail.layout.S, vtail.layout.S

    end

    # Initialize previous weight iterations and convergence flag
    # 1st-previous to 3rd previous iteration weight for convergence criterion
    WMTO1 = 0.0
    WMTO2 = 0.0
    WMTO3 = 0.0 
    Lconv = false

    # Initialize wing layout parameters for first iteration
    wing.layout.span= wing.layout.S = 0.0

    # Initialize choke flags for all mission points
    ichoke5 = ichoke7 = zeros(Int, iptotal)

    # -------------------------------------------------------    
    #                   Weight loop
    # -------------------------------------------------------    

    @inbounds for iterw = 1:itermax
        if iterw == itermax
            @warn "Reached max iterations in weight sizing loop!"
        end

        # Set relaxation factor
        rlx = if iterw <= (initwgt == 0 ? 5 : 2)
            wrlx1  # Under-relax first few iterations
        elseif iterw >= 0.75 * itermax
            wrlx3  # Under-relax final iterations
        else
            wrlx2  # Default relaxation
        end

        # Max tail lifts at maneuver qne
        Lhmax = qne * Sh * htail.CL_max
        Lvmax = qne * Sv * vtail.CL_max / vtail.ntails

        # Max Δp (fuselage pressure) at end of cruise-climb
        wcd = para[iafracW, ipcruisen] / para[iafracW, ipcruise1]
        Δp = parg[igpcabin] - pare[iep0, ipcruise1] * wcd
        parg[igdeltap] = Δp

       # Engine weight mounted on tailcone, if any
        if compare_strings(options.opt_engine_location, "wing") # Eng on "wing" or aft "fuselage"
            Wengtail = 0.0
            Waftfuel = 0.0
        elseif compare_strings(options.opt_engine_location, "fuselage")
            Wengtail = parg[igWeng]
        else
            error("Engine location provided is \"$options.opt_engine_location\". Engine position can only be:
                        > \"wing\" - engines under wing
                        > \"fuselage\" - engines on aft fuselage")
        end

        # Extract relevant weights and positions
        Whtail, Wvtail = htail.weight, vtail.weight
        xhtail, xvtail, xwing = htail.layout.x, vtail.layout.x, wing.layout.x
        xeng = parg[igxeng]
        Wtesys = parg[igWtesys]
        nftanks = fuse_tank.tank_count
        
        if !(options.has_wing_fuel) #fuselage fuel store
            tank_placement = fuse_tank.placement
            Wftank_single = parg[igWftank] / nftanks #Weight of a single tank
            ltank = parg[iglftank] #length of a fuel tank

            #Calculate the weight of the fuel near the tail depending on the tank location
            if tank_placement == "rear"
                Waftfuel = parg[igWfuel]
                xftank_fuse = parg[igxftankaft] #assumed location of tank for fuselage sizing
            elseif tank_placement == "both"
                Waftfuel = parg[igWfuel] / 2.0
                xftank_fuse = parg[igxftankaft]

            elseif tank_placement == "front" #The case when the fuel is at the front is treated specially
                #The code assumes that the fuel is located at the back for the purpose of sizing of the symmetric fuselage
                Waftfuel = parg[igWfuel]
                xftank_fuse = fuse.layout.x_end - parg[igxftank]
            end
        else
            tank_placement = ""
            xftank_fuse = 0.0
            Wftank_single = 0.0
        end
        #Note that fuselage is sized for a maximum payload weight in off-design missions
        parg[igcabVol] = fusew!(fuse, Nland, Wpaymax, Wengtail, 
             nftanks,
            Waftfuel,  Wftank_single, ltank, xftank_fuse, tank_placement,
             Δp,
            Whtail, Wvtail, rMh, rMv, Lhmax, Lvmax,
            bv, vtail.outboard.λ, vtail.ntails,
            xhtail, xvtail,
            xwing, wing.layout.box_x, cbox,
            xeng)

        # Use cabin volume to get actual buoyancy weight
        ρcab = max(parg[igpcabin], pare[iep0, ipcruise1]) / (RSL * TSL)
        WbuoyCR = (ρcab - pare[ierho0, ipcruise1]) * gee * parg[igcabVol]

        if (iterw == 1 && initwgt == 0)

            feng = 0.08
            fsum = feng + ffuel + fuse.HPE_sys.W + flgnose + flgmain
            WMTO = (Wpay + fuse.weight + Wwing + Wstrut + Whtail + Wvtail) / (1.0 - fsum)

            Weng, Wfuel = WMTO .* [feng, ffuel]
            parg[igWMTO] = WMTO
            parg[igWeng] = Weng
            parg[igWfuel] = Wfuel
            if printiter
                println("Wfuel initial = ", (ffuel * WMTO))
            end

        else
            # Call a better update_weights! function
            update_WMTO!(ac, rlx)

            parm[imWTO] = parg[igWMTO]
            parm[imWfuel] = parg[igWfuel]

        end
        # this calculated WMTO is the design-mission WTO
        parm[imWTO] = parg[igWMTO]
        # Convergence tests
        WMTO = parg[igWMTO]
        errw1 = (WMTO - WMTO1) / WMTO
        errw2 = (WMTO - WMTO2) / WMTO
        errw3 = (WMTO - WMTO3) / WMTO

        errw = max(abs(errw1), abs(errw2), abs(errw3))

        # Print weight/ convergnce started
        if (printiter && iterw == 1)
            @printf("%5s %15s %15s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s \n",
                "iterw", "errW", "errW1", "WMTO", "Wpay", "Wfuel", "Weng", "Wfuse", "Wwing", "span", "area", "HTarea", "xwbox")
        end
        if printiter
            @printf("%5d %+13.8e %+13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",
                iterw, errw, errw1, parm[imWTO], parg[igWpaymax], parg[igWfuel], parg[igWeng],
                fuse.weight, wing.weight, wing.layout.span, wing.layout.S,
                htail.layout.S, wing.layout.box_x)
        end
        if (errw <= tolerW)
            Lconv = true
            break
        end

        #--------------------------------
        ##  Wing sizing section
        #--------------------------------
        WMTO = parg[igWMTO]

        # Size wing area and chords at start-of-cruise
        ip = ipcruise1
        We = WMTO * para[iafracW, ip]
        CL = para[iaCL, ip]
        ρ0 = pare[ierho0, ip]
        u0 = pare[ieu0, ip]
        qinf = 0.5 * ρ0 * u0^2
        BW = We + WbuoyCR # Weight including buoyancy

        # Size the wing area and chords
        set_wing_geometry!(BW, CL, qinf, wing)

        # Update wing box chord for fuseW in next iteration
        cbox = wing.layout.root_chord * wing.inboard.cross_section.width_to_chord

        # Calculate wing centroid and mean aerodynamic chord
        calculate_centroid_offset!(wing,calc_cma=true)
        xwing = wing.layout.x
        
        # Update wing pitching moment constants, uses the first mission point 
        # in the list as the reference for the others in the range.
        update_wing_pitching_moments!(para, ipstatic:ipclimb1, wing)
        update_wing_pitching_moments!(para, ipclimb1+1:ipdescentn-1, wing)
        update_wing_pitching_moments!(para, ipdescentn:ipdescentn, wing)

        # Calculate wing center load
        ip = ipcruise1
        γt, γs = wing.outboard.λ * para[iarclt, ip], wing.inboard.λ * para[iarcls, ip]
        Lhtail = WMTO * htail.CL_CLmax * htail.layout.S / wing.layout.S

        po = wing_loading(wing, para[iarclt, ip], para[iarcls, ip], Nlift, BW, Lhtail)

        # Calculate wing engine weight
        if compare_strings(options.opt_engine_location,"wing")
            if compare_strings(options.opt_prop_sys_arch,"te")
                @error "Turboelectric architectures are not currently supported. Their reintroduction with `struct`s is on the roadmap."
            elseif compare_strings(options.opt_prop_sys_arch,"tf")  || compare_strings(options.opt_prop_sys_arch,"constant_tsfc")
                Weng1 = parg[igWeng] / parg[igneng]
            end
        else
            Weng1 = 0.0
        end

        # Set up parameters for wing_weights function
        Winn, Wout = wing.inboard.weight, wing.outboard.weight
        dyWinn, dyWout = wing.inboard.dyW, wing.outboard.dyW
        rhofuel = !(options.has_wing_fuel) ? 0.0 : parg[igrhofuel]

        # Call wing_weights function
        Wwing,Wsinn,Wsout,
        dyWsinn,dyWsout,
        Wfcen,Wfinn,Wfout,
        dxWfinn,dxWfout,
        dyWfinn,dyWfout,lstrutp = wing_weights!(wing, po, γt, γs,
                                            Nlift, Weng1, 0, 0.0, 1, wing.layout.ηs,
                                            parg[igsigfac], rhofuel)

        # Calculate fuel weight if stored in wings
        Wfmax, dxWfmax, rfmax = 0.0, 0.0, 0.0
        if (options.has_wing_fuel)
            Wfmax = 2.0 * ((options.has_centerbox_fuel ? Wfcen : 0.0) + Wfinn + Wfout)
            dxWfmax = 2.0 * (dxWfinn + dxWfout)
            Wfuelmp = Wpay - Wpaymax + parg[igWfuel]
            rfmax = Wfuelmp / Wfmax
        end

        # Update wing properties
        wing.weight = Wwing * rlx + wing.weight * (1.0 - rlx)
        parg[igWfmax] = Wfmax
        # wing.dxW = dxWwing
        parg[igdxWfuel] = dxWfmax * rfmax

        wing.outboard.webs.weight = wing.inboard.webs.weight
        wing.outboard.caps.weight = wing.inboard.caps.weight

        # Calculate strut properties
        if wing.has_strut
            cstrut = sqrt(0.5 * wing.strut.axial_force / (tohstrut * wing.strut.thickness_to_chord))
            wing.strut.chord = cstrut
            wing.strut.S = 2.0 * cstrut * lstrutp
        end

        # Update wing panel weights
        wing.inboard.weight = Wsinn * (1.0 + fwadd) + rfmax * Wfinn
        wing.outboard.weight = Wsout * (1.0 + fwadd) + rfmax * Wfout
        wing.inboard.dyW = dyWsinn * (1.0 + fwadd) + rfmax * dyWfinn
        wing.outboard.dyW = dyWsout * (1.0 + fwadd) + rfmax * dyWfout

        #TODO: No reason why above lines shouldnt be inside wing_weights
        # -------------------------------
        #      Tail sizing section
        # -------------------------------

        # Set tail CL derivative
        dϵdα = htail.downwash_factor
        tanL = tand(wing.layout.sweep)
        tanLh = tand(htail.layout.sweep)

        ip = ipcruise1
        Mach = para[iaMach, ip]
        β = sqrt(1.0 - Mach^2) # Prandtl-Glauert factor

        # Calculate tail lift-curve slope
        dCLhdCL = (β + 2.0 / wing.layout.AR) / (β + 2.0 / htail.layout.AR) * 
                  sqrt(β^2 + tanL^2) / sqrt(β^2 + tanLh^2) * (1.0 - dϵdα)
        parg[igdCLhdCL] = dCLhdCL

        # Set Nacelle CL derivative fraction
        dCLnda = parg[igdCLnda]
        dCLndCL = dCLnda * (β + 2.0 / wing.layout.AR) * sqrt(β^2 + tanL^2) / 
                  (2.0 * π * (1.0 + 0.5 * wing.inboard.cross_section.thickness_to_chord))
        parg[igdCLndCL] = dCLndCL

        # Fuselage pitching moment
        #=       Use this with caution - slender body theory can be used here to estimate the fuselage 
        #       pitching moment - this ofc isn't true if the aircraft fuselage isn't "slender"
        #       Drela used a 3D panel method to actually calculate the CMVf1 and CMV0  for the aircraft studied in the N+3 work
        #       If sizes are roughly that of the 737/ 777 or D8 perhaps best to use those values and and specify 
        #       `calculates_fuse_moment_vol_deriv = false`, which would better capture the "thick body" effects. 
        =#
        if fuse.calculates_pitching_moment_volume
            cosL = cos(wing.layout.sweep * π / 180.0)
            Mperp = Mach * cosL
            βn = sqrt(1 - Mperp^2) # PG correction factor with M⟂ 

            # Estimate finite wing ∂CL/∂α from thin airfoil lift-slope 2π and 
            #  corrections for sweep and compressibility:
            CLα = 2π * cosL / (sqrt(βn^2 + (2 * cosL / wing.layout.AR)^2) + 2 * cosL / wing.layout.AR)
            # Estimate CMVf1 via slender body theory: dM/dα = 𝒱 ⟹ dM/dCL = dM/dα × dα/dCL = 𝒱/(dCL/dα)
            parg[igCMVf1] = fuse.volume/CLα
        #else, the offline specified value is kept
        end 

        # Size HT
        if (iterw <= 2 && initwgt == 0)
            lhtail = xhtail - xwing
            Vh = htail.volume
            Sh = Vh * wing.layout.S * wing.mean_aero_chord / lhtail
            htail.layout.S = Sh
        else
            size_htail(ac, view(para, :, ipdescentn), view(para, :, ipcruise1), view(para, :, ipcruise1);
                    Ldebug=Ldebug)
            wing.layout.box_x, xwing = wing.layout.box_x, wing.layout.x
            lhtail = xhtail - xwing
            Sh = htail.layout.S
            htail.volume = Sh * lhtail / (wing.layout.S * wing.mean_aero_chord)
        end

        # Vertical tail sizing 
        ip = iprotate
        qstall = 0.5 * pare[ierho0, ip] * (pare[ieu0, ip] / 1.2)^2
        dfan = parg[igdfan]
        CDAe = parg[igcdefan] * 0.25π * dfan^2
        De = qstall * CDAe
        Fe = pare[ieFe, ip]

        # Calculate max eng out moment
        Me = (Fe + De) * yeng

        #Size vertical tail ("size_vtail()")
        if compare_strings(vtail.opt_sizing,"fixed_Vv")
            lvtail = xvtail - xwing
            Vv = vtail.volume
            Sv = Vv * wing.layout.S * wing.layout.span/ lvtail
            vtail.layout.S = Sv
            parg[igCLveout] = Me / (qstall * Sv * lvtail)
        elseif compare_strings(vtail.opt_sizing,"OEI")
            lvtail = xvtail - xwing
            CLveout = parg[igCLveout]
            Sv = Me / (qstall * CLveout * lvtail)
            vtail.layout.S = Sv
            vtail.volume = Sv * lvtail / (wing.layout.S * wing.layout.span)
        end

        # Set HT max loading magnitude
        poh,htail.layout.span = tail_loading!(htail,Sh, qne)
        htail.layout.ηs = htail.layout.ηo
        
        # Set VT max loading magnitude, based on single tail + its bottom image
        pov,bv2 = tail_loading!(vtail,2.0 * Sv / vtail.ntails, qne; t_fac=2.0)
        bv = bv2 / 2
        vtail.layout.span = bv2
        vtail.layout.ηs = vtail.layout.ηo

        # HT weight
        htail.weight, _ = wing_weights!(htail, poh, htail.outboard.λ, htail.inboard.λ,
            0.0, 0.0, 0, 0.0, 0, 0.0,
            parg[igsigfac], rhofuel)
        
        # HT centroid x-offset
        calculate_centroid_offset!(htail, htail.layout.span, λhs)
        # HT pitching moment coeff
        fLoh, fLth = 0.0, fLt
        CMh0, CMh1 = wing_CM(htail.layout.span, htail.layout.root_span, htail.layout.root_span, htail.layout.sweep, wing.layout.spar_box_x_c, htail.outboard.λ, 1.0, htail.outboard.λ, 1.0,
            htail.layout.AR, fLoh, fLth, 0.0, 0.0, 0.0)
        para[iaCMh0, :] .= CMh0
        para[iaCMh1, :] .= CMh1

        # VT weight
        vtail.weight, _ = wing_weights!(vtail, pov, vtail.outboard.λ, vtail.inboard.λ,
            0.0, 0.0, 0, 0.0, 0, 0.0,
            parg[igsigfac], rhofuel; n_wings=vtail.ntails)
        # Set VT span
        vtail.layout.span = vtail.layout.span/2.0
        
        # VT centroid x-offset
        calculate_centroid_offset!(vtail, bv2, λhs)
        
        # ----------------------
        #     Fuselage Fuel Tank weight
        # ----------------------
        if !(options.has_wing_fuel) #If fuel is stored in the fuselage
            
            #Size fuel tank and calculate weight
            tanksize!(ac)

            # Update fuselage according to tank requirements
            update_fuse!(ac) #update fuselage length to accommodate tank; boundary layer also recalculated
            
            #Use homogeneous tank model to calculate required venting
            _, ps, _, _, _, _, _, Mvents, _, _ = CryoTank.analyze_TASOPT_tank(ac, fuse_tank.t_hold_orig, fuse_tank.t_hold_dest)
            parm[imWfvent] = Mvents[end] * gee #Store total fuel weight that is vented
            parg[igWfvent] = parm[imWfvent] #Store vented weight as parg parameter too
            fuse_tank.pmin = minimum(ps) #Store minimum tank pressure across mission

            if iterw > 2 #Calculate takeoff engine state and time
                #This is needed to know the TO duration to integrate the tank state
                # set static thrust for takeoff routine
                ip = ipstatic
                case = "off_design"
                eng.enginecalc!(ac, case, imission, ip, initializes_engine)

                # set rotation thrust for takeoff routine
                # (already available from cooling calculations)
                ip = iprotate
                case = "off_design"
                eng.enginecalc!(ac, case, imission, ip, initializes_engine)

                takeoff!(ac; printTO = false)
            end
        end

        # -----------------------------
        # Heat exchanger design and operation
        # ------------------------------
        ipdes = ipcruise1 #Design point: start of cruise

        if iterw > 2 #Only include heat exchangers after second iteration
            if eng.model.model_name == "fuel_cell_with_ducted_fan"
                ipdes = iprotate #Design point: takeoff rotation
                pare[ieRadiatorCoolantT,:] = eng.data.FC_temperature[:,imission]
                pare[ieRadiatorCoolantP,:] = eng.data.FC_pressure[:,imission]
                pare[ieRadiatorHeat,:] = eng.data.FC_heat[:,imission]
            end
            eng.heat_exchangers = hxdesign!(ac, ipdes, imission, rlx = 0.5) #design and off-design HX performance

            for HX in eng.heat_exchangers
                if HX.type == "Radiator"
                    TASOPT.engine.VerifyRadiatorHeat(engine, imission)
                end
            end
            #Note that engine state at takeoff should be calculated every iteration for correct balance-field. 
            #With fuel storage in tanks, this is done in the block above.           
        end

        # -----------------------------
        # Landing gear sizing
        # ------------------------------
        size_landing_gear!(ac)

        # -----------------------------
        # Drag and engine calculations
        # ------------------------------
        # Total Drag

        WMTO = parg[igWMTO]
        #calculate for start-of-cruise point
        # ip = ipclimbn
        ip = ipcruise1

        #Calculate fuel weight moment for balance
        if (options.has_wing_fuel) #If fuel is stored in the wings
            xfuel = wing.layout.box_x + parg[igdxWfuel] / parg[igWfuel]
            parg[igxWfuel] = parg[igWfuel] * wing.layout.box_x + parg[igdxWfuel] #Store fuel weight moment
        end

        # Pitch trim by adjusting Clh or by moving wing
        Wzero = WMTO - parg[igWfuel] #Zero fuel weight
        Wf = para[iafracW, ip] * WMTO - Wzero
        rfuel = Wf / parg[igWfuel]
        rpay = 1.0
        ξpay = 0.0
        opt_trim_var = "CL_htail"
        balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

        # Set N.P. at cruise
        parg[igxNP] = para[iaxNP, ip]

        para[iaalt, ipclimbn] = para[iaalt, ipcruise1]

        # Drag buildup aircraft_drag!()
        computes_wing_direct = true
        aircraft_drag!(ac, imission, ip, computes_wing_direct)

        # L/D and Design point thrust
        # println("CD = ", para[iaCD,ip])
        LoD = para[iaCL, ip] / para[iaCD, ip]
        gamV = para[iagamV, ip]
        We = para[iafracW, ip] * WMTO
        BW = We + WbuoyCR
        Fdes = BW * (1 / LoD + gamV)

        pare[ieFe, ip] = Fdes / neng

        # Size engine for TOC
        case = "design" #Design the engine for this mission point
        eng.enginecalc!(ac, case, imission, ip, initializes_engine, iterw)

        #Calculate engine mass properties
        eng.engineweight!(ac)

        _mission_iteration!(ac, imission, Ldebug)

        # this calculated fuel is the design-mission fuel 
        parg[igWfuel] = parm[imWfuel]
        
        # Store all OPRs for diagnostics
        pare[ieOPR, :] .= pare[iepilc, :] .* pare[iepihc, :]
        # size cooling mass flow at takeoff rotation condition (at Vstall)
        ip = iprotate

        # must define CDwing for this point in case there's wing BLI
        cdfw = para[iacdfw, ip] * para[iafexcdw, ip]
        cdpw = para[iacdpw, ip] * para[iafexcdw, ip]
        cosL = cosd(wing.layout.sweep)
        para[iaCDwing, ip] = cdfw + cdpw * cosL^3

        case = "cooling_sizing"
        eng.enginecalc!(ac, case, imission, ip, initializes_engine, iterw)

        # Recalculate weight update_weights!()
        ip = ipcruise1
        update_weights!(ac, rlx)

        parm[imWTO] = parg[igWMTO]
        parm[imWfuel] = parg[igWfuel]

        # Set previous iteration weights 
        WMTO3 = WMTO2
        WMTO2 = WMTO1
        WMTO1 = parg[igWMTO]

        ifirst = false

        # Get mission fuel burn (check if fuel capacity is sufficent)

        # Recalculate weight update_weights!()
        ip = ipcruise1
        update_weights!(ac, rlx)

        parm[imWTO] = parg[igWMTO]
        parm[imWfuel] = parg[igWfuel]
        # printstyled("Wfuel = $(parg[igWfuel]) \n", color=:blue)

        # Set previous iteration weights 
        WMTO3 = WMTO2
        WMTO2 = WMTO1
        WMTO1 = parg[igWMTO]

        # BFL calculations/ Noise? / Engine perf 

    end

    # normal takeoff and balanced-field takeoff calculations
    # set static thrust for takeoff routine
    ip = ipstatic
    case = "off_design"
    eng.enginecalc!(ac, case, imission, ip, initializes_engine)

    # set rotation thrust for takeoff routine
    # (already available from cooling calculations)
    ip = iprotate
    case = "off_design"
    eng.enginecalc!(ac, case, imission, ip, initializes_engine)

    if options.calculate_takeoff == true #If the engine can model the takeoff performance
        # calculate takeoff and balanced-field lengths
        takeoff!(ac, printTO = printiter)
    else
        @warn "Engine model does not allow takeoff calculations"
    end

    # calculate CG limits from worst-case payload fractions and packings
    rfuel0, rfuel1, rpay0, rpay1, xCG0, xCG1 = CG_limits(ac; Ldebug = Ldebug)
    parg[igxCGfwd] = xCG0
    parg[igxCGaft] = xCG1
    parg[igrpayfwd] = rpay0
    parg[igrpayaft] = rpay1

    # set neutral point at cruise
    ip = ipcruise1
    Wzero = WMTO - parg[igWfuel]
    Wf = para[iafracW, ip] * WMTO - Wzero
    rfuel = Wf / parg[igWfuel]
    rpay = 1.0
    ξpay = 0.0
    opt_trim_var = "none"
    balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, opt_trim_var; 
                        Ldebug = Ldebug)

    #Check if all engine points have converged, warn if not
    check_engine_convergence_failure(pare)
    
    #Warn user if HX effectiveness is overwritten
    check_HX_overwriting(eng.heat_exchangers) 
end

#TODO: update_WMTO! and update_weights! docstrings need full description
"""
    update_WMTO!(ac, rlx)

update_WMTO! updates the max takeoff weight of the aircraft (WMTO). Uses relaxation factor rlx.
Formerly, `Wupdate0!()`.
"""
function update_WMTO!(ac, rlx)
    parg, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear = unpack_ac_components(ac)

    WMTO = parg[igWMTO]
    
    ftotadd = fuse.HPE_sys.W #TODO this should be stored as a weight fraction, not a weight

    Wsum = parg[igWpay] +
           fuse.weight +
           wing.weight +
           wing.strut.weight +
           htail.weight +
           vtail.weight +
           parg[igWeng] +
           parg[igWfuel] +
           parg[igWtesys] +
           parg[igWftank] +
           landing_gear.nose_gear.weight.W + 
           landing_gear.main_gear.weight.W

    #update WMTO
    WMTO = rlx * Wsum / (1.0 - ftotadd) + (1.0 - rlx) * WMTO
    parg[igWMTO] = WMTO

    #check that fsum <= 1.0
    fwing = wing.weight / WMTO
    fstrut = wing.strut.weight / WMTO
    fhtail = htail.weight / WMTO
    fvtail = vtail.weight / WMTO
    feng = parg[igWeng] / WMTO
    ffuel = parg[igWfuel] / WMTO
    flgnose = landing_gear.nose_gear.weight.W / WMTO
    flgmain = landing_gear.main_gear.weight.W / WMTO
    ftesys = parg[igWtesys] / WMTO
    ftank = parg[igWftank] / WMTO

    fsum = fwing + fstrut + fhtail + fvtail + feng + ffuel + fuse.HPE_sys.W +
           flgnose + flgmain + ftank + ftesys

    if (fsum >= 1.0)
        @error "Something is wrong!! fsum ≥ 1.0"
    end

end


"""
    update_weights!(ac, rlx)

Adjusts the aircraft's maximum takeoff weight (WMTO) and other component weights using a relaxation factor (`rlx`). 
Formerly, `Wupdate!()`.

"""
function update_weights!(ac, rlx)
    parg, options, fuse, fuse_tank, wing, htail, vtail, engine, landing_gear = unpack_ac_components(ac)

    WMTO = parg[igWMTO]

    fwing = wing.weight / WMTO
    fstrut = wing.strut.weight / WMTO
    fhtail = htail.weight / WMTO
    fvtail = vtail.weight / WMTO
    feng = parg[igWeng] / WMTO
    ffuel = parg[igWfuel] / WMTO
    flgnose = landing_gear.nose_gear.weight.W / WMTO
    flgmain = landing_gear.main_gear.weight.W / WMTO

    ftesys = parg[igWtesys] / WMTO

    Wtesys = parg[igWtesys]
    Wftank = parg[igWftank]
    Wpay = parg[igWpay]
    Wfuse = fuse.weight

    ftank = parg[igWftank] / WMTO

    fsum = fwing + fstrut + fhtail + fvtail + feng + ffuel + fuse.HPE_sys.W +
           flgnose + flgmain + ftank + ftesys

    if (fsum ≥ 1.0)
        @error "SOMETHING IS WRONG fsum ≥ 1"
    end

    # WMTO = rlx*(Wpay + Wfuse + Wtesys + Wftank)/(1.0-fsum) + (1.0-rlx)*WMTO
    WMTO = rlx * (Wpay + fuse.weight) / (1.0 - fsum) + (1.0 - rlx) * WMTO

    parg[igWMTO] = WMTO
    wing.weight = WMTO * fwing
    wing.strut.weight = WMTO * fstrut
    htail.weight = WMTO * fhtail
    vtail.weight = WMTO * fvtail
    parg[igWeng] = WMTO * feng
    parg[igWfuel] = WMTO * ffuel
   
    parg[igWftank] = WMTO * ftank 

    parg[igWtesys] = WMTO * ftesys
end

"""
    set_ambient_conditions!(ac, ip, Mach=NaN; im = 1)

Sets ambient condition at the given mission point `ip` and mission `im` (default is 1).

Creates a `FlightCondition` from altitude and Mach, stores it in `ac.flight_conditions`,
and updates legacy `para`/`pare` arrays for backward compatibility.
"""
function set_ambient_conditions!(ac, ip, Mach=NaN; im = 1)
    ΔTatmos = ac.parm[imDeltaTatm, im]
    alt = ac.para[iaalt, ip, im]
    if Mach === NaN
        Mach = ac.para[iaMach, ip, im]
    end
    γ = ac.para[iagamV, ip, im]

    # Create FlightCondition and store it
    fc = FlightCondition(alt, Mach; ΔT=ΔTatmos, climb_angle=γ)
    set_flight_condition!(ac, ip, fc; im=im)
end  # function set_ambient_conditions

"""
update_wing_pitching_moments!(para, ip_range, wing)

Updates wing pitching moments and calls wing_CM for mission points
"""
function update_wing_pitching_moments!(para, ip_range, wing)
    ip = ip_range[1]
    cmpo, cmps, cmpt = para[iacmpo, ip], para[iacmps, ip], para[iacmpt, ip]
    γt = wing.outboard.λ * para[iarclt, ip]
    γs = wing.inboard.λ * para[iarcls, ip]
    
    CMw0, CMw1 = wing_CM(
        wing.layout.span, wing.layout.break_span, wing.layout.root_span, 
        wing.layout.sweep, wing.layout.spar_box_x_c, wing.outboard.λ, wing.inboard.λ, 
        γt, γs, wing.layout.AR, wing.fuse_lift_carryover, wing.tip_lift_loss, cmpo, cmps, cmpt
    )
    for ip in ip_range
        para[iaCMw0, ip] = CMw0
        para[iaCMw1, ip] = CMw1
    end
end

"""
    setup_fuel_storage!(options, fuse, fuse_tank, parg, pare)

Initializes fuel storage parameters based on whether fuel is stored in the
wings or fuselage. For fuselage storage, computes fuel properties from the
saturated mixture model and sets tank positions. #TODO this implicitly assumes that 
if the fuel is stored in the tanks then it is cryogenic fuel... which is not 
necessarily true.

Also resets heat exchanger values for the engine.
"""
function setup_fuel_storage!(options, fuse, fuse_tank, parg, pare)
    if options.has_wing_fuel
        xftank = xftankaft = 0.0
    else
        xftank = fuse.layout.x_start_cylinder + 1.0 * ft_to_m
        xftankaft = fuse.layout.x_end_cylinder

        # Calculate fuel properties from saturated mixture model
        β0 = 1 - fuse_tank.ullage_frac
        fuel_mix = SaturatedMixture(fuse_tank.fueltype, fuse_tank.pvent, β0)
        Tfuel = fuel_mix.liquid.T
        ρliq = fuel_mix.liquid.ρ
        ρgas = fuel_mix.gas.ρ
        hvap = fuel_mix.hvap

        # Set fuel properties in parameter arrays and tank struct
        pare[ieTft, :] .= Tfuel
        pare[ieTfuel, :] .= Tfuel
        parg[igrhofuel] = fuel_mix.ρ
        fuse_tank.rhofuel = ρliq
        fuse_tank.Tfuel = Tfuel
        fuse_tank.hvap = hvap
        fuse_tank.rhofuelgas = ρgas
    end

    # Update fuel tank positions
    parg[igxftank] = xftank
    parg[igxftankaft] = xftankaft

    # Reset engine values for heat exchangers
    resetHXs(pare)

    return nothing
end