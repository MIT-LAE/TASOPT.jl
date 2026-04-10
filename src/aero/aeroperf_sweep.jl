using Printf
export calc_wing_aoa, calc_ac_AoA, calc_ac_AoA!, calc_ac_Theta!, calc_mission_attitude!

"""
    aeroperf_sweep(ac_orig, CL_range; Mach_ac=nothing, imission=1, ip=ipcruise1, rfuel=1, rpay=1, ξpay=0.5)

Performs a sweep over a range of lift coefficients (`CL_range`) for a given aircraft model, evaluating aerodynamic and performance metrics at each point.

This function deep-copies the input aircraft model, sets the lift coefficient, balances the aircraft, and computes drag and related quantities for each value in `CL_range`. 
    The results are collected in arrays for further analysis or plotting and returned in a NamedTuple.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac_orig`: Aircraft model object (deep-copied internally).
    - `CL_range`: Range of lift coefficients to sweep.
    - `Mach`: Mach number at which points are evaluated. If none specified, defaults to value at specified ip, imission (optl.).
    - `imission::Integer`: Mission index (default: 1, optl.).
    - `ip::Integer`: Flight point index, determining altitude (default: `ipcruise1`, optl.).
    - `rfuel::Float64`: Fuel fraction (default: 1, optl.).
    - `rpay::Float64`: Payload fraction (default: 1, optl.).
    - `ξpay::Float64`: Payload distribution factor (0.0 = front-loaded, 1.0 = rear-loaded; default: 0.5, optl.). 
                    Doesn't matter if rpay = 1.

    **Outputs:**
    - Returns a named tuple of vectors containing:
        - `Mach`: Mach number at which points are evaluated.
        - `CLs`: Lift coefficients at each point.
        - `CDs`: Drag coefficients at each point.
        - `LDs`: Lift-to-drag ratios.
        - `CLhs`: Horizontal tail lift coefficients.
        - `CDis`: Induced drag coefficients.
        - `CDwings`, `CDfuses`, `CDhtails`, `CDvtails`, `CDothers`: Component drag breakdowns.
        - `clpos`, `clpss`, `clpts`: airfoil section lift coefficients (root, spanbreak, and tip).

    Sample usage:

        ```julia
        CL_range = 0.2:0.05:0.8
        results_nb = aeroperf_sweep(ac, CL_range, print_results=true)
        ```

See also: [`TASOPT.plot_drag_polar`](@ref), [`TASOPT.balance_aircraft!`](@ref), [`TASOPT.aerodynamics.aircraft_drag!`](@ref).
"""
function aeroperf_sweep(ac_orig, CL_range; Mach_ac=nothing, imission=1, ip=ipcruise1, rfuel=1, rpay=1, ξpay=0.5,
                        print_results = false)

    #confirm aircraft is sized
    if !ac_orig.is_sized[1]
        error("Aircraft must be sized via `size_aircraft!()` before performing an aeroperformance sweep.")
    end

    #deepcopy to avoid overwrites to original aircraft
    ac = deepcopy(ac_orig)

    #substitute Mach if specified
    if !isnothing(Mach_ac) && Mach_ac isa Number
        ac.para[iaMach, ip, imission] = Mach_ac
    end

    #initalize results tuple
    n = length(CL_range)
    results = (Mach = ac.para[iaMach, ip, imission],
               CLs = Vector{Float64}(undef, n),
               CDs = Vector{Float64}(undef, n),
               LDs = Vector{Float64}(undef, n),
               CLhs = Vector{Float64}(undef, n),
               CDis = Vector{Float64}(undef, n),
               CDwings = Vector{Float64}(undef, n),
               CDfuses = Vector{Float64}(undef, n),
               CDhtails = Vector{Float64}(undef, n),
               CDvtails = Vector{Float64}(undef, n),
               CDothers = Vector{Float64}(undef, n),
               
               #spanbreak quantities
               clpos = Vector{Float64}(undef, n),
               clpss = Vector{Float64}(undef, n),
               clpts = Vector{Float64}(undef, n),

               cdfss = Vector{Float64}(undef, n),
               cdpss = Vector{Float64}(undef, n),
               cdwss = Vector{Float64}(undef, n),
               cdss =  Vector{Float64}(undef, n),

               #aoa quantities
               aoaps = Vector{Float64}(undef, n),
               aoaws = Vector{Float64}(undef, n),
               AoA = Vector{Float64}(undef, n),

               #pitching moment quantities
               cmss    = Vector{Float64}(undef, n), #section cm at spanbreak (from airfun)
               CMwings = Vector{Float64}(undef, n), #wing CM about CG
               CMtails = Vector{Float64}(undef, n), #tail CM about CG
               CMfuses = Vector{Float64}(undef, n), #fuselage CM about CG
               CMs     = Vector{Float64}(undef, n), #total aircraft CM (≈ 0 in trimmed flight)
                )

    # geometry constants for component CM calculations (invariant over CL sweep with "CL_htail" trim)
    co    = ac.wing.layout.root_chord
    cma   = ac.wing.mean_aero_chord
    xwbox = ac.wing.layout.box_x
    S_ref = ac.wing.layout.S
    coh   = ac.htail.layout.root_chord
    xhbox = ac.htail.layout.box_x
    Sh    = ac.htail.layout.S
    CMVf1 = ac.parg[igCMVf1]
    CLMf0 = ac.parg[igCLMf0]

    for (i, CL) in enumerate(CL_range)
        ac.para[iaCL, ip, imission] = CL
        balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, TrimVar.CLHtail)
        aircraft_drag!(ac, imission, ip, true)

        #get basic performance
        results.CLs[i] = ac.para[iaCL, ip, imission]
        results.CDs[i] = ac.para[iaCD, ip, imission]
        results.CLhs[i] = ac.para[iaCLh, ip, imission]
        results.CDis[i] = ac.para[iaCDi, ip, imission]
        results.CDwings[i] = ac.para[iaCDwing, ip, imission]
        results.CDfuses[i] = ac.para[iaCDfuse, ip, imission] + ac.para[iaCDover, ip, imission] #"fuselage carryover CD", SET TO 0
        results.CDhtails[i] = ac.para[iaCDhtail, ip, imission]
        results.CDvtails[i] = ac.para[iaCDvtail, ip, imission]
        results.CDothers[i] = ac.para[iaCDstrut, ip, imission] + #strut drag
                                ac.para[iaCDnace, ip, imission] + #nacelle drag
                                ac.para[iadCDBLIf, ip, imission] + #fuse bli "drag" effect
                                ac.para[iadCDBLIw, ip, imission] #wing bli "drag" effect

        #get max clp for airfoil bounds
        results.clpos[i] = ac.para[iaclpo, ip, imission]
        results.clpss[i] = ac.para[iaclps, ip, imission]
        results.clpts[i] = ac.para[iaclpt, ip, imission]

        #get airfoil section drag via airfun (only for spanbreak section)
        toc = ac.wing.outboard.cross_section.thickness_to_chord
        sweep = ac.wing.sweep
        Mach_ac = ac.para[iaMach,ip,imission]
        Mach_perp = Mach_ac*cosd(sweep)  # Perpendicular Mach number
        airf = ac.wing.airsection

        cdfss, cdpss, cdwss, cms, aoaps = airfun(results.clpss[i], toc, Mach_perp, airf)
        results.cdfss[i] = cdfss #skin friction
        results.cdpss[i] = cdpss #pressure
        results.cdwss[i] = cdwss #wave (= 0, with models at time of writing)
        results.cdss[i] = cdfss + cdpss + cdwss

        #compute aoa for section and AoA for aircraft (all in radians)
        wing_mounting_angle = ac.wing.mounting_angle
        results.aoaps[i] = aoaps   #perpendicular AoA for spanbreak section [rad]
        results.aoaws[i] = atan(cosd(sweep)*tan(aoaps)) #body-aligned AoA for spanbreak section [rad]
        results.AoA[i] = results.aoaws[i] - wing_mounting_angle #aircraft AoA [rad]

        #store section cm from airfun
        results.cmss[i] = cms

        #compute aircraft component pitching moments about CG (formulas from balance.jl)
        CL_i  = results.CLs[i]
        CLh_i = results.CLhs[i]
        xCG   = ac.para[iaxCG, ip, imission]
        CMw0  = ac.para[iaCMw0, ip, imission]
        CMw1  = ac.para[iaCMw1, ip, imission]
        CMh0  = ac.para[iaCMh0, ip, imission]
        CMh1  = ac.para[iaCMh1, ip, imission]

        #TODO: confirm these calculations are correct (see Drela TASOPT PDF 2.11.2)
        #they were borrowed from commented-out sections of balance.jl
        #REMEMBER TO CHECK IF THE AIRCRAFT IS TRIMMED, AND IF THE REFERENCE POINTS ARE CORRECT/CONSISTENT
        results.CMwings[i] = (co * CMw0 + (co * CMw1 - xwbox - xCG) * (CL_i - CLh_i * Sh / S_ref)) / cma
        results.CMtails[i] = (coh * CMh0 * Sh / S_ref + (coh * CMh1 - xhbox - xCG) * CLh_i * Sh / S_ref) / cma
        results.CMfuses[i] = CMVf1 * (CL_i - CLMf0) / (S_ref * cma)
        results.CMs[i]     = results.CMwings[i] + results.CMtails[i] + results.CMfuses[i]
    end
    
    # Calculate L/D for each point
    for i in 1:n
        results.LDs[i] = results.CLs[i] / results.CDs[i]
    end

    if print_results
        # get headers
        fields = propertynames(results)
        headers = [string(f)[1:end-1] for f in fields] #Make headers by stripping the trailing "s" (e.g. :CLs → "CL")
        # Collect the columns
        table_data = hcat([getfield(results, f) for f in fields]...)
        
        #table header row
        println(join(headers, "\t"))

        # Print each row with each value formatted
        for (_, row) in enumerate(eachrow(table_data))
            formatted_row = [@sprintf("%5.4f", v) for v in row]
            println(join(formatted_row, "\t"))
        end
    end
        
    return results
end 

"""
    calc_wing_aoa(ac; ip=ipcruise1, imission=1)

Return the body-frame angle of attack (radians) of the wing spanbreak section at mission
point `ip` of mission `imission`.

Uses `airfun` to look up the perpendicular AoA (`aoaps`, radians) for the spanbreak airfoil
at the local section lift coefficient (`iaclps`) and perpendicular Mach number, then projects
back to the body-aligned frame: `aoaws = atan(cosd(sweep) * tan(aoaps))`.
Wing `sweep` is in degrees; `cosd` is used to stay consistent with the `WingLayout` convention.
"""
function calc_wing_aoa(ac; ip=ipcruise1, imission=1)
    # look up perpendicular AoA for spanbreak section (aoaps in radians from airfun)
    toc = ac.wing.outboard.cross_section.thickness_to_chord
    sweep = ac.wing.sweep
    Mach = ac.para[iaMach, ip, imission]
    Mach_perp = Mach * cosd(sweep)
    airf = ac.wing.airsection
    clps = ac.para[iaclps, ip, imission]
    _, _, _, _, aoaps = airfun(clps, toc, Mach_perp, airf)

    # project perpendicular AoA to body-aligned frame [rad]
    aoaws = atan(cosd(sweep) * tan(aoaps))
    return aoaws
end

"""
    calc_ac_AoA(ac; ip=ipcruise1, imission=1)

Return the aircraft angle of attack (radians) at mission point `ip` of mission `imission`.

Aircraft AoA is the wing spanbreak body-frame AoA minus the wing mounting angle:
`AoA = aoaws - wing.mounting_angle`, both in radians. Non-mutating; use [`calc_ac_AoA!`](@ref)
to store the result in `ac.para[iaAoA, ip, imission]`.
"""
function calc_ac_AoA(ac; ip=ipcruise1, imission=1)
    aoaws = calc_wing_aoa(ac, ip=ip, imission=imission)
    wing_mounting_angle = ac.wing.mounting_angle
    AoA = aoaws - wing_mounting_angle
    return AoA
end

"""
    calc_ac_AoA!(ac; ip=ipcruise1, imission=1)

Compute the aircraft angle of attack at mission point `ip` of mission `imission` and store
the result in `ac.para[iaAoA, ip, imission]`. Returns the AoA in radians.

See also [`calc_ac_AoA`](@ref) for the non-mutating version.
"""
function calc_ac_AoA!(ac; ip=ipcruise1, imission=1)
    AoA = calc_ac_AoA(ac, ip=ip, imission=imission)
    ac.para[iaAoA,ip,imission] = AoA
    return AoA
end

"""
    calc_ac_Theta!(ac; ip=ipcruise1, imission=1)

Compute the aircraft pitch attitude Θ at mission point `ip` of mission `imission` and store
the result in `ac.para[iaTheta, ip, imission]`. Returns Θ in radians.

Θ is the sum of the aircraft AoA and the flight-path angle γ (`iagamV`), both in radians:
`Theta = AoA + gamma`. Also stores AoA via [`calc_ac_AoA!`](@ref) as a side effect.
"""
function calc_ac_Theta!(ac; ip=ipcruise1, imission=1)
    AoA = calc_ac_AoA!(ac, ip=ip, imission=imission)
    gamma = ac.para[iagamV,ip,imission] #flight path angle
    Theta = AoA + gamma
    ac.para[iaTheta,ip,imission] = Theta
    return Theta
end

"""
    calc_mission_attitude!(ac; imission=1, ips=ipclimb1:ipdescentn)

Compute AoA and Theta for each mission point in `ips`.

Skips ground/takeoff phases (`ipstatic`–`ipcutback`) by default, as `iagamV`,
`iaclps`, and `iaMach` are not in steady trimmed-flight states for those points.
`ips` can be overridden to restrict or extend the range as needed.
"""
function calc_mission_attitude!(ac; imission=1, ips=ipclimb1:ipdescentn)
    for ip in ips
        calc_ac_Theta!(ac, ip=ip, imission=imission)
        #calc_ac_AoA!() is called by calc_ac_Theta!()
    end
end