using Printf

"""
    aeroperf_sweep(ac_orig, CL_range; Mach=nothing, imission=1, ip=ipcruise1, rfuel=1, rpay=1, ξpay=0.5)

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

See also: [`TASOPT.DragPolar`](@ref), [`TASOPT.balance_aircraft!`](@ref), [`TASOPT.aerodynamics.aircraft_drag!`](@ref).
"""
function aeroperf_sweep(ac_orig, CL_range; Mach=nothing, imission=1, ip=ipcruise1, rfuel=1, rpay=1, ξpay=0.5,
                        print_results = false)

    #confirm aircraft is sized
    if !ac_orig.is_sized[1]
        error("Aircraft must be sized via `size_aircraft!()` before performing an aeroperformance sweep.")
    end

    #deepcopy to avoid overwrites to original aircraft
    ac = deepcopy(ac_orig)

    #substitute Mach if specified
    if !isnothing(Mach) && Mach isa Number
        ac.para[iaMach, ip, imission] = Mach
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
               
               clpos = Vector{Float64}(undef, n),
               clpss = Vector{Float64}(undef, n),
               clpts = Vector{Float64}(undef, n),

               cdfss = Vector{Float64}(undef, n),
               cdpss = Vector{Float64}(undef, n),
               cdwss = Vector{Float64}(undef, n),
               cdss =  Vector{Float64}(undef, n),
                )

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
        tau = ac.wing.outboard.cross_section.thickness_to_chord
        sweep = ac.wing.sweep
        Mach = ac.para[iaMach,ip,imission]
        Mach_perp = Mach*cosd(sweep)  # Perpendicular Mach number
        airf = ac.wing.airsection

        cdfss, cdpss, cdwss, _ = airfun(results.clpss[i], tau, Mach_perp, airf)
        results.cdfss[i] = cdfss #skin friction
        results.cdpss[i] = cdpss #pressure
        results.cdwss[i] = cdwss #wave (= 0, with models at time of writing)
        results.cdss[i] = cdfss + cdpss + cdwss
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
