import ..TASOPT: balance_aircraft!
using Printf

"""
    aeroperf_sweep(ac_orig, CL_vec; imission=1, ip=ipcruise1, rfuel=1, rpay=1, ξpay=0.5)

Performs a sweep over a vector of lift coefficients (`CL_vec`) for a given aircraft model, evaluating aerodynamic and performance metrics at each point.

This function deep-copies the input aircraft model, sets the target lift coefficient, balances the aircraft, and computes drag and related quantities for each value in `CL_vec`. 
    The results are collected in arrays for further analysis or plotting and returned in a NamedTuple.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac_orig`: Aircraft model object (deep-copied internally).
    - `CL_vec::AbstractVector{Float64}`: Vector of target lift coefficients to sweep.
    - `imission::Integer`: Mission index (default: 1).
    - `ip::Integer`: Flight point index (default: `ipcruise1`).
    - `rfuel::Float64`: Fuel fraction (default: 1).
    - `rpay::Float64`: Payload fraction (default: 1).
    - `ξpay::Float64`: Payload distribution factor (0.0 = front-loaded, 1.0 = rear-loaded; default: 0.5). Doesn't matter if rpay = 1.

    **Outputs:**
    - Returns a named tuple of vectors containing:
        - `CLs`: Lift coefficients at each point.
        - `CDs`: Drag coefficients at each point.
        - `LDs`: Lift-to-drag ratios.
        - `CLhs`: Horizontal tail lift coefficients.
        - `CDis`: Induced drag coefficients.
        - `CDwings`, `CDfuses`, `CDhtails`, `CDvtails`, `CDothers`: Component drag breakdowns.
        - `clpos`, `clpss`, `clpts`: airfoil section lift coefficients (root, spanbreak, and tip).

    Sample usage:

        ```julia
        CL_vec = [0.2:0.05:0.8...]
        results_nb = aeroperf_sweep(ac, CL_vec)
        ```

See also: [`TASOPT.DragPolar`](@ref), [`TASOPT.balance_aircraft!`](@ref), [`TASOPT.aerodynamics.aircraft_drag!`](@ref).
"""
function aeroperf_sweep(ac_orig, CL_vec; imission=1, ip=ipcruise1, rfuel=1, rpay=1, ξpay=0.5,
                        print_results = False)
    #initalize results tuple
    n = length(CL_vec)
    results = (CLs = Vector{Float64}(undef, n),
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
                )

    ac = deepcopy(ac_orig)

    for (i, CL) in enumerate(CL_vec)
        ac.para[iaCL, ip, imission] = CL
        balance_aircraft!(ac, imission, ip, rfuel, rpay, ξpay, "CL_htail")
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

        #get max clp for airfoil bounds if specified
        results.clpos[i] = ac.para[iaclpo, ip, imission]
        results.clpss[i] = ac.para[iaclps, ip, imission]
        results.clpts[i] = ac.para[iaclpt, ip, imission]

    end
    
    # Calculate L/D for each point
    for i in 1:n
        results.LDs[i] = results.CLs[i] / results.CDs[i]
    end

        # add a line to the first plot sampling the airfoil database
    airf = ac.wing.airsection
    tau = ac.wing.inboard.cross_section.thickness_to_chord
    Mach = ac.para[iaMach,ip,imission]
    sweep = ac.wing.sweep
    Mach_perp = Mach*cosd(sweep)  # Perpendicular Mach number

    # # Sample the airfoil database and plot as a reference curve
    # for (i, CL) in enumerate(CL_vec)
    #     CL_perp = CL / cosd(sweep)^2  # Perpendicular lift coefficient
    #     cdf, cdp, cdw, cm = TASOPT.aerodynamics.airfun(CL_perp, tau, Mach_perp, airf)
    #     # Save airfoil performance to results
    #     results.CLperps[i] = CL_perp  # airfoil CL
    #     results.CDperps[i] = cdf + cdp + cdw  # Total drag from airfoil
    #     # TODO: add sweep corrections
    #     results.CDs_inf_swept_wing[i] = cdf + cdp + cdw  # Total drag from airfoil
    #     results.CLs_inf_swept_wing[i] = CL  # Keep the same CL for comparison
    # end

    if print_results
        # get headers
        fields = propertynames(results)
        headers = [string(f)[1:end-1] for f in fields] #Make headers by stripping the trailing "s" (e.g. :CLs → "CL")
        # Collect the columns
        table_data = hcat([getfield(results, f) for f in fields]...)
        
        #table header row
        println(join(headers, "\t"))

        # Print each row with each value formatted
        for i in 1:size(table_data, 1)
            row = table_data[i, :]
            formatted_row = [@sprintf("%5.4f", v) for v in row]
            println(join(formatted_row, "\t"))
        end
    end
        
    return results
end 
