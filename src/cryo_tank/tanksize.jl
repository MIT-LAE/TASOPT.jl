"""Minimum insulation layer thickness (m). Small but positive, since a zero-thickness layer gives a singular thermal circuit."""
const T_INSUL_MIN = 1e-4

"""
        check_insulation_thickness(t)

Throws a `DomainError` if any insulation layer thickness is negative or non-finite. Used inside the insulation solve so that 
infeasible trial points are turned into large residuals by the solver wrapper.
"""
function check_insulation_thickness(t::AbstractVector{<:Real})
        for ti in t
                (isfinite(ti) && ti >= 0) || throw(DomainError(ti, "insulation layer thickness must be non-negative and finite"))
        end
        return nothing
end

"""
        clamp_insulation_thickness!(t, iinsuldes = eachindex(t))

Raises any layer in `t` (restricted to indices `iinsuldes`) that is below `T_INSUL_MIN` to `T_INSUL_MIN`, with a warning.
Returns `true` if any layer was clamped.
"""
function clamp_insulation_thickness!(t::AbstractVector{<:Real}, iinsuldes = eachindex(t))
        clamped = false
        for ind in iinsuldes
                if t[ind] < T_INSUL_MIN
                        clamped |= t[ind] < T_INSUL_MIN * (1 - 1e-6)
                        t[ind] = T_INSUL_MIN
                end
        end
        clamped && @warn "tanksize!: insulation layer thickness fell below the minimum; clamped to $(T_INSUL_MIN) m. " *
                         "The boil-off target may not be met by this insulation stack." maxlog=3
        return clamped
end

"""
        boiloff_heat_rate(fuse_tank)

Heat rate (W) into one tank that corresponds to the design boil-off rate `fuse_tank.boiloff_rate` (%/h), after removing the 
extra heat leak factor `qfac`.
"""
function boiloff_heat_rate(fuse_tank::fuselage_tank)
        mdot_boiloff = fuse_tank.boiloff_rate * fuse_tank.Wfuelintank / (gee * 100) / 3600 #Boil-off mass flow rate (kg/s)
        Q_net = mdot_boiloff * fuse_tank.hvap #Heat rate from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
        return Q_net / fuse_tank.qfac
end

"""
        min_insulation_shift(t_cond, iinsuldes)

Smallest uniform shift `Δt` that can be added to the sizeable layers `iinsuldes` of `t_cond` while keeping every one of them at 
or above `T_INSUL_MIN`. The insulation solve cannot go below this without producing a non-physical stack.
"""
min_insulation_shift(t_cond::AbstractVector{<:Real}, iinsuldes) = maximum(T_INSUL_MIN - t_cond[ind] for ind in iinsuldes)

"""
        max_heat_rate(fuse, fuse_tank, Δt, z, TSL, Mair, xftank, ifuel)

Heat rate (W, before the `qfac` factor) into the tank when the sizeable layers are shifted by `Δt`, found by solving the thermal 
network with the heat rate as an unknown. Called with the smallest admissible shift (see [`min_insulation_shift`](@ref)), it gives the 
largest heat rate (i.e. boil-off) that the insulation stack can reach. Returns `NaN` if the thermal solve fails.
"""
function max_heat_rate(fuse::Fuselage, fuse_tank::fuselage_tank, Δt::Real, z::Float64, TSL::Float64, Mair::Float64, 
                       xftank::Float64, ifuel::Int64)
        t_all = copy(fuse_tank.t_insul)
        for ind in fuse_tank.iinsuldes
                t_all[ind] += Δt
        end
        inner_tank = size_inner_tank(fuse, fuse_tank, t_all)

        p = thermal_params{typeof(fuse.layout.cross_section), Float64}()
        p.l_cyl = inner_tank.l_cyl
        p.l_tank = inner_tank.l_tank
        p.r_tank = inner_tank.Rtank_outer
        p.Shead = inner_tank.Shead_insul
        p.t_cond = t_all
        p.material = fuse_tank.material_insul
        p.Tfuel = fuse_tank.Tfuel
        p.z = z
        p.TSL = TSL
        p.Mair = Mair
        p.xftank = xftank
        p.ifuel = ifuel
        p.fuse_cs = fuse.layout.cross_section

        function residual(x)
                try
                        return residuals_Q(x, p, "Q_unknown")
                catch e
                        (e isa MethodError && eltype(x) <: ForwardDiff.Dual) && rethrow()
                        return ones(length(x)) * 1e3
                end
        end

        _, _, Taw = freestream_heat_coeff(z, TSL, Mair, xftank)
        ΔT = Taw - fuse_tank.Tfuel
        thickness = sum(t_all)
        #Same cold guess as in tankWthermal, with a guess for the total resistance that is tried over a wide range. The thinnest
        #admissible stack has a very small resistance (order 1e-4 K/W for a non-vacuum stack, where a guess of 0.01 K/W or more
        #stagnates at a non-physical point), while a vacuum layer or a small tank has a large one.
        for Rguess in (1e-3, 1e-2, 1e-4, 0.1, 1.0, 10.0)
                g = zeros(length(t_all) + 2)
                g[1] = ΔT / Rguess
                g[2] = fuse_tank.Tfuel + 1.0
                for i in eachindex(t_all)
                        g[i + 2] = fuse_tank.Tfuel + ΔT * sum(t_all[1:i]) / thickness
                end
                g[end] -= 1.0
                sol = nlsolve_thermal(residual, g, xtol = 1e-5, ftol = 1e-5, iterations = 500)
                #Only a true root is accepted (small residual, and temperatures between the fuel and the hot air). Convergence on step
                #size alone can stagnate at non-physical points with a spuriously high Q (e.g. a 0 K wall).
                T_phys = all(T -> fuse_tank.Tfuel - 1e-3 <= T <= Taw + 1e-3, @view sol.zero[2:end])
                if sol.f_converged && T_phys && isfinite(sol.zero[1]) && sol.zero[1] > 0
                        return sol.zero[1]
                end
        end
        return NaN
end

"""
        tanksize!(ac, imission)

`tanksize` sizes a cryogenic fuel tank for a cryogenic-fuel aircraft

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `ac::aircraft`: TASOPT aircraft object.
        - `imission::Int64`: design mission index (default is 1).
        
        **Outputs:**
        No direct outputs. Aircraft object gets modified with tank parameters.

See [here](@ref fueltanks).
"""
function tanksize!(ac, imission::Int64 = 1)
        #Unpack storage arrays
        fuse = ac.fuselage
        fuse_tank = ac.fuse_tank
        para = view(ac.para, :, :, imission)
        parg = ac.parg

        #Unpack variables in ac
        tank_placement = fuse_tank.placement
        nftanks = fuse_tank.tank_count
       
        #Convective cooling
        if tank_placement == "rear"
                xftank_heat = parg[igxftankaft]
        else
                xftank_heat = parg[igxftank]
        end
        ifuel = ac.options.ifuel
        Mair = para[iaMach, ipcruise1]
        z = para[iaalt, ipcruise1]
            
        fuse_tank.Wfuelintank = parg[igWfuel] / nftanks #Each fuel tank carries 1/nftanks of the fuel

        #Unpack variables in fuse_tank
        rhofuel = fuse_tank.rhofuel
        t_cond = fuse_tank.t_insul #thickness of each insulation layer
        Wfuelintank = fuse_tank.Wfuelintank #weight of fuel in tank
        Tfuel = fuse_tank.Tfuel #fuel temperature
        sizes_insulation = fuse_tank.sizes_insulation #Boolean for whether to size for a boiloff rate
        TSL = fuse_tank.TSLtank[imission] #sea-level temperature for tank design

        #------Size insulation, if requested------
        if sizes_insulation #If insulation is sized for a given boiloff rate
                boiloff_percent = fuse_tank.boiloff_rate
                iinsuldes = fuse_tank.iinsuldes

                #Thermal calculations
                _, _, Taw = freestream_heat_coeff(z, TSL, Mair, xftank_heat) #Find adiabatic wall temperature 

                ΔT = Taw - Tfuel

                #Create inner function with residuals as a function of x
                function residual(x)
                        try
                                return res_MLI_thick(x, fuse, fuse_tank, z, TSL, Mair, xftank_heat, ifuel) #Residual in boiloff rate as a function of Δt
                        catch e
                                #A MethodError with ForwardDiff duals means the residual is not AD-compatible: rethrow so that 
                                #`nlsolve_autodiff` can fall back to finite differences instead of silently masking it.
                                (e isa MethodError && eltype(x) <: ForwardDiff.Dual) && rethrow()
                                return ones(length(x))*1e3
                        end
                end

                #Assemble guess for non linear solver
                #x[1] = Δt; x[2] = T_tank; x[3:(end-1)]: T at edge of insulation layer; x[end] = T at fuselage wall
                function cold_guess()
                        g = zeros(length(t_cond) + 2)
                        g[1] = 0.0
                        g[2] = Tfuel + 1.0

                        n_layers = length(t_cond)
                        for i = 1:n_layers
                                g[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ sum(t_cond)
                        end
                        g[end] = g[end] - 1.0 #fuselage wall temperature
                        return g
                end

                #Feasibility check before running the solver. The heat rate (boil-off) is largest for the most-reduced admissible stack
                #(every sizeable layer shifted by `Δt_min`). If even that stack leaks less than the target, no thickness shift can match the
                #target, so the solver is skipped and the stack is left in that state instead of letting the solver chase a non-physical root.
                #The reverse case is always feasible, since a thicker stack leaks less.
                Δt_min = min_insulation_shift(t_cond, iinsuldes)
                Q_max = max_heat_rate(fuse, fuse_tank, Δt_min, z, TSL, Mair, xftank_heat, ifuel)
                Q_target = boiloff_heat_rate(fuse_tank)
                infeasible = isfinite(Q_max) && Q_max < Q_target #A failed thermal solve (NaN) is inconclusive, so the solver still runs

                if infeasible
                        boiloff_max = boiloff_percent * Q_max / Q_target #Boil-off rate (%/h) attained by the thinnest stack
                        @warn "tanksize!: the boil-off target of $(boiloff_percent) %/h cannot be reached by this insulation stack " *
                              "(highest attainable: $(round(boiloff_max, sigdigits = 3)) %/h). The sizeable layers are reduced as far as the " *
                              "uniform thickness shift allows (the thinnest layer is set to $(T_INSUL_MIN) m). Check the insulation materials " *
                              "and base thicknesses, or the target." maxlog=3
                        Δt = Δt_min
                else
                        warm = length(fuse_tank.last_MLI_sol) == length(t_cond) + 2
                        guess = warm ? copy(fuse_tank.last_MLI_sol) : cold_guess()
                        #The previous Δt has already been added to t_cond, so the remaining correction starts at zero;
                        #only the temperature states are worth warm-starting.
                        warm && (guess[1] = 0.0)

                        #Forward-mode AD for the Jacobian, except when a vacuum layer is present: in that case the residual has multiple roots 
                        #(including unphysical ones) and finite differences and AD can converge to different ones, so the
                        #original finite-difference Jacobian is kept to preserve results.
                        use_ad = !check_vacuum(fuse_tank.material_insul)

                        #Solve non-linear problem with NLsolve.jl with looser tolerances and iteration cap.
                        sol = nlsolve_autodiff(MLI_AUTODIFF, "MLI-thickness", residual, guess, use_ad = use_ad, ftol = 1e-4, iterations = 500)

                        if sol.x_converged || sol.f_converged
                                Δt = sol.zero[1]
                                # cache
                                fuse_tank.last_MLI_sol = copy(sol.zero)
                        else
                                # retry once from the generic guess (stale warm start may be the cause) with tighter tolerances
                                sol_retry = nlsolve_autodiff(MLI_AUTODIFF, "MLI-thickness", residual, cold_guess(), use_ad = use_ad, ftol = 1e-5, iterations = 1000)
                                if sol_retry.x_converged || sol_retry.f_converged
                                        Δt = sol_retry.zero[1]
                                        fuse_tank.last_MLI_sol = copy(sol_retry.zero)
                                elseif warm
                                        # reuse last converged
                                        Δt = 0.0
                                else
                                        error("tanksize!: MLI thickness nonlinear solve failed to converge and no warm-start available")
                                end
                        end
                end

                for ind in iinsuldes #For every segment whose thickness can be changed
                        t_cond[ind] = t_cond[ind] + Δt #This will modify fuse_tank.t_insul
                end
                clamp_insulation_thickness!(t_cond, iinsuldes) #Keep the layers physical
                
        end
        
        #------Calculate tank weight------
        thickness_insul = sum(t_cond)
        
        #Evaluate tank weight
        inner_tank        = size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        Winnertank = inner_tank.Wtank
        Winsul_sum = inner_tank.Winsul_sum
        Vfuel      = inner_tank.Vfuel
        Rinnertank = inner_tank.Rtank_outer
        l_inner    = inner_tank.l_tank
        lcyl1      = inner_tank.l_cyl
        AR_head    = inner_tank.AR_head
        #Store in fuse_tank
        fuse_tank.Rinnertank  = Rinnertank
        fuse_tank.l_inner     = l_inner
        fuse_tank.l_cyl_inner = lcyl1
        fuse_tank.Shead_insul = inner_tank.Shead_insul

        has_vacuum = check_vacuum(fuse_tank.material_insul) #check if there is a vacuum layer

        if has_vacuum #If tank is double-walled
                #The outer vessel encloses the inner vessel and its insulation. When the inner vessel has its maximum size this is 
                #Rfuse - clearance; when it is still smaller (spherical regime) the outer vessel shrinks with it.
                Routertank = Rinnertank + thickness_insul
                lcyl2 = lcyl1 * Routertank / Rinnertank #Scale outer vessel length for geometric similarity
                Winner_tot = Winnertank + Wfuelintank #weight of inner vessel + fuel

                Ninterm = optimize_outer_tank(fuse, fuse_tank, Winner_tot, lcyl2; Router = Routertank, AR_head = AR_head) #Find optimal number of intermediate stiffeners
                
                fuse_tank.Ninterm = Ninterm #Store in fuse_tank to use as guess in next _size_aircraft! iteration

                outer_tank = size_outer_tank(fuse, fuse_tank, Winner_tot, lcyl2, Ninterm; Router = Routertank, AR_head = AR_head)

                Wtank = Winnertank + outer_tank.Wtank
                ltank = outer_tank.l_outer #If there is an outer vessel, the total length is the length of this tank
                Rtank = Routertank

        else
                Wtank = Winnertank
                ltank = l_inner #Tank length when there is only an inner vessel
                Rtank = Rinnertank
        end

        #------Store outputs in arrays------
        parg[igWfmax] = Vfuel * rhofuel * gee * nftanks #If more than one tank, max fuel capacity is nftanks times that of one tank
        parg[igWftank] = nftanks * Wtank #total weight of fuel tanks (including insulation)
        parg[iglftank] = ltank
        parg[igRftank] = Rtank
        parg[igWinsftank] = nftanks * Winsul_sum #total weight of insulation in fuel tanks

        #Tank placement and weight moment
        lcabin = fuse.layout.l_cabin_cylinder
        if tank_placement == "front"
                flag_front = 1
                flag_aft = 0
                xftank = fuse.layout.x_start_cylinder + 1.0*ft_to_m + ltank/2.0
                xftankaft = 0.0
        elseif tank_placement == "rear"
                flag_front = 0
                flag_aft = 1
                xftank = 0.0
                xftankaft = fuse.layout.x_start_cylinder + lcabin + 1.0*ft_to_m + ltank/2.0
        elseif tank_placement == "both"
                flag_front = 1
                flag_aft = 1
                xftank = fuse.layout.x_start_cylinder + 1.0*ft_to_m + ltank/2.0
                xftankaft = fuse.layout.x_start_cylinder + 1.0*ft_to_m + ltank + 1.0*ft_to_m + lcabin + 1.0*ft_to_m + ltank/2.0
        end
        
        parg[igxftank] = xftank
        parg[igxftankaft] = xftankaft
        parg[igxWftank] = Wtank * (flag_front * xftank + flag_aft * xftankaft) 
        xfuel = (flag_front * xftank + flag_aft * xftankaft) / (flag_front + flag_aft)
        parg[igxWfuel] = parg[igWfuel] * xfuel
        
end

"""
        res_MLI_thick(x::Vector{Float64}, fuse::Fuselage, fuse_tank::fuselage_tank, z::Float64, TSL::Float64, Mair::Float64, xftank::Float64, ifuel::Int64)

This function evaluates the residual vector for a given state containing change in wall thickness, heat transfer rate and 
insulation interface temperatures.

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `x::Float64`: vector with states
        - `fuse::Fuselage`: fuselage object.
        - `fuse_tank::fuselage_tank`: fuselage tank object.
        - `z::Float64`: flight altitude (m)
        - `TSL::Float64`: sea-level temperature (K)
        - `Mair::Float64`: external air Mach number
        - `xftank::Float64`: longitudinal coordinate of fuel tank centroid from nose (m)
        - `ifuel::Int64`: fuel index.

        **Outputs:**
        - `res::Vector{Float64}`: residuals vector.
"""
function res_MLI_thick(x::AbstractVector{<:Real}, fuse::Fuselage, fuse_tank::fuselage_tank, z::Float64, TSL::Float64, Mair::Float64, xftank::Float64, ifuel::Int64)

        #Extract parameters from fuse_tank
        boiloff_percent = fuse_tank.boiloff_rate
        t_cond = fuse_tank.t_insul
        qfac = fuse_tank.qfac
        iinsuldes = fuse_tank.iinsuldes
        Tfuel = fuse_tank.Tfuel
        Wfuel = fuse_tank.Wfuelintank
        h_v = fuse_tank.hvap #heat of vaporization

        # Extract states
        Δt = x[1]
        x_thermal = x[2:end]

        #Prepare to call residual_Q for thermal-related residuals
        t_all = Vector{eltype(x)}(t_cond) #copy (with the element type of x, for ForwardDiff) to avoid modifying input

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_all[ind] = t_all[ind] + Δt
        end
        check_insulation_thickness(t_all) #Reject infeasible (negative) trial thicknesses

        inner_tank = size_inner_tank(fuse, fuse_tank, t_all)

        Q = boiloff_heat_rate(fuse_tank) #Heat rate corresponding to the boil-off target (boil-off equals heat rate over heat of vaporization)

        #Assemble struct with parameters for residual_Q
        p = thermal_params{typeof(fuse.layout.cross_section), eltype(x)}()
        p.Q = Q #Store heat rate as it is known
        p.l_cyl = inner_tank.l_cyl
        p.l_tank = inner_tank.l_tank
        p.r_tank = inner_tank.Rtank_outer
        p.Shead = inner_tank.Shead_insul
        p.t_cond = t_all
        p.material = fuse_tank.material_insul
        p.Tfuel = Tfuel
        p.z = z
        p.TSL = TSL
        p.Mair = Mair
        p.xftank = xftank
        p.ifuel = ifuel
        p.fuse_cs = fuse.layout.cross_section

        res = residuals_Q(x_thermal, p, "Q_known") #Find thermal-related residuals
        return res
end 

"""
        check_vacuum(material_insul::Vector{String})

This function checks if any of the insulation layers requires a vacuum.

!!! details "🔃 Inputs and Outputs"
        **Inputs:**
        - `material_insul::Vector{ThermalInsulator}`: vector with layer materials

        **Outputs:**
        - `has_vacuum::Bool`: flag for vacuum, true if a vacuum is needed.
"""
function check_vacuum(material_insul::Vector{ThermalInsulator})
        vacuum_materials = ["vacuum", "microspheres"] #currently supported options are vacuum and microspheres

        has_vacuum = false
        for material in material_insul
                if lowercase(material.name) in vacuum_materials
                        has_vacuum = true
                end
        end

        return has_vacuum
end