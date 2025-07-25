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
                        catch #Return some high residual if it fails
                                return ones(length(x))*1e3
                        end
                end
                #Assemble guess for non linear solver
                #x[1] = Δt; x[2] = T_tank; x[3:(end-1)]: T at edge of insulation layer; x[end] = T at fuselage wall
                guess = zeros(length(t_cond) + 2) 
                guess[1] = 0.0
                guess[2] = Tfuel + 1.0

                n_layers = length(t_cond)
                for i = 1:n_layers
                        guess[i + 2] = Tfuel + ΔT * sum(t_cond[1:i])/ sum(t_cond)
                end
                guess[end] = guess[end] - 1.0 #fuselage wall temperature
                #Solve non-linear problem with NLsolve.jl
                sol = nlsolve(residual, guess, ftol = 1e-7) 
                Δt = sol.zero[1] 

                for ind in iinsuldes #For every segment whose thickness can be changed
                        t_cond[ind] = t_cond[ind] + Δt #This will modify fuse_tank.t_insul
                end
                
        end
        
        #------Calculate tank weight------
        thickness_insul = sum(t_cond)
        
        #Evaluate tank weight
        Winnertank, Winsul_sum, Vfuel, Shead_insul, Rinnertank, l_inner, lcyl1 = size_inner_tank(fuse, fuse_tank, fuse_tank.t_insul)
        #Store in fuse_tank
        fuse_tank.Rinnertank = Rinnertank 
        fuse_tank.l_inner = l_inner 
        fuse_tank.l_cyl_inner = lcyl1 
        fuse_tank.Shead_insul = Shead_insul 

        has_vacuum = check_vacuum(fuse_tank.material_insul) #check if there is a vacuum layer

        if has_vacuum #If tank is double-walled
                Routertank = fuse.layout.radius - fuse_tank.clearance_fuse
                lcyl2 = lcyl1 * Routertank / Rinnertank #Scale outer vessel length for geometric similarity
                Winner_tot = Winnertank + Wfuelintank #weight of inner vessel + fuel

                Ninterm = optimize_outer_tank(fuse, fuse_tank, Winner_tot, lcyl2) #Find optimal number of intermediate stiffeners
                
                fuse_tank.Ninterm = Ninterm #Store in fuse_tank to use as guess in next _size_aircraft! iteration

                Wtank2, Wcyl2, Whead2, Wstiff2, Souter, Shead2, Scyl2, 
                t_cyl2, t_head2, l_outer = size_outer_tank(fuse, fuse_tank, Winner_tot, lcyl2, Ninterm)

                Wtank = Winnertank + Wtank2
                ltank = l_outer #If there is an outer vessel, the total length is the length of this tank
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
function res_MLI_thick(x::Vector{Float64}, fuse::Fuselage, fuse_tank::fuselage_tank, z::Float64, TSL::Float64, Mair::Float64, xftank::Float64, ifuel::Int64)

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
        t_all = deepcopy(t_cond) #copy to avoid modifying input

        for ind in iinsuldes #For every segment whose thickness can be changed
                t_all[ind] = t_all[ind] + Δt
        end

        Wtank, Winsul_sum, Vfuel, Shead, r_tank, l_tank, l_cyl = size_inner_tank(fuse, fuse_tank, t_all)

        mdot_boiloff = boiloff_percent *  Wfuel / (gee * 100) / 3600  
        # Boil-off rate equals the heat rate divided by heat of vaporization
        Q_net = mdot_boiloff * h_v  # Heat rate from ambient to cryo fuel, including extra heat leak from valves etc as in eq 3.20 by Verstraete
        Q = Q_net / qfac

        #Assemble struct with parameters for residual_Q
        p = thermal_params{typeof(fuse.layout.cross_section)}()
        p.Q = Q #Store heat rate as it is known
        p.l_cyl = l_cyl
        p.l_tank = l_tank
        p.r_tank = r_tank
        p.Shead = Shead
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