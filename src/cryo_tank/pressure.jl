"""
Structure with parameters for tank pressure calculation.
"""
struct tank_params
    """Initial mixture in the tank"""
    mixture_init::SaturatedMixture
    """Tank volume (m^3)"""
    V::Float64
    """Tank maximum pressure (Pa)"""
    pmax::Float64
    """Quality of mixture being drawn from tank"""
    xout::Float64
    """Quality of mixture being vented"""
    xvent::Float64
    """Pressure fudge factor"""
    α::Float64
end

"""
Structure with input time-dependent functions used for tank pressure calculation.
"""
struct tank_inputs{F1, F2, F3}
    """Heat rate as a function of time (W)"""
    Q_calc::F1
    """Work rate as a function of time (W)"""
    W_calc::F2
    """Liquid mass flow rate as a function of time (kg/s)"""
    mdot_calc::F3
end

"""
    dpdt(mixture, Q, W, mdot, xout, V)

This function returns the time derivatives for the pressure in a cryogenic tank.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mixture::SaturatedMixture`: liquid/vapor mixture in tank
    - `Q::Float64`: heat rate (W)
    - `W::Float64`: work rate (W)
    - `mdot::Float64`: imposed mass flow rate out of the tank, e.g., to burn in combustor (kg/s)
    - `xout::Float64`: quality of mass flow rate out of tank
    - `mdot_vent::Float64`: required venting mass flow rate to keep tank pressure constant (kg/s)
    - `xvent::Float64`: quality of mixture that is vented
    - `V::Float64`: tank volume (m^3) 
    - `α::Float64`: fudge factor to account for effect of stratification
    
    **Outputs:**
    - `dp_dt::Float64`: pressure derivative in time (Pa/s)
"""
function dpdt(mixture::SaturatedMixture, Q::Float64, W::Float64, mdot::Float64, xout::Float64, mdot_vent::Float64, xvent::Float64, V::Float64, α::Float64 = 1.0)

    dp_dt = α * mixture.ϕ / V * (Q + W - mdot * mixture.hvap * (xout + mixture.ρ_star) - mdot_vent * mixture.hvap * (xvent + mixture.ρ_star))

    return dp_dt
end

"""
    dβdt(mixture, dp_dt, mdot, V)

This function returns the time derivatives for the liquid fill volume fraction in a cryogenic tank.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mixture::SaturatedMixture`: liquid/vapor mixture in tank
    - `dp_dt::Float64`: pressure derivative (Pa/s)
    - `mdot_tot::Float64`: total mass flow rate out of the tank, including venting (kg/s)
    - `V::Float64`: tank volume (m^3) 
    
    **Outputs:**
    - `dβ_dt::Float64`: fill fraction derivative in time (1/s)
"""
function dβdt(mixture::SaturatedMixture, dp_dt::Float64, mdot_tot::Float64, V::Float64)
    dρ_dt = -mdot_tot / V #Density derivative in tank

    gas = mixture.gas
    liquid = mixture.liquid
    dβ_dt = (dρ_dt - mixture.β * liquid.ρ_p * dp_dt - (1 - mixture.β) * gas.ρ_p * dp_dt) / (liquid.ρ - gas.ρ)
    return dβ_dt
end

"""
    venting_mass_flow(mixture, Q, W, mdot, xout, xvent)

This function returns the mass flow rate that needs to be vented to keep the tank pressure constant.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mixture::SaturatedMixture`: liquid/vapor mixture in tank
    - `Q::Float64`: heat rate (W)
    - `W::Float64`: work rate (W)
    - `mdot::Float64`: mass flow rate out of the tank (kg/s)
    - `xout::Float64`: quality of mass flow rate out of tank
    - `xvent::Float64`: quality of the mixture that is vented
    
    **Outputs:**
    - `mdot_vent::Float64`: venting mass flow rate (kg/s)
"""
function venting_mass_flow(mixture::SaturatedMixture, Q::Float64, W::Float64, mdot::Float64, xout::Float64, xvent::Float64)
    mdot_vent = (Q + W - mdot * mixture.hvap * (xout + mixture.ρ_star)) / (mixture.hvap * (xvent + mixture.ρ_star))
    return max(mdot_vent, 0.0)
end

"""
    mdot_boiloff(β, dβ_dt, dp_dt, ρl, ρl_p, mdot, V)

This function returns the rate of mass boiloff in the tank.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `mixture::SaturatedMixture`: liquid/vapor mixture in tank
    - `dβ_dt::Float64`: fill fraction derivative in time (1/s)
    - `dp_dt::Float64`: pressure derivative (Pa/s)
    - `mdot_liq::Float64`: liquid-phase mass flow rate out of the tank (kg/s)
    - `V::Float64`: tank volume (m^3) 
    
    **Outputs:**
    - `mdot_boiloff::Float64`: rate of liquid mass boiloff (kg/s)
"""
function mdot_boiloff(mixture::SaturatedMixture, dβ_dt::Float64, dp_dt::Float64, mdot_liq::Float64, V::Float64)
    β = mixture.β
    ρl = mixture.liquid.ρ
    ρl_p = mixture.liquid.ρ_p

    mdot_boiloff = -(dβ_dt * V * ρl + β * V * dp_dt * ρl_p + mdot_liq)
    return mdot_boiloff
end

"""
    TankDerivatives(t, y, u, params)

This function returns the time derivatives for pressure and liquid volume fill fraction in a cryogenic tank.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `t::Float64`: time (s)
    - `y::Vector{Float64}`: vector with the states; `y[1]` is pressure, `y[2]` is fill fraction, y[3] is tank mass, 
    y[4] is total liquid mass extracted, y[5] is the total mass that is vented and y[6] is the total mass that is boiled off.
    - `u::Struct`: structure with inputs; functions for heat rate, work rate and mass flow rate 
    - `params::Struct`: structure with parameters; including initial tank mixture, tank max pressure and tank volume 
    
    **Outputs:**
    - `dydt::Vector{Float64}`: vector with the state derivatives in time
"""
function TankDerivatives(t::Float64, y::SVector{6, Float64}, u::tank_inputs, params::tank_params)
    #Extract states 
    p = y[1]
    β = y[2]
    
    #Extract inputs
    Q = u.Q_calc(t)
    W = u.W_calc(t)
    mdot = u.mdot_calc(t)

    #Extract parameters
    mixture_init = params.mixture_init
    V = params.V
    pmax = params.pmax #Maximum tank pressure
    xout = params.xout #Only liquid is being drawn from the tank
    xvent = params.xvent #Only vent gas
    α = params.α #homogenous pressure fudge factor

    #Update mixture with current state
    mix_current = SaturatedMixture(mixture_init.species, p, β)

    #Make sure tank has liquid mass inside before taking mass out
    if β <= 0
        mdot = 0.0
    end

    #Calculate derivatives
    if p < pmax
        mdot_vent = 0.0
    else
        mdot_vent = venting_mass_flow(mix_current, Q, W, mdot, xout, xvent)
    end

    #Calculate liquid mass flow rate for boiloff
    mdot_liq = (1 - xout) * mdot + (1 - xvent) * mdot_vent

    #Total mass flow rate for beta
    mdot_tot = mdot + mdot_vent

    dp_dt = dpdt(mix_current, Q, W, mdot, xout, mdot_vent, xvent, V, α)
    dβ_dt = dβdt(mix_current, dp_dt, mdot_tot, V) #dβ/dt
    #Calculate derivatives and store them
    dydt = @SVector Float64[
        dp_dt, 
        dβ_dt,
        -mdot_tot, #dM_tank/dt
        mdot, #dMburn/dt
        mdot_vent, #dMvent/dt
        mdot_boiloff(mix_current, dβ_dt, dp_dt, mdot_liq, V) #dMboiloff/dt
    ]
    return dydt
end

"""
    calculate_venting_mass(tank_mass_vec, liquid_out_vec)

This function calculates how much mass is being vented out of a cryogenic tank. It returns both the cumulative mass 
vented out in time and the total vented mass at the end.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `tank_mass_vec::Vector{Float64}`: vector with the mass inside of the tank at different points in time (kg).
    - `liquid_out_vec::Vector{Float64}`: vector with the cumulative liquid mass that has been taken out of the tank (kg).

    **Outputs:**
    - `vented_mass_vec::Vector{Float64}`: vector with the cumulative mass that has been vented at different points in time (kg)
    - `vented_mass_tot:Float64`: total mass vented out of the tank (kg)
"""
function calculate_venting_mass(tank_mass_vec::Vector{Float64}, liquid_out_vec::Vector{Float64})
    vented_mass_vec =  tank_mass_vec[1] .- tank_mass_vec .- liquid_out_vec #vented mass is total mass out of tank minus liquid mass
    vented_mass_tot = vented_mass_vec[end]
    return vented_mass_vec, vented_mass_tot
end

"""
    calculate_boiloff_rate(times, Mboil_vec)

This function uses finite differences to calculate the boiloff rate from the vector with boil mass.
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `times::Vector{Float64}`: vector with time (s)
    - `Mboil_vec::Vector{Float64}`: vector with the cumulative liquid mass that has boiled off (kg).

    **Outputs:**
    - `mdot_boil::Vector{Float64}`: vector with the boiloff rate (kg/s)
"""
function calculate_boiloff_rate(times, Mboil_vec)
    mdot_boil = zeros(length(times))
    mdot_boil[1:(end - 1)] = diff(Mboil_vec) ./ diff(times)
    if length(mdot_boil) > 1 #Only do this if the vector has more than one element
        mdot_boil[end] = mdot_boil[end - 1]
    end
    return mdot_boil
end