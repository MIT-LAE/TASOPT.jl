module ElectricMachine

abstract type AbstractElectricMachine end
abstract type AbstractMagnets end

const μ₀ = 1.25663706127e-6 #N⋅A⁻² https://physics.nist.gov/cgi-bin/cuu/Value?mu0 

@kwdef mutable struct Motor
    rotor::RotorGeometry = RotorGeometry()
    stator::StatorGeometry = StatorGeometry()
    teeth::Teeth = Teeth()
    magnet::PermanentMagnet = PermanentMagnet()
    windings::Windings = Windings()

    """Design shaft power [W]"""
    P_design::Float64 = 1e6
    """Angular velocity [rad/s]"""
    Ω::Float64 = 10e3
    """Frequency [Hz]"""
    f::Float64 = 12e3
    """Slots per pole [-]"""
    Nsp::Int = 3
    """Phases [-]"""
    phases::Int = 3
    """Phase resistance [Ω]"""
    phase_resistance::Float64 = 0.0
    """Design current [A]"""
    Id::Float64 = 0.0
    """Design Voltage [V]"""
    Vd::Float64 = 0.0
    """Current [A]"""
    I::Float64 = 0.0
    """Voltage [V]"""
    V::Float64 = 0.0

    """Pole pairs [-]"""
    N_pole_pairs::Int = 8
    """Max rotor tip speed [m/s]"""
    U_max::Float64 = 200.0
    """Saturation Flux [Wb/m²]"""
    B_sat::Float64 = 1.8
    """Maximum Current Density [A/m²]"""
    J_max::Float64 = 1e6
    """Max ratio to saturation"""
    rB_sat::Float64 = 0.98
    """Stack length [m]"""
    l::Float64 = 1.0
    """Mean internal temperature [K],
     used to calculate air density and viscosity, and winding resistance"""
    T::Float64 = 273.15 + 90
    """Thickness of air gap [m]"""
    airgap_thickness::Float64 = 2e-3
    """Pole pairs [-]"""
    p::Int = 8
    """Mass of machine [kg]"""
    mass::Float64 = 0.0
end

@kwdef struct PermanentMagnet <: AbstractMagnets
    """Magnet thickness [m]"""
    thickness::Float64 = 20e-3
    """Magnet Density [kg/m⁻³]"""
    density::Float64 = 7501.0 #Neodymium magnets
    """Magnetization constant [A/m]"""
    M::Float64 = 8.604e5
    """Mass of magnets [kg]"""
    mass::Float64 = 0.0
end

@kwdef struct Windings
    """Winding conductor material"""
    conductor::Conductor
    """Winding insulator material"""
    insulator::Insulator
    """Number of turns [-]"""
    N_turns::Int = 1
    """Wire area [m²]"""
    A_wire::Float64 = 0.518e-6 #20AWG
    """Winding temperature [K]"""
    T::Float64 = 273.15 + 90
    """Packing factor [-]"""
    kpf::Float64 = 0.35
    """Mass of windings [kg]"""
    mass::Float64 = 0.0
end

@kwdef mutable struct Teeth
    """Tooth width [m]"""
    width::Float64 = 0.005
    """Tooth thickness [m]"""
    thickness::Float64 = 0.02
    """Inner radius [m]"""
    Ri::Float64 = 0.01
    """Flux density [Wb/m² = T]"""
    B::Float64 = 1.0
    """Material"""
    material::ElectricSteel = ElectricSteel("M19")
    """Mass [kg]"""
    mass::Float64 = 0.0
end

@kwdef mutable struct RotorGeometry
    """Rotor thickness [m]"""
    thickness::Float64 = 0.01
    """Inner radius [m]"""
    Ri::Float64 = 0.01
    """Outer radius [m]"""
    Ro::Float64 = 0.02
    """Flux density [Wb/m² = T]"""
    B::Float64 = 1.0
    """Material"""
    material::ElectricSteel = ElectricSteel("M19")
    """Mass [kg]"""
    mass::Float64 = 0.0
end

@kwdef mutable struct StatorGeometry
    """Stator thickness [m]"""
    thickness::Float64 = 0.01
    """Inner radius [m]"""
    Ri::Float64 = 0.01
    """Outer radius [m]"""
    Ro::Float64 = 0.02
    """Flux density [Wb/m² = T]"""
    B::Float64 = 1.0
    """Material"""
    material::ElectricSteel = ElectricSteel("M19")
    """Mass [kg]"""
    mass::Float64 = 0.0
end

@kwdef mutable struct ShaftGeometry
    """Inner Radius [m]"""
    Ri::Float64 = 0.01
    """Outer Radius [m]"""
    Ro::Float64 = 0.01
    """Overhang fraction of shaft [-]"""
    l_extra::Float64 = 1.2
    """Material"""
    material::AbstractMaterial
    """Mass [kg]"""
    mass::Float64 = 0.0
end
# Layout:
# 		
#          +-------------------------------+     <----------------+
#          |        Stator back Iron       | }ts                  |
#          |   +--+   +--+   +--+   +--+   |     <------------+   |
#   Teeth  |   |  |   |  |   |  |   |  | wt| }tt              |   |
#          *---*  *---*  *---*  *---*  *---*     <--------+   |   |
# air gap {                                               |   |   |
#   (g)    +---------------+---------------+   ^ <----+   |   |   |
#          | Magnets       |               |   tm     |   |   |   |
#          |               |               |   v      |   |   |   |
#          +---------------+---------------+     <-+  |   |   |   |
#          |                               |   ^   |  |   |   |   |
#          |         Rotor back Iron       |   tr  |  |   |   |   |
#          |                               |   v   |  |   |   |   |
#          +-------------------------------+       |  |   |   |   |
#                                          ^       |  |   |   |   |
#                                          |       |  |   |   |   |
#                                          +       +  +   +   +   +
#                                         Rri    Rro  Rg  Rt  Rsi Rso




"""
$TYPEDEF

Simplified permanent magnet synchronous machine (PMSM) sizing.
"""
function size_PMSM!(motor::Motor, shaft_speed::AbstractFloat, design_power::AbstractFloat)
    rotor = motor.rotor
    stator = motor.stator
    magnet = motor.magnet
    teeth = motor.teeth
    windings = motor.windings
    shaft = motor.shaft

    motor.Ω = shaft_speed * (2 * π / 60)
    motor.f = motor.N_pole_pairs * shaft_speed / 60

    #Outer radius of the magnet - subject to the highest rotational speeds
    radius_gap = motor.U_max / motor.Ω

    # Calculate maximum flux through the stator and rotor back iron
    # Really, this needs to be a constraint Bsbi ≤ Bsat (at the top level) rather 
    # than a prescribed setting...
    stator.B = motor.rB_sat * motor.B_sat
    rotor.B = motor.rB_sat * motor.B_sat

    B_gap = airgap_flux(magnet.M, magnet.thickness, motor.airgap_thickness)

    stator.thickness = B_gap / stator.B * π * radius_gap / (2 * motor.N_pole_pairs)
    rotor.thickness = B_gap / rotor.B * π * radius_gap / (2 * motor.N_pole_pairs)

    rotor.Ro = radius_gap - magnet.thickness
    rotor.Ri = rotor.Ro - rotor.thickness

    teeth.Ri = radius_gap + motor.airgap_thickness

    stator.Ri = teeth.Ri + teeth.thickness
    stator.Ro = stator.Ri + stator.thickness

    N_teeth = motor.Nsp * 2 * motor.N_pole_pairs

    A_teeth = N_teeth * teeth.thickness * teeth.width
    A_slots = π * (teeth.Ri + stator.Ri) * teeth.thickness - A_teeth
    A_slot = A_slots / N_teeth
    λ = A_slots / (A_slots + A_teeth)
    l_end_turns = π / (2 * motor.N_pole_pairs) * teeth.Ri * (λ / sqrt(1 - λ^2))

    #Armature reaction
    B_windings = μ₀ * motor.J_max * teeth.thickness

    teeth.B = sqrt((B_gap / λ)^2 + B_windings^2)
    if teeth.B > motor.B_sat
        @warn "Flux density in teeth exceeds saturation flux density"
    end

    torque = design_power / motor.Ω
    slot_current = J_max * windings.kpf * A_slot #Peak slot current
    windings.N_turns = floor(windings.kpf * A_slot / windings.A_wire)

    #Assume 2 phases excited at any given time
    motor.I = motor.Id = (2 / motor.phases * motor.Nsp * 2 * motor.N_pole_pairs) * slot_current
    force_length = motor.Id * B_gap 
    motor.l = torque / (force_length * radius_gap)

    #Size shaft
    shaft.Ro = rotor.Ri
    if shaft.Ro^4 > torque / shaft.material.τmax * 2 * shaft.Ro / π
        shaft.Ri = (shaft.Ro^4 - torque / shaft.material.τmax * 2 * shaft.Ro / π)^0.25
    else
        # shaft can be solid as well
        shaft.Ri = 0.0
    end

    #Masses
    rotor.mass = cross_sectional_area(rotor) * motor.l * rotor.material.density
    stator.mass = cross_sectional_area(stator) * motor.l * stator.material.density
    magnet.mass = cross_sectional_area(radius_gap, rotor.Ro) * motor.l * magnet.material.density
    teeth.mass = A_teeth * motor.l * teeth.material.density

    total_winding_volume = A_slots * (motor.l + 2 * l_end_turns)
    windings.mass = windings.kpf * total_winding_volume * windings.conductor.material.density +
                    (1 - windings.kpf) * total_winding_volume * windings.insulator.material.density

    shaft.mass = cross_sectional_area(shaft) * (motor.l * shaft.l_extra) * shaft.material.density

    motor.mass = rotor.mass + stator.mass + magnet.mass + teeth.mass + windings.mass + shaft.mass
    #size_rotor!(rotor, R_gap)
    #size_stator!()
    #size_windings!()

    ## Calculate losses
    windings.T = motor.T # Set temperature
    motor.phase_resistance = motor.Nsp / motor.phases * 2 * motor.N_pole_pairs *
                             slot_resistance(motor.windings,
                                 windings.kpf * A_slot,
                                 motor.l + 2 * l_end_turns)


end  # function size_PMSM!

"""
"""
function calculate_losses(EM::Motor)
    Q_ohmic = ohmic_loss(EM.I, EM.phase_resistance)
    Q_core = core_loss(EM)
    Q_wind = windage_loss(EM)
    return Q_ohmic + Q_core + Q_wind
end  # function calculate_losses

"""

Returns the power dissipated due to resistive heating assuming a certain 
number of phases are energized at any given time. 
"""
function ohmic_loss(I::AbstractFloat, phase_resistance::AbstractFloat;
    energized_phases::Int=2)
    return I^2 * (energized_phases * phase_resistance)
end  # function ohmic_loss

core_loss(motor) = hysteresis_loss(motor) + eddy_loss(motor)

"""
"""
function hysteresis_loss(motor::Motor)
    return hysteresis_loss(motor.stator, motor.f) + hysteresis_loss(motor.teeth, motor.f)
end  # function hysteresis_loss

function hysteresis_loss(steel, f)
    return steel.mass * steel.material.kₕ * f * steel.B^steel.material.α
end

"""
"""
function eddy_loss(motor::Motor)
    return eddy_loss(motor.stator, motor.f) + eddy_loss(motor.teeth, motor.f)
end #function eddy_loss

function eddy_loss(steel, f)
    return steel.mass * steel.material.kₑ * f^2 * steel.B^2    
end  # function eddy_loss
    
"""
"""
function windage_loss(motor::Motor, air)
    #air from IdealGases?
    Re = motor.Ω * motor.radius_gap * motor.airgap_thickness * air.ρ/air.μ
    Cf = 0.0725*Re^-0.2 #avg. skin friction approx. from fully turbulent flat plate
    return Cf*π*air.ρ*motor.Ω^3 * motor.radius_gap^4 * motor.l
end  # function windage_loss

function slot_resistance(wind::Windings, A, l)
    resist = resistivity(wind.conductor, wind.T)
    return resist * l / A
end  # function slot_resistance

"""
"""
function cross_sectional_area(Ro, Ri)
    return π * (Ro^2 - Ri^2)
end  # function cross_section

# Might not need to define each of this and just let anything go
# cross_sectional_area(X) = cross_sectional_area(X.Ro, X.Ri)
# if it doesn't have the right fields it'll throw an error.
cross_sectional_area(R::RotorGeometry) = cross_sectional_area(R.Ro, R.Ri)
cross_sectional_area(R::StatorGeometry) = cross_sectional_area(R.Ro, R.Ri)
cross_sectional_area(R::ShaftGeometry) = cross_sectional_area(R.Ro, R.Ri)


"""
"""
function size_rotor!(rotor::RotorGeometry, R_gap::AbstractFloat)
    Ro = R_gap - magnet.thickness
    rotor.thickness = (B_gap / rotor.B) * π * R_gap / (2p)
    rotor.Ri = Ro

end  # function size_rotor!

"""
    airgap_flux(M, thickness, airgap)

Given the magnetization constant (M), thickness of the magnet, and the 
air gap thickness, returns the flux density in the air gap. 
Simplified expression from a magnetic circuit analysis. Assumes no MMF drop 
across the steel. 
"""
function airgap_flux(M, thickness, airgap)
    B_gap = μ₀ * M * thickness / (thickness + airgap)
end  # function airgap_flux

end