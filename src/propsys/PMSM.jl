module ElectricMachine

abstract type AbstractElectricMachine end
abstract type AbstractMagnets end

struct Motor <: AbstractElectricMachine end
struct Generator <: AbstractElectricMachine end

import ..propsys: __TASOPTindices__, __TASOPTroot__

using ..materials
using ..engine
using DocStringExtensions
using Roots

include(joinpath(__TASOPTroot__,"utils/constants.jl"))
include("inverter.jl")
include("cable.jl")

@kwdef mutable struct PermanentMagnet <: AbstractMagnets
    """Magnet thickness [m]"""
    thickness::Float64 = 20e-3
    """Magnet Density [kg/m⁻³]"""
    ρ::Float64 = 7501.0 #Neodymium magnets
    """Magnetization constant [A/m]"""
    M::Float64 = 8.604e5
    """Mass of magnets [kg]"""
    mass::Float64 = 0.0
end

@kwdef mutable struct Windings
    """Winding conductor material"""
    conductor::Conductor = Conductor("Cu")
    """Winding insulator material"""
    insulator::Insulator = Insulator("PTFE")
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
    material::StructuralAlloy = StructuralAlloy("AISI-4340")
    """Mass [kg]"""
    mass::Float64 = 0.0
end

@kwdef mutable struct Air
    """Gas name"""
    gas::String = "air"
    """Temperature [K]"""
    T::Float64 = 363.15
    """Pressure [Pa]"""
    p::Float64 = 101325
end

#Override air to return other gas properties
function Base.getproperty(obj::Air, sym::Symbol)
    if sym === :ρ #Density
        R, _, _, _, _, _ = gasPr(obj.gas, obj.T)
        return obj.p /(R * obj.T)
    elseif sym === :μ
        _, _, _, _, μ, _  = gasPr(obj.gas, obj.T)
        return μ
    elseif sym === :ν
        return obj.μ / obj.ρ
    else
        return getfield(obj, sym)
    end
end  # function Base.getproperty

"""
$TYPEDEF

$TYPEDFIELDS

Structure that defines a permanent magnet synchronous machine (PMSM).
"""
@kwdef mutable struct PMSM{T<:AbstractElectricMachine}
    type::T = T()
    
    rotor::RotorGeometry = RotorGeometry()
    stator::StatorGeometry = StatorGeometry()
    teeth::Teeth = Teeth()
    magnet::PermanentMagnet = PermanentMagnet()
    windings::Windings = Windings()
    shaft::ShaftGeometry = ShaftGeometry()
    air::Air = Air()

    """Design shaft power [W]"""
    P_design::Float64 = 250e3
    """Shaft power [W]"""
    P::Float64 = 250e3
    """Angular velocity [rad/s]"""
    Ω::Float64 = 10e3
    """Slots per pole [-]"""
    Nsp::Int64 = 3
    """Phases [-]"""
    phases::Int64 = 3
    """Energized phases [-]"""
    energized_phases::Int64 = 2
    """Phase resistance [Ω]"""
    phase_resistance::Float64 = 0.0
    """Design current in a slot [A]"""
    Id::Float64 = 0.0
    """Design Voltage [V]"""
    Vd::Float64 = 0.0
    """Current in a slot [A]"""
    I::Float64 = 0.0
    """Voltage or back emf [V]"""
    V::Float64 = 0.0

    """Pole pairs [-]"""
    N_pole_pairs::Int = 8
    """Max rotor tip speed [m/s]"""
    U_max::Float64 = 200.0
    """Saturation Flux [Wb/m²]"""
    B_sat::Float64 = 1.8
    """Maximum Current Density [A/m²]"""
    J_max::Float64 = 5e6
    """Max ratio to saturation"""
    rB_sat::Float64 = 0.98 #Default assumption is that the metal is almost saturated at design
    """Stack length [m]"""
    l::Float64 = 0.0
    """Area of a single slot [m]"""
    A_slot::Float64 = 0.0
    """Radius of gap start (top of magnet) [m]"""
    radius_gap::Float64 = 0.0
    """Thickness of air gap [m]"""
    airgap_thickness::Float64 = 2e-3
    """Mass of machine [kg]"""
    mass::Float64 = 0.0

    """Number of multi-phase inverters [-]"""
    N_inverters::Int64 = 1
end

#Override PMSMBase to return other properties
function Base.getproperty(obj::PMSM, sym::Symbol)
    if sym === :f #Electric frequency
        return obj.Ω * obj.N_pole_pairs/ (2 * pi)
    elseif sym === :torque #Torque
        return obj.P / obj.Ω
    elseif sym === :N_slots #Number of slots
        return obj.Nsp * 2 * obj.N_pole_pairs
    elseif sym === :N_slots_per_phase #Number of slots divided by number of phases
        return (obj.Nsp * 2 * obj.N_pole_pairs / obj.phases) 
    elseif sym === :N_energized_slots #Energized slots
        return (obj.energized_phases * obj.N_slots_per_phase) 
    elseif sym === :λ #Geometric factor
        return obj.N_slots * obj.A_slot / cross_sectional_area(obj.teeth.Ri + obj.teeth.thickness, obj.teeth.Ri) 
    elseif sym === :B_gap #Air gap flux density
        return airgap_flux(obj.magnet.M, obj.magnet.thickness, obj.airgap_thickness)
    elseif (sym === :P_input && typeof(obj) == PMSM{Motor}) # Input electrical power
        #Assumes that multiple inverters are powering the motor
        return obj.phases * 2/pi * obj.I * obj.V * obj.N_inverters 
    elseif (sym === :P_output && typeof(obj) == PMSM{Generator}) # Output electrical power
        return obj.phases * 2/pi * obj.I * obj.V
    else
        return getfield(obj, sym)
    end
end  # function Base.getproperty

# Define convenience constructors
PMSM_generator() = PMSM{Generator}()
PMSM_motor() = PMSM{Motor}()

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
    size_PMSM!(PMSM::AbstractElectricMachine, shaft_speed::AbstractFloat, design_power::AbstractFloat)

Simplified permanent magnet synchronous machine (PMSM) sizing.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `PMSM::AbstractElectricMachine`: motor or generator object
    - `shaft_speed::Float64`: shaft rotational speed [rpm]
    - `design_power::Float64`: design shaft power [W]
    
    **Outputs:**
    No direct outputs. `PMSM` object gets modified with the input power and mass.
"""
function size_PMSM!(PMSM::PMSM, shaft_speed::AbstractFloat, design_power::AbstractFloat)
    rotor = PMSM.rotor
    stator = PMSM.stator
    magnet = PMSM.magnet
    teeth = PMSM.teeth
    windings = PMSM.windings
    shaft = PMSM.shaft

    PMSM.Ω = shaft_speed * (2 * π / 60)
    PMSM.P_design = design_power
    PMSM.P = design_power

    #Outer radius of the magnet - subject to the highest rotational speeds
    PMSM.radius_gap = PMSM.U_max / PMSM.Ω

    #-------Rotor and stator sizing-------
    # Calculate maximum flux through the stator and rotor back iron
    # Really, this needs to be a constraint Bsbi ≤ Bsat (at the top level) rather 
    # than a prescribed setting. Or rB_sat can be a global optimization variable.
    # The higher Bsbi, the thinner the sbi, but higher Bsbi ⟹ higher core losses
    stator.B = PMSM.rB_sat * PMSM.B_sat
    rotor.B = PMSM.rB_sat * PMSM.B_sat

    B_gap = PMSM.B_gap #Air gap flux density

    stator.thickness = B_gap / stator.B * π * PMSM.radius_gap / (2 * PMSM.N_pole_pairs)
    rotor.thickness = B_gap / rotor.B * π * PMSM.radius_gap / (2 * PMSM.N_pole_pairs)

    rotor.Ro = PMSM.radius_gap - magnet.thickness
    rotor.Ri = rotor.Ro - rotor.thickness
    if rotor.Ri < 0.0
        error("Rotor radius too small")
    end

    teeth.Ri = PMSM.radius_gap + PMSM.airgap_thickness

    stator.Ri = teeth.Ri + teeth.thickness
    stator.Ro = stator.Ri + stator.thickness

    #-------Teeth sizing-------
    #Armature reaction
    B_windings = μ₀ * PMSM.J_max * windings.kpf * teeth.thickness
    #TODO: Dowdle uses confusing notation for the maximum current density; it is possible
    # that the kpf factor is not needed here
    
    #Teeth B-field
    teeth.B = PMSM.rB_sat * PMSM.B_sat

    #Size teeth
    A_teeth_annulus = cross_sectional_area(teeth.Ri + teeth.thickness, teeth.Ri) #Annulus area where teeth lie
    A_teeth = A_teeth_annulus - B_gap * A_teeth_annulus / sqrt(teeth.B^2 - B_windings^2) #Analytical solution
    N_teeth = PMSM.N_slots #Number of teeth
    
    teeth.width = A_teeth /(N_teeth * teeth.thickness)

    #-------Slot sizing-------
    A_slots = A_teeth_annulus - A_teeth
    PMSM.A_slot = A_slots / PMSM.N_slots
    λ = A_slots / (A_slots + A_teeth)
    l_end_turns = π / (2 * PMSM.N_pole_pairs) * teeth.Ri * (λ / sqrt(1 - λ^2))

    #-------Calculate length-------
    torque = design_power / PMSM.Ω
    slot_current = PMSM.J_max * windings.kpf * PMSM.A_slot #Peak slot current
    windings.N_turns = floor(windings.kpf * PMSM.A_slot / windings.A_wire)

    PMSM.I =
        PMSM.Id = slot_current
    force_length = PMSM.Id * PMSM.N_energized_slots * B_gap 
    PMSM.l = torque / (force_length * PMSM.radius_gap)

    #-------Size shaft-------
    shaft.Ro = rotor.Ri
    if shaft.Ro^4 > torque / shaft.material.τmax * 2 * shaft.Ro / π
        shaft.Ri = (shaft.Ro^4 - torque / shaft.material.τmax * 2 * shaft.Ro / π)^0.25
    else
        error("Shaft radius too small")
    end

    if PMSM.Ω^2*shaft.Ri^2*shaft.material.ρ > shaft.material.YTS
        error("Stress due to rotational speed exceeds shaft yield strength")
    end

    #-------Compute masses-------
    #TODO this could go into the structure definitions by overriding Base.getproperty
    rotor.mass = cross_sectional_area(rotor) * PMSM.l * rotor.material.ρ
    stator.mass = cross_sectional_area(stator) * PMSM.l * stator.material.ρ
    magnet.mass =
        cross_sectional_area(PMSM.radius_gap, rotor.Ro) * PMSM.l * magnet.ρ
    teeth.mass = A_teeth * PMSM.l * teeth.material.ρ

    total_winding_volume = A_slots * (PMSM.l + 2 * l_end_turns)
    windings.mass =
        windings.kpf * total_winding_volume * windings.conductor.ρ +
        (1 - windings.kpf) * total_winding_volume * windings.insulator.ρ

    shaft.mass =
        cross_sectional_area(shaft) * (PMSM.l * shaft.l_extra) * shaft.material.ρ

    PMSM.mass =
        rotor.mass + stator.mass + magnet.mass + teeth.mass + windings.mass + shaft.mass

    #-------Calculate phase resistance-------
    windings.T = PMSM.air.T # Set temperature
    PMSM.phase_resistance =
        PMSM.N_slots_per_phase *
        slot_resistance(PMSM.windings, windings.kpf * PMSM.A_slot, PMSM.l + 2 * l_end_turns)

    #-------Calculate design voltage-------
    operate_PMSM!(PMSM, shaft_speed, design_power)
    PMSM.Vd = PMSM.V #Store design voltage

end  # function size_PMSM!

"""
    operate_PMSM!(motor::Motor, shaft_speed::AbstractFloat, shaft_power::AbstractFloat)

Runs a motor with a given shaft speed and power and calculates the back emf voltage. 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `motor::Motor`: motor object
    - `shaft_speed::Float64`: shaft rotational speed [rpm]
    - `design_power::Float64`: design shaft power [W]
    
    **Outputs:**
    No direct outputs. `motor` object gets modified with the input power.
"""
function operate_PMSM!(motor::PMSM{Motor}, shaft_speed::AbstractFloat, shaft_power::AbstractFloat)
    motor.Ω = shaft_speed * (2 * π / 60)
    B_gap = motor.B_gap
    motor.P = shaft_power

    torque = shaft_power/motor.Ω
    motor.I = torque / (motor.radius_gap * B_gap * motor.l * motor.N_energized_slots)

    Q_loss = calculate_losses(motor)

    P_elec = shaft_power + Q_loss #Total required power
    #TODO the phase calculations in the source paper are sketchy. Verify this or do it properly without "energized" phases
    motor.V = P_elec / (motor.phases * 2/pi * motor.I) / motor.N_inverters #Compute required back emf voltage
    
end  # function operate_PMSM!

"""
    operate_PMSM!(generator::Generator, shaft_speed::AbstractFloat, shaft_power::AbstractFloat)

Runs a generator with a given shaft speed and power and calculates the back emf voltage. 

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `generator::Generator`: motor object
    - `shaft_speed::Float64`: shaft rotational speed [rpm]
    - `design_power::Float64`: design shaft power [W]
    
    **Outputs:**
    No direct outputs. `generator` object gets modified with the input power.
"""
function operate_PMSM!(generator::PMSM{Generator}, shaft_speed::AbstractFloat, shaft_power::AbstractFloat)
    generator.Ω = shaft_speed * (2 * π / 60)
    B_gap = generator.B_gap
    generator.P = shaft_power

    torque = shaft_power/generator.Ω
    generator.I = torque / (generator.radius_gap * B_gap * generator.l * generator.N_energized_slots)

    Q_loss = calculate_losses(generator)

    P_elec = shaft_power - Q_loss #Total required power
    #TODO the phase calculations in the source paper are sketchy. Verify this or do it properly without "energized" phases
    generator.V = P_elec / (generator.phases * 2/pi * generator.I) #Compute produced voltage
    
end  # function operate_PMSM!

"""
Calculates the power losses in a PMSM, including ohmic, core and windage losses.
"""
function calculate_losses(EM::PMSM)
    Q_ohmic = ohmic_loss(EM.I, EM.phase_resistance, EM.energized_phases)
    Q_core = core_loss(EM)
    Q_wind = windage_loss(EM)

    return Q_ohmic + Q_core + Q_wind #Power loss
end  # function calculate_losses

"""

Returns the power dissipated due to resistive heating assuming a certain 
number of phases are energized at any given time. By default assumes that 2 
phases are energized at any given time.
"""
function ohmic_loss(
    I::AbstractFloat,
    phase_resistance::AbstractFloat,
    energized_phases::Int = 2,
)
    return I^2 * (energized_phases * phase_resistance)
end  # function ohmic_loss

core_loss(motor) = hysteresis_loss(motor) + eddy_loss(motor)

"""
$TYPEDSIGNATURES

Calculates the hysteresis loss in the motor (i.e., in the rotor + the stator)
"""
function hysteresis_loss(motor::PMSM)
    return hysteresis_loss(motor.stator, motor.f) + hysteresis_loss(motor.teeth, motor.f)
end  # function hysteresis_loss

"""
Calculates the hysteresis losses in steel. The steel needs to be of type
[`ElectricSteel`](@ref).
"""
function hysteresis_loss(steel, f)
    return steel.mass * steel.material.kₕ * f * steel.B^steel.material.α
end

"""
$TYPEDSIGNATURES

Calculates the eddy losses in the motor (i.e., due to the rotor + the stator)
"""
function eddy_loss(motor::PMSM)
    return eddy_loss(motor.stator, motor.f) + eddy_loss(motor.teeth, motor.f)
end #function eddy_loss

"""
Calculates the eddy losses in an electrical steel. The steel needs to be of type
[`ElectricSteel`](@ref).
"""
function eddy_loss(steel, f)
    return steel.mass * steel.material.kₑ * f^2 * steel.B^2
end  # function eddy_loss

"""
Calculates the windage (i.e., air friction losses) for the rotor.
From Vrancik (1968) - Prediction of windage power loss in alternators
"""
function windage_loss(motor::PMSM)
    #air from IdealGases?
    Re = motor.Ω * motor.radius_gap * motor.airgap_thickness / motor.air.ν
    res(Cf) = 1/sqrt(Cf) - 2.04 - 1.768*log(Re*sqrt(Cf)) #From Vrancik, residual to be zeroed
    Cf = find_zero(res, 1e-2) #Solve for skin friction coefficient
    return Cf * π * motor.air.ρ * motor.Ω^3 * motor.radius_gap^4 * motor.l
end  # function windage_loss

"""
Returns the slot resistance given a winding temperature and material.
The winding material needs to be of type [`Conductor`](@ref).
"""
function slot_resistance(wind::Windings, A, l)
    resist = resistivity(wind.conductor, wind.T)
    return resist * l / A
end  # function slot_resistance

"""
$TYPEDSIGNATURES

Returns the cross sectional area of an annulus with inner radius Rᵢ and outer radius Rₒ.
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
    airgap_flux(M, thickness, airgap)

Given the magnetization constant (M), thickness of the magnet, and the 
air gap thickness, returns the flux density in the air gap. 
Simplified expression from a magnetic circuit analysis. Assumes no MMF drop 
across the steel. 
"""
function airgap_flux(M, thickness, airgap)
    B_gap = μ₀ * M * thickness / (thickness + airgap)
    return B_gap
end  # function airgap_flux

"""
Linear model for magnet remanent flux as a function of temperature.
This will be important if we model the TMS system along with the PMSM.
"""
function remanent_flux(magnet::AbstractMagnets, T)
    magnet.remanent_flux * (1 - magnet.α * (T - (273.15 + magnet.Tbase)) / 100.0)
end  # function remanent_flux

end