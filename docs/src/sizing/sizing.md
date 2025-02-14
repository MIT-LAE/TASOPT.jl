# Design and evaluation

## [Sizing the aircraft] (@id sizing)

The aircraft is sized via a fixed point iteration for the design mission ([`TASOPT.wsize`](@ref)). The performance of the design can be evaluated for an off-design mission ([`TASOPT.fly_off_design!`](@ref)).

[`TASOPT.wsize`](@ref) is typically the driving script in an analysis, as is the case in the `size_aircraft!` call (as demonstrated in the [first example] (@ref firstexample)). The sizing analysis calls the various performance routines (e.g., `fusebl`, `get_wing_weights`, `cdsum`, `mission`, etc.) as shown in the [TASOPT flowchart](@ref flowchart).

!!! details "🖥️ Code structure - Aircraft sizing" 
    The aircraft-sizing function requires an `aircraft` object as input. This object is unpacked into storage arrays and other component objects, such as wing, fuselage or engine. The eventual aim is to eliminate all data storage array and replace them by component objects but this is still work in progress. The main major function called within `wsize` is [`TASOPT.fusebl!`](@ref), which calculates the fuselage boundary layer properties and drag coefficients for start-of-cruise; these are then used in other mission points. 

    `wsize` then uses simplified methods to initialize the relevant aircraft weights and parameters, unless the user specifies otherwise with an optional input ('init_weight=true'). The bulk of the computational cost and time is spent in the weight sizing loop.

    After the weight sizing loop is completed, the aircraft takeoff performance and field lengths are calculated using [`TASOPT.takeoff!`](@ref).

    ### Weight sizing loop
    `TASOPT.wsize` performs a fixed point iteration by sequentially running weight and performance models for the different aircraft components. This is done via a `for` loop that gets terminated once the maximum aircraft weight has converged within a desired tolerance. Unfortunately, there is no assurance that this process will result in a converged aircraft; most failures correspond to infeasible aircraft but it is also possible that the iteration fails to size feasible ones if the starting guess is poor.

    The first weight that gets calculated inside the sizing loop is the fuselage weight, through [`TASOPT.fusew!`](@ref). After this is calculated, the total maximum takeoff weight gets recomputed and there is a check for whether the sizing loop is terminated. If weight has not converged, the loop continues.

    The wing geometry is set by running [`TASOPT.set_wing_geometry!`](@ref) and the wing pitching moments are computed through [`TASOPT.surfcm`](@ref). The horizontal and vertical tail geometry is computed through [`TASOPT.tailpo!`](@ref). Finally, the weights of the three aerodynamic surfaces are calculated by running [`TASOPT.get_wing_weights!`](@ref).

    If the aircraft requires an insulated fuel tank in the fuselage, for example, if the fuel is cryogenic, the tank is sized using [`TASOPT.tanksize!`](@ref); this function calculates the structural weight and sizes the thermal insulation. For details on how the fuel tank is sized, see [Fuel tanks](@ref fueltanks). The sized tank dimensions are then use to recalculate the fuselage geometry to accommodate the tank in [`TASOPT.update_fuse!`](@ref). 

    The weight and balance of the aircraft at start-of-cruise is adjusted using [`TASOPT.balance`](@ref). This can move the wing, resize the horizontal tail, or change the tail trim to achieve a desired metric for longitudinal stability (e.g., a set static margin).  

    The total drag at start-of-cruise, which is the engine design point, is calculated using [`aerodynamics.cdsum!`](@ref). This function calls a combination of models for the drag of aerodynamic surfaces, engine nacelle, and induced drag at the Trefftz plane. 

    The engines are sized at the start-of-cruise to produce a total thrust force equal to the aircraft drag. The generalized function `ac.engine.enginecalc!` is used for this purpose. This field stores a user defined function for the engine performance. The only engine option with explicit support in TASOPT is a two-spool turbofan engine, although the user is free to use alternative models by modifying the engine object. The engine functions get called via a wrapper, [`TASOPT.tfwrap!`](@ref), which in turns calls the engine calculation function `tfcalc`.

    Once the engines are sized, the fuel demand at every point in the mission is calculated using [`TASOPT.mission!`](@ref). This function in turn recalculates the balance, drag, and engine performance at every point. Further details on mission are provided below. `mission!` is usually the greatest time sink in an aircraft sizing. The weight gets updated after running `mission` and the loop restarts.

## [Mission evaluation] (@id mission)
The function [`TASOPT.mission!`](@ref) contains the fuel burn calculation for the entire mission. It can be used both in sizing, as part of the iteration to obtain a converged aircraft, or in off-design, to calculate the performance of an already-designed airplane. 

!!! details "🖥️ Code structure - Mission" 
    The `mission!` function simulates the entire mission of an aircraft, calculating fuel burn and other mission variables. 

    The function sets the initial conditions for takeoff, including altitude, temperature, pressure, and density. It then interpolates the lift coefficient over the climb and descent points to ensure smooth transitions between different phases of the mission. 

    Next, the function estimates the takeoff speed and sets the velocity and Reynolds number over the climb and descent points. This involves calculating the takeoff lift coefficient and using it to estimate the takeoff speed. The Reynolds number is then calculated based on the takeoff speed and other atmospheric conditions.

    The function proceeds to calculate the climb and descent parameters using aerodynamic and engine performance models. It integrates the climb and descent trajectories using a predictor-corrector scheme to update the range, time, and weight fractions. 

    Once the climb and descent parameters are set, the function sets the conditions for the cruise phase, including altitude, speed, and fuel consumption. It calculates the fuel burn and weight fractions for the entire mission via calls to `engine.enginecalc!`, and adds any vented fuel. This involves adjusting the aircraft's balance and trim settings via calls to `balance` to ensure stability throughout the mission, and recalculating the drag via `cdsum!`.

    The function also sets up the climb points at equal altitude intervals and calculates the available thrust assuming max throttle climb. It initializes the climb integrands and integrates the trajectory over the climb phase. The function calculates the cruise-climb angle based on available thrust and atmospheric conditions.

    The descent phase is then set up, with the function interpolating the descent points and integrating the time and weight over the descent. It calculates the velocity, Mach number, and Reynolds number for each descent point and adjusts the pitch trim by adjusting the horizontal tail lift coefficient.

    Finally, the function calculates the mission fuel fractions and weights, including the weight of any vented fuel. It updates the mission parameters, such as the takeoff weight and fuel weight, and calculates the mission payload-fuel energy intensity (PFEI).

## [Off-design performance] (@id missionexec)
A sized aircraft's mission performance can be obtained using `mission`, along with operation constraints via a pitch trim calculation (`balance`) and balanced field length calculation (`takeoff!`). The function [`TASOPT.fly_off_design!`](@ref) performs this calculation: it runs the aircraft through a mission with different range, payload, and conditions to the design mission.

## Function documentation
```@docs
TASOPT.wsize

TASOPT.fly_off_design!

TASOPT.size_aircraft!

TASOPT.fusebl!

TASOPT.set_wing_geometry!

TASOPT.surfcm

TASOPT.tailpo!

TASOPT.get_wing_weights!

TASOPT.update_fuse!

TASOPT.tfwrap!

TASOPT.mission!(ac, imission, Ldebug)

TASOPT.takeoff!(ac; printTO)

```

