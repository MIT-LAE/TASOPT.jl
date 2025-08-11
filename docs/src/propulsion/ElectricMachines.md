# Electric Machines

## Permanent-magnet synchronous machines

Permanent-magnet synchronous machines (PMSMs), particularly motors and generators, can be modeled to estimate their performance and weight. The function [`size_PMSM!()`](@ref propsys.ElectricMachine.size_PMSM!) can be used to find the thickness and lengths of the PMSM components, their masses, as well as to calculate the phase electrical resistance. Once a PMSM has been sized, its off-design performance can be computed using [`operate_PMSM!()`](@ref propsys.ElectricMachine.operate_PMSM!), which calculates power losses and computes the input electrical power required (motor) or output shaft power (generator).

!!! details "📖 Theory - Permanent-magnet synchronous machines"
    The model for PMSMs is based on that in Dowdle et al.[^1], with some modifications for increased fidelity. The PMSM consists of:
    - A **rotor back iron** connected to a shaft.
    - A **magnet** with a number of pole pairs (`p`) attached to the rotor back iron.
    - A thin **air gap** separating the rotor from the stator.
    - A **stator** consisting of metallic teeth supported by a back iron. Slots between the teeth contain Litz wires.

    This implementation supports an arbitrary number of electrical phases, although the source model was designed for three-phase power.
    ### PMSM sizing
    The following parameters are required for PMSM sizing:
    - Number of pole pairs (``p``).
    - Magnet thickness (``t_M``).
    - Air gap thickness (``t_\mathrm{gap}``).
    - Maximum rotor linear velocity (``U_\mathrm{max}``).
    - Ratio of metal flux density to saturation flux density (``r_\mathrm{sat}``).
    - Number of slots per pole (``N_{sp}``).
    - Teeth thickness (``t_\mathrm{teeth}``).
    - Shaft speed (``\Omega``).
    - Shaft power (``P_\mathrm{shaft}``).
    - Maximum wire current density (``j_\mathrm{max}``).
    - Slot packing factor (``k_{pf}``)
    - Component materials.

    The radius at the gap start ``R_\mathrm{gap}`` is calculated from the maximum rotational speed as
    ```math
        R_\mathrm{gap} = \frac{U_\mathrm{max}}{\Omega}.
    ```
    The air gap magnetic field flux density, ``B_\mathrm{gap}`` can be computed from the magnet properties as[^1]
    ```math
        B_\mathrm{gap} = \frac{\mu_0 M t_M}{t_M + t_\mathrm{gap}},
    ```
    where ``\mu_0 `` is the vacuum permeability and ``M`` is the magnet's magnetization constant. The magnetic flux densities in the rotor and stator back irons can be computed by
    ```math
        B_x = r_\mathrm{sat} B_\mathrm{sat},
    ```
    where ``x`` refers to the rotor or stator back irons and ``B_\mathrm{sat}`` is the back iron material's saturation flux density. The stator and rotor back iron thicknesses can then be computed by 
    ```math
        t_x = \frac{B_\mathrm{gap} \pi R_\mathrm{gap}}{B_x 2 p}.
    ```
    From the geometry, the rotor back iron outer radius is simply ``R_{r,o}=R_\mathrm{gap}-t_M`` and its inner radius is ``R_{r,i}=R_{r,o}-t_{r}``, where ``t_{r}`` is the rotor back iron thickness.

    The magnetic flux density produced by the windings on the teeth is given by[^1]
    ```math
        B_\mathrm{wind} = \mu_0 j_\mathrm{max} k_{pf} t_\mathrm{teeth},
    ```
    where ``k_{pf}`` is the slot packing factor. This adds vectorially to the teeth flux due to the permanent magnets,
    ```math
        B_\mathrm{teeth}^2 = B_\mathrm{wind}^2 + \left(B_\mathrm{gap}\frac{A_\mathrm{ann}}{A_\mathrm{slots}}\right)^2,
    ```
    where ``A_\mathrm{ann} = \pi\left((R_{s,i}+t_\mathrm{teeth})^2-R_{s,i}^2\right)`` is the area of the annulus where the teeth and slots are located and ``A_\mathrm{slots}`` is the total area of the slots. If the flux density on the teeth is given by ``B_\mathrm{teeth} = r_\mathrm{sat} B_\mathrm{sat}``, the slot total area is 
    ```math
        A_\mathrm{teeth} = A_\mathrm{ann} - A_\mathrm{ann}\frac{B_\mathrm{gap}}{\sqrt{B_\mathrm{teeth}^2 - B_\mathrm{wind}^2}}  
    ```
    Since the number of teeth is equal to the number of slots and is ``N_\mathrm{teeth}=2pN_{sp}``, the width of a single tooth is ``w_\mathrm{tooth} = \frac{A_\mathrm{teeth}}{t_\mathrm{teeth} N_\mathrm{teeth}}``. The area of a single slot, ``A_\mathrm{slot}`` is given by 
    ```math
        A_\mathrm{slot} = \frac{A_\mathrm{ann} - A_\mathrm{teeth}}{N_\mathrm{teeth}}. 
    ```
    The current in a slot, ``I``, is given by
    ```math
        I = j_\mathrm{max} A_\mathrm{slot} k_{pf}. 
    ```

    The length of the PMSM, ``l``, can be computed from knowledge of the torque ``T=P_\mathrm{shaft}/\Omega``,
    ```math
        l = \frac{T}{R_\mathrm{gap} I B_\mathrm{gap} N_{es}}, 
    ```
    where ``N_{es}`` is the number of energized slots. The masses of the different components can then be estimated from their respective volumes and densities.

    ### PMSM operation

    #### Motor Operation
    For motors, the input power is calculated as:
    ```math
    P_\mathrm{input} = P_\mathrm{shaft} + \dot{Q}_\mathrm{loss},
    ```
    where ``\dot{Q}_\mathrm{loss}`` includes all relevant power losses.

    #### Generator Operation
    For generators, the output power is calculated as:
    ```math
    P_\mathrm{output} = P_\mathrm{shaft} - \dot{Q}_\mathrm{loss}.
    ```
    ### Loss Calculations
    The total power loss in a PMSM is the sum of ohmic, core, and windage losses,
    ```math
    \dot{Q}_\mathrm{loss} = \dot{Q}_\mathrm{ohmic} + \dot{Q}_\mathrm{core} + \dot{Q}_\mathrm{wind}.
    ```
    
    **Ohmic Losses**:
    ```math
    Q_\mathrm{ohmic} = I^2 R_\mathrm{phase} N_\mathrm{ep},
    ```
    where ``I`` is the current, ``R_\mathrm{phase}`` is the phase resistance, and ``N_\mathrm{ep}`` is the number of energized phases.

    **Core Losses**:
    Core losses include hysteresis and eddy current losses:
    ```math
    \dot{Q}_\mathrm{core} = \dot{Q}_\mathrm{hysteresis} + \dot{Q}_\mathrm{eddy}.
    ```

    **Hysteresis Losses**:
    ```math
    \dot{Q}_\mathrm{hysteresis} = k_h f B^α,
    ```
    where ``k_h`` is the hysteresis coefficient, ``f`` is the frequency, ``B`` is the flux density, and ``α`` is a material constant.

    **Eddy Current Losses**:
    ```math
    \dot{Q}_\mathrm{eddy} = k_e f^2 B^2,
    ```
    where ``k_e`` is the eddy current coefficient.

    **Windage Losses**:
    Windage losses are caused by air friction and are calculated using Vrancik's model[^2],
    ```math
    \dot{Q}_\mathrm{wind} = C_f \pi \rho \Omega^3 R_\mathrm{gap}^4 l,
    ```
    where ``C_f`` is the skin friction coefficient, ``\rho`` is the air density, ``\Omega`` is the angular velocity, ``R_\mathrm{gap}`` is the gap radius, and ``l`` is the stack length.

    #### Back e.m.f.
    For motors, the voltage ``V`` corresponding to the power demand of the motor is calculated as
    ```math
        V = \frac{P_\mathrm{input}}{N_\mathrm{phases} \frac{2}{\pi} I N_\mathrm{inverters}},
    ```
    where ``N_\mathrm{phases}`` is the number of phases and ``N_\mathrm{inverters}`` is the number of inverters used to power the motor. A similar expression is used to calculate the output voltage of the generator, without the number of inverters.

## Inverters
An inverter is an electronic device which converts a direct current into an alternating current. This usually occurs through the use of switches which alternate the current flow. A simplified sizing model for an inverter, [`size_inverter!()`](@ref propsys.ElectricMachine.size_inverter!), calculates the mass of an inverter from its output power and the user-specified specific power (power per unit mass). The efficiency of the inverter can be calculated with [`operate_inverter!()`](@ref propsys.ElectricMachine.operate_inverter!), which uses a simplified efficiency model based on those in Enders[^3] and Faranda et al.[^4] to compute the input power to the inverter.

## Functions
```@docs
propsys.ElectricMachine.size_PMSM!
```
```@docs
propsys.ElectricMachine.operate_PMSM!
```
```@docs
propsys.ElectricMachine.size_inverter!
```
```@docs
propsys.ElectricMachine.operate_inverter!
```


[^1]: Dowdle, Aidan P., David K. Hall, and Jeffrey H. Lang. "Electric propulsion architecture assessment via signomial programming." 2018 AIAA/IEEE Electric Aircraft Technologies Symposium (EATS). IEEE, 2018.
[^2]: Vrancik, James E. Prediction of windage power loss in alternators. No. NASA-TN-D-4849. 1968.
[^3]: Enders, Wilhelm. Development of electrical powertrain simulation methods for hybrid- and turboelectric commercial aircraft design. Master's thesis, 2020.
[^4]: Faranda, Roberto S., et al. "The optimum PV plant for a given solar DC/AC converter." Energies 8.6 (2015): 4853-4870.