# [Homogeneous model](@id CryoTank)

In a tank containing a boiling liquid, such as a cryogenic liquid that has reached its saturation point, two fluid states exist: a liquid phase and a gas phase (or vapor). In a real tank, there would be some level of stratification leading to temperature differences between the two phases, which exist at the same pressure. This section describes the models used to simulate the evolution in time of the conditions inside a cryogenic tank. It is based on the homogeneous tank model[^1][^2][^3], which treats the fluid inside the tank as a well-mixed saturated mixture at homogeneous temperature and pressure.

## Theory

!!! details "📖 Theory - Homogeneous model" 
    ## Ordinary-differential equation system
    Let us consider a tank that contains a mixture of saturated liquid and gas phases at a saturation pressure ``p``. The tank has a heat rate ``\dot{Q}`` into it and a work rate ``\dot{W}`` is performed in the tank. The tank also has a mass flow rate ``\dot{m}`` out of the tank of a fluid with quality ``x_{out}``. In general, the tank may also be vented with a venting mass flow rate ``\dot{m}_{vent}`` of fluid with a quality ``x_{vent}``. From conservation of energy, it can be shown that the derivative in time of the saturation pressure in this homogeneous case is[^2]
    ```math
        \left( \frac{dp}{dt}\right)_h = \frac{\phi}{V}\left[\dot{Q} + \dot{W} - \dot{m} h_{vap} (x_{out} + \rho ^*) - \dot{m}_{vent} h_{vap} (x_{vent} + \rho ^*)\right],
    ```
    where ``\phi`` is the energy derivative of the fluid mixture (to be defined later), ``V`` is the tank volume, ``h_{vap}`` is the enthalpy of vaporization of the liquid, and ``\rho ^*=\frac{\rho_g}{\rho_l-\rho_g}`` is the density ratio (``l`` refers to the properties of the liquid phase and ``g`` to those of the gaseous phase). The enthalpy of vaporization is simply the enthalpy difference between the gaseous and liquid phases.

    Experiments suggest that the pressure derivative in a real tank can be significantly greater than that of the homogeneous model due to stratification. Authors have accounted for this with a scaling factor ``\alpha``, such that[^2][^3]
    ```math
        \frac{dp}{dt} = \alpha \left( \frac{dp}{dt}\right)_h. 
    ```

    If the liquid volume fraction (fill fraction) if the tank is ``\beta``, the density of the mixture can be expressed as 
    ```math
        \rho = \beta \rho_l + (1-\beta) \rho_g.
    ```
    Differentiating this in time,
    ```math
        \frac{d\rho}{dt} = \frac{d \beta}{dt} (\rho_l - \rho_g) + \beta \frac{d\rho_l}{dp}\frac{dp}{dt} + (1-\beta) \frac{d\rho_g}{dp} \frac{dp}{dt},
    ```
    where the liquid- and gas-phase densities are functions of pressure only since the fluid is in saturated conditions. The derivative in time of density in the tank is equal to 
    ```math
        \frac{d\rho}{dt} = -\frac{\dot{m} + \dot{m}_{vent}}{V}.
    ```
    Solving for the rate of change of the fill volume in the tank,
    ```math
        \frac{d\beta}{dt} = \frac{1}{\rho_l - \rho_g}\left[\frac{d\rho}{dt} - \frac{dp}{dt} \left(\beta \frac{d\rho_l}{dp}+(1-\beta) \frac{d\rho_g}{dp}\right)\right].
    ```
    The evolution in time of properties inside the tank can be modeled by solving the system of ordinary differential equations (ODEs) given by
    ```math
        \frac{d}{dt} \mathbf{x} = \mathbf{f}(t, \mathbf{x}, \mathbf{u}, \mathbf{p}),
    ```
    where ``\mathbf{x} = [p, \beta]^T``, ``\mathbf{u} = [\dot{Q}(t), \dot{W}(t), \dot{m}(t)]^T`` and ``\mathbf{p} = [V, \alpha, p_{vent}]^T``, where ``p_{vent}`` is the pressure at which the tank is vented. These ODEs can be integrated in time from knowledge of the initial conditions, namely the initial pressure and initial tank fill fraction.

    ### Energy derivative
    The energy derivative ``\phi`` is given by[^1]
    ```math
        \phi = \frac{1}{\rho \left(\frac{\partial u}{\partial p}\right)_\rho},
    ```
    where ``u`` represents the internal energy of the mixture. The partial derivative ``\left(\frac{\partial u}{\partial p}\right)_\rho`` is taken at constant mixture density. The internal energy of the mixture is given by
    ```math
        u = x u_g + (1-x) u_l,
    ```
    where ``x`` represents the quality of the mixture (the mass fraction of the gaseous phase). Differentiating this in pressure at constant mixture density,
    ```math
        \left(\frac{\partial u}{\partial p}\right)_\rho = x \frac{du_g}{dp} + (1-x)\frac{du_l}{dp} + (u_g - u_l) \left(\frac{\partial x}{\partial p}\right)_\rho.
    ```
    The density of the mixture can be related to the quality by[^1]
    ```math
       \frac{1}{\rho} = \frac{x}{\rho_g} + \frac{1-x}{\rho_l}.
    ```
    Differentiating in pressure at constant mixture density yields
    ```math
       \left(\frac{\partial x}{\partial p}\right)_\rho = \left(\frac{1}{1/\rho_g - 1/\rho_l}\right)\left(\frac{x}{\rho_g^2}\frac{d\rho_g}{dp} + \frac{1-x}{\rho_l^2}\frac{d\rho_l}{dp}\right),
    ```
    which completes the expression for the energy derivative.

    ### Venting
    Venting may be needed when the tank pressure has reached some maximum pressure. The purpose of mass venting from the tank is to prevent the tank pressure from increasing any further. To find the required venting mass flow rate, the equation for ``\frac{dp}{dt}`` can be set to zero whenever the venting pressure ``p_{vent}`` is reached or exceeded, yielding,
    ```math
       \dot{m}_{vent} = \begin{cases} 
       0, & p<p_{vent}\\ 
       \frac{\dot{Q} + \dot{W}}{h_{vap} (x_{vent} + \rho^*)} - \frac{x_{out} + \rho^*}{x_{vent} + \rho ^*}\dot{m}, &p\geq p_{vent}
       \end{cases}.
    ```

## Functions
### Thermodynamic properties
```@docs
CryoTank.gas_properties
```
```@docs
CryoTank.liquid_properties
```
### Saturated mixtures
```@docs
CryoTank.SaturatedGas
```
```@docs
CryoTank.SaturatedLiquid
```
```@docs
CryoTank.SaturatedMixture
```
```@docs
CryoTank.update_pβ!
```
### Pressure evolution
```@docs
CryoTank.dpdt
```
```@docs
CryoTank.dβdt
```
```@docs
CryoTank.venting_mass_flow
```
```@docs
CryoTank.mdot_boiloff
```
```@docs
CryoTank.TankDerivatives
```
### TASOPT interfacing
```@docs
CryoTank.find_mdot_time
```
```@docs
CryoTank.calc_Q_points
```
```@docs
CryoTank.find_Q_time_interp
```
```@docs
CryoTank.analyze_TASOPT_tank
```

[^1]: Lin, Chin S., Neil T. Van Dresar, and Mohammad M. Hasan. "Pressure control analysis of cryogenic storage systems." Journal of propulsion and power 20.3 (2004): 480-485.
[^2]: Verstraete, Dries, et al. "Hydrogen fuel tanks for subsonic transport aircraft." International journal of hydrogen energy 35.20 (2010): 11085-11098.
[^3]: Winnefeld, Christopher, et al. "Modelling and designing cryogenic hydrogen tanks for future aircraft applications." Energies 11.1 (2018): 105.