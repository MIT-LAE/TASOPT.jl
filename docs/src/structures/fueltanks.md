# [Fuel tanks](@id fueltanks)

Liquid hydrocarbon fuel is assumed to be stored in the interior of the wings and no additional tanks are needed. The weight of the fuel is accounted for while sizing the wing structure. See [`structures.surfw`](@ref).

However, alternate fuels such as cryogenic liquid hydrogen require additional storage tanks that are insulated pressure vessels.

## Theory

!!! details "📖 Theory - Thermal and structural sizing of cryogenic fuel tanks" 
    ### Thermal design
    The fuel tanks in TASOPT are assumed to consist of circular cylinders with two hemiellipsoidal caps. There are two distinct layers: the tank itself, assumed to be made of an isotropic material with high thermal conductivity; and an insulation layer, which does not carry structural loads and has very low thermal conductivity. The insulation layer itself may consist of additional sublayers of different materials, forming a multi-layer insulation (MLI).

    As the insulation layer consists of two different geometries across which heat can be transferred (the cylinder and the hemiellipsoids), two slightly different models for thermal resistance must be used. In the case of heat transfer across a layer between two concentric cylinders, it can be shown from Fourier's equation that the thermal resistance, ``R_{cyl}``, is given by 
    ```math
        R_{cyl} = \frac{\ln\left( \frac{R_f}{R_0}\right)} {2\pi l_{cyl} k},
    ``` 
    where ``R_0`` is the layer's inner radius, ``R_f = R_0 + t`` is the outer radius (with ``t`` being the layer thickness), ``l_{cyl}`` is the cylinder length, and ``k`` is the thermal conductivity of the layer. For the hemiellipsoids, an approximate solution can be used, given by 
    ```math
        R_{ell} = \frac{t} {k \left(S_f + S_0 - \frac{S_{he}}{R^2} t^2\right)},
    ``` 
    where ``S_f`` is the surface area of the hemiellipsoid at the final radius, ``S_0`` is the surface area at the initial radius, and ``\frac{S_{he}}{R^2}`` is the ratio of hemiellipsoid surface area to radius squared; for example, this ratio is ``\frac{S_{he}}{R^2}=2\pi`` for a hemisphere. The equation above makes use of the fact that there are two hemiellipsoids at both ends of the tank and represents their total resistance. By combining these resistances, the total resistance of a layer of MLI can be found using parallel resistance addition,
    ```math
        R_{l} = \frac{R_{cyl}R_{ell}} {R_{cyl} + R_{ell}}.
    ``` 
    In an MLI, the total resitance across the insulation is the serial addition of the resitances across each layer,
    ```math
        R_{MLI} = \sum_i R_l^i.
    ``` 

    In addition to the insulation resistance, the convective (from fuel to tank wall and from exterior wall to freestream) and radiative heat transfers have to be taken into account. The heat transfer to the freestream can be modeled as having two components: radiation and convection. The radiative component has an equivalent heat transfer coefficient
    ```math
        h_{rad} = σ ε (T_a^2 + T_f^2) (T_a + T_f),
    ``` 
    where ``σ`` is the Stefan-Boltzmann constant, ``ε`` is the emissivity of the surface, ``T_a`` is the freestream temperature, and ``T_f`` is the temperaure of the fuel. Similarly, the heat transfer coefficient from froced convection from the external wall to the freestream can be modeled using the Chilton-Colburn analogy,
    ```math
        h_{convair} = \frac{C_f}{2 Pr^{2/3}}  ρ u c_p,
    ``` 
    where ``C_f`` is the skin-friction coefficient, ``Pr`` is the Prandtl number (``Pr\approx 0.85`` for turbulent air), ``ρ`` is the freestream air density, ``u`` is the freestream velocity, and ``c_p`` is the specific heat of the freestream air at constant pressure. The skin-friction coefficient can be modeled using a flat-plate solution,
    ```math
        C_f = \frac{0.455}{\log_{10}(\mathrm{Re}_x)^{2.58}  (1 + 0.144 M^2)^{0.65}},
    ```
    where ``M`` is the freestream Mach number and ``\mathrm{Re}_x`` is the Reynolds number at the location fo the fuel tank in the fuselage. The equivalent heat transfer coefficient to the freestream air is ``h_{air} = h_{convair}+h_{rad} ``, such that the equivalent resistance is
    ```math
        R_{air} = \frac{1}{h_{air} (2\pi l_{cyl} R_{fuse} +2 S_{he})},
    ```
    where ``l_{cyl}`` is the length of the cylindrical portion of the tank, ``R_{fuse}`` is the fuselage radius and ``S_{he}`` is the outer area of the hemiellipsoidal caps. 

    Inside the tank, there is a heat transfer from the bulk of the liquid fluid to the tank via natural convection. The Nusselt number for this heat transfer process can be modeled as
    ```math
        \mathrm{Nu}_l = 0.0605 \mathrm{Ra}_l^{1/3},
    ```
    where ``\mathrm{Ra}_l`` is the tank-length based Rayleigh number and is given by
    ```math
        \mathrm{Ra}_l = \frac{g \beta (T_w-T_f) l^3 Pr}{\nu^2},
    ```
    where ``g`` is the gravitational acceleration, ``\beta`` is the fuel's coefficient of thermal expansion, ``T_w`` is the temperature at the wall, ``l`` is the tank length, ``Pr`` is the Prandlt number of the liquid fuel, and ``\nu`` is the kinematic viscosity of the fuel. The thermal resistance due to natural convection is then
    ```math
        R_{liq} = \frac{l}{\mathrm{Nu}_l k S_{int}},
    ```
    where ``k`` is the thermal conductivity of the liquid fuel and ``S_{int}`` is the internal surface area of the tank. 

    The combined thermal resistance is ``R_{eq} = R_{liq} + R_{MLI} + R_{air}``, such that the total heat transfer rate is ``\dot{Q} = \frac{T_a - T_f}{R_{eq}}``. Once the heat transfer rate is known, the boiloff rate is simply ``\dot{m}_{boil}=\frac{\dot{Q}}{h_v}``, where ``h_v`` is the heat of vaporization of the fuel. 

    #### Notes on implementation
    In the current version of TASOPT, the desired boiloff rate (in percentage per hour) is an input and the thicknesses of some desired layers of the MLI insulation are changed until the desired boiloff rate is met. The non-linear solver in NLsolve.jl is used to find the change in layer thickness needed to meet this requirement. 

    In general, the thermal conductivity of insulation materials is a function of temperature. In the code, the average of the conductivities at the interface of each layer is used to find the resistance of the layer; if the conductivity changes linearly with temperature, this method is exact. However, this calculation requires knowing the temperatures at every interface in the MLI layer. The nonlinear solver in NLsolve is used to find these temperatures and the overall heat transfer rate.

    ### Structural design
    The structural part of the tank is sized for a given pressure difference between the high-pressure interior and the exterior. As the boiling temperature of a liquid is a function of temperature, it is preferable to keep the interior pressure constant. If the desired pressure difference is ``\Delta p``, the tank is actually designed for a pressure difference ``\Delta p_{des} = \alpha \Delta p``, where ``\alpha>1`` is a safety factor. The outer radius of the strcutrual portion of the tank is
    ```math
        R_{t,o} = R_{fuse} - d_{fclear} - t_{MLI},
    ```
    where ``d_{fclear} `` is the fuselage clearance distance and ``t_{MLI}`` is the total thickness of the MLI. The skin thickness of the cylindrical portion of the tank wall is sized using
    ```math
        t_{s,cyl} = \frac{2 \Delta p_{des} R_{t,o}}{2 \sigma_a f_{weld} + 0.8 \Delta p_{des}},
    ```
    where ``\sigma_a`` is the maximum allowable stress for the wall material and ``f_{weld}<1`` is a factor that accounts for structural weakening due to welding. The wall thickness of the hemiellipsoidal caps is given by
    ```math
        t_{s,cap} = \frac{2 \Delta p_{des} R_{t,o} K }{2 \sigma_a f_{weld} + 2 \Delta p_{des} (K- 0.1)},
    ```
    where ``K=\frac{1}{6}(AR^2+2)`` is a factor that accounts for the ellipsoidal aspect ratio (``AR``). Once the wall thicknesses have been determined, the internal tank radius is given by ``R_{t,i}=R_{t,o}-t_{s,cyl} ``. The volume required for the fuel is 
    ```math
        V_{fuel} = (1+f_{ull})\frac{m_{fuel} + m_{boil}}{\rho},
    ```
    where ``m_{boil}`` is the total fuel mass that boils off during flight and ``f_{ull}>0`` is a factor to account for the fact that the tank must contain some empty volume to account for ullage. The internal volume of a hemiellipsoidal cap is given by 
    ```math
        V_{cap} =  \frac{2πR_{t,i}^3}{3AR}.
    ```
    Therefore, the internal volume of the cylindrical portion of the tank is ``V_{cyl}=V_{fuel}-2V_{cap}``, and the length of the cylindrical portion can then be found to be
    ```math
        l_{cyl} =  \frac{V_{cyl}}{\pi R_{t,i}^2}.
    ```
    Once this length is known, the masses of the tank and insulation layers can be found from their respective volumes and densities.

```@docs
structures.tanksize
structures.tankWmech
structures.tankWthermal

```