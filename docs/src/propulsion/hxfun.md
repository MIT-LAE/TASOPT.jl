# [Heat exchangers](@id hxers)

To study novel aircraft designs that may leverage cryogenic hydrogen, the following heat exchanger models were integrated, along with [cryogenic tank](@ref fueltanks) and [fuel cell](@ref fuelcells) models.

## Theory

!!! details "📖 Theory - Tubular heat exchangers" 
    ### Effectiveness–NTU method
    There are many different heat exchanger (HX) topologies. The HX currently implemented in TASOPT consists of staggered tubes in a cross flow; this geometry was selected because it is simple to integrate into the space between two concentric cylinders. The HX is designed using the effectiveness–NTU method, described below [^1].  
    
    Consider a heat exchanger in which a coolant gas is used to cool or heat gas. The minimum and maximum heat capacity rates are defined as
    ```math
    C_{\mathrm{min}} = \mathrm{min} (\dot{m}_c c_{p,c}, \dot{m}_p c_{p,p})
    ```
    ```math
    C_{\mathrm{max}} = \mathrm{max} (\dot{m}_c c_{p,c}, \dot{m}_p c_{p,p}),
    ```
    where ``\dot{m}`` is the mass flow rate, ``c_p`` is the specific heat at constant pressure, and the subscripts ``p`` and ``c`` refer to the process-side and coolant streams. The capacity ratio is ``C_r = \frac{C_{\mathrm{min}}}{C_{\mathrm{max}}}``. 
    
    A maximum heat transfer rate is defined as
    ```math
    \dot{Q}_{max} = C_{\mathrm{min}} (T_{i,p} - T_{i,c}),
    ```
    where ``T`` is the absolute temperature. A measure of the HX performance is the effectiveness, ``\varepsilon``, defined as
    
    ```math
    \varepsilon = \frac{\dot{Q}}{\dot{Q}_{max}},
    ```
    where ``\dot{Q}`` is the actual heat transfer rate. The effectiveness can range between 0 and 1. A related quantity known as the number of transfer units (NTU) is defined as
    ```math
    \mathrm{NTU} = \frac{1}{C_{\mathrm{min}} R_o},
    ```
    where ``R_o`` is the overall thermal resistance. For any heat exchanger geometry, it can be shown that the effectiveness is a function of the NTU and the ratio ``\frac{C_{\mathrm{min}}}{C_{\mathrm{max}}}``. 
    
    In the case of a cross flow heat exchanger, this functional relationship between ``\varepsilon`` and NTU depends on whether there is internal mixing within the stream and on which stream has the minimum capacity rate. In a HX within a jet engine, it is reasonable to assume that the coolant stream has the minimum capacity rate as the fuel-air ratio is small. For tubular HX, the process stream is mixed but the coolant stream is unmixed as the tubes are independent. In this case, the relationship between ``\varepsilon`` and NTU is[^1]
    ```math
    \varepsilon = \frac{1}{C_{r}} [1 - \exp(-C_{r} (1 - \exp(-\mathrm{NTU})))],
    ```
    or conversely
    ```math
    \mathrm{NTU} = -\ln\left(1 + \frac{\ln(1 - C_r\varepsilon)}{C_r}\right).
    ```
    On the other hand, if the coolant stream has the maximum heat capacity rate, the effectiveness is given by
    ```math
    \varepsilon = 1 - \exp\left[-\frac{1}{C_{r}} (1 - \exp(-C_r \mathrm{NTU}))\right],
    ```
    and the corresponding NTU is
    ```math
    \mathrm{NTU} = -\frac{1}{C_{r}} \ln\left[ 1 + C_r \ln(1 - \varepsilon)\right].
    ```
    A notable property of these expressions is that there is a maximum effectiveness ``\varepsilon_\mathrm{max} < 1`` as the NTU tends to infinity. In the code, there is a check to see if the desired effectiveness exceeds the maximum possible one, in which case there is no solution. The effectiveness is limited to 95% of ``\varepsilon_\mathrm{max}`` to prevent a very large NTU.

    Once the effectiveness is known, the outlet specific enthalpies can be computed using
    ```math
    h_{o,p} = h_{i,p} - \frac{\dot{Q}}{\dot{m}_p}
    ```
    ```math
    h_{o,c} = h_{i,c} + \frac{\dot{Q}}{\dot{m}_c},
    ```
    where ``h`` represents the specific enthalpy, and the outlet temperatures can be determined from these.

    ### Recirculation
    If the coolant is originally in liquid form, it needs to undergo a phase change at some point in the heat exchanger. For cryogenic liquids, such as hydrogen or methane, it is unadvisable to expose air to these cryogenic temperatures as it can result in freezing or liquefaction of some species in air. A possible approach to overcome this is to introduce recirculation in the heat exchanger: this increases the coolant mass flow rate and allows for a higher coolant temperature while still having the same heat transfer. 

    The way recirculation is currently modeled in TASOPT is via a virtual "mixing chamber", where the hot recirculating mass flow that comes of of the HX is mixed with the colder coolant (which may be liquid in general) and heats it up to a desired HX inlet temperature,``T_{i,c}``. Neglecting the kinetic energy in the fluids, conservation of energy requires that
    ```math
    \dot{m}_{c,\infty} (h_{c,\infty}-h_{lat}) + \dot{m}_{r} h_{o,c} = (\dot{m}_{c,\infty}+\dot{m}_{r})h_{i,c},
    ```
    where ``\dot{m}_r`` is the recirculating mass flow rate and ``\dot{m}_{c,\infty}`` is the coolant mass flow rate before mixing with recirculation and, by mass conservation, the mass flow rate that leaves the system. The term ``h_{lat}`` represents a latent heat and may account for vaporization or, in the case of hydrogen, the ortho- to parahydrogen conversion. The specific enthalpies into and out of the HX are related by 
    ```math
    h_{o,c} = h_{i,c} + \frac{\dot{Q}}{\dot{m}_{c,\infty}+\dot{m}_{r}}.
    ```
    As the heat transfer rate is given by ``\dot{Q} = \varepsilon C_{\mathrm{min}} (T_{i,p} - T_{i,c})``, we can distinguish two cases depending on whether the coolant has the minimum or maximum heat capacity rate. If ``C_{\mathrm{min}}=(\dot{m}_{c,\infty}+\dot{m}_{r}) c_{pi,c}``,
    ```math
    \dot{m}_r = \dot{m}_{c,\infty} \frac{h_{i,c} - h_{c,\infty} + h_{lat}}{\varepsilon c_{pi,c}(T_{i,p} - T_{i,c})},
    ```
    and if ``C_{\mathrm{max}}=(\dot{m}_{c,\infty}+\dot{m}_{r}) c_{pi,c}``,
    ```math
    \dot{m}_r = \dot{m}_{c,\infty}\frac{A}{1-A},
    ```
    where 
    ```math
    A = \dot{m}_{c,\infty}\frac{h_{i,c} -h_{c,\infty} + h_{lat}}{\varepsilon C_h (T_{i,p} - T_{i,c})}.
    ```
    A failure case exists if ``A\geq 1``: in this case, the hot stream does not have enough heat capacity to provide the heat needed to get the coolant to the desired ``T_{i,c}``.

    ### Heat exchanger geometry
    ![HXfig](../assets/tubular_HX.svg)

    Two HX cross-sections are currently supported: rectangular and concentric cylinders (e.g., for a jet engine core). If the geometry is rectangular, the length of the tubes dictates the width of the rectangle, but if it is concentric, the tubes could have any length greater or equal to the distance between the cylinders as involute tubes could be used. The heat exchanger is assumed to have ``N_\mathrm{stages}`` different serial stages, with each stage containing ``N_\mathrm{passes}`` tube rows, each corresponding to a coolant pass. For example, the HX in the figure above with 2 stages and 3 coolant passes has a total of ``N_L=N_\mathrm{stages}N_\mathrm{passes}=6`` tube rows. The number of tubes at every row is ``N_t``; this parameter can be calculated from the mass flow rate through the cold side or from the geometry of the stage
    ```math
    N_t = \frac{4 \dot{m}_c}{\rho_{c,i}V_{c,i} \pi D_{t,i}^2 N_\mathrm{stages}} = \frac{b}{\frac{x_t}{D} D_{t,o}},
    ```
    where ``b`` is length across which the tubes are distributed. If the cross-section is concentric, ``b=\pi D_{c,i}`` with ``D_{c,i}`` being the inner cylinder diameter; if it is rectangular, ``b = \frac{A_{cs}}{l}``. In this expressions, ``A_{cs}=\frac{\dot{m}_p}{\rho_{p,i}V_{p,i}}`` is the freestream cross-sectional area. Since the tube inner diameter can be expressed as ``D_{t,i} = D_{t,o} - 2t``, this equation can be solved for the tube outer diameter
    ```math
    D_{t,o} = \frac{4  K  t + \sqrt{8  K t + 1} + 1}{2 K},
    ```
    with ``K = \frac{\pi b N_\mathrm{stages} \rho_{c,i}V_{c,i}}{4 \frac{x_t}{D} \dot{m}_c}``.
    
    The total length of the HX is simply ``L = N_L \frac{x_l}{D} D_{t,o}``.

    Some calculations rely on knowing the tangential pitch between tubes. This pitch may vary in the radial direction as the circumference changes but the tube diameter remains unchanged. In addition to this, the tubes are generally involute, which makes calculating this pitch even more challenging. In the code, a mean tangential pitch ``x_{t,m}`` is used to compute parameters such as the pressure drop and the Nusselt number. This mean pitch is calculated as
    ```math
    x_{t,m} = \frac{A_{cs}}{N_t l},
    ```
    where ``l`` is the length of each involute tube. The mass flow rate per unit area at the minimum free flow area is 
    ```math
    G = \frac{\dot{m}_p}{A_{cs} - N_t l D_{t,o}}.
    ```
    Note that for this expression to be valid, it is sufficient that $x_l/D\geq 1$. If the general geometry and total hot-side heat transfer area are known (e.g., from the NTU), but the number of coolant passes has not been determined yet, this can be calculated as
    ```math
    N_\mathrm{passes} = \frac{A_h}{N_t N_\mathrm{stages} \pi D_{t,o} l}.
    ```

    ### Heat transfer coefficients
    The above analysis relies on being able to determine the overall thermal resistance. In general, the thermal resistance has five components: hot- and cold-side gas resistances, wall resistance, and hot- and cold-side fouling resistances. The gas resistances are the aerodynamic resistances due to the thermal boundary layers, the wall resistance depends on the material conductivity and thickness, and the fouling resistances account for buildup of dirt layers during operation. The product of thermal resistance and heat transfer area (thermal insulance) is in practice easier to compute
    ```math
    R_o A_p = \frac{1}{h_p} + \frac{1}{h_c \frac{A_c}{A_p}} + \frac{t}{k_w} + R_{f,p}A_p + R_{f,c}A_c \frac{A_p}{A_c},
    ```
    where ``h`` is the the aerodynamic heat transfer coefficient, ``A`` is the heat transfer area, ``t`` is the wall thickness, ``k`` is the thermal conductivity, ``w`` denotes the wall, and ``R_fA`` is the fouling factor. A list of design fouling factors can be found in [^2] and [^3].
    
    The heat transfer coefficients depend on the gas temperature, which changes as heat is added to the flows. A mean gas temperature, taken as the average between the inlet and outlet temperatures, is used to calculate the gas properties,
    ```math
    T_{p,m} = \frac{T_{p,o} - T_{p,i}}{2}
    ```
    ```math
    T_{c,m} = \frac{T_{c,o} - T_{c,i}}{2}.
    ```
    
    #### Coolant-side heat transfer coefficient
    The flow inside the tubes can be modeled by assuming that it is fully-developed turbulent flow in a smooth pipe. In this case, the 1913 Blasius correlation provides a method to calculate the skin-friction coefficient, ``C_f``
    ```math
    C_f = \frac{\tau_w}{\frac{1}{2}\rho_{c,m} V_{c,m}^2} = 0.0791 \mathrm{Re}_{D,c}^{-1/4},
    ```
    where ``\tau_w`` is the wall shear stress, ``\rho`` is the mass density, ``V`` is the velocity and the cold-side diameter-based Reynolds number is defined as ``\mathrm{Re}_{D,c}=\frac{V_{c,m}\rho_{c,m} D_{t,i}}{\mu_{c,m}}``, with ``D_{t,i}`` being the tube inner diameter and ``\mu`` being the viscosity.
    
    When the skin-friction coefficient is known, the Colburn j-factor can be calculated using the Reynolds analogy
    ```math
    j = \mathrm{St} \mathrm{Pr}^{2/3} = \frac{C_f}{2},
    ```
    where ``\mathrm{St} = \frac{h}{\rho V c_p}`` is the Stanton number and ``\mathrm{Pr} = \frac{c_p \mu}{k}`` is the Prandtl number. Once ``j`` is determined, the heat transfer coefficient ``h_c`` can be computed from the cold gas properties.
    
    #### Process-side heat transfer coefficient
    The flow past a set of staggered tubes is complex. Žkauskas[^4] provides simplified correlations that can be used to model the heat transfer properties of these tubes. Generally, the Nusselt number can be expressed as
    ```math
    \mathrm{Nu} = C_1 C_2 \mathrm{Re}^m \mathrm{Pr}^n,
    ```
    where the Reynolds number is defined as ``\mathrm{Re}= \frac{G D_{t,o}}{\mu_{p,m}}``, ``D_{t,o}`` is the tube outer diameter, and ``G`` is the hot-side mass flow rate per unit area at the minimum free-flow area. Hence, this Reynolds number accounts for blockage effects due to the presence of the tubes.
    
    The following table shows the value of the parameters ``C_1``, ``m`` and ``n`` as a function of Reynolds number.
    
    | Re       | ``C_1`` | ``m``   | ``n`` |
    | -------- | ------- |-------- | ------|
    | 0–40     | 1.04    |0.4      |0.36   |
    | 40–1000  | 0.71     |0.5     |0.36   |
    | 1000–``2\times 10^5`` & ``x_t/x_l<2``| ``0.35  (x_t / x_l) ^ {0.2}``|0.6|0.36|
    | 1000–``2\times 10^5`` & ``x_t/x_l\geq 2``| 0.4    |0.6|0.36|
    | ``>2\times 10^5``| ``0.031  (x_t / x_l) ^ {0.2}``|0.8|0.4|
    
    The parameters in the table can be affected by the distances ``x_t`` and ``x_l``, which are the distances between tubes in the tangential and longitudinal directions. Note that the distance ``x_{t}`` used in the calculations is ``x_{t,m}`` as the pitch varies in the radial direction. The ratios of this distances to the tube outer diameter, ``\frac{x_t}{D}`` and ``\frac{x_t}{D}``, are design parameters. 
    
    The parameter ``C_2`` is a correction that accounts for the number of rows, ``N_L``, and tends to 1 as the number of rows goes to infinity. It can be approximated as 
    ```math
    C_2 = 1-\exp(-N_L^{1 / \sqrt{3}})
    ```
    if ``Re>1000`` and otherwise as
    ```math
    C_2 = 1-\exp(-\sqrt{3 N_L^{1 / \sqrt{2}}}).
    ```
    
    Once the Nusselt number is known, the hot-side heat transfer coefficient can be computed as ``h_p = \frac{\mathrm{Nu} k_p}{D_{t,o}}``.
    
    ### Pressure drops
    The pressure drop in the hot-side (across the staggered tube bank) can be estimated using the method of Gunter and Shaw[^5]. The first necessary parameter is the volumetric hydraulic diameter, defined as
    ```math
    D_v = \frac{4(\mathrm{Net\,free\,volume})}{\mathrm{Friction\,surface}} = \frac{4 L \frac{\dot{m}_p}{V_{p,i}\rho_{p,i}}-N_t N_\mathrm{passes} N_\mathrm{stages} \pi D_{t,o}^2 l}{A_p}.
    ```
    From this, the pressure drop across the process side can be computed as
    ```math
    \Delta p_p = \frac{G^2 L }{D_v \rho_{p,m}} \frac{f}{2} \left(\frac{D_v}{x_t}\right)^{0.4}\left(\frac{x_l}{x_t}\right)^{0.6} \left(\frac{\mu}{\mu_w}\right)^{-0.14},
    ```
    where ``\frac{f}{2}`` is a friction factor that can be related to the Reynolds number, ``Re_{D_v} = \frac{G D_v}{\mu_{p,m}}``, as ``\frac{f}{2}= 90 / Re_{D_v}`` for ``Re_{D_v}<200`` and ``\frac{f}{2}= 0.96  Re_{D_v}^{-0.145}`` otherwise. As in the heat transfer coefficient case, note that the distance ``x_{t}`` used in the calculations is ``x_{t,m}`` since the pitch varies in the radial direction.
    
    The cold-side pressure drop can be calculated from the skin-friction coefficient, ignoring the minor losses due to flow turning at each pass,
    ```math
    \Delta p_c = \frac{4 \tau_w N_\mathrm{passes} \pi D_{t,i} l}{\pi D_{t,i}^2}= \frac{4 \tau_w N_\mathrm{passes} l}{D_{t,i}},
    ```
    with ``\tau_w = C_f \frac{1}{2}\rho_{c,m} V_{c,m}^2``.

    ### Correction for temperature changes across thermal boundary layer
    Temperature changes across the thermal boundary layer in both the coolant and process streams. Identifying the temperature at which to evaluate the properties is not trivial. In the code, the method in Kays and London[^1] is used: all properties are evaluated at the freestream temperatures (``T_{p,m}`` and ``T_{c,m}``) and a correction is applied to the Nusselt number and the friction coefficient (``f`` in the process stream and ``C_f`` for the coolant). These take the form
    ```math
    \mathrm{Nu} = \mathrm{Nu}_m\left(\frac{T_w}{T_m}\right)^n
    ```
    ```math
    f = f_m\left(\frac{T_w}{T_m}\right)^m
    ```
    where the coefficients are ``n=0.0`` and ``m=0.0`` for the process-side properties and ``n=-0.5`` and ``m=-0.1`` for the coolant side[^1]. The wall temperature is calculated from the overall thermal resistance and the thermal resistance of the coolant side,
    ```math
    T_w = T_{c,m} + \frac{T_{p,m}-T_{c,m}}{R_o A_p}\left(\frac{1}{h_c \frac{A_c}{A_p}} + R_{f,c}A_c \frac{A_p}{A_c}\right)
    ```

## Structures
```@docs
engine.HX_gas
```
```@docs
engine.HX_tubular
```
```@docs
engine.HeatExchanger
```
## Functions
### Heat exchanger sizing and off-design operations
```@docs
engine.hxsize!
```
```@docs
engine.hxoper!
```
### Optimization
```@docs
engine.hxoptim!
```
```@docs
engine.hxobjf
```

### Overall design and analysis
```@docs
engine.hxdesign!
```
### Heating and pressure calculations
```@docs
engine.jcalc_pipe
```
```@docs
engine.Nu_calc_staggered_cyl
```
```@docs
engine.Δp_calc_staggered_cyl
```
### Weight estimation
```@docs
engine.hxweight
```

[^1]: Kays, W. M., & London, A. L. (1984). Compact heat exchangers.
[^2]: Standards of the Tubular Exchanger Manufacturers Association
[^3]: [Powder Process](https://powderprocess.net/Tools_html/Data_Diagrams/Heat_Exchanger_Fouling_Factor.html)
[^4]: Žkauskas, A. (1987). Heat transfer from tubes in crossflow. In Advances in heat transfer (Vol. 18, pp. 87-159). Elsevier.
[^5]: Gunter, A. Y., & Shaw, W. A. (1945). A general correlation of friction factors for various types of surfaces in crossflow. Transactions of the American Society of Mechanical Engineers, 67(8), 643-656.