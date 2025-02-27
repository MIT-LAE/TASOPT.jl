# [Landing gear](@id landing_gear)

The landing gear in TASOPT.jl can be modeled in two ways depending on the user preference in the input file. The simplest model assumes that the nose and main landing gears are simply constant fractions of the maximum takeoff weight. A more elaborate model calculates the landing gear length and uses this to find its mass via a correlation to historical data from Raymer[^1]. The landing gear is sized in [`size_landing_gear!()`](@ref structures.size_landing_gear!).

!!! details "📖 Theory - Landing gear sizing via historical correlations" 
    The landing gear length can be sized by (a) the need to avoid a tailstrike at rotation during takeoff, or (b) the need to provide a minimum engine ground clearance, as illustrated in the figure below. In TASOPT.jl, the nose gear length is calculated first and the main gear length is subsequently calculated based on the wing geometry. 

    ![LGlength](../assets/landing_gear_diags.svg)

    In the tailstrike case, the length of the landing gear (``l_{ts}``) can be determined from the geometry as
    ```math
        l_{ts}= (x_{end} - x_{mg}) \tan(\theta_{ts}) - D_f,
    ```
    where ``x_{end}`` is the distance from the front of the aircraft to the end of the fuselage, ``x_{mg}`` is the main gear location, ``\theta_{ts}`` is the desired tailstrike angle, and ``D_f`` is the fuselage diameter.

    In the engine ground clearance case, the length of the landing gear (``l_c``) is given by
    ```math
        l_{c}= d_c + D_{fan} - y_{eng}\tan(\Gamma),
    ```
    where ``d_c`` is the desired ground clearance, ``D_{fan}`` is the fan diameter, ``y_{eng}`` is the distance from the fuselage centerline to the engine centerline, and ``\Gamma`` is the wing dihedral angle.

    The selected nose gear length is ``l_{ng}=\max(l_{ts},l_c)``. The main gear length is
    ```math
        l_{mg} = l_{ng}+y_{mg}\tan(\Gamma),
    ```
    where ``y_{mg}`` is the spanwise location of the main gear.
    
    Once this has been determined, the historical-data correlations in Raymer[^1] can be used to calculate the mass of the nose and main landing gears. For the main gear, the mass (in lbs) is given by
    ```math
        m_{mg} = 0.0106\, \mathrm{MTOW}^{0.888}\, N_l^{0.25}\, l_{mg}^{0.4}\, N_{w,m}^{0.321}\, N_{s,m}^{-0.5}\, V_{stall}^{0.1},
    ```
    where ``\mathrm{MTOW}`` is the maximum takeoff weight in lbs, ``N_l=4.5`` is the ultimate load factor that the landing gear is sized for, ``l_mg`` is the main gear length in inches, ``N_{w,m}`` is the total number of wheels in the main gear, ``N_{s,m}`` is the number of main gear shock struts, and ``V_{stall}`` is the stall speed in knots.

    Similarly, the mass of the nose gear (in lbs) is determined by[^1]
    ```math
        m_{ng} = 0.032\,\mathrm{MTOW}^{0.646}\, N_l^{0.2}\, l_{ng}^{0.5}\, N_{w,n}^{0.45},
    ```
    where ``l_{ng}`` is the length of the nose gear in inches and ``N_{w,n}`` is the number of wheels in the nose gear.

```@docs
TASOPT.structures.size_landing_gear!
```

[^1]: Raymer, Daniel. Aircraft design: a conceptual approach. American Institute of Aeronautics and Astronautics, Inc., 2012.