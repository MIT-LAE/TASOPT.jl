# Aircraft Stability

This module provides functions to size the tail surfaces and balance the aircraft to achieve pitch trim throughout flight when called by [`_size_aircraft!()`](@ref TASOPT._size_aircraft!).

- [`balance_aircraft!()`](@ref TASOPT.balance_aircraft!) makes adjustments as described below to achieve pitch trim. It calculates the aircraft's resulting **center of gravity (`xCG`)**, **center of pressure (`xCP`)**, and **neutral point (`xNP`)** at a specific flight point. To meet the pitch trim requirement, the routine adjusts one of (i) the horizontal tail's lift coefficient, (ii) its area, or (iii) the axial location of the wing box. In *almost all* the use cases, only option (i) adjusting the horizontal tail's lift coefficient is relevant for trim calculations.

!!! warning "Some `balance_aircraft!()` options change the aircraft design"
    Exercise caution if options (ii) or (iii) are being used as these options *alter the design of the aircraft* (i.e., they should not be used in off-design cases).

- [`size_htail()`](@ref TASOPT.size_htail) performs a more involved stability analysis: it determines the **horizontal tail area (`Sh`)** and **wing box position (`xwbox`)**, solving for pitch trim and static stability across all flight conditions.

Both stability analyses call the following helper functions.

- [`CG_limits()`](@ref TASOPT.CG_limits) computes **forward (`xcgF`) and aft (`xcgB`) CG limits** based on payload and fuel distribution. 
- Lastly, [`cabin_centroid()`](@ref TASOPT.cabin_centroid) calculates the **cabin centroid (`xcabin`) and length (`lcabin`)**, accounting for fuel tank placement and fuselage geometry.

!!! details "📖 Theory - Pitch trim and stability requirements"
    Every operating point must meet the requirement of pitch trim, which is equivalent to the centers of weight and pressure coinciding. This is enforced by requiring that the following total-moment residual is zero.

    ```math
    \mathcal{R}_M\left( x_{\text{wbox}}, S_h, C_{Lh}, C_L, r_{\text{fuel}}, r_{\text{pay}}, \xi_{\text{pay}} \right) \equiv x_{\text{CG}} - x_{\text{CP}} = \frac{xW}{W} + \frac{c_o \, C_m}{C_L} = 0
    ```
    
    An aircraft must also have some minimum amount of static pitch stability, which means that the rearmost center of gravity must be ahead of the neutral point by the static margin fraction $f_{SM}$ of the mean aerodynamic chord. This is met when the following stability residual is zero:

    ```math
    \mathcal{R}_{S_h}\left( x_{\text{wbox}}, S_h, r_{\text{fuel}}, r_{\text{pay}}, \xi_{\text{pay}} \right) \equiv x_{\text{CG}} - x_{\text{NP}} + f_{\text{SM}} \, c_{\text{MA}} = 0
    ```

```@docs
TASOPT.balance_aircraft!
TASOPT.CG_limits
TASOPT.size_htail
TASOPT.cabin_centroid
```
