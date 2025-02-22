# Aircraft Stability

This module provides key functions for aircraft **stability and center of gravity (CG) analysis**, ensuring proper trim and balance throughout flight. The `balance` routine calculates the aircraft's **center of gravity (`xCG`)**, **center of pressure (`xCP`)**, and **neutral point (`xNP`)**, adjusting horizontal tail lift, area, or wing box location to maintain stability. The `htsize` routine determines the **horizontal tail area (`Sh`)** and **wing box position (`xwbox`)**, solving for pitch trim and static stability across different flight conditions. The `cglpay` routine computes **forward (`xcgF`) and aft (`xcgB`) CG limits** based on payload and fuel distribution. Lastly, the `cabin_centroid` routine calculates the **cabin centroid (`xcabin`) and length (`lcabin`)**, accounting for fuel tank placement and fuselage geometry. These functions collectively ensure proper **aircraft weight distribution, stability, and trim performance** across various loading scenarios.


```@docs
TASOPT.balance
TASOPT.htsize
TASOPT.cglpay
TASOPT.cabin_centroid
```