!!! details "📖 Theory - Airfoil database lookup"
    
    ## The drag database
    The airfoil drag database (e.g., `C_airfoil.csv`) stores the 2D airfoil section performance characteristics as a function of:
    - airfoil thickness/chord ratio, $\frac{t}{c} =$ `toc`,
    - section lift coefficient, $c_l$, and
    - Mach number, $M$.

    At each point, the raw database includes the section's angle-of-attack ($AoA$), drag ($c_d$), pitching moment ($c_m$), viscous, wave, pressure, and skin friction drag components ($c_d = c_{d,v} + c_{d,w} = c_{d,p} + c_{d,f}$), and data source information (source file, status, interpolation status and axis). Since the processing requires a rectilinear grid of points, the status column ("Missing") indicates where data points are absent, usually due to performance infeasibility.

    The [`airtable`](@ref TASOPT.aerodynamics.airtable) function builds an `airfoil` structure that stores the database $C_D, C_{D,p}, C_{D,v}, C_{D,w}, C_M, AoA = f(\frac{t}{c}, c_l, M)$ as well as the derivatives at the *knots* (as described below) to use for tricubic interpolation.

    ## Tricubic interpolation
    Below is a terse description of the cubic interpolation as relevant to the drag interpolation. There are several online resources available for a more in-depth description of the topic.
    
    The cubic polynomials are represented as:
    ```math
    \newcommand{\im}{{i-1}}
    \newcommand{\ip}{{i+1}}
    \newcommand{\p}[1]{{#1 + 1}}
    \newcommand{\m}[1]{{#1 - 1}}
    q_i(t) = (1-t)y_\m{i} + ty_i + t(1-t)\left((1-t)a_i - tb_i\right),\\
    ```
    where,
    ```math
    \newcommand{\im}{{i-1}}
    \newcommand{\ip}{{i+1}}
    \newcommand{\p}[1]{{#1 + 1}}
    \newcommand{\m}[1]{{#1 - 1}}
    \begin{aligned}
    
    t &= \frac{x- x_\im}{x_i - x_\im}\\
    a_i &= k_\m{i}(x_i - x_\m{i}) - (y_i - y_\m{i}),\\
    b_i &= k_i(x_i - x_\m{i})  - (y_i - y_\m{i})

    \end{aligned}
    ```
    and 
    ```math
    k_i = \left. \frac{\mathrm{d} y}{\mathrm{d} x} \right|_i
    ```
    is the slope at the $i^{th}$ point in the interpolation array.

    These slopes (in each direction) are calculated and stored in the `airfoil.A` structure by the [`airtable`](@ref TASOPT.aerodynamics.airtable) function.

    For brevity only an example of how these are stored is depicted below, the interested user may directly read `airtable` to see the matrix structure where the knots are stored. Here is an example of taking a slice of the 3-D array storing $c_{d_p}$ at a fixed $c_l = c_{l_i}$ (the data is really a 4-D array with each of the performance characteristics being stored along the 4$^\mathrm{th}$ dimension),

    ```math
    \newcommand{\pder}[3]{\frac{\partial^{#1} #2}{\partial #3^{#1}} }
    \begin{aligned}
    \mathbf{A}[c_{l_i}, :, :] &= 
    \begin{bmatrix}
    	\left. c_{d_p}\right._{\frac{t}{c}_1, M_1}& \left. c_{d_p}\right._    {\frac{t}{c}_1, M_2} & \cdots & \left. c_{d_p}\right._{\frac{t}{c}_1,     M_N} \\
    	\left. c_{d_p}\right._{\frac{t}{c}_1, M_2}& \cdots& \cdots &     \left. c_{d_p}\right._{\frac{t}{c}_2, M_N}\\
    	\vdots & \ddots & & \vdots \\
    	\left. c_{d_p}\right._{\frac{t}{c}_N, M_1}& \cdots& \cdots &     \left. c_{d_p}\right._{\frac{t}{c}_N, M_N}\\
    	
    & \\
    \end{bmatrix}\\ 
    \\
    
    \pder{}{\mathbf{A}[c_{l_i}, :, :] }{M} & =
    \begin{bmatrix}
    	\left. \pder{}{c_{d_p}}{M} \right|_{1,1}& \cdots & \left.     \pder{}{c_{d_p}}{M} \right|_{1,N} \\
    	\vdots & \ddots &  \vdots \\
    	\left. \pder{}{c_{d_p}}{M} \right|_{N,1}& \cdots & \left.     \pder{}{c_{d_p}}{M} \right|_{N,N} \\
    & \\
    \end{bmatrix}\\	
    
    \end{aligned}
    ```
    Similarly, the slopes in the $\frac{t}{c}$ direction are stored in `A_toc` and the cross-derivatives are stored in `A_M_toc` and so on. 

    ## Evaluating the drag calculations

    The stored slopes are used to evaluate the drag components by first locating the interval to use based on the $c_l, \frac{t}{c}, M$ values in the call to [`airfun`](@ref TASOPT.aerodynamics.airfun), and then performing the interpolation in each dimension sequentially. For any values outside of the non-missing data in the database, a quadratic drag penalty is added.  


