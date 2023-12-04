!!! details "📖 Theory - Chord distribution"
    The wing or tail surface is assumed to have a two-piece linear planform
    with constant sweep $\Lambda$, shown in
    the figure below. The
    inner and outer surface planforms are defined in terms of the center
    chord $c_o$ and the inner and outer taper ratios. 

    $$\begin{aligned}
    \lambda_s & = & c_s/c_o
    \\
    \lambda_t & = & c_t/c_o
    \end{aligned}$$ 

    Similarly, the spanwise dimensions are defined in terms
    of the span $b$ and the normalized spanwise coordinate $\eta$.

    $$\begin{aligned}
    \eta   & = & 2y / b  
    \\
    \eta_o & = & b_o / b
    \\
    \eta_s & = & b_s/b 
    \end{aligned}$$ 

    For generality, the wing center box width $b_o$ is
    assumed to be different from the fuselage width to allow possibly
    strongly non-circular fuselage cross-sections. It will also be different
    for the tail surfaces. A planform break inner span $b_s$ is defined,
    where possibly also a strut or engine is attached. Setting
    $b_s \!=\! b_o$ and $c_s \!=\! c_o$ will recover a single-taper surface.

    It's convenient to define the piecewise-linear normalized chord function

    $$\begin{aligned}
    \frac{c {\scriptstyle (\eta)}}{c_o} \; \equiv \; 
    C {\scriptstyle (\eta \, ; \, \eta_o,\eta_s, \lambda_s,\lambda_t)}
      & = & 
    \left\{
    \begin{array}{lcl}
    \; 1 & , & 0 < \eta < \eta_o 
    \\[0.5em]
    \displaystyle
    \; 1 \, + (\lambda_s\!-\!1 \,) \frac{\eta   \!-\! \eta_o}{\eta_s\!-\!\eta_o}
    & , & \eta_o < \eta < \eta_s 
    \\[0.25em]
    \displaystyle
    \lambda_s + (\lambda_t\!-\!\lambda_s) \frac{\eta   \!-\! \eta_s}{1\!-\!\eta_s}
    & , & \eta_s < \eta < 1
    \end{array}
    \right.
    %\label{ceta}
    \end{aligned}$$ 

    The following integrals will be useful for area, volume,
    shear, and moment calculations. 

    $$\begin{aligned}
    \int_0^{\eta_o} C \:\: {\rm d}\eta 
    & =  \eta_o
    \\
    \int_{\eta_o}^{\eta_s} C \:\: {\rm d}\eta
    & = \frac{1}{2} ( 1\!+\!\lambda_s)(\eta_s\!-\!\eta_o)
    \\
    \int_{\eta_s}^1 C \:\: {\rm d}\eta 
    & =  \frac{1}{2} (\lambda_s\!+\!\lambda_t)(1\!-\!\eta_s)
    \\[0.25em]
    \int_0^{\eta_o} C^2\:\: {\rm d}\eta 
    & =  \eta_o
    \\
    \int_{\eta_o}^{\eta_s} C^2\:\: {\rm d}\eta
    & =  \frac{1}{3} ( 1\!+\!\lambda_s\!+\!\lambda_s^2)(\eta_s\!-\!\eta_o)
    \\
    \int_{\eta_s}^1 C^2\:\: {\rm d}\eta 
    & =  \frac{1}{3} (\lambda_s^2\!+\!\lambda_s \lambda_t\!+\!\lambda_t^2)(1\!-\!\eta_s)
    \\
    %
    \int_{\eta_o}^{\eta_s} C \: (\eta\!-\!\eta_o) \:\: {\rm d}\eta
    & =  \frac{1}{6} ( 1\!+\!2\lambda_s)(\eta_s\!-\!\eta_o)^2
    \\
    \int_{\eta_s}^1 C \: (\eta\!-\!\eta_s) \:\: {\rm d}\eta
    & =  \frac{1}{6} (\lambda_s\!+\!2\lambda_t)(1\!-\!\eta_s)^2
    \\
    \int_{\eta_o}^{\eta_s} C^2 \: (\eta\!-\!\eta_o) \:\: {\rm d}\eta
    & =  \frac{1}{12} ( 1\!+\!2\lambda_s\!+\!3\lambda_s^2)(\eta_s\!-\!\eta_o)^2
    \\
    \int_{\eta_s}^1 C^2 \: (\eta\!-\!\eta_s) \:\: {\rm d}\eta
    & =  \frac{1}{12} (\lambda_s^2\!+\!2\lambda_s\lambda_t\!+\!3\lambda_t^2)(1\!-\!\eta_s)^2
    %\label{Cint2}
    \end{aligned}$$

!!! details "📖 Theory - Surface area and aspect ratio"

    The surface area $S$ is defined as the exposed surface area plus the
    fuselage carryover area. 

    $$\begin{aligned}
    S & =  2 \int_0^{b/2} \! c \:\: {\rm d}y \;=\; c_o \, b \, K_c
    %\label{Sdef}
    \\
    \mathrm{where} \hspace{3em}
    K_c & =   \int_0^1  C \:\: {\rm d}\eta \;=\; 
    \eta_o + {\textstyle \frac{1}{2}}( 1       \!+\!\lambda_s)(\eta_s\!-\!\eta_o)
          + {\textstyle \frac{1}{2}}(\lambda_s\!+\!\lambda_t  )(1     \!-\!\eta_s)
    \end{aligned}$$ 

    The aspect ratio is then defined in the usual way. This
    will also allow relating the root chord to the span and the taper
    ratios. 

    $$\begin{aligned}
    {A\hspace{-0.5ex}R}& = & \frac{b^2}{S}
    %\label{ARdef}
    \end{aligned}$$ 

    It is also useful to define the wing's mean aerodynamic
    chord $c_{\rm ma}$ and area-centroid offset
    ${\scriptstyle \Delta}x_{\rm wing}$ from the center axis.

    $$\begin{aligned}
    \frac{c_{\rm ma}}{c_o} &\!=
    \frac{1}{c_o} \frac{2}{S} \int_0^{b/2} \! c^2 \; {\rm d}y 
    \;=\; \frac{K_{cc}}{K_c}
    \\
    {\scriptstyle \Delta}x_{\rm wing}
    &\!=
    \frac{2}{S} \int_{b_o/2}^{b/2} \! c \, (y \!-\! y_o) \tan\Lambda \; {\rm d}y 
    \;=\; \frac{K_{cx}}{K_c} \, b \: \tan\Lambda
    \\
    x_{\rm wing}& =  x_{\rm wbox}\,+\, {\scriptstyle \Delta}x_{\rm wing}
    \\[0.25em]
    \mathrm{where} \hspace{1em}
    K_{cc} &\!= \int_0^1  C^2 \:\: {\rm d}\eta
    \nonumber \\
    &\!= \eta_o
    + \frac{1}{3} ( 1         \!+\!\lambda_s\!+\!\lambda_s^2)(\eta_s\!-\!\eta_o)
    + \frac{1}{3} (\lambda_s^2\!+\!\lambda_s \lambda_t\!+\!\lambda_t^2)(1\!-\!\eta_s)
    \hspace{3em}
    \\[0.25em]
    K_{cx} &\!=  \int_{\eta_o}^1  C \: (\eta\!-\!\eta_o) \:\: {\rm d}\eta 
    \nonumber \\
    &\!=
    \frac{1}{12} ( 1       \!+\!2\lambda_s)(\eta_s\!-\!\eta_o)^2 +
    \frac{1}{12} (\lambda_s\!+\!2\lambda_t)(1\!-\!\eta_s)^2      +
    \frac{1}{4}  (\lambda_s\!+\!\lambda_t)(1\!-\!\eta_s)(\eta_s\!-\!\eta_o)
    \hspace{3em}
    \end{aligned}$$

    The wing area centroid is used in the fuselage bending
    load calculations as described earlier.

## Reference quantities

The aircraft reference quantities are chosen to be simply the values for
the wing. 

$$\begin{aligned}
b_{\rm ref}   &\!= (b)_{\rm wing}\\
S_{\rm ref}   &\!= (S)_{\rm wing}\\
{A\hspace{-0.5ex}R}_{\rm ref} &\!= ({A\hspace{-0.5ex}R})_{\rm wing}\\
c_{\rm ref}   &\!= (c_{\rm ma})_{\rm wing}
\end{aligned}$$

