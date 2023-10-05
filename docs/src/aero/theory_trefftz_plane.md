!!! details "📖 Theory - induced drag $C_{d,i}$" 

    The induced drag is calculated using a discrete vortex Trefftz-Plane
    analysis. The circulation of the wing wake immediately behind the
    trailing edge is 
    
    $$\begin{aligned}
    \Gamma_{\!{\rm wing}} {\scriptstyle (\eta)}
    & = \frac{\tilde{p} {\scriptstyle (\eta)}}{\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}}
    \; \simeq \; \frac{p {\scriptstyle (\eta)}}{\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}} \, \sqrt{1 \!-\! \eta^{k_t}}
    \\
    k_t & \simeq 16
    \end{aligned}$$ 
    
    where the approximation realistically represents the tip
    lift rolloff for typical taper ratios, and is consistent with the
    assumed $f_{L_{\scriptstyle t}}\simeq -0.05$ value for the tip lift loss
    factor. This circulation is convected into the wake along streamlines
    which will typically constrict behind the fuselage by continuity.
    The Figure above shows two possible aft fuselage taper shapes, giving two different wake
    constrictions.
    
    An annular streamtube at the wing contracts to another annular
    streamtube in the wake with the same cross-sectional area. The $y$ and
    $y'$ locations on the wing and wake which are connected by a streamline
    are therefore related by the correspondence function. 
    
    $$\begin{aligned}
    y'{\scriptstyle (y)}& = & \sqrt{ y^2 - y_o^2 + {y'_o}^2}
    \end{aligned}$$ 
    
    The Trefftz Plane circulation
    $\Gamma {\scriptstyle (y')}$ is then given by the coordinate shift. The
    mapping function $y'{\scriptstyle (y)}$ is not defined for
    $y \! < \! y_o$, so the circulation there is simply set from the $y_o$
    value. 
    
    $$\begin{aligned}
    \Gamma_{\rm wake} {\scriptstyle (y')}& = & 
    \left\{
    \begin{array}{lcl}
    \Gamma_{\!{\rm wing}} \left( y {\scriptstyle (y')}\right)
    & , & y \!>\! y_o \\
    \Gamma_{\!{\rm wing}} \left( y_o \right)
    \end{array}
    \right.
    \end{aligned}$$
    
    The Trefftz Plane analysis uses point vortices. The circulation
    (\[Gamwake\]) is
    evaluated at the midpoints of $n$ intervals along the wake trace, spaced
    more or less evenly in the Glauert angle to give a cosine distribution
    in physical space. The wake's vertical $z$ positions are simply taken
    directly from the wing. 

    $$\begin{aligned}
    \theta_{i+1/2} & = & \frac{\pi}{2} \, \frac{i-1/2}{n}
    \hspace{2ex} , \hspace{2ex} i = 1 \ldots n
    \\
    y_{i+1/2} & = & \frac{b}{2} \, \cos \theta_{i+1/2}
    \\
    y'_{i+1/2} & = & \sqrt{y_{i+1/2}^2 - y_o^2 + {y'_o}^2}
    \\
    z'_{i+1/2} & = & z_{i+1/2}
    \\
    \Gamma_{\!i+1/2} & = & \Gamma_{\!{\rm wing}} (y_{i+1/2})
    \end{aligned}$$ 
    
    The locations of $n+1$ trailing vortices are computed
    similarly. 
    
    $$\begin{aligned}
    \theta_i & = & \frac{\pi}{2} \, \frac{i-1}{n}
    \hspace{2ex} , \hspace{2ex} i = 1 \ldots n\!+\!1
    \\
    y_i & = & \frac{b}{2} \, \cos \theta_i
    \\
    y'_i & = & \sqrt{y_i^2 - y_o^2 + {y'_o}^2}
    \\
    z'_i & = & z_i
    \end{aligned}$$ 

    The circulations of these trailing vortices are the
    differences of the adjacent bound circulations, with the circulation
    beyond the tips effectively zero. 
    
    $$\begin{aligned}
    \bar{\Gamma}_{\!i} & = & 
    \left\{
    \begin{array}{lcl}
    \hspace{6.5ex}
    -\Gamma_{\!i-1/2} & , & i = 1 \hspace{3.5em} \mathrm{(left tip)}
    \\
    \Gamma_{\!i+1/2}-\Gamma_{\!i-1/2} & , & i = 2 \ldots n
    \\
    \Gamma_{\!i+1/2} & , & i = n\!+\!1 \hspace{2em} \mathrm{(right tip)}
    \end{array}
    \right.
    \end{aligned}$$ 

    The above definitions are also applied to the horizontal
    tail, with its discrete points simply appended to the list and $n$
    increased accordingly.
    
    The Trefftz plane calculation proceeds by first calculating the $y$-$z$
    wake velocity components at the $y'_{i+1/2},z'_{i+1/2}$ interval
    midpoints, induced by all the trailing vortices and their left-side
    images. 

    $$\begin{aligned}
    v_{i+1/2} & = & \sum_{j=1}^{n+1} \frac{\bar{\Gamma}_{\!j}}{2\pi} 
    \left[
    \frac{-(z'_{i+1/2}\!-\!z'_j)}
          {(y'_{i+1/2}\!-\!y'_j)^2
         + (z'_{i+1/2}\!-\!z'_j)^2}
    \,-\,
    \frac{-(z'_{i+1/2}\!-\!z'_j)}
          {(y'_{i+1/2}\!+\!y'_j)^2
         + (z'_{i+1/2}\!-\!z'_j)^2}
    \right]
    \hspace{2em}
    \\
    w_{i+1/2} & = & \sum_{j=1}^{n+1} \frac{\bar{\Gamma}_{\!j}}{2\pi}
    \left[
    \frac{ y'_{i+1/2}\!-\!y'_j}
         {(y'_{i+1/2}\!-\!y'_j)^2
        + (z'_{i+1/2}\!-\!z'_j)^2}
    \,-\,
    \frac{ y'_{i+1/2}\!+\!y'_j}
         {(y'_{i+1/2}\!+\!y'_j)^2
        + (z'_{i+1/2}\!-\!z'_j)^2}
    \right]
    \end{aligned}$$ 

    The overall lift and induced drag are then computed
    using the Trefftz Plane vertical impulse and kinetic energy. The sums
    are doubled to account for the left side image. 

    $$\begin{aligned}
    C_{\!L_{\scriptscriptstyle T\!P}} & = & \frac{2}{{\textstyle \frac{1}{2}}\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}^2 S} \: \sum_{i=1}^{n} 
    \: \rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}\, \Gamma_{\!i+1/2} \: (y'_{i+1}-y'_i)
    \\
    C_{\!D_{\scriptscriptstyle T\!P}} & = & \frac{2}{{\textstyle \frac{1}{2}}\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}^2 S} \: \sum_{i=1}^{n} 
    -\frac{\rho_{\scriptscriptstyle \infty}}{2} \, \Gamma_{\!i+1/2} \left[
          w_{i+1/2} \,(y'_{i+1}\!-\!y'_i)
     \:-\:v_{i+1/2} \,(z'_{i+1}\!-\!z'_i) \right]
    \end{aligned}$$
    
    To minimize any modeling and numerical errors incurred in the wake
    contraction model and the point-vortex summations, the final induced
    drag value is scaled by the square of the surface-integral and
    Trefftz-Plane drag values. 

    $$\begin{aligned}
    C_{\!D_i} & = & C_{\!D_{\scriptscriptstyle T\!P}} \left(\frac{C_{\!L}}{C_{\!L_{\scriptscriptstyle T\!P}}} \right)^{\!2}
    \end{aligned}$$ 
    
    This is equivalent to using the Trefftz Plane analysis
    to calculate the span efficiency rather than the actual induced drag
    coefficient.