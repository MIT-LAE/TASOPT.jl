!!! details "Theory" 
    The viscous calculation
    produces displacement, momentum, and kinetic energy areas
    $\Delta^*, \Theta, \Theta^* {\scriptstyle (x)}$.

    The cross-sectional area over the center cylindrical portion is
    $A_{\rm fuse}$, which has already been defined by

    $$\begin{aligned}
    A_{\rm fuse} =  \left[ \pi + n_{\rm fweb}\left( 2\theta_{\rm fb}+ \sin 2 \theta_{\rm fb}\right) 
     \right] R_{\rm fuse}^2 \;+\; 2 \left[ R_{\rm fuse}+ n_{\rm fweb}w_{\rm fb}\right] \Delta R_{\rm fuse}.
    \end{aligned}$$ 

    This also defines the radius of the equivalent round cylinder. 
    $$\begin{aligned}
    R_{\rm cyl}& = & \sqrt{\frac{A_{\rm fuse}}{\pi}}
    \end{aligned}.$$ 
    
    The equivalent radii over the tapering nose and radius
    are then defined via the following convenient functions.
    
    $$\begin{aligned}
    R {\scriptstyle (x)}& = &
    \left\{
    \begin{array}{lcl}
    \displaystyle
    R_{\rm cyl}
    \left[ \: 
    1 - \left( \frac{x_{{\rm blend}_1} \!-\! x}{x_{{\rm blend}_1} \!-\! x_{\rm nose}} 
          \right)^{\!\! a} \;
          \right]^{\! 1/a}
    & , & x_{\rm nose}< x < x_{{\rm blend}_1} 
    \\[1.0em]
    \displaystyle
    R_{\rm cyl}
    & , & x_{{\rm blend}_1} < x < x_{{\rm blend}_{\,2}}
    \\[0.5em]
    \displaystyle
    R_{\rm cyl}
    \left[ \:
    1 - \left( \frac{x \!-\! x_{{\rm blend}_{\,2}}}{x_{\rm end} \!-\! x_{{\rm blend}_{\,    2}}} 
          \right)^{\!\! b} \;
          \right]
    & , & x_{{\rm blend}_{\,2}} < x < x_{\rm tail}
    \end{array}
    \right.
    \\
    a & \simeq 1.6 &
    \\
    b & \simeq 2.0 &
    \end{aligned}$$

    The $x_{{\rm blend}_1}$ and $x_{{\rm blend}_{\,2}}$
    locations are the nose and tailcone blend points, and do not necessarily
    have to be exactly the same as the $x_{{\rm shell}_1}$ and
    $x_{{\rm shell}_{\,2}}$ locations which define the loaded pressure
    shell. Likewise, $x_{\rm end}$ is the aerodynamic endpoint of the
    tailcone, and is distinct from its structural endpoint $x_{\rm conend}$.
    The $a$ and $b$ constant values above give reasonable typical fuselage
    shapes.
    
    If the fuselage is nearly round, the necessary area and perimeter
    distributions follow immediately. 
    $$\begin{aligned}
    A {\scriptstyle (x)}& = & \pi \, {R{\scriptstyle (x)}}^2
    \\
    b_0 {\scriptstyle (x)}& = & 2 \pi R {\scriptstyle (x)}
    \end{aligned}$$ This would be suitably modified for non-circular
    cross-sections.
    
    With this geometry definition, the viscous/inviscid calculation
    procedure provides the momentum and kinetic energy area distributions
    along the body and wake, 
    $$\begin{aligned}
    \left\{ \Theta {\scriptstyle (s)}\, , \: \Theta^* {\scriptstyle (s)}\right\} & = &
    f_{\rm f_{excr}} \:
    {\cal F}(M_{{\scriptscriptstyle \infty}}, Re_\ell\, ; \, A {\scriptstyle (x)}, b_0 {\scriptstyle (x)})
    \end{aligned}$$ where ${\cal F}$ denotes the overall viscous/inviscid
    calculation procedure, and $f_{\rm f_{excr}} \geq 1$ is an empirical
    factor to allow for fuselage excrescence drag sources.
    
    Specific values of interest are the momentum area $\Theta_{\rm wake}$ at
    the wake numerical endpoint $s_{\rm wake}$, the far-downstream momentum
    area $\Theta_{\scriptscriptstyle \infty}$, and the kinetic energy area
    $\Theta_{\scriptscriptstyle T\!E}$ at the body endpoint or trailing
    edge. 

    $$\begin{aligned}
    \Theta_{\rm wake}& = & \Theta (s_{\rm wake}) 
    \\
    H_{\rm avg}  & = & 
    {\textstyle \frac{1}{2}}\left[ H(s_{\rm wake}) + 1 + (\gamma\!-\!1) M_{{\scriptscriptstyle \infty}}^2 \right]
    \\
    \Theta_{\scriptscriptstyle \infty}& = & \Theta_{\rm wake}
    \left( \frac{u_e (s_{\rm wake})}{V_{\!{\scriptscriptstyle \infty}}} \right)^{H_{\rm avg}} 
    \\
    \Theta^*_{\scriptscriptstyle T\!E}& = & \Theta^*(s_{\scriptscriptstyle T\!E}) 
    \end{aligned}$$ 
    
    The equation above is the Squire-Young formula, with $H_{\rm avg}$
    being the average shape parameter between the end of the wake and far
    downstream.
    
    The fuselage surface + wake dissipated power in the absence of BLI is
    then evaluated as follows, consistent with the usual wake momentum
    defect relations. 
    $$\begin{aligned}
    C'_{\!D_{\rm fuse}} & \equiv &
    \frac{\Phi_{\rm surf}-P_{V_{\rm surf}} + \Phi_{\rm wake}-P_{V_{\rm wake}}}
               {{\textstyle \frac{1}{2}}\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}^3 S} 
    \hspace{6ex}
    \rm{(without BLI)}
    \hspace{-9.0ex}
    \\
    C'_{\!D_{\rm fuse}}
    & = & \frac{D_{\rm fuse}}{{\textstyle \frac{1}{2}}\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}^2 S}
    \;=\; \frac{2 \Theta_{\scriptscriptstyle \infty}}{S} 
    \hspace{17ex}
    \rm{(without BLI)}
    \hspace{-9.0ex}
    \end{aligned}$$
    
    If BLI is present at or near the trailing edge, the upstream boundary
    layer and corresponding surface dissipation $\Phi_{\rm surf}$ will be
    mostly unaffected. But the viscous fluid flowing into the wake is now
    reduced by the ingestion fraction
    ${f_{\rm {\scriptscriptstyle BLI}_{\scriptstyle \,f}}}$, so that the
    wake dissipation $\Phi_{\rm wake}$ will be reduced by the same fraction.
    This then gives the following overall fuselage dissipation coefficient
    for the BLI case. 
    $$\begin{aligned}
    C_{\!D_{\rm fuse}} &\!=\! &
    \frac{\Phi_{\rm surf}\!-\!P_{V_{\rm surf}} \,+\, 
           (\Phi_{\rm wake}\!-\!P_{V_{\rm wake}})(1\!-\!{f_{\rm {\scriptscriptstyle BLI}_{\scriptstyle \,f}}})}
         {{\textstyle \frac{1}{2}}\rho_{\scriptscriptstyle \infty}V_{\!{\scriptscriptstyle \infty}}^3 S} 
    \hspace{3ex}
    \rm{(with BLI)}
    \hspace{-2ex}
    \\
    C_{\!D_{\rm fuse}} 
    & \!\simeq\! &
    C_{\!D_{\rm fuse}} \,-\, C_{\Phi_{\rm wake}} {f_{\rm {\scriptscriptstyle BLI}_{\scriptstyle \,f}}}
    \hspace{23ex}
    \rm{(with BLI)}
    \hspace{-2ex}
    \\[0.5em]
    \rm{where} \hspace{3ex}
    C_{\Phi_{\rm wake}} & \!=\! & 
    \frac{2 \Theta_{\scriptscriptstyle \infty}}{S} \:-\: \frac{\Theta^*_{\scriptscriptstyle T\!E}}{S} 
    \end{aligned}$$