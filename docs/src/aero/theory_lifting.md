!!! details "📖 Theory - Wing lift distribution"

    The lift distribution "taper ratios" are defined using local section $c_\ell$
    factors $r_{c_{\ell s}}$ and $r_{c_{\ell t}}$. 

    $$\begin{aligned}
    \gamma_s & =  r_{c_{\ell s}} \, \lambda_s 
    \\
    \gamma_t & =  r_{c_{\ell t}} \, \lambda_t 
    \\
    \frac{p {\scriptstyle (\eta)}}{p_o} \; \equiv \; 
    P {\scriptstyle (\eta \, ; \, \eta_o,\eta_s, \gamma_s,\gamma_t)}
    & =  
    \left\{
    \begin{array}{lcl}
    \; 1 & , & 0 < \eta < \eta_o 
    \\[0.5em]
    \displaystyle
        \: 1 + (\gamma_s\!- 1\,)  \frac{\eta \!-\! \eta_o}{\eta_s\!-\!\eta_o} 
    & , & \eta_o < \eta < \eta_s 
    \\[0.25em]
    \displaystyle
    \gamma_s + (\gamma_t\!-\!\gamma_s) \frac{\eta \!-\! \eta_s}{1\!-\!\eta_s} 
    & , & \eta_s < \eta < 1
    \end{array}
    \right.
    %\label{peta}
    \end{aligned}$$

    To get the actual aerodynamic load ``\tilde{p}``, lift corrections
    $\Delta L_o$ and $\Delta L_t$ are applied to account for the fuselage
    carryover and tip lift rolloff, as sketched in the figure above. The detailed
    shapes of these modifications are not specified, but instead only their
    integrated loads are defined by the following integral relation.

    $$\begin{aligned}
    \frac{L_{\rm wing}}{2} \:\: = \;
    \int_0^{b/2} \!\! \tilde{p} \:\: {\rm d}y & = & \int_0^{b/2} \!\! p \:\: {\rm d}y
    \:+\: \Delta L_o
    \:+\: \Delta L_t
    \end{aligned}$$ 

    The corrections are specified in terms of the center
    load magnitude ``p_o`` and the
    $f_{L_{\scriptstyle o}}, f_{L_{\scriptstyle t}}$ adjustment factors.

    $$\begin{aligned}
    \Delta L_o & =  f_{L_{\scriptstyle o}}\, p_o \, \frac{b_o}{2} 
    \;=\; f_{L_{\scriptstyle o}}\, p_o \, \frac{b}{2} \, \eta_o
    %\label{DLo}
    \\[0.25em]
    \Delta L_t & =  f_{L_{\scriptstyle t}}\, p_t \, c_t 
    \;=\; f_{L_{\scriptstyle t}}\, p_o \, c_o \, \gamma_t\, \lambda_t
    %\label{DLt}
    \\[0.5em]
    f_{L_{\scriptstyle o}}& \simeq  -0.5
    \\
    f_{L_{\scriptstyle t}}& \simeq  -0.05
    \end{aligned}$$

    ### Lift load magnitude (Wing only)

    The wing's $p_o$ center loading magnitude is determined by requiring
    that the aerodynamic loading integrated over the whole span is equal to
    the total weight times the load factor, minus the tail lift.

    $$\begin{aligned}
    L_{\rm total} = 2 \int_0^{b/2} \!\tilde{p}{\scriptstyle (\eta)}\:\: {\rm d}y \;=\; 
    p_o \, b \!\int_0^1 \!P{\scriptstyle (\eta)}\:\: {\rm d}\eta \;+\; 2 \Delta L_o
    \,+\, 2 \Delta L_t 
    &\!=\!& N W \,-\, (L_{\rm htail})_N
    \hspace{2em}
    %\label{totlift}
    \end{aligned}$$ 

    For structural sizing calculations
    $N \!=\! N_{\rm lift}$ is chosen, and the appropriate value of
    $(L_{\rm htail})_N$ is the worst-case (most negative) tail lift expected
    in the critical sizing case. One possible choice is the trimmed tail
    load at dive speed, where $N_{\rm lift}$ is most likely to occur.

    The wing area $S_{def}$ and aspect ratio $AR_{def}$ definitions allow the root chord and the tip lift
    drop $\Delta L_t$ to be
    expressed as 

    $$\begin{aligned}
    c_o & =  b \, K_o
    \\
    \Delta L_t & =  
    f_{L_{\scriptstyle t}}\, p_o \, b \, K_o \: \gamma_t\, \lambda_t
    %\label{DLt2}
    \\
    \mathrm{where} \hspace{2em}
    K_o &\!\equiv\! \frac{1}{K_c \, {A\hspace{-0.5ex}R}} 
    \end{aligned}$$ 

    so that $L_{\rm total}$ can be evaluated to the following. The
    $P{\scriptstyle (\eta)}$ integrals have the form as for
    $C{\scriptstyle (\eta)}$, given by $\int_0^{\eta_o} C \:\: {\rm d}\eta$
    -- $\int_{\eta_s}^1 C^2 \: (\eta\!-\!\eta_s) \:\: {\rm d}\eta$, but with the $\lambda$'s replaced by $\gamma$'s.

    $$\begin{aligned}
    p_o \, b \, K_p & \!=\!  N W - (L_{\rm htail})_N \hspace{2em}
    \\
    \mathrm{where} \hspace{3em}
    K_p & \!=\! 
    \eta_o + {\textstyle \frac{1}{2}}(1 \!+\! \gamma_s) (\eta_s \!\!-\! \eta_o)
        + {\textstyle \frac{1}{2}}(\gamma_s \!+\! \gamma_t) (1 \!-\! \eta_s)
    \nonumber \\
    & \!+\!  f_{L_{\scriptstyle o}}\eta_o \,+\, 2 f_{L_{\scriptstyle t}}K_o \gamma_t\lambda_t 
    \hspace{2em}
    \end{aligned}$$ 

    The root and planform-break loadings can then be
    explicitly determined. 

    $$\begin{aligned}
    p_o & \!=\!  \frac{N W - (L_{\rm htail})_N}{K_p \, b} 
    %\label{podef}
    \\
    p_s & \!=\! p_o \, \gamma_s
    \\
    p_t & \!=\! p_o \, \gamma_t
    \end{aligned}$$