# [Wings and tails](@id wingtail)

## Wing or Tail Planform

![Piecewise-linear wing or tail surface planform, with break at
$\eta_s$. ](../assets/wingplan.png)

The surface geometry relations derived below correspond to the wing.
Most of these apply equally to the tails if the wing parameters are
simply replaced with the tail counterparts. The exceptions which pertain
to only the wing will be indicated with "(Wing only)" in the subsection
title.

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

### Reference quantities

The aircraft reference quantities are chosen to be simply the values for
the wing. 

$$\begin{aligned}
b_{\rm ref}   &\!= (b)_{\rm wing}\\
S_{\rm ref}   &\!= (S)_{\rm wing}\\
{A\hspace{-0.5ex}R}_{\rm ref} &\!= ({A\hspace{-0.5ex}R})_{\rm wing}\\
c_{\rm ref}   &\!= (c_{\rm ma})_{\rm wing}
\end{aligned}$$

## Surface Airloads

### Lift distribution

The surface lift distribution $\tilde{p}$ is defined in terms of a
baseline piecewise-linear distribution $p{\scriptstyle (\eta)}$ defined
like the chord planform, but with its own taper ratios $\gamma_s$ and
$\gamma_t$. These are actually defined using local section $c_\ell$
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

![Piecewise-linear aerodynamic load $\tilde{p}{\scriptstyle (\eta)}$,
with modifications at center and tip.](../assets/pload.png)

To get the actual aerodynamic load $\tilde{p}$, lift corrections
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
load magnitude $p_o$ and the
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

!!! details "📖 Theory - Surface pitching moment"

    The surface's reference axis is at some specified chordwise fractional
    location $\xi_{\rm ax}$, as shown in the first figure. The
    profile pitching moment acts along the span-axis coordinate
    $y_{\scriptscriptstyle \perp}$, and scales with the normal-plane chord
    $c_{\scriptscriptstyle \perp}$. These are shown in
    the first figure, and
    related to the spanwise and streamwise quantities via the sweep angle.

    $$\begin{aligned}
    y_{\scriptscriptstyle \perp}& = & y /  \cos \Lambda 
    %\label{yperp} \\
    c_{\scriptscriptstyle \perp}& = & c \: \cos \Lambda 
    %\label{cperp} \\
    V_{\scriptscriptstyle \perp}& = & V_{\!{\scriptscriptstyle \infty}}\: \cos \Lambda 
    %\label{Vperp} 
    \end{aligned}$$ 

    The airfoil's pitching moment contribution shown in the figure below is

    $$\begin{aligned}
    dM_{y_{\scriptscriptstyle \perp}} & =  
    \frac{1}{2} \rho V_{\scriptscriptstyle \perp}^2 \: c_{\scriptscriptstyle \perp}^2 \: c_m  \:\: {\rm d}y_{\scriptscriptstyle \perp}
    \\ \\
    c_m {\scriptstyle (\eta)}& =  \left\{
    \begin{array}{lcl}
    \; c_{m_o} & , & 0 < \eta < \eta_o
    \\[0.5em]
    \displaystyle
    \; c_{m_o} + (c_{m_s} \!-\! c_{m_o}\,) 
                \frac{\eta \!-\! \eta_o}{\eta_s\!-\!\eta_o}
    & , & \eta_o < \eta < \eta_s
    \\[0.25em]
    \displaystyle
    \; c_{m_s} + (c_{m_t} \!-\! c_{m_s}\,) 
                \frac{\eta \!-\! \eta_s}{1\!-\!\eta_s}
    & , & \eta_s < \eta < 1
    \end{array}
    \right.
    %\label{cmeta}
    \end{aligned}$$ 

    and including the contribution of the lift load
    $\tilde{p}$ with its moment arm gives the following overall wing
    pitching moment $\Delta M_{\rm wing}$ increment about the axis center
    location. 
    $$\begin{aligned}
    {\rm d}\Delta M_{\rm wing}& \!=\! & 
    \tilde{p}
    \left[ c_{\scriptscriptstyle \perp}\! \left(\xi_{\rm ax}\!-\!{\textstyle \frac{1}{4}}\right) \! \cos\Lambda
          \,-\, (y\!-\!y_o) \tan\Lambda \right]  {\rm d}y 
    \:+\: {\rm d}M_{y_{\scriptscriptstyle \perp}} \, \cos\Lambda \hspace{3em}
    \end{aligned}$$ 

    Integrating this along the whole span then gives the
    total surface pitching moment about its root axis. 

    $$\begin{aligned}
    \Delta M_{\rm wing}& \!=\!  
    (p_o \, b_o + 2\Delta L_o) \, c_o \! \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) 
    \nonumber \\
    & \!+ \;\! 
    \cos^2 \! \Lambda \;
    b \int_{\eta_o}^1 p{\scriptstyle (\eta)}\: c{\scriptstyle (\eta)}\left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) 
    \: {\rm d}\eta 
    \nonumber \\
    & \!-\!\; 
    \frac{b}{2} \, \tan\!\Lambda \: b \int_{\eta_o}^1 p{\scriptstyle (\eta)}(\eta\!-\!\eta_o) 
    \; {\rm d}\eta 
    \nonumber \\
    & \!+\!\;  2 \Delta L_t \left[
    c_o \lambda_t \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) \cos^2\!\Lambda 
    \;-\;
    \frac{b}{2} (1\!-\!\eta_o) \, \tan\!\Lambda \right]
    \nonumber \\
    & \!+\!\; \frac{1}{2} \rho V_{\!{\scriptscriptstyle \infty}}^2 \, \cos^4 \! \Lambda \; b \!
    \int_{\eta_o}^1 \!\! c_m {\scriptstyle (\eta)}\: {c{\scriptstyle (\eta)}}^2 \:\: {\rm d}\eta 
    \\
    %%
    \Delta M_{\rm wing}& \!=\; 
    p_o \, b \, c_o \: \eta_o \,
    (1 \!+\! f_{L_{\scriptstyle o}}) \! \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) 
    \nonumber \\
    & \!+\!\;
    p_o \, b \, c_o \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) \cos^2 \! \Lambda \;
    \frac{1}{3}
    \left[ \,
    \left( 1 + {\textstyle \frac{1}{2}}\left(\lambda_s \!+\!\gamma_s\right) + \lambda_s \gamma_s
          \right)(\eta_s \!-\! \eta_o) 
    \right.
    \nonumber \\
    & \hspace{9em}
    \left.
    \,+\,
    \left( \lambda_s \gamma_s + 
          {\textstyle \frac{1}{2}}\left(\lambda_s \gamma_t\!+\!\gamma_s\lambda_t\right)
        + \lambda_t \gamma_t
          \right)(1 \!-\! \eta_s) \, \right]
    \nonumber \\
    & \!-\!\;
    p_o \, b \, c_o \: \frac{\tan\Lambda}{K_o} \, 
    \frac{1}{12} \left[ \,
          \left( 1 \!+\! 2\gamma_s \right) (\eta_s \!-\! \eta_o)^2 
    \,+\,  \left( \gamma_s \!+\! 2\gamma_t\right) (1 \!-\! \eta_s)^2
    \,+\,  3\left( \gamma_s \!+\! \gamma_t \right) (\eta_s \!-\! \eta_o)(1 \!-\! \eta_s)
            \, \right]
    \nonumber \\
    & \!+\!\;
    2 \, p_o \, b \, c_o \, f_{L_{\scriptstyle t}}\, \lambda_t \, \gamma_t \left[
    K_o \lambda_t \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) \cos^2\!\Lambda 
    \;-\; {\textstyle \frac{1}{2}}(1\!-\!\eta_o) \, \tan\!\Lambda \right]
    \nonumber \\
    & \!+\!\; 
    \frac{1}{2} \rho V_{\!{\scriptscriptstyle \infty}}^2 \, S \, c_o \, \frac{\cos^4\!\Lambda}{K_c} 
    \frac{1}{12}
    \left[ \left(
          c_{m_o} (3             \!+\! 2           \lambda_s \!+\! \lambda_s^2)
        + c_{m_s} (3 \lambda_s^2 \!+\! 2           \lambda_s \!+\!   1        )
          \right)(\eta_s\!-\!\eta_o)  \right.
    \nonumber \\
    &
    \left.
    +     \left(
          c_{m_s} (3 \lambda_s^2 \!+\! 2 \lambda_s \lambda_t \!+\! \lambda_t^2)
        + c_{m_t} (3 \lambda_t^2 \!+\! 2 \lambda_s \lambda_t \!+\! \lambda_s^2)
          \right)(1\!-\!\eta_s)  
    \right]
    \hspace{3em}
    %\label{DMwing}
    \end{aligned}$$ 

    By using the relation 

    $$\begin{aligned}
    p_o \, b & = & \frac{1}{2} \rho V_{\!{\scriptscriptstyle \infty}}^2 \, S \:
    \frac{1}{K_p} \left( C_{\!L} \!-\! \frac{S_{\rm h}}{S} C_{\!L_{\rm h}} \right)
    \end{aligned}$$ 

    and the equation above it gives the equivalent pitching moment coefficient
    constant and $C_{\!L}$ derivative. 

    $$\begin{aligned}
    \Delta C_{\!M_{\rm wing}} \;\equiv\;
    \frac{\Delta M_{\rm wing}}{{\textstyle \frac{1}{2}}\rho V_{\!{\scriptscriptstyle \infty}}^2 S c_o}
    & \!=\;  \Delta C_{m_0} 
    \,+\, \frac{{\rm d}C_m}{{\rm d}C_{\!L}} 
    \left( C_{\!L} \!-\! \frac{S_{\rm h}}{S} C_{\!L_{\rm h}} \right)
    \\
    \frac{{\rm d}C_m}{{\rm d}C_{\!L}} 
    & \!=\! 
    \frac{1}{K_p} \left\{ \rule[-1.25ex]{0ex}{4.5ex}
    \eta_o \,
    (1 \!+\! f_{L_{\scriptstyle o}}) \! \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) 
    \right.
    \nonumber \\
    & 
    + \, \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) \cos^2 \! \Lambda \;
    \frac{1}{3}
    \left[ \,
    \left( 1 + {\textstyle \frac{1}{2}}\left(\lambda_s \!+\!\gamma_s\right) + \lambda_s \gamma_s
          \right)(\eta_s \!-\! \eta_o) 
    \right.
    \nonumber \\
    &
    \left.
    \,+\,
    \left( \lambda_s \gamma_s + 
          {\textstyle \frac{1}{2}}\left(\lambda_s \gamma_t\!+\!\gamma_s\lambda_t\right)
        + \lambda_t \gamma_t
          \right)(1 \!-\! \eta_s) \, \right]
    \nonumber \\
    & 
    - \; \frac{\tan\!\Lambda}{K_o} \,
    \frac{1}{12} \left[ \,
          \left( 1 \!+\! 2\gamma_s \right) (\eta_s \!-\! \eta_o)^2 
    \,+\,  \left( \gamma_s \!+\! 2\gamma_t\right) (1 \!-\! \eta_s)^2 
            \right.
    \nonumber \\
    & 
    \left.
    \,+\, 3\left( \gamma_s \!+\! \gamma_t \right) (\eta_s \!-\! \eta_o)(1 \!-\! \eta_s)
          \rule[-0.5ex]{0ex}{3ex}    \, \right]
    \nonumber \\
    & 
    \left.
    + \; 2 \, f_{L_{\scriptstyle t}}\, \lambda_t \, \gamma_t \left[
    K_o \lambda_t \left( \xi_{\rm ax} \!-\! {\textstyle \frac{1}{4}}\right) \cos^2\!\Lambda 
    \;-\; {\textstyle \frac{1}{2}}(1\!-\!\eta_o) \, \tan\Lambda \right]
    \rule[-1.25ex]{0ex}{4.5ex}  \right\}

    \\
    %
    \Delta C_{m_0} 
    & \!=\! 
    \frac{\cos^4 \!\Lambda}{K_c} \,
    \frac{1}{12}
    \left[ \left(
          c_{m_o} (3             \!+\! 2           \lambda_s \!+\! \lambda_s^2)
        + c_{m_s} (3 \lambda_s^2 \!+\! 2           \lambda_s \!+\!   1        )
          \right)(\eta_s\!-\!\eta_o)  \right.
    \nonumber \\
    & \hspace{1em}
    \left.
    +     \left(
          c_{m_s} (3 \lambda_s^2 \!+\! 2 \lambda_s \lambda_t \!+\! \lambda_t^2)
        + c_{m_t} (3 \lambda_t^2 \!+\! 2 \lambda_s \lambda_t \!+\! \lambda_s^2)
          \right)(1\!-\!\eta_s)  
    \right]
    \hspace{3em}
    %\label{CM0wing}
    \hspace{2em}
    \end{aligned}$$


![Wing pitching moment quantities.](../assets/dMwing.png)

## Wing or Tail Structural Loads

The figure below shows
the airload $\tilde{p}$ again, partly offset by weight load
distributions of the structure and fuel, producing shear and bending
moment distributions.

![Aerodynamic load $\tilde{p}{\scriptstyle (\eta)}$ and weight load
$w{\scriptstyle (\eta)}$, with resulting shear and bending moments. An
optional strut modifies the shear and bending moment as
indicated.](../assets/wingload.png)

## Wing or Tail Stresses

### Normal-plane quantities

The wing and tail surface stress and weight analyses are performed in
the cross-sectional plane, normal to the spanwise axis
$y_{\scriptscriptstyle \perp}$ running along the wing box sketched in the first and third figures. Together
with the normal-plane coordinate and chord relations, the shear
and bending moment are related to the corresponding airplane-axes
quantities and to the sweep angle $\Lambda$ as follows.

$$\begin{aligned}
{\cal S}_{\scriptscriptstyle \perp}& = & {\cal S}
%\label{Sperp} \\
{\cal M}_{\scriptscriptstyle \perp}& = & {\cal M}/ \cos \Lambda
%\label{Mperp}
\end{aligned}$$

### Wing or tail section

The assumed wing or tail airfoil and structural box cross-section is
shown in the figure below. The box is assumed to be the only
structurally-significant element, with the slats, flaps, and spoilers
(if any), represented only by added weight. It is convenient to define
all dimensions as ratios with the local normal-plane chord
$c_{\scriptscriptstyle \perp}$.

$$\begin{aligned}
\bar{h} &\!\equiv\!& \frac{h_{\rm wbox}}{c_{\scriptscriptstyle \perp}} \\
\bar{w} &\!\equiv\!& \frac{w_{\rm wbox}}{c_{\scriptscriptstyle \perp}} \\
\bar{t}_{\rm cap}&\!\equiv\!& \frac{t_{\rm cap}}{c_{\scriptscriptstyle \perp}} \\
\bar{t}_{\rm web}&\!\equiv\!& \frac{t_{\rm web}}{c_{\scriptscriptstyle \perp}} 
\end{aligned}$$

![Wing or tail airfoil and structure cross-section, shown perpendicular
to spar axis. Leading edges, fairings, slats, flaps, and spoilers
contribute to weight but not to the primary
structure.](../assets/wingbox.png)

The maximum height $h_{\rm wbox}$ at the box center corresponds to the
airfoil thickness, so that $\bar{h}$ is the usual "$t/c$" airfoil
thickness ratio. The height is assumed to taper off quadratically to a
fraction $r_h$ at the webs, so that the local height
$h {\scriptstyle (\xi)}$ is 
$$\begin{aligned}
h {\scriptstyle (\xi)}& = & h_{\rm wbox}\left[ \: 1 - (1\!-\!r_h) \xi^2 \: \right]
\end{aligned}$$ 

where $\xi = -1 \ldots 1$ runs chordwise over the
sparbox extent. Typical metal wings and airfoils have
$\bar{w} \simeq 0.5$, $r_h \simeq 0.75$, although these are left as
input parameters. For evaluating areas and approximating the bending
inertia, it's useful to define the simple average and r.m.s. average
normalized box heights. 

$$\begin{aligned}
\bar{h}_{\rm avg}& = & \frac{1}{c_{\scriptscriptstyle \perp}} \int_0^1 h {\scriptstyle (\xi)}\; {\rm d}\xi 
\;=\; \bar{h} \left[ \: 1 - \frac{1}{3}(1\!-\!r_h) \, \right]
\\
\bar{h}_{\rm rms}^2 & = & \frac{1}{c_{\scriptscriptstyle \perp}^2} \int_0^1 h^2 {\scriptstyle (\xi)}\; {\rm d}\xi 
\;=\; \bar{h}^2 \left[ \: 
  1 - \frac{2}{3}(1\!-\!r_h) + \frac{1}{5} (1\!-\!r_h)^2 \, \right]
\end{aligned}$$

The areas and the bending and torsion inertias, all normalized by the
normal chord, can now be determined. 

$$\begin{aligned}
\bar{A}_{\rm fuel}&\!\equiv\!& \frac{A_{\rm fuel}}{c_{\scriptscriptstyle \perp}^2} \;=\;
(\bar{w} - 2 \bar{t}_{\rm web})(\bar{h}_{\rm avg}- 2 \bar{t}_{\rm cap}) 
\\
\bar{A}_{\rm cap}&\!\equiv\!& \frac{A_{\rm cap}}{c_{\scriptscriptstyle \perp}^2} \;=\; 2 \, \bar{t}_{\rm cap}\bar{w} 
\\
\bar{A}_{\rm web}&\!\equiv\!& \frac{A_{\rm web}}{c_{\scriptscriptstyle \perp}^2} \;=\; 2 \, \bar{t}_{\rm web}\, r_h \, \bar{h} 
\\
\bar{I}_{\rm cap}& \simeq & \frac{I_{\rm cap}}{c_{\scriptscriptstyle \perp}^4} \;=\;
  \frac{\bar{w}}{12} 
 \left[ \bar{h}_{\rm rms}^3 - (\bar{h}_{\rm rms}\!\!-\!2\bar{t}_{\rm cap})^3 \right]
\\
\bar{I}_{\rm web}&\!\equiv\!& \frac{I_{\rm web}}{c_{\scriptscriptstyle \perp}^4} \;=\;
 \frac{\bar{t}_{\rm web}\, r_h^3 \, \bar{h}^3}{6}  \; \ll \; \bar{I}_{\rm cap}
\hspace{2em} \mathrm{(typically)}
\\
G\bar{J} &\!\equiv\!& 
\frac{4 (\bar{w} - \bar{t}_{\rm web})^2 (\bar{h}_{\rm avg}- \bar{t}_{\rm cap})^2}
{ \displaystyle 
  2 \frac{ r_h \bar{h} \!-\! \bar{t}_{\rm cap}}{G_{\rm web}\bar{t}_{\rm web}} \:+\:
  2 \frac{     \bar{w} \!-\! \bar{t}_{\rm web}}{G_{\rm cap}\bar{t}_{\rm cap}} }
\end{aligned}$$

```@docs

structures.surfw(po,b,bs,bo,co,zs,
	lambdat,lambdas,gammat,gammas,
	Nload,iwplan,We,neout, dyeout, neinn, dyeinn,
	Winn,Wout,dyWinn,dyWout,
	sweep,wbox,hboxo,hboxs,rh, fLt,
	tauweb,sigcap,sigstrut,Ecap,Eweb,Gcap,Gweb,
	rhoweb,rhocap,rhostrut,rhofuel)

structures.tailpo(S, AR, λa, qne, CLmax)

structures.surfdx(b, bs, bo, λt, λs, sweep)
```
