"""
    gas_tset(alpha, n, hspec, tguess)
    
Calculates temperature for a specified enthalpy.
The constant-cp equivalent is

              t = (hspec - hf) /cp  

     where hf is the heat of formation included in h[t]

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of constituents present
    - `hspec`: specified enthalpy
    - `tguess`: first guess for temperature

    **Output:**
    - `t`: temperature
"""
function gas_tset(alpha, n, hspec, tguess)

      itmax = 10
      ttol = 0.000001

      t = tguess

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)
            res = h - hspec
            res_t = h_t

            dt = -res / res_t

            if (abs(dt) < ttol)
                  return t
            end

            t = t + dt
      end
      println("gas_tset: convergence failed.  dT =", dt)

end # gas_tset

"""
    gas_tsetd(alpha, n, hspec, tguess) 
    
Same as gas_tset, but also returns derivative

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of constituents present
    - `hspec`: specified enthalpy
    - `tguess`: first guess for temperature

    **Outputs:**
    - `t`: temperature
    - `t_hspec`: ?
    - `t_al`: ?
"""
function gas_tsetd(alpha, n, hspec, tguess)
      #-------------------------------------------------------------------
      #     Same as gas_tset, but also returns derivative
      #-------------------------------------------------------------------

      itmax = 10
      ttol = 0.000001

      t = tguess

      t_al = zeros(n)

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)
            res = h - hspec
            res_t = h_t

            dt = -res / res_t

            if (abs(dt) < ttol)
                  t_hspec = 1.0 / h_t

                  for i = 1:n
                        si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
                        t_al[i] = -hi / h_t
                  end

                  return t, t_hspec, t_al
            end

            t = t + dt
      end
      println("gas_tsetd: convergence failed.  dT =", dt)

end # gas_tsetd

"""
    gassum(alpha, n, t)
    
Calculates all gas-mixture properties at specified temperature T, and mixing fractions alpha(.)

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of gas constituents
    - `t`: temperature T,  Kelvin

    **Outputs:**
    - `s`: entropy-complement function s[T]
    - `s_t`: ds/dT
    - `h`: complete enthalpy function h[T]
    - `h_t`: dh/dT
    - `cp`: specific heat cp[T]
    - `r`: ideal-gas constant R
"""
function gassum(alpha, n, t)
      s = 0.0
      h = 0.0
      cp = 0.0
      r = 0.0
      s_t = 0.0
      h_t = 0.0
      for i = 1:n
            si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
            s = s + si * alpha[i]
            h = h + hi * alpha[i]
            cp = cp + cpi * alpha[i]
            r = r + ri * alpha[i]
            s_t = s_t + s_ti * alpha[i]
            h_t = h_t + h_ti * alpha[i]
      end
      return s, s_t, h, h_t, cp, r
end # gassum

"""
    gassumd(alpha, n, t)
    
Same as gassum, but also returns cp_t

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of gas constituents
    - `t`: temperature T,  Kelvin

    **Outputs:**
    - `s`: entropy-complement function s[T]
    - `s_t`: ds/dT
    - `h`: complete enthalpy function h[T]
    - `h_t`: dh/dT
    - `cp`: specific heat cp[T]
    - `r`: ideal-gas constant R
    - `cp_t`: dcp / dT
"""
function gassumd(alpha, n, t)

      s = 0.0
      h = 0.0
      cp = 0.0
      r = 0.0
      s_t = 0.0
      h_t = 0.0
      for i = 1:n
            si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
            s = s + si * alpha[i]
            h = h + hi * alpha[i]
            cp = cp + cpi * alpha[i]
            r = r + ri * alpha[i]
            s_t = s_t + s_ti * alpha[i]
            h_t = h_t + h_ti * alpha[i]
      end

      #---- set dcp/dt via finite difference
      dt = 0.01
      cp_t = 0.0
      for i = 1:n
            si, s_ti, hi, h_ti, cpi1, ri = gasfun(i, t - dt)
            si, s_ti, hi, h_ti, cpi2, ri = gasfun(i, t + dt)
            cp_ti = (cpi2 - cpi1) * 0.5 / dt
            cp_t = cp_t + cp_ti * alpha[i]
      end

      return s, s_t, h, h_t, cp, cp_t, r
end # gassumd

"""
    gas_prat(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)
    
Calculates state change for a specified pressure ratio.
The constant-cp equivalent is the usual isentropic relations, but with epol included.

    g = cp/(cp-r)
    gexp = (g-1)/(g*epol)
    tau = pratio^gexp
    p = po * pratio
    t = to * tau
    (h-hf) = (ho-hf) * tau

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of gas constituents
    - `po`: starting pressure
    - `to`: starting temperature
    - `ho`: starting enthalpy
    - `so`: starting entropy-complement
    - `cpo`: starting specific heat
    - `ro`: starting gas constant
    - `pratio`: pressure ratio
    - `epol`: polytropic efficiency of change process  , if compression
    - `epol`: 1/(polytropic efficiency of change process) , if expansion

    **Outputs:**
    - `p`: ending pressure
    - `t`: ending temperature
    - `h`: ending enthalpy
    - `s`: ending entropy-complement
    - `cp`: ending specific heat
    - `r`: ending gas constant (this will be the same as starting ro)
"""
function gas_prat(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)

      itmax = 10
      ttol = 0.000001

      gexp = ro / (cpo * epol)
      t = to * pratio^gexp

      pile = log(pratio) / epol

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)
            res = (s - so) / r - pile
            res_t = s_t / r

            dt = -res / res_t


            if (abs(dt) < ttol)
                  p = po * pratio
                  return p, t, h, s, cp, r
            end

            t = t + dt
      end

      println("gas_prat: convergence failed.  dT =", dt)
      println("To pratio ='", To, pratio)

end # gas_prat


"""
    gas_pratd(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)
    
Same as gas_prat, but also returns Jacobians w.r.t. po,to,pratio,epol

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of gas constituents
    - `po`: starting pressure
    - `to`: starting temperature
    - `ho`: starting enthalpy
    - `so`: starting entropy-complement
    - `cpo`: starting specific heat
    - `ro`: starting gas constant
    - `pratio`: pressure ratio
    - `epol`: polytropic efficiency of change process  , if compression
    - `epol`: 1/(polytropic efficiency of change process) , if expansion

    **Outputs:**
    - `p`: ending pressure
    - `t`: ending temperature
    - `h`: ending enthalpy
    - `s`: ending entropy-complement
    - `cp`: ending specific heat
    - `r`: ending gas constant (this will be the same as starting ro)
"""
function gas_pratd(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)


      itmax = 10
      ttol = 0.000001

      gexp = ro / (cpo * epol)
      t = to * pratio^gexp

      pile = log(pratio) / epol

      t_al = zeros(n)
      h_al = zeros(n)
      s_al = zeros(n)
      p_al = zeros(n)
      cp_al = zeros(n)
      r_al = zeros(n)

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)

            #------ res( s(t,al), r[al], so, pratio )
            res = (s - so) / r - pile
            res_t = s_t / r

            dt = -res / res_t

            #c        write(*,*) iter, t, dt

            t = t + dt

            if (abs(dt) < ttol)

                  s, s_t, h, h_t, cp, cp_t, r = gassumd(alpha, n, t)
                  res_t = s_t / r

                  #---- res( so, pratio, al )
                  res_so = -1.0 / r
                  res_pi = -1.0 / (epol * pratio)
                  res_ep = pile / epol
                  #c    res_al =  si/r - (s-so)/r^2 * ri

                  t_so = -res_so / res_t
                  t_pi = -res_pi / res_t
                  t_ep = -res_ep / res_t
                  #c    t_al = -res_al/res_t 

                  #c          t(so,pi,ep,al)
                  #c       h( t(so,pi,ep,al),al )
                  #c       s( t(so,pi,ep,al),al )
                  #c      cp( t(so,pi,ep,al),al )
                  #c       r( al )
                  #c    

                  #---- p( po, so, pi )
                  p = po * pratio
                  p_so = 0.0
                  p_po = pratio
                  p_pi = po
                  p_ep = 0.0

                  #---- h( t(so,pi,ep,al),al )
                  h_so = h_t * t_so
                  h_pi = h_t * t_pi
                  h_ep = h_t * t_ep
                  #c    h_al = h_t*t_al + hi

                  #---- s( t(so,pi,ep,al),al )
                  s_so = s_t * t_so
                  s_pi = s_t * t_pi
                  s_ep = s_t * t_ep
                  #c    s_al = s_t*t_al + si

                  for i = 1:n
                        si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
                        res_al = si / r - (s - so) / r^2 * ri
                        t_al[i] = -res_al / res_t
                        h_al[i] = h_t * (-res_al / res_t) + hi
                        s_al[i] = s_t * (-res_al / res_t) + si
                        p_al[i] = 0.0
                        cp_al[i] = cp_t * (-res_al / res_t) + cpi
                        r_al[i] = ri
                  end

                  return p, t, h, s, cp, r, p_po, p_so, t_so, h_so, s_so, p_pi, t_pi, h_pi, s_pi,
                  p_ep, t_ep, h_ep, s_ep,
                  p_al, t_al, h_al, s_al, cp_al, r_al
            end

      end
      println("gas_pratd: convergence failed.  dT =", dt)
      println("To pratio =", to, pratio)

end # gas_pratd

"""
    gas_delh(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

Calculates state change for a specified enthalpy change.
The constant-cp equivalent is the usual isentropic relations, but with epol included.
      
      t - to = delh/cp
      g = cp/(cp-r)
      gexp = (g-1)/(g*epol)
      tau = t/to
      pi = tau^(1/gexp)
      p = po * pi
      (h-hf) = (ho-hf) * tau
      
!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of constituents present
    - `po`: starting pressure
    - `to`: starting temperature
    - `ho`: starting enthalpy
    - `so`: starting entropy-complement
    - `cpo`: starting specific heat
    - `ro`: starting gas constant
    - `delh`: enthalpy change
    - `epol`: polytropic efficiency of change process, if compression and 1/(polytropic efficiency of change process) , if expansion
     
    **Output:**
    - `p`: ending pressure
    - `t`: ending temperature
    - `h`: ending enthalpy
    - `s`: ending entropy-complement
    - `cp`: ending specific heat
    - `r`: ending gas constant (this will be the same as starting ro)

"""
function gas_delh(alpha, n, po, to, ho, so, cpo, ro, delh, epol)
   itmax = 10
      ttol = 0.000001


      t = to + delh / cpo

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)
            res = h - ho - delh
            res_t = h_t

            dt = -res / res_t

            if (abs(dt) < ttol)
                  p = po * exp(epol * (s - so) / r)
                  return p, t, h, s, cp, r
            end

            t = t + dt
      end
      println("gas_delh: convergence failed.  dT =", dt)

end # gas_delh

"""
    gas_delhd(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

Same as gas_delh, but also returns Jacobians w.r.t. po,to,delh

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `alpha(.)`: mass fractions for gas constituents i = 1..n
    - `n`: number of constituents present
    - `po`: starting pressure
    - `to`: starting temperature
    - `ho`: starting enthalpy
    - `so`: starting entropy-complement
    - `cpo`: starting specific heat
    - `ro`: starting gas constant
    - `delh`: enthalpy change
    - `epol`: polytropic efficiency of change process, if compression and 1/(polytropic efficiency of change process) , if expansion
     
    **Output:**
    - `p`: ending pressure
    - `t`: ending temperature
    - `h`: ending enthalpy
    - `s`: ending entropy-complement
    - `cp`: ending specific heat
    - `r`: ending gas constant (this will be the same as starting ro)
    - `p_so`: 
    - `p_po`: 
    - `p_ep`: 
    - `p_ho`: 
    - `t_ho`: 
    - `h_ho`: 
    - `s_ho`: 
    - `p_dh`: 
    - `t_dh`: 
    - `h_dh`: 
    - `s_dh`: 
    - `p_al`: 
    - `t_al`: 
    - `h_al`: 
    - `s_al`: 
    - `cp_al`: 
    - `r_al`: 

"""
function gas_delhd(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

      itmax = 15
      ttol = 0.000001


      t = to + delh / cpo
      t = max(t, 0.2 * to)

      t_al = zeros(n)
      h_al = zeros(n)
      s_al = zeros(n)
      p_al = zeros(n)
      cp_al = zeros(n)
      r_al = zeros(n)

      dt = 0.0
      for iter = 1:itmax

            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)

            res = h - ho - delh
            res_t = h_t

            dt = -res / res_t

            rlx = 1.0
            if (rlx * dt < -0.4 * t)
                  rlx = -0.4 * t / dt
            end
            if (rlx * dt > 0.6 * t)
                  rlx = 0.6 * t / dt
            end

            t = t + rlx * dt

            if (abs(dt) < ttol)
                  s, s_t, h, h_t, cp, cp_t, r = gassumd(alpha, n, t)
                  res_t = h_t

                  res_ho = -1.0
                  res_dh = -1.0

                  t_ho = -res_ho / res_t
                  t_dh = -res_dh / res_t

                  p = po * exp(epol * (s - so) / r)
                  p_r = -p * epol * (s - so) / r^2
                  p_s = p * epol / r
                  p_so = -p * epol / r
                  p_po = exp(epol * (s - so) / r)
                  p_ep = p * (s - so) / r

                  p_t = p_s * s_t
                  p_ho = p_t * t_ho
                  p_dh = p_t * t_dh

                  h_ho = h_t * t_ho
                  h_dh = h_t * t_dh

                  s_ho = s_t * t_ho
                  s_dh = s_t * t_dh

                  for i = 1:n
                        si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
                        res_al = hi
                        t_al[i] = -res_al / res_t
                        p_al[i] = p_t * (-res_al / res_t) + p_s * si + p_r * ri
                        h_al[i] = h_t * (-res_al / res_t) + hi
                        s_al[i] = s_t * (-res_al / res_t) + si
                        cp_al[i] = cp_t * (-res_al / res_t) + cpi
                        r_al[i] = ri
                  end

                  return p, t, h, s, cp, r, p_so, p_po, p_ep,
                  p_ho, t_ho, h_ho, s_ho, p_dh,
                  t_dh, h_dh, s_dh, p_al, t_al,
                  h_al, s_al, cp_al, r_al
            end

      end
      println("gas_delh: convergence failed.  dT =", dt)

end # gas_delhd

"""
    gas_burn(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

Calculates fuel/air mass fraction in combustion with specified start and end temperatures to,t .
Calculates mass fractions of post-combustion constituents
!!! details "🔃 Inputs and Outputs"   
    **Input:**
    `alpha(.)`: mass fractions for air  constituents i = 1..n
    - `beta(.)`: mass fractions for fuel constituents i = 1..n
    - `gamma(.)`: mass fraction changes of air constituents due to combustion
    - `n`: number of constituents present, air is 1..n-1, fuel is n
    - `ifuel`: index specifying fuel molecule
    - `to`: starting air temperatur
    - `tf`: starting fuel temperature
    - `t`: temperature of combustion products
    - `hvap`: fuel enthalpy of vaporization
      
    **Output:**
    `f`: fuel/air mass fraction
    `lambda(.)`: mass fractions for combustion product constituents

"""
function gas_burn(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

      nm = n - 1
      so, s_t, ho, h_t, cpo, ro = gassum(alpha, nm, to)
      sa, s_t, ha, h_t, cpa, ra = gassum(alpha, nm, t)
      sc, s_t, hc, h_t, cpc, rb = gassum(gamma, nm, t)
      sf, s_t, hf, h_t, cpf, rf = gassum(beta, nm, tf)

      #---- add on fuel contribution to hf and hc, which gassum cannot index via ifuel
      si, s_t, hi, h_t, cpi, ri = gasfun(ifuel, tf)
      hf = hf + (hi - hvap) * beta[n]

      #sb, s_t, hb, h_t, cpb, rb = gasfun(ifuel, t)
      #hc = hc + hb * gamma[n] #Add contribution from unburnt fuel when etab < 1

      f = (ha - ho) / (hf - hc)

      lambda = zeros(n)
      buf = Zygote.Buffer(lambda, length(lambda))
      for i = 1:n
            buf[i] = (alpha[i] + f * gamma[i]) / (1.0 + f)
      end
      lambda = copy(buf)

      return f, lambda
end # gas_burn

"""
    gas_burnd(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)
    
Same as gas_burn, but also returns derivatives.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    `alpha(.)`: mass fractions for air  constituents i = 1..n
    - `beta(.)`: mass fractions for fuel constituents i = 1..n
    - `gamma(.)`: mass fraction changes of air constituents due to combustion
    - `n`: number of constituents present, air is 1..n-1, fuel is n
    - `ifuel`: index specifying fuel molecule
    - `to`: starting air temperatur
    - `tf`: starting fuel temperature
    - `t`: temperature of combustion products
    - `hvap`: fuel enthalpy of vaporization
      
    **Output:**
    `f`: fuel/air mass fraction
    `lambda(.)`: mass fractions for combustion product constituents
    `f_t`: 
    `l_to`: 
    `l_tf`: 
    `l_t`: 

"""
function gas_burnd(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

      nm = n - 1
      so, s_t, ho, ho_to, cpo, ro = gassum(alpha, nm, to)
      sa, s_t, ha, ha_t, cpa, ra = gassum(alpha, nm, t)
      sc, s_t, hc, hc_t, cpc, rb = gassum(gamma, nm, t)
      sf, s_t, hf, hf_tf, cpf, rf = gassum(beta, nm, tf)

      #---- add on fuel contribution to hf and hc, which gassum cannot index via ifuel
      si, s_t, hi, hi_tf, cpi, ri = gasfun(ifuel, tf)
      hf = hf + (hi - hvap) * beta[n]
      hf_tf = hf_tf + hi_tf * beta[n]

      #sb, s_t, hb, hb_t, cpb, rb = gasfun(ifuel, t)
      #hc = hc + hb * gamma[n] #Add contribution from unburnt fuel when etab < 1
      #hc_t = hc_t + hb_t * gamma[n] 

      f = (ha - ho) / (hf - hc)
      f_ho = -1.0 / (hf - hc)
      f_ha = 1.0 / (hf - hc)
      f_hc = f / (hf - hc)
      f_hf = -f / (hf - hc)

      f_to = f_ho * ho_to
      f_tf = f_hf * hf_tf
      f_t = f_ha * ha_t +
            f_hc * hc_t

      lambda = zeros(n)
      l_to = zeros(n)
      l_tf = zeros(n)
      l_t = zeros(n)
      for i = 1:n
            lambda[i] = (alpha[i] + f * gamma[i]) / (1.0 + f)
            l_f = (gamma[i] - lambda[i]) / (1.0 + f)
            l_to[i] = l_f * f_to
            l_tf[i] = l_f * f_tf
            l_t[i] = l_f * f_t
      end

      return f, lambda, f_to, f_tf, f_t, l_to, l_tf, l_t
end # gas_burnd

"""
    gas_mach(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

Calculates state change for a specified Mach number change.
The constant-cp equivalent is the usual isentropic relations, but with epol included.
      
    g = cp/(cp-r)
    gexp = (g-1)/(g*epol)
    tau = (1 + 0.5*(g-1)*mo^2) / (1 + 0.5*(g-1)*m^2)
    pi = tau^(1/gexp)
    p = po * pi
    t = to * tau
    (h-hf) = (ho-hf) * tau
      
    !!! details "🔃 Inputs and Outputs"     
    **Input:**
    `alpha(.)`: mass fractions for gas constituents i = 1..n
    `n`: number of constituents present
    `po`: starting pressure
    `to`: starting temperature
    `ho`: starting enthalpy
    `so`: starting entropy-complement
    `cpo`: starting specific heat
    `ro`: starting gas constant
    `mo`: starting Mach number
    `m`: ending Mach number
    `epol`: polytropic efficiency of change process  , if compression
    `epol`: 1/(polytropic efficiency of change process) , if expansion
      
    **Output:**
    `p`: ending pressure
    `t`: ending temperature
    `h`: ending enthalpy
    `s`: ending entropy-complement
    `cp`: ending specific heat
    `r`: ending gas constant (this will be the same as starting ro)

"""
function gas_mach(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

      itmax = 10
      ttol = 0.000001

      uosq = mo^2 * cpo * ro / (cpo - ro) * to

      #---- initial guess for temperature, using constant-gamma relation
      t = to * (1.0 + 0.5 * ro / (cpo - ro) * mo^2) /
          (1.0 + 0.5 * ro / (cpo - ro) * m^2)


      dt = 0.0
      #---- Newton iteration for actual temperature
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)
            cp_t = 0.0    # could evaluate this from the cp[T] splines (later)
            #
            #------ usq( t, cp(t,al), r[al] , m )
            usq = m^2 * cp * r / (cp - r) * t
            usq_t = m^2 * cp * r / (cp - r)
            usq_cp = m^2 * r / (cp - r) * t - usq / (cp - r)

            #------ res( h[t], usq( t, cp[t] ) )
            res = h + 0.5 * usq - ho - 0.5 * uosq
            res_t = h_t + 0.5 * (usq_t + usq_cp * cp_t)
            #
            dt = -res / res_t

            #        write(*,*) iter, t, dt

            if (abs(dt) < ttol)
                  p = po * exp(epol * (s - so) / ro)
                  return p, t, h, s, cp, r

            end

            t = t + dt
      end
      println("gas_mach: convergence failed.  dT =", dt)


end # gas_mach

"""
    gas_machd(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

Same as gas_mach, but also returns derivatives

"""
function gas_machd(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

      itmax = 10
      ttol = 0.000001

      uosq = mo^2 * cpo * ro / (cpo - ro) * to
      uosq_mo = mo * 2.0 * cpo * ro / (cpo - ro) * to
      uosq_cpo = mo^2 * ro / (cpo - ro) * to - uosq / (cpo - ro)
      uosq_ro = mo^2 * cpo / (cpo - ro) * to + uosq / (cpo - ro)
      uosq_to = mo^2 * cpo * ro / (cpo - ro)

      #---- initial guess for temperature, using constant-gamma relation
      t = to * (1.0 + 0.5 * ro / (cpo - ro) * mo^2) /
          (1.0 + 0.5 * ro / (cpo - ro) * m^2)

      t_al = zeros(n)
      p_al = zeros(n)
      h_al = zeros(n)
      s_al = zeros(n)
      cp_al = zeros(n)
      r_al = zeros(n)

      dt = 0.0
      #---- Newton iteration for actual temperature
      for iter = 1:itmax
            #c       gassum(alpha,n, t, s,s_t, h,h_t, cp,r)
            s, s_t, h, h_t, cp, cp_t, r = gassumd(alpha, n, t)
            #c      cp_t = 0.    # could evaluate this from the cp[T] splines (later)
            #
            #------ usq( m, cp(t,al), r[al] , t )
            usq = m^2 * cp * r / (cp - r) * t
            usq_m = m * 2.0 * cp * r / (cp - r) * t
            usq_cp = m^2 * r / (cp - r) * t - usq / (cp - r)
            usq_r = m^2 * cp / (cp - r) * t + usq / (cp - r)
            usq_t = m^2 * cp * r / (cp - r)

            #------ res( h(t,al), usq( m, cp(t,al), r[al], t ), ho, uosq(mo,cpo,ro,to) )
            res = h + 0.5 * usq - ho - 0.5 * uosq
            res_t = h_t + 0.5 * (usq_t + usq_cp * cp_t)
            #
            dt = -res / res_t

            if (abs(dt) < ttol)

                  #---- res( h(t,al), usq( m, cp(t,al), r[al], t ), ho, uosq(mo,cpo,ro,to) )
                  res_m = 0.5 * usq_m
                  res_ho = -1.0
                  res_mo = -0.5 * uosq_mo
                  res_cpo = -0.5 * uosq_cpo
                  res_ro = -0.5 * uosq_ro
                  res_to = -0.5 * uosq_to
                  #c    res_al = hi + 0.5*(usq_cp*cpi + usq_r*ri)

                  #c    t(to,ho,cpo,ro, mo, m, al)
                  t_to = -res_to / res_t
                  t_ho = -res_ho / res_t
                  t_cpo = -res_cpo / res_t
                  t_ro = -res_ro / res_t
                  t_mo = -res_mo / res_t
                  t_m = -res_m / res_t
                  #c    t_al  = -res_al /res_t

                  #---- p( r[al], s( t(to,ho,cpo,ro,mo,m,al),al ), so, po )
                  p = po * exp(epol * (s - so) / r)
                  p_r = -p * epol * (s - so) / r^2
                  p_s = p * epol / r
                  p_so = -p * epol / r
                  p_po = exp(epol * (s - so) / r)
                  p_ep = p * (s - so) / r

                  #---- p( t(to,ho,cpo,ro,mo,m,al), so, po, al )
                  p_t = p_s * s_t
                  #c    p_al = p_s*si + p_r*ri

                  #---- p( t(to,ho,cpo,ro,mo,m), so, po, al )
                  #c    p_al = p_s*si + p_r*ri + p_t*t_al

                  #---- p( so, po, to,ho,cpo,ro,mo,m, al )
                  #cc   p_so =            p_so
                  #cc   p_po =            p_po
                  p_to = p_t * t_to
                  p_ho = p_t * t_ho
                  p_cpo = p_t * t_cpo
                  p_ro = p_t * t_ro
                  p_mo = p_t * t_mo
                  p_m = p_t * t_m
                  #cc   p_al  =           p_al

                  #---- h( t(to,ho,so,cpo,ro,mo,m,al),al )
                  h_to = h_t * t_to
                  h_ho = h_t * t_ho
                  h_cpo = h_t * t_cpo
                  h_ro = h_t * t_ro
                  h_mo = h_t * t_mo
                  h_m = h_t * t_m
                  #c    h_al  = h_t*t_al + hi

                  #---- s( t(to,ho,so,cpo,ro,mo,m,al),al )
                  s_to = s_t * t_to
                  s_ho = s_t * t_ho
                  s_cpo = s_t * t_cpo
                  s_ro = s_t * t_ro
                  s_mo = s_t * t_mo
                  s_m = s_t * t_m
                  #c    s_al  = s_t*t_al + si

                  for i = 1:n
                        si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
                        res_al = hi + 0.5 * (usq_cp * cpi + usq_r * ri)
                        t_al[i] = -res_al / res_t
                        p_al[i] = p_t * (-res_al / res_t) + p_s * si + p_r * ri
                        h_al[i] = h_t * (-res_al / res_t) + hi
                        s_al[i] = s_t * (-res_al / res_t) + si
                        cp_al[i] = cp_t * (-res_al / res_t) + cpi
                        r_al[i] = ri
                  end

                  return p, t, h, s, cp, r,
                  p_so, p_po, p_ep,
                  p_to, t_to, h_to, s_to,
                  p_ho, t_ho, h_ho, s_ho,
                  p_m, t_m, h_m, s_m,
                  p_al, t_al, h_al, s_al, cp_al, r_al


            end

            t = t + dt
      end

      println("gas_mach: convergence failed.  dT =", dt)

end # gas_machd

"""
    gas_mass(alpha, n, po, to, ho, so, cpo, ro, mflux, Mguess)

Calculates state a specified mass flux.
Mguess specifies the initial guess, and also selects either the subsonic or the supersonic branch.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    `alpha(.)  mass fractions for gas constituents i = 1..n
    `n         number of constituents present
    `po        total pressure
    `to        total temperature
    `ho        total enthalpy
    `so        total entropy-complement
    `cpo       total specific heat
    `ro        total gas constant
    `mflux     specified mass flux = rho u = mdot/A
    `Mguess    specifies the initial guess for the static quantities
      
    **Output:**      
    `p       static pressure
    `t       static temperature
    `h       static enthalpy
    `s       static entropy-complement
    `cp      static specific heat
    `r       static gas constant (this will be the same as total ro)
"""
function gas_mass(alpha, n, po, to, ho, so, cpo, ro, mflux, Mguess)

      itmax = 25
      ttol = 0.000001

      t = to / (1 + 0.5 * ro / (cpo - ro) * Mguess^2)

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)

            p = po * exp((s - so) / r)
            p_t = p * s_t / r

            res = 2.0 * p^2 / (t * r)^2 * (ho - h) - mflux^2
            res_t = 4.0 * p * p_t / (t * r)^2 * (ho - h) -
                    4.0 * p^2 / (t * r)^2 * (ho - h) / t +
                    2.0 * p^2 / (t * r)^2 * (-h_t)

            dt = -res / res_t

            rlx = 1.0
            if (rlx * dt > 0.8 * (to - t))

                  rlx = 0.8 * (to - t) / dt

            end

            if (abs(dt) < ttol)

                  return p, t, h, s, cp, r
            end

            t = t + rlx * dt
      end
      println("gas_mass: convergence failed.  dT =", dt)

end # gas_mass

"""
    gasfuel(ifuel, n)

Returns mass fraction of constituent changes as a result of combustion with atmospheric oxygen

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ifuel   index of fuel  (see function gasfun)
    - `n       number of constituents in reaction

    **Output:**
    - `gamma(.)  mass fraction changes due to reaction for i = 1..n
"""
function gasfuel(ifuel, n)

      #---- molar weights of C, H, O, N atoms
      wchon = [12.01078, 1.00795, 15.99943, 14.00672]
      in2 = 1
      io2 = 2
      ico2 = 3
      ih2o = 4
      #
      #---- get number of atoms in fuel molecule
      nchon = gaschem(ifuel)
      #
      #---- molar weight of fuel
      wfuel = wchon[1] * float(nchon[1]) +
              wchon[2] * float(nchon[2]) +
              wchon[3] * float(nchon[3]) +
              wchon[4] * float(nchon[4])
      #
      #---- molar weights of combustion gases, to balance atoms in reaction equation
      #-     Wfuel + W_O2  =  W_N2 + W_CO2 + W_H2O + W_N2
      wn2 = (wchon[4] * 2.0) * float(nchon[4]) * 0.5
      wco2 = (wchon[1] + wchon[3] * 2.0) * float(nchon[1])
      wh2o = (wchon[3] + wchon[2] * 2.0) * float(nchon[2]) * 0.5
      wo2 = (wchon[3] * 2.0) * (float(nchon[1]) +
                                float(nchon[2]) * 0.25 -
                                float(nchon[3]) * 0.5)
      #
      #---- store weight fractions of reaction gases, relative to fuel weight
      #-    (left side is negative, right side is positive)
      # Zygote cannot handle this...
      # gamma = zeros(n)
      # for i = 1:n
      #       gamma[i] = 0.0
      # end
      # gamma[in2] = wn2 / wfuel
      # gamma[io2] = -wo2 / wfuel
      # gamma[ico2] = wco2 / wfuel
      # gamma[ih2o] = wh2o / wfuel

      gamma = zeros(n)
      buf = Zygote.Buffer(gamma, length(gamma))
      for i = 1:n
            buf[i] = 0.0
      end
      buf[in2] = wn2 / wfuel
      buf[io2] = -wo2 / wfuel
      buf[ico2] = wco2 / wfuel
      buf[ih2o] = wh2o / wfuel

      gamma = copy(buf)

      return gamma
end # gasfuel

"""
      fuelLHV(ifuel)

Calculates the lower heating value of a fuel.

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ifuel::Int64`: index of fuel  (see function gasfun)

    **Output:**
    - `LHV::Float64``: lower heating value (J/kg)
"""
function fuelLHV(ifuel)
      # air fractions  
      #        N2      O2      CO2    H2O      Ar       fuel
      alpha = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0] #Air composition does not affect LHV

      # fuel fractions
      beta = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

      n = 6

      #Reactant and product temperatures
      to = 298.15
      tf = 298.15
      t = 298.15

      gamma = gasfuel(ifuel, n) #Change in mass fractions after combustion

      #Complete combustion
      f = -alpha[2]/gamma[2]
      lambda = (alpha + f * gamma) / (1.0 + f) #Product mass fractions after combustion

      nm = n - 1
      _, _, ho, _, _, ro = gassum(alpha, nm, to) #Enthalpy of oxidizer before combustion
      _, _, hp, _, _, rp = gassum(lambda, nm, t) #Enthalpy of products
      _, _, hf, _, _, rf = gassum(beta, nm, tf)

      #---- add on fuel contribution to hf, which gassum cannot index via ifuel
      si, s_t, hi, h_t, cpi, ri = gasfun(ifuel, tf)
      hf = hf + hi * beta[n] #Enthalpy of fuel reactant

      Δh = ho + f*hf - (1 + f) * hp #Heat released in cooling products back to reference temperature
      LHV = Δh / f

      return LHV
end

"""
    gasPr(gas, T)

This function calculates some gas thermodynamic properties of different species, including
viscosity, thermal conductivity, specific heat, and Prandtl number.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `gas::char`: gas name
    - `T::Float64`: temperature (K)
 
    **Outputs:**
    - `R::Float64`: gas constant (J/kg/K)
    - `Pr::Float64`: Prandtl number
    - `γ::Float64`: ratio of specific heats
    - `cp::Float64`: specific heat at constant pressure (J/kg/K)
    - `μ::Float64`: dynamic viscosity (Pa s)
    - `k::Float64`: thermal conductivity (W/m/K) 
"""
function gasPr(gas, T)
      #TODO: replace with new gas model
      if (gas == "air") || (gas == "air_simple")
            μ0 = 1.716e-5
            S_μ = 111
            K0 = 0.0241
            S_k = 194
            T0 = 273

            if gas == "air" 
                  alpha = [0.7532, 0.2315, 0.0006, 0.0020, 0.0127, 0.0]
                  nair = 5
                  s, dsdt, ht, dhdt, cp, R = gassum(alpha, nair, T)
                  
             #The simple model should be used when results are not very sensitive to cp but a speedup is desired
            else #"air_simple"
                  R = 287.1 
                  cp = 1005.0 #Just return constant cp
            end

      elseif (gas == "co2")
            μ0 = 1.370e-5
            S_μ = 222
            K0 = 0.0146
            S_k = 1800
            T0 = 273

            igas = 3
            s, s_t, h, h_t, cp, R = gasfun(igas, T)

      elseif (gas == "n2")
            μ0 = 1.663e-5
            S_μ = 107
            K0 = 0.0242
            S_k = 150
            T0 = 273

            igas = 1
            s, s_t, h, h_t, cp, R = gasfun(igas, T)
      elseif (gas == "o2")
            μ0 = 1.919e-5
            S_μ = 139
            K0 = 0.0244
            S_k = 240  
            T0 = 273   

            igas = 2
            s, s_t, h, h_t, cp, R = gasfun(igas, T)

      elseif (gas == "h2o")
            #parameters obtained by a fit to NIST data
            μ0 = 1.147e-5
            S_μ = 1010
            K0 = 0.02133
            S_k = 8331   
            T0 = 350 
            
            igas = 4
            s, s_t, h, h_t, cp, R = gasfun(igas, T)

      elseif (gas == "ch4")
            #parameters obtained by a fit to NIST data
            μ0 = 1.024e-05 
            S_μ = 179.1
            K0 = 0.02977
            S_k = 2859   
            T0 = 273 
            
            igas = 11
            s, s_t, h, h_t, cp, R = gasfun(igas, T)

      elseif (gas == "h2")
            #parameters obtained by a fit to NIST data
            μ0 = 8.485e-06
            S_μ = 99.95
            K0 = 0.1701
            S_k = 161.1    
            T0 = 273
            
            igas = 40
            s, s_t, h, h_t, cp, R = gasfun(igas, T)
      end
      

      # Apply Sutherland`s laws for viscosity and thermal conductivity
      μ = μ0 * (T / T0)^(3 / 2) * ( (T0 + S_μ) / (T + S_μ) )
      k = K0 * (T / T0)^(3 / 2) * ( (T0 + S_k) / (T + S_k) )

      Pr = cp * μ / k
      γ = cp / (cp - R)

      return R, Pr, γ, cp, μ, k 
end

"""
     Calculates temperature for a specified enthalpy for a single gas species.
     The constant-cp equivalent is

              t = (hspec - hf) /cp  

     where hf is the heat of formation included in h[t]

  # Input:
     igas      gas index
     hspec     specified enthalpy
     tguess    first guess for temperature

  # Output:
     t        temperature
"""
function gas_tset_single(igas, hspec, tguess)

      itmax = 15
      ttol = 0.000001

      t = tguess

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gasfun(igas, t)
            res = h - hspec
            res_t = h_t

            dt = -res / res_t

            if (abs(dt) < ttol)
                  return t
            end

            t = t + dt
      end
      println("gas_tset_single: convergence failed.  dT =", dt)

end # gas_tset