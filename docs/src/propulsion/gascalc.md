# Gas Calculations

Gas calculations used with TASOPT.jl propulsion systems.

## Ideal gas with variable $c_p$
```@docs

engine.gas_tset(alpha, n, hspec, tguess)

engine.gas_tsetd(alpha, n, hspec, tguess)

engine.gasfun

engine.gaschem(igas)

engine.gassum(alpha, n, t)

engine.gassumd(alpha, n, t)

engine.gas_prat(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)

engine.gas_pratd(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)

engine.gas_delh(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

engine.gas_delhd(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

engine.gas_burn(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

engine.gas_burnd(alpha, beta, gamma, n, ifuel, to, tf, t, hvap)

engine.gas_mach(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

engine.gas_machd(alpha, n, po, to, ho, so, cpo, ro, mo, m, epol)

engine.gas_mass(alpha, n, po, to, ho, so, cpo, ro, mflux, Mguess)

engine.gasfuel(ifuel, n)

engine.gasPr(gas, T)



```
