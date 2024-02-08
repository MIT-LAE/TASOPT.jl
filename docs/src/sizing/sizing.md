# Design and evaluation

## [Sizing the aircraft] (@id sizing)

The aircraft is sized via a fixed point iteration for the design mission (`wsize`). The performance of the design can be evaluated for an off-design mission (`woper`).

`wsize` is typically the driving script in an analysis, as is the case in the `size_aircraft!` call (as demonstrated in the [first example] (@ref firstexample)). The sizing analysis calls the various performance routines (e.g., `fusebl`, `surfw`, `cdsum`, `mission`, etc.) as shown in the [TASOPT flowchart](@ref flowchart).

```@docs
<<<<<<< HEAD
TASOPT.wsize(ac, imis,
            itermax, wrlx1, wrlx2, wrlx3,
            initwgt, initeng, iairf, Ldebug, printiter, saveODperf)
=======
TASOPT.wsize(ac)
>>>>>>> af22732b798c848f8abdd63011006b99a0856978

TASOPT.woper(pari, parg, parm, para, pare, 
          parad, pared, itermax, initeng, NPSS_PT, NPSS)
```
---

## [Mission evaluation] (@id missionexec)
A sized aircraft's mission performance can be obtained (`mission!`), along with operation constraints via a pitch trim calculation (`balance`) and balanced field length calculation (`takeoff!`).

```@docs
TASOPT.mission!(pari, parg, parm, para, pare, Ldebug)

TASOPT.takeoff!(pari, parg, parm, para, pare,
    initeng,
    ichoke5, ichoke7)

TASOPT.balance(pari, parg, para, rfuel, rpay, ξpay, itrim)

```

