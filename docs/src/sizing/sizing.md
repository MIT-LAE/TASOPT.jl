# Design and evaluation

## [Sizing the aircraft] (@id sizing)

The aircraft is sized via a fixed point iteration for the design mission (`wsize`). The performance of the design can be evaluated for an off-design mission (`woper`).

`wsize` is typically the driving script in an analysis, as is the case in the `size_aircraft!` call (as demonstrated in the [first example] (@ref firstexample)). The sizing analysis calls the various performance routines (e.g., `fusebl`, `surfw`, `cdsum`, `mission`, etc.) as shown in the [TASOPT flowchart](@ref flowchart).

```@docs
TASOPT.wsize(ac; imission = 1, itermax=35,
    wrlx1=0.5, wrlx2=0.9, wrlx3=0.5, initwgt=false, initeng=0, 
    iairf=1, Ldebug=false, printiter=true, saveODperf=false)

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

