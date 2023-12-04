# Propulsion system

A turbofan model is provided in `TASOPT.jl`

## Turbofan model

```@docs
engine.tfcalc!

engine.tfsize!

engine.tfweight

```
## Turbofan Maps

```@docs
engine.Ncmap(pratio, mb, piD, mbD, NbD, Cmap)

engine.ecmap(pratio, mb, piD, mbD, Cmap, effo, piK, effK)

engine.Ncmap1(pratio, m, piD, mbD, NbD, ABCDm, iabcd, Tr, pr)

engine.ecmap1(pratio, m, piD, mbD, ABCDm, iabcd, effo, Tr, pr)

engine.etmap(dh, mb, Nb, piD, mbD, NbD, ept0, Tmap, Tt, cpt, Rt)

engine.Pimap(mb, Nb, piD, mbD, NbD, Cmap)

engine.tfoper!

```

## Turbofan Cooling

```@docs

engine.mcool(ncrowx, Tmrow, Tt3, Tt4, dTstreak, Trrat, efilm, tfilm, StA)

engine.Tmcalc(ncrowx, ncrow, Tt3, Tt4, dTstreak, Trrat, efilm, tfilm, StA, epsrow)


```
