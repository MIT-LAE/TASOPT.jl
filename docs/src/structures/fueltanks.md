# [Fuel tanks](@id fueltanks)

Liquid hydrocarbon fuel is assumed to be stored in the interior of the wings and no additional tanks are needed. The weight of the fuel is accounted for while sizing the wing structure. See [`structures.surfw`](@ref).

However, alternate fuels such as cryogenic liquid hydrogen require additional storage tanks that are insulated pressure vessels.

```@docs
structures.tanksize
structures.tankWmech
structures.tankWthermal

```