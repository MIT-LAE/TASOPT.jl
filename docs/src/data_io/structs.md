# [Data structure details](@id datastructs)
We're incorporating `struct`s as part of modernizing TASOPT from Fortran. Note that as we represent more parts of an `aircraft` as `struct`s, these may change and will hopefully become a natural and intuitive data structure. 

## Primary `struct`s

Here are the main `struct`s that comprise an `aircraft` object.

```@docs
aircraft

TASOPT.Options

TASOPT.Fuselage

TASOPT.fuselage_tank

TASOPT.Wing

TASOPT.Tail

TASOPT.Engine

```

## Subordinate `struct`s

The above `struct`s are in turn partially composed of subordinate `struct`s, including some to represent materials, geometric layouts, and airfoil aerodynamic performance. 

*Users shouldn't need to mess with these*, but a few for reference:

```@docs

TASOPT.structures.FuselageLayout

TASOPT.structures.Cabin

TASOPT.structures.StructuralMember

TASOPT.structures.WingLayout

TASOPT.structures.WingSection

TASOPT.structures.WingCrossSection

aerodynamics.airfoil

TASOPT.materials.StructuralAlloy

```
