# [Cabin sizing](@id cabin_sizing)
The TASOPT input file requires the user to specify a radius and general fuselage layout, including the start and end of the pressure vessel, and the start and end of the cylindrical portions. Whether this provided layout is used or not depends on the flag `calculate_cabin_length` in the `Fuselage.Geometry` input field: this flag is `False` by default, implying that the fuselage length and layout does not get recalculated for the input radius. However, if the flag is set to `True`, the function `update_fuse_for_pax!()` is used to update the fuselage length and layout to accommodate the desired number of passengers.

```@docs
structures.update_fuse_for_pax!
```

To resize the fuselage, information about the number of passengers and seat properties are needed. The number of passengers used the resize the fuselage is set in the input file as `exit_limit` in `Mission`; this is assumed to be the maximum number of people that could fit in the cabin in an all-economy layout. The goal then is to determine the minimum cabin length corresponding to the specified radius and seat properties: seat pitch, width and height. Two key functions are used for this purpose:

```@docs
structures.find_floor_angles
```

This function aims to determine the angular placement of the cabin floor inside the fuselage cross-section. If the cabin has a single deck, this is the angle that maximizes the cabin width. If the aircraft is a double-decker, the main (lower) cabin could be at any angle; the main-cabin angle must then be specified so that the upper-cabin angle can be calculated from the cabin height. The double-decker geometry problem is more involved and is described in the subsequent section.

If the fuselage has a single deck, the angular position of the cabin floor is then used by the following function to calculate the maximum width corresponding to that section. For example, if the seat height is 0, the angular placement would be 0 and the maximum width of a single bubble fuselage is ``2R_{fuse}``. When the seat height is greater than 0, the effective cabin width is the distance between the top or bottom of the seats in the cabin, whichever is smallest.

```@docs
structures.find_cabin_width
```

Once the geometric properties of the floor and cabin width have been determined, the seats are placed in the cabin so as to determine the total cabin length. A 10 ft offset from the front of the cabin is assumed. First, the number of seats abreast is determined by 

```@docs
structures.findSeatsAbreast
```

This function assumes an aisle width of 20 inches and a fuselage offset of 6 inches, as seats do not start at the edge of the fuselage. Once the number of seats abreast is known, the number of rows is calculated from the exit-limit number of passengers and the cabin length is calculated from the seat pitch, taking emergency exits into account.

When the cabin length has been determined, the function `update_fuse_for_pax!` proceeds to update the position of all the structural elements that have to be specified by the user, such as the cylinder and pressure vessel ends, horizontal and vertical tail placements, etc. To do this, the relative distances between these elements specified by the user are maintained. For example, the distance from the end of the fuselage to the horizontal and vertical tails is the same in the resized aircraft as in the user-defined parameters.

## Double-decker aircraft
When the cabin has two decks, determining the best layout is more involved as the position of the lower floor and the passenger split between the lower and upper cabins is not known in advance. In TASOPT, this is solved as an optimization problem: the optimization variables are the number of passengers in the lower cabin and the lower-cabin floor angle. The geometry of the cabin is shown in the figure below.

![cabin](../assets/cabin_layout.svg)

The objective function to be minimized is the cabin length, which is the maximum of the lengths of the upper and lower cabins. An optimal design is expected to be one in which both cabins have similar lengths. Constraints are applied on the minimum cargo bay height, ``d_{cargo}`` (by default 1.626 m), and on the minimum height of the upper cabin, ``d_{cabin}`` (by default 2 m). The objective function is 

```@docs
structures.find_double_decker_cabin_length
```
This function returns the maximum of the lengths of either cabins (which is used as the objective function), and the number of seats abreast in the main cabin for a later check (ignored in the optimization process). The optimization is done using the NLopt.jl optimization suite, with the function 

```@docs
structures.optimize_double_decker_cabin
```

Since there are only two optimization variables, this function uses the `GN_AGS` global optimizer. Once the minimum cabin length has been determined, the fuselage and component layouts get updated in `update_fuse_for_pax!()` as in the single deck case.

## Determining radius from seats abreast
In some circumstances, it may be of interest to calculate the minimum fuselage radius that results in a given number of seats abreast. If multiple fuselage radii are being tested, it is likely that the ones providing best performance are those with the narrowest cabin for a given number of seats abreast (and therefore cabin length).

In TASOPT, this problem is solved as an optimization problem: the goal is to find the minimum fuselage radius such that the seats abreast in the main cabin are exactly equal to the desired number. This is done in the function

```@docs
structures.find_minimum_radius_for_seats_per_row
```

In this function, two optimizers are used. First, the global optimizer `GN_DIRECT` is used to find the approximate location of the solution with up to 500 function evaluations. This optimizer cannot handle equality constraints directly; instead, the constraint is introduced with a penalty function. The objective function is simply 
```math
f = R_{fuse}+10^3 \Delta s,
```
where ``\Delta s`` is the difference between the desired and actual seats per row for a given radius. The function that computes this difference for a given radius is
```@docs
structures.check_seats_per_row_diff
```

Once the approximate solution has been found, a local optimizer `LN_NELDERMEAD` is used to polish off the solution at a lower computational cost. Finally, a check is done to verify that the solution meets the equality constraint.
