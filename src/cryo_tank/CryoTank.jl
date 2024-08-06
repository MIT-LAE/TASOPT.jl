"""
`CryoTank` is a module that contains the structural, thermal and energy models of a cryogenic tank.
"""
module CryoTank

using ..engine
using ..atmosphere
using NLsolve
using Roots
using NLopt
using DifferentialEquations

import ..TASOPT: __TASOPTindices__, __TASOPTroot__

export SaturatedMixture, tanksize!

include(__TASOPTindices__)
include("../misc/constants.jl")
include("../utils/bubble_geom.jl")

include("tankWmech.jl")
include("tankWthermal.jl")
include("tanksize.jl")

include("mixture.jl")
include("fuel_thermo.jl")
include("pressure.jl")
include("tanktools.jl")



end