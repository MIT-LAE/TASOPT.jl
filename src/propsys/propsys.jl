"""
Module containing all propulsion system related code.

This includes alternate propulsion systems such as turbo-electric architectures.
Will eventually cover alternate engine models - NPSS vs Drela's orig. model vs pyCycle, 
replacing the turbofan code such as `tfoper`, `tfsize`etc.
"""
module propsys

export NPSS_run, startNPSS, endNPSS
import ..TASOPT: __TASOPTindices__, __TASOPTroot__
using ..materials
using DocStringExtensions

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"misc/constants.jl"))
include("cable.jl")
include("PMSM.jl")
include("NPSS_functions.jl")

end