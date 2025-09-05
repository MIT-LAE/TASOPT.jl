"""
Module containing all propulsion system related code.

This includes alternate propulsion systems such as turbo-electric architectures.
Will eventually cover alternate engine models - NPSS vs Drela's orig. model vs pyCycle, 
replacing the turbofan code such as `tfoper`, `tfsize`etc.
"""
module propsys

import ..TASOPT: __TASOPTindices__, __TASOPTroot__
using ..materials
using ..engine
using DocStringExtensions

include(__TASOPTindices__)
include("cable.jl")
include("PMSM.jl")

end