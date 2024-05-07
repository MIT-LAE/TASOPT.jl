"""
Module containing all propulsion system related code.

Will cover alternate engine models - NPSS vs Drela's orig. model vs pyCycle
"""
module propsys

export NPSS_run, startNPSS, endNPSS
import ..TASOPT: __TASOPTindices__, __TASOPTroot__

include(__TASOPTindices__)
include(joinpath(__TASOPTroot__,"misc/constants.jl"))

include("NPSS_functions.jl")

end