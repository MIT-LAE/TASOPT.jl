"""
Module containing all propulsion system related code.

Will cover alternate engine models - NPSS vs Drela's orig. model vs pyCycle
"""
module propsys

export NPSS_run, startNPSS, endNPSS, NPSS_TEsys, NPSS_TEsysOD
include("../misc/index.inc")
include("../misc/constants.jl")

include("PMSM.jl")
include("NPSS_functions.jl")

end