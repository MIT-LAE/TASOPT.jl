"""
Module containing all propulsion system related code.

Will cover alternate engine models - NPSS vs Drela's orig. model vs pyCycle
"""
module propsys

export NPSS_run, startNPSS, endNPSS, NPSS_TEsys, NPSS_TEsysOD, 
NPSS_TFsys, NPSS_TFsysOD, NPSS_TFsysOD2

include("../misc/index.inc")
include("../misc/constants.jl")

include("PMSM.jl")
include("NPSS_functions.jl")
include("NPSS_TF_funcs.jl")

end