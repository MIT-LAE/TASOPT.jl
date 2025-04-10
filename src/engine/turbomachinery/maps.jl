# TASOPT turbomachinery maps
include("map_functions.jl")
include("FanMap.jl")
include("LPCMap.jl")
include("HPCMap.jl")

# These maps are used in exponential-form equations for the polytropic efficiency

# The parameters come from "tfmap.inc" in Mark Drela's original TASOPT 

#----------Compressor maps----------
# Compressor map parameters optimized to best match the HBTF compressor maps in
# pyCycle (https://github.com/OpenMDAO/pyCycle)

const Cmapf = [3.31140687,  0.77839352,   0.03086818,  0.57042461,  -0.81725615, 6.24179886, 15.42860808,  2.95705214,  0.61792148]
const Cmapl = [2.7602901 ,  1.0172024 ,   0.01628814,  0.70333439,  -0.11845742, 4.18337024, 10.73075544,  2.27229342,  7.38958654]
const Cmaph = [1.61159761,  3.01716623, 0.0115221233,  1.11644367, -0.248346891, 3.27456201, 10.2924591,   12.6465383, 0.717378851]

#The optimized epol0 parameters were Fan: 0.91281925; LPC: 0.9259578; HPC: 0.905517572

#----------Alternative compressor maps----------
# Based on one of the maps in Drela's TASOPT (labeled with the comment 'ccc'), with the a, b, k 
# parameters in his documentation used for the HPC. These were the values that appeared
# to best match the E3 compressor data. See the journal paper (https://arc.aiaa.org/doi/10.2514/3.23024)
# and the NASA report (https://ntrs.nasa.gov/citations/19840021807) for more details on the E3 initiative.
#               a     b     k     mo     da    c    d     C    D
# const Cmapf = [3.50, 0.80, 0.03, 0.75, -0.50, 3.0, 6.0,  2.5, 15.0]
# const Cmapl = [2.50, 1.00, 0.03, 0.75, -0.20, 3.0, 5.5,  4.0,  6.0]
# const Cmaph = [1.5,  5.00, 0.03, 0.75, -0.35, 3.0, 5.0,  10.5, 3.0]

# Constant efficiency maps
# const Cmapf = [3.50, 0.80, 0.03, 0.95, -0.50, 3.0, 6.0, 0.0, 0.0]
# const Cmapl = [1.90, 1.00, 0.03, 0.95, -0.20, 3.0, 5.5, 0.0, 0.0]
# const Cmaph = [1.75, 2.00, 0.03, 0.95, -0.35, 3.0, 5.0, 0.0, 0.0]

#----------Turbine maps----------
const Tmapl = [0.15, 0.15]
const Tmaph = [0.15, 0.15]