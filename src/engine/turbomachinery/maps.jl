# TASOPT turbomachinery maps
# These maps are used in log-form equations for the polytropic efficiency

# The parameters come from "tfmap.inc" in Mark Drela's original TASOPT 

# Compressor maps
# Based on the `ccc` option in Drela's TASOPT, with the a, b, k parameters
# in his documentation used for the HPC
#               a     b     k     mo     da    c    d     C    D
const Cmapf = [3.50, 0.80, 0.03, 0.75, -0.50, 3.0, 6.0,  2.5, 15.0]
const Cmapl = [2.50, 1.00, 0.03, 0.75, -0.20, 3.0, 5.5,  4.0,  6.0]
const Cmaph = [1.5,  5.00, 0.03, 0.75, -0.35, 3.0, 5.0,  10.5, 3.0]

#              Pcon   Ncon
const Tmapl = [0.15, 0.15]
const Tmaph = [0.15, 0.15]

# Constant efficiency maps
# const Cmapf = [3.50, 0.80, 0.03, 0.95, -0.50, 3.0, 6.0, 0.0, 0.0]
# const Cmapl = [1.90, 1.00, 0.03, 0.95, -0.20, 3.0, 5.5, 0.0, 0.0]
# const Cmaph = [1.75, 2.00, 0.03, 0.95, -0.35, 3.0, 5.0, 0.0, 0.0]
