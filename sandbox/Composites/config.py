# config.py (Centralized Imports)
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import stiffness as stiff
import solve as sol
import mesh as msh
import boundary as bound
import postProcess as pproc
