
import scipy.sparse.linalg as spla

def solve_system(K, F):
    K = K.tocsc()
    U = spla.spsolve(K, F)
    return U