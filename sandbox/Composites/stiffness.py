import numpy as np
import scipy.sparse as sp


def compute_ABD_matrices(E1, E2, G12, nu12, ply_angles, t):
    num_plies = len(ply_angles)
    Q = np.array([[E1 / (1 - nu12**2 * E2 / E1), nu12 * E2 / (1 - nu12**2 * E2 / E1), 0],
                  [nu12 * E2 / (1 - nu12**2 * E2 / E1), E2 / (1 - nu12**2 * E2 / E1), 0],
                  [0, 0, G12]])
    ABD = np.zeros((6, 6))
    z_k = np.linspace(-t / 2, t / 2, num_plies + 1)
    
    for k in range(num_plies):
        theta = np.radians(ply_angles[k])
        c, s = np.cos(theta), np.sin(theta)
        T = np.array([[c**2, s**2, 2*c*s], [s**2, c**2, -2*c*s], [-c*s, c*s, c**2 - s**2]])
        Q_bar = np.linalg.inv(T) @ Q @ T.T
        z_upper, z_lower = z_k[k + 1], z_k[k]
        ABD[:3, :3] += Q_bar * (z_upper - z_lower)
        ABD[3:, 3:] += Q_bar * ((z_upper**3 - z_lower**3) / 3)
        ABD[:3, 3:] += Q_bar * ((z_upper**2 - z_lower**2) / 2)
        ABD[3:, :3] = ABD[:3, 3:].T
    return ABD


