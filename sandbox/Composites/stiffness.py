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
        #Q_bar = np.linalg.inv(T) @ Q @ T.T   << Incorrect transofrmation
        T_inv = np.linalg.inv(T)
        Q_bar = T_inv @ Q @ T_inv.T
        z_upper, z_lower = z_k[k + 1], z_k[k]
        ABD[:3, :3] += Q_bar * (z_upper - z_lower)
        ABD[3:, 3:] += Q_bar * ((z_upper**3 - z_lower**3) / 3)
        ABD[:3, 3:] += Q_bar * ((z_upper**2 - z_lower**2) / 2)
        ABD[3:, :3] = ABD[:3, 3:].T
    return ABD


def compute_stiffness_matrix(E, nu, thickness, element_nodes, element_connectivity):
    """
    Computes the global stiffness matrix for a 2D linear elastic homogeneous material.
    
    Parameters:
    E : float
        Young's modulus of the material
    nu : float
        Poisson's ratio of the material
    thickness : float
        Thickness of the 2D plane stress/strain model
    element_nodes : ndarray
        Nodal coordinates of the elements
    element_connectivity : ndarray
        Element connectivity matrix
    
    Returns:
    K_global : sparse matrix
        The assembled global stiffness matrix
    """
    num_nodes = element_nodes.shape[0]
    K_global = sp.lil_matrix((2 * num_nodes, 2 * num_nodes))
    
    C = (E / (1 - nu**2)) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu) / 2]
    ])
    
    for elem in element_connectivity:
        node_indices = elem
        coords = element_nodes[node_indices]
        
        # Compute element stiffness matrix (assuming linear triangular elements)
        x = coords[:, 0]
        y = coords[:, 1]
        A = 0.5 * abs(x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]))
        
        B = np.array([
            [y[1] - y[2], 0, y[2] - y[0], 0, y[0] - y[1], 0],
            [0, x[2] - x[1], 0, x[0] - x[2], 0, x[1] - x[0]],
            [x[2] - x[1], y[1] - y[2], x[0] - x[2], y[2] - y[0], x[1] - x[0], y[0] - y[1]]
        ]) / (2 * A)
        
        k_elem = thickness * A * (B.T @ C @ B)
        
        # Assemble into global stiffness matrix
        dof_indices = np.ravel([[2 * i, 2 * i + 1] for i in node_indices])
        for i in range(6):
            for j in range(6):
                K_global[dof_indices[i], dof_indices[j]] += k_elem[i, j]
    
    return K_global.tocsr()