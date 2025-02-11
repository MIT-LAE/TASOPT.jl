import numpy as np

def apply_cantilever_bc(K, F, nodes):
    constrained_dofs = []
    for i, (x, _) in enumerate(nodes):
        if np.isclose(x, 0.0):  # Fix all nodes along x = 0 (cantilevered edge)
            constrained_dofs.extend([3 * i, 3 * i + 1, 3 * i + 2])
    
    for dof in constrained_dofs:
        K[dof, :] = 0
        K[:, dof] = 0
        K[dof, dof] = 1
        F[dof] = 0
    
    return K, F