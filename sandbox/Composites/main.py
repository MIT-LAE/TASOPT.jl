# Execute the provided code to validate results

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import stiffness as stiff
import mesh as msh
import solve as sol
import boundary as bound

def compute_stress_fields(nodes, U, nx, ny, ABD):
    w_disp = U[0::3].reshape(ny + 1, nx + 1)
    sigma_xx = ABD[0, 0] * w_disp
    sigma_yy = ABD[1, 1] * w_disp
    sigma_res = np.sqrt(sigma_xx**2 + sigma_yy**2)
    return sigma_xx, sigma_yy, sigma_res

def plot_stress_fields(nodes, sigma_xx, sigma_yy, sigma_res, nx, ny):
    X, Y = nodes[:, 0].reshape(ny + 1, nx + 1), nodes[:, 1].reshape(ny + 1, nx + 1)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=300)
    
    cmap = "coolwarm"
    for ax, sigma, title in zip(axes, [sigma_xx, sigma_yy, sigma_res],
                                ["Sigma_xx", "Sigma_yy", "Resultant Stress"]):
        c = ax.contourf(X, Y, sigma, cmap=cmap, levels=20)
        fig.colorbar(c, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
    
    plt.tight_layout()
    plt.show()


def plot_displacement(U, nodes, nx, ny):
    w_disp = U[0::3].reshape(ny + 1, nx + 1) * 1e6  # Scale displacement for visibility
    X, Y = nodes[:, 0].reshape(ny + 1, nx + 1), nodes[:, 1].reshape(ny + 1, nx + 1)
    plt.figure(figsize=(8, 6))
    plt.contourf(X, Y, w_disp, cmap="coolwarm", levels=20)
    plt.colorbar(label="Displacement (µm)")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.title("Displacement Contour (Scaled in µm)")
    plt.show()

# Run main simulation
Lx, Ly = 5.0, 1.0
t = 0.005
nx, ny = 10, 10
q0 = 1e6
E1, E2, G12, nu12 = 150e8, 10e8, 5e8, 0.3
ply_angles = [0, 90, 0, 90]

ABD = stiff.compute_ABD_matrices(E1, E2, G12, nu12, ply_angles, t)
nodes, elements, dx, dy = msh.generate_mesh(Lx, Ly, nx, ny)
ndof = 3 * len(nodes)
K = sp.lil_matrix((ndof, ndof))
F = np.zeros(ndof)

# Assemble Stiffness Matrix and Apply Force
for el in elements:
    node_indices = np.array(el)
    dofs = np.ravel([[3 * n, 3 * n + 1, 3 * n + 2] for n in node_indices])
    ke = np.kron(np.eye(4), ABD[:3, :3])
    K[np.ix_(dofs, dofs)] += ke
    force_per_node = q0 * dx * dy / 4  
    for i in range(4):
        F[dofs[i * 3]] += force_per_node

# Apply Cantilever Boundary Conditions
K, F = bound.apply_cantilever_bc(K, F, nodes)

# Solve for Displacements
U = sol.solve_system(K, F)

# Display Results
print("Max displacement:", np.max(np.abs(U[0::3])))
plot_displacement(U, nodes, nx, ny)

sigma_xx, sigma_yy, sigma_res = compute_stress_fields(nodes, U, nx, ny, ABD)
plot_stress_fields(nodes, sigma_xx, sigma_yy, sigma_res, nx, ny)


