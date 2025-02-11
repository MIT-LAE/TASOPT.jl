from config import plt

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