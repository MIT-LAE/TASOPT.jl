import numpy as np

def generate_mesh(Lx, Ly, nx, ny):
    dx, dy = Lx / nx, Ly / ny
    nodes = np.array([(i * dx, j * dy) for j in range(ny + 1) for i in range(nx + 1)])
    elements = [(i + j * (nx + 1), i + 1 + j * (nx + 1), i + 1 + (j + 1) * (nx + 1), i + (j + 1) * (nx + 1)) 
                for j in range(ny) for i in range(nx)]
    return nodes, elements, dx, dy