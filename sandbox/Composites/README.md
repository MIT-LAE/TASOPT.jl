# Finite Element Analysis (FEA) Stiffness Computation

## Overview
This project implements stiffness matrix computations for both composite laminates and homogeneous linear elastic materials in a finite element framework. The code is modular and can be used for structural analysis using the Finite Element Method (FEM).

## Features
- **Composite Laminate Analysis**: Computes the **ABD stiffness matrix** using Classical Laminate Theory (CLT).
- **Homogeneous Material Analysis**: Computes the **global stiffness matrix** for a 2D linear elastic material.
- **Sparse Matrix Representation**: Efficient handling of large-scale problems using SciPy sparse matrices.

## Classical Laminate Theory (CLT)
The stress-strain relationship for a lamina in the principal material coordinates is given by:

\[
 \begin{bmatrix} \sigma_1 \\ \sigma_2 \\ \tau_{12} \end{bmatrix} = 
 \begin{bmatrix} Q_{11} & Q_{12} & 0 \\ Q_{12} & Q_{22} & 0 \\ 0 & 0 & Q_{66} \end{bmatrix} 
 \begin{bmatrix} \varepsilon_1 \\ \varepsilon_2 \\ \gamma_{12} \end{bmatrix}
\]

where \( Q_{ij} \) are the transformed reduced stiffness coefficients.

The ABD matrix is computed as:

\[
 \begin{bmatrix} N_x \\ N_y \\ N_{xy} \\ M_x \\ M_y \\ M_{xy} \end{bmatrix} =
 \begin{bmatrix} A_{11} & A_{12} & A_{16} & B_{11} & B_{12} & B_{16} \\
                 A_{12} & A_{22} & A_{26} & B_{12} & B_{22} & B_{26} \\
                 A_{16} & A_{26} & A_{66} & B_{16} & B_{26} & B_{66} \\
                 B_{11} & B_{12} & B_{16} & D_{11} & D_{12} & D_{16} \\
                 B_{12} & B_{22} & B_{26} & D_{12} & D_{22} & D_{26} \\
                 B_{16} & B_{26} & B_{66} & D_{16} & D_{26} & D_{66} 
 \end{bmatrix} 
 \begin{bmatrix} \varepsilon_x^0 \\ \varepsilon_y^0 \\ \gamma_{xy}^0 \\ \kappa_x \\ \kappa_y \\ \kappa_{xy} \end{bmatrix}
\]

where:
- \( A_{ij} \) are the extensional stiffness coefficients,
- \( B_{ij} \) are the coupling stiffness coefficients,
- \( D_{ij} \) are the bending stiffness coefficients.

## Installation
Ensure you have Python 3.x installed, then install dependencies:
```sh
pip install numpy scipy
```

## Usage
### Compute ABD Matrix for Laminates
```python
from stiffness import compute_ABD_matrices

E1 = 135e9   # Longitudinal Young's modulus (Pa)
E2 = 10e9    # Transverse Young's modulus (Pa)
G12 = 5e9    # Shear modulus (Pa)
nu12 = 0.3   # Poisson's ratio
ply_angles = [0, 45, -45, 90]  # Ply orientations (degrees)
t = 0.005    # Total laminate thickness (m)

ABD_matrix = compute_ABD_matrices(E1, E2, G12, nu12, ply_angles, t)
print(ABD_matrix)
```

### Compute Stiffness Matrix for Homogeneous Material
```python
from stiffness import compute_stiffness_matrix
import numpy as np

E = 200e9   # Young's modulus (Pa)
nu = 0.3    # Poisson's ratio
t = 0.01    # Thickness (m)

# Define nodes and elements
nodes = np.array([[0, 0], [1, 0], [0, 1]])  # Nodal coordinates
connectivity = np.array([[0, 1, 2]])  # Single triangular element

K = compute_stiffness_matrix(E, nu, t, nodes, connectivity)
print(K.toarray())
```

## Code Structure
```
/project-directory/
 ├── stiffness.py      # Stiffness matrix computation functions
 ├── boundary.py       # Boundary conditions handling
 ├── mesh.py           # Mesh generation functions
 ├── solve.py          # Finite element solver
 ├── postProcess.py    # Post-processing and visualization
 ├── config.py         # Centralized configuration and imports
 ├── main.py           # Main execution script
 ├── README.md         # Documentation (this file)
```

## Future Improvements
- Extend to **3D elasticity** for solid elements.
- Implement **quadrilateral elements**.
- Parallelize global stiffness assembly.

## License
This project is licensed under the MIT License.
