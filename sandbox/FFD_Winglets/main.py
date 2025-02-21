import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def plot_linearly_spaced_points(ax, ffd_vertics, num_points  =10):
    """ Plots linearly spaced points starting and ending at the edges of the FFD box"""
    min_x, max_x = np.min(ffd_vertices[:, 0]), np.max(ffd_vertices[:, 0])
    min_y, max_y = np.min(ffd_vertices[:, 1]), np.max(ffd_vertices[:, 1])
    min_z, max_z = np.min(ffd_vertices[:, 2]), np.max(ffd_vertices[:, 2])
    
    # First vertex
    start_point = np.array([ffd_vertices[0,0], ffd_vertices[0,1], ffd_vertices[0,2]])
    end_point = np.array([ffd_vertices[3,0], ffd_vertices[3,1], ffd_vertices[3,2]])
    
    points_1 = np.linspace(start_point, end_point, num_points)
    
    # Second vertex
    start_point_2 = np.array([ffd_vertices[1,0], ffd_vertices[1,1], ffd_vertices[1,2]])
    end_point_2 = np.array([ffd_vertices[2,0], ffd_vertices[2,1], ffd_vertices[2,2]])
    
    # Third vertex
    start_point_3 = np.array([ffd_vertices[4,0], ffd_vertices[4,1], ffd_vertices[4,2]])
    end_point_3 = np.array([ffd_vertices[7,0], ffd_vertices[7,1], ffd_vertices[7,2]])
    
    # Fourth vertex
    start_point_4 = np.array([ffd_vertices[5,0], ffd_vertices[5,1], ffd_vertices[5,2]])
    end_point_4 = np.array([ffd_vertices[6,0], ffd_vertices[6,1], ffd_vertices[6,2]])
    
    points_2 = np.linspace(start_point_2, end_point_2, num_points)
    points_3 = np.linspace(start_point_3, end_point_3, num_points)
    points_4 = np.linspace(start_point_4, end_point_4, num_points)
    
    ax.plot(points_1[:, 0], points_1[:, 1], points_1[:, 2], color='C0', ms = 12, marker = "o")
    ax.plot(points_2[:, 0], points_2[:, 1], points_2[:, 2], color='C0', ms = 12, marker = "o")
    ax.plot(points_3[:, 0], points_3[:, 1], points_3[:, 2], color='C0', ms = 12, marker = "o")
    ax.plot(points_4[:, 0], points_4[:, 1], points_4[:, 2], color='C0', ms = 12, marker = "o")
    

def compute_end_face_centroids(ffd_vertices):
    """Computes the centroids of the two end faces of an FFD box."""
    min_y = np.min(ffd_vertices[:, 1])
    max_y = np.max(ffd_vertices[:, 1])

    min_y_face = np.array([v for v in ffd_vertices if np.isclose(v[1], min_y)])
    max_y_face = np.array([v for v in ffd_vertices if np.isclose(v[1], max_y)])

    centroid_min_y = np.mean(min_y_face, axis=0) if len(min_y_face) > 0 else None
    centroid_max_y = np.mean(max_y_face, axis=0) if len(max_y_face) > 0 else None

    return centroid_min_y, centroid_max_y

def compute_intersection(new_point_1, new_point_2):
    """Computes the intersection of the line segment with the arbitrary plane at Y = new_point_1[1]."""
    t_intersection = (new_point_1[1] - new_point_2[1]) / (new_point_1[1] - new_point_2[1]) if (new_point_1[1] - new_point_2[1]) != 0 else None

    if t_intersection is not None and 0 <= t_intersection <= 1:
        intersection_x = new_point_2[0] + t_intersection * (new_point_1[0] - new_point_2[0])
        intersection_z = new_point_2[2] + t_intersection * (new_point_1[2] - new_point_2[2])
        return np.array([intersection_x, new_point_1[1], intersection_z])
    
    return None

def plot_ffd_box(ffd_vertices, new_point_1, new_point_2, intersection_point):
    """Plots the FFD box, the line segment, and the intersection point."""
    faces = [
        [ffd_vertices[i] for i in [0, 1, 5, 4]],
        [ffd_vertices[i] for i in [1, 2, 6, 5]],
        [ffd_vertices[i] for i in [2, 3, 7, 6]],
        [ffd_vertices[i] for i in [3, 0, 4, 7]],
        [ffd_vertices[i] for i in [4, 5, 6, 7]],
        [ffd_vertices[i] for i in [0, 1, 2, 3]]
    ]

    arbitrary_plane_y = new_point_1[1]
    arbitrary_plane_vertices = np.array([
        [min(ffd_vertices[:, 0]), arbitrary_plane_y, min(ffd_vertices[:, 2])],
        [max(ffd_vertices[:, 0]), arbitrary_plane_y, min(ffd_vertices[:, 2])],
        [max(ffd_vertices[:, 0]), arbitrary_plane_y, max(ffd_vertices[:, 2])],
        [min(ffd_vertices[:, 0]), arbitrary_plane_y, max(ffd_vertices[:, 2])]
    ])
    arbitrary_plane_faces = [[arbitrary_plane_vertices[i] for i in [0, 1, 2, 3]]]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    ax.add_collection3d(Poly3DCollection(faces, facecolors='whitesmoke', linewidths=1, edgecolors='k', alpha=0.5))
    ax.add_collection3d(Poly3DCollection(arbitrary_plane_faces, facecolors='yellow', linewidths=1, edgecolors='k', alpha=0.3))

    ax.plot([new_point_1[0], new_point_2[0]], [new_point_1[1], new_point_2[1]], [new_point_1[2], new_point_2[2]], 
            'r-', lw=2, label="Line Segment")

    ax.scatter([new_point_1[0], new_point_2[0]], [new_point_1[1], new_point_2[1]], [new_point_1[2], new_point_2[2]], 
               color='red', s=100, label="Segment Endpoints")

    if intersection_point is not None:
        ax.scatter(intersection_point[0], intersection_point[1], intersection_point[2], 
                   color='blue', s=150, edgecolors='k', label="Intersection Point")
    
    plot_linearly_spaced_points(ax, ffd_vertices, num_points=10)

    ax.set_xlim(min(ffd_vertices[:, 0]), max(ffd_vertices[:, 0]))
    ax.set_ylim(min(ffd_vertices[:, 1]), max(ffd_vertices[:, 1]) + 5)
    ax.set_zlim(min(ffd_vertices[:, 2]), max(ffd_vertices[:, 2]))
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    ax.set_title("FFD Box with Line Segment and Intersection")
    ax.legend()
    #plt.show()
    plt.savefig("FFD_box.png")

# Main Execution
ffd_vertices = np.array([
    [24.5, 2.8818, 2.5],
    [37.5, 2.8818, 2.5],
    [48.1, 29.7, 5],
    [44.5, 29.2639, 5],
    [24.5, 2.8818, 5.9150],
    [37.5, 2.8818, 5.9150],
    [48.1, 29.7, 7.5],
    [44.5, 29.5923, 7.5]
])

new_point_1, new_point_2 = compute_end_face_centroids(ffd_vertices)
intersection_point = compute_intersection(new_point_1, new_point_2)

print(f"new_point_1 (Min Y face centroid): {new_point_1}")
print(f"new_point_2 (Max Y face centroid): {new_point_2}")
if intersection_point is not None:
    print(f"Intersection Point: {intersection_point}")
else:
    print("No intersection point found.")

plot_ffd_box(ffd_vertices, new_point_1, new_point_2, intersection_point)
