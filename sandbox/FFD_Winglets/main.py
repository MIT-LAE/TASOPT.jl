import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
plt.style.use('seaborn-v0_8-deep')

def plot_intersections(ax, intersections):
    """Plots intersection points of the line segment with the planes."""
    for point in intersections:
        ax.scatter(point[0], point[1], point[2], marker='P', color='C2', s=100, edgecolors='k')
        
        ax.quiver(point[0], point[1], point[2], 0, 1, 0, color='black', length=3, normalize=False)

def plot_linearly_spaced_points(new_point_1, new_point_2, ax, ffd_vertices, num_points  =10):
    """ Plots linearly spaced points starting and ending at the edges of the FFD box"""
    # First vertex
    start_point = np.array([ffd_vertices[0,0], ffd_vertices[0,1], ffd_vertices[0,2]])
    end_point = np.array([ffd_vertices[3,0], ffd_vertices[3,1], ffd_vertices[3,2]])
    points_1 = np.linspace(start_point, end_point, num_points)
        
    # Second vertex
    start_point_2 = np.array([ffd_vertices[1,0], ffd_vertices[1,1], ffd_vertices[1,2]])
    end_point_2 = np.array([ffd_vertices[2,0], ffd_vertices[2,1], ffd_vertices[2,2]])
    points_2 = np.linspace(start_point_2, end_point_2, num_points)
    
    # Third vertex
    start_point_3 = np.array([ffd_vertices[4,0], ffd_vertices[4,1], ffd_vertices[4,2]])
    end_point_3 = np.array([ffd_vertices[7,0], ffd_vertices[7,1], ffd_vertices[7,2]])
    points_3 = np.linspace(start_point_3, end_point_3, num_points)
    
    # Fourth vertex
    start_point_4 = np.array([ffd_vertices[5,0], ffd_vertices[5,1], ffd_vertices[5,2]])
    end_point_4 = np.array([ffd_vertices[6,0], ffd_vertices[6,1], ffd_vertices[6,2]])
    points_4 = np.linspace(start_point_4, end_point_4, num_points)
    
    ax.plot(points_1[:, 0], points_1[:, 1], points_1[:, 2], color='C0', ms = 12, marker = "o")
    ax.plot(points_2[:, 0], points_2[:, 1], points_2[:, 2], color='C0', ms = 12, marker = "o")
    ax.plot(points_3[:, 0], points_3[:, 1], points_3[:, 2], color='C0', ms = 12, marker = "o")
    ax.plot(points_4[:, 0], points_4[:, 1], points_4[:, 2], color='C0', ms = 12, marker = "o")
    
    # Define arbitrary plane locations based on Y values of the points
    arbitrary_planes_y = [points_1[i, 1] for i in range(len(points_1[:,1]))]
    
    # Compute the intestions with all planes
    intersections = compute_intersection_with_planes(new_point_1, new_point_2, arbitrary_planes_y)

    
    
    for y_plane in arbitrary_planes_y:
        plane_vertices = np.array([
            [np.min(ffd_vertices[:, 0]), y_plane, np.min(ffd_vertices[:, 2])],
            [np.max(ffd_vertices[:, 0]), y_plane, np.min(ffd_vertices[:, 2])],
            [np.max(ffd_vertices[:, 0]), y_plane, np.max(ffd_vertices[:, 2])],
            [np.min(ffd_vertices[:, 0]), y_plane, np.max(ffd_vertices[:, 2])]
        ])
        ax.add_collection3d(Poly3DCollection([plane_vertices], facecolors='whitesmoke', linewidths=0.2, linestyle = "--", edgecolors='k', alpha = 0.0))
    
    return intersections

def compute_end_face_centroids(ffd_vertices):
    """Computes the centroids of the two end faces of an FFD box."""
    min_y = np.min(ffd_vertices[:, 1])
    max_y = np.max(ffd_vertices[:, 1])

    min_y_face = np.array([v for v in ffd_vertices if np.isclose(v[1], min_y)])
    max_y_face = np.array([v for v in ffd_vertices if np.isclose(v[1], max_y)])

    centroid_min_y = np.mean(min_y_face, axis=0) if len(min_y_face) > 0 else None
    centroid_max_y = np.mean(max_y_face, axis=0) if len(max_y_face) > 0 else None

    return centroid_min_y, centroid_max_y

def compute_intersection_with_planes(new_point_1, new_point_2, plane_ys):
    """Computes the intersection of the line segment with multiple planes at given Y values."""
    intersections = []
    for y_plane in plane_ys:
        if (new_point_1[1] - new_point_2[1]) != 0:
            t_intersection = (y_plane - new_point_2[1]) / (new_point_1[1] - new_point_2[1])
            if 0 <= t_intersection <= 1:
                intersection_x = new_point_2[0] + t_intersection * (new_point_1[0] - new_point_2[0])
                intersection_z = new_point_2[2] + t_intersection * (new_point_1[2] - new_point_2[2])
                intersections.append(np.array([intersection_x, y_plane, intersection_z]))
    return intersections

def plot_ffd_box(ffd_vertices, new_point_1, new_point_2):
    """Plots the FFD box, the line segment"""
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

     # Create figure and axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.add_collection3d(Poly3DCollection(faces, facecolors='whitesmoke', linewidths=1, edgecolors='k', alpha=0.0))
    #ax.add_collection3d(Poly3DCollection(arbitrary_plane_faces, facecolors='whitesmoke', linewidth=1, linestyle = ':', edgecolors='k', alpha=0.5))

    ax.plot([new_point_1[0], new_point_2[0]], [new_point_1[1], new_point_2[1]], [new_point_1[2], new_point_2[2]], 
            'r', linestyle = '--', lw=2, label="Rotation axis")

    #ax.scatter([new_point_1[0], new_point_2[0]], [new_point_1[1], new_point_2[1]], [new_point_1[2], new_point_2[2]], 
    #           color='red', s=100, label="Segment Endpoints")

    
    intersections = plot_linearly_spaced_points(new_point_1, new_point_2, ax, ffd_vertices, num_points=10)

    plot_intersections(ax, intersections)

    
    ax.set_xlim(min(ffd_vertices[:, 0]), max(ffd_vertices[:, 0]))
    ax.set_ylim(min(ffd_vertices[:, 1]), max(ffd_vertices[:, 1]) + 5)
    ax.set_zlim(min(ffd_vertices[:, 2]), max(ffd_vertices[:, 2]))
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    #ax.set_title("FFD Box with Line Segment and Intersection")
    ax.set_axis_off()

    plt.legend(frameon=False, loc='upper right', prop={'size': 16, 'family': 'Times New Roman'}, ncol=3)
    plt.tight_layout()
    # Adjust figure size
    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 2.5, Size[1] * 1.5, forward=True)
    # High resolution settings
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    
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

#new_point_1, new_point_2 = compute_end_face_centroids(ffd_vertices)

new_point_1 = (28.17598676450, 3.05447084050, 4.24026012500)
new_point_2 = (45.91271297500, 29.38152874000, 6.74589198825)

#plot_ffd_box(ffd_vertices, new_point_1, new_point_2)

