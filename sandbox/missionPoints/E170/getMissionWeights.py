import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the dataset
file_path = "data/E170_histogram_data.csv"
df = pd.read_csv(file_path)

# Define the bounding box limits
mach_min, mach_max = 0.70, 0.78
cl_min, cl_max = 0.49, 0.54

# Compute the center of the bounding box
mach_center = (mach_min + mach_max) / 2
cl_center = (cl_min + cl_max) / 2

# Define 5 integration points (center + 4 surrounding points)
mach_points_5 = np.array([mach_min, mach_max, mach_center, mach_min, mach_max])
cl_points_5 = np.array([cl_min, cl_min, cl_center, cl_max, cl_max])

# Function to extract the average frequency from the 5 surrounding bins
def get_average_frequency(mach, cl, df):
    """Find the average frequency from the 5 surrounding bins for a given Mach and CL."""
    surrounding_bins = df[(df["Mach_bin"] >= mach - 0.01) & (df["Mach_bin"] <= mach + 0.01) &
                          (df["CL_bin"] >= cl - 0.005) & (df["CL_bin"] <= cl + 0.005)]
    return surrounding_bins["Frequency"].mean() if not surrounding_bins.empty else 0

# Get frequencies by averaging over the 5 surrounding bins
frequencies_5 = np.array([get_average_frequency(m, cl, df) for m, cl in zip(mach_points_5, cl_points_5)])

# Normalize weights so they sum to 1
weights_5 = frequencies_5 / np.sum(frequencies_5)

# Plot the 2D histogram with 5 integration points
plt.figure(figsize=(8, 6))
histogram_data = df.pivot(index="CL_bin", columns="Mach_bin", values="Frequency")
log_histogram_data = np.log1p(histogram_data)

# Set color range limits for consistency
vmin, vmax = log_histogram_data.min().min(), log_histogram_data.max().max()

plt.imshow(log_histogram_data, aspect="auto", origin="lower",
           extent=[df["Mach_bin"].min(), df["Mach_bin"].max(), df["CL_bin"].min(), df["CL_bin"].max()],
           cmap="viridis", vmin=vmin, vmax=vmax)

# Add the 5 integration points
plt.scatter(mach_points_5, cl_points_5, color='white', marker='o', edgecolor='black', s=100, label="Integration Points")

# Labels and title
plt.colorbar(label="Log(Frequency)")
plt.xlabel("Mach Number")
plt.ylabel("Lift Coefficient (C_L)")
plt.title("2D Histogram with 5 Integration Points (Log Scale)")
plt.legend()
plt.show()

# Display the numerical weights
weights_5
