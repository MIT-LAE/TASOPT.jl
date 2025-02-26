import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_histogram_data(file_path):
    """Load the CSV file and return the dataframe."""
    return pd.read_csv(file_path)


def plot_2d_histogram(df, mach_min, mach_max, cl_min, cl_max, mach_points_5, cl_points_5):
    """Plot the 2D histogram with integration points using exact frequency values."""
    
    # Use Seaborn rocket_r colormap
    cmap = sns.color_palette("rocket_r", as_cmap=True)
    
    MMo = 0.82
    
    fig, ax1 = plt.subplots()
    histogram_data = df.pivot(index="CL_bin", columns="Mach_bin", values="Frequency")
    mach_bins = df["Mach_bin"].unique()
    cl_bins = df["CL_bin"].unique()
    
    pcm = ax1.pcolor(mach_bins, cl_bins, histogram_data, cmap=cmap, shading='auto', vmin = 0, vmax = 200, alpha = 0.8)
    
    # Configure the colorbar
    cbar = fig.colorbar(pcm, ax=ax1, label='Frequency')
    
    # Define linear ticks for the colorbar
    ticks = [0, 25, 50, 75, 100, 125, 150, 175, 200]  # Linear tick positions
    cbar.set_ticks(ticks)
    
    # Customize the colorbar appearance
    cbar.ax.tick_params(labelsize=20, labelrotation=0, labelcolor='black')
    cbar.set_label('Frequency', fontsize=20, fontname="Times New Roman")
    
    # Ensure tick labels use the desired font
    for label in cbar.ax.get_yticklabels():
        label.set_fontname("Times New Roman")
        
    # Add the bounding box
    bbox_x = [mach_min, mach_max, mach_max, mach_min, mach_min]
    bbox_y = [cl_min, cl_min, cl_max, cl_max, cl_min]
    plt.plot(bbox_x, bbox_y, color='red', linewidth=3, linestyle='dashed')
    
    # Add the 5 integration points
    plt.scatter(mach_points_5, cl_points_5, color='C0', marker='o', edgecolor='black', s=100, label="Integration Points")
    
    # Add a vertical dashed line at Mach = 0.80 for all Reynolds numbers
    ax1.axvline(x=MMo, color='C0', linewidth=2, linestyle=':')    

    # Add labels and title
    ax1.set_xlabel('Mach', fontsize=22, fontname="Times New Roman")
    ax1.set_ylabel(r'C$_L$', fontsize=22, fontname="Times New Roman")
    
    # Adjust axis ticks and spines
    ax1.tick_params(bottom=True, top=False, left=True, right=True)
    ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
    ax1.tick_params(which='major', length=10, width=1.2, direction='in')
    ax1.tick_params(which='minor', length=5, width=1.2, direction='in')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(1.5)
    
    plt.xticks(fontname="Times New Roman", fontsize=20)
    plt.yticks(fontname="Times New Roman", fontsize=20)
    
    F = plt.gcf()
    Size = F.get_size_inches()
    F.set_size_inches(Size[0] * 1.5, Size[1] * 1.5, forward=True)
    
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    
    ax1.set_xlim(0.59, 0.90)
    ax1.set_ylim(0.44, 0.58)
    plt.tight_layout()

    plt.show()


def compute_frequencies_and_weights(mach_points, cl_points, df, mach_bin_size, cl_bin_size):
    """Compute average frequencies from surrounding bins and normalize weights. Ensure sum of weights is 1."""
    frequencies = []
    surrounding_bins_list = []
    for mach, cl in zip(mach_points, cl_points):
        surrounding_bins = df[(df["Mach_bin"] >= mach - 4 * mach_bin_size) & (df["Mach_bin"] <= mach + 4 * mach_bin_size) &
                              (df["CL_bin"] >= cl - 4 * cl_bin_size) & (df["CL_bin"] <= cl + 4 * cl_bin_size)]
        surrounding_bins_list.append(surrounding_bins)
        avg_freq = surrounding_bins["Frequency"].mean() if not surrounding_bins.empty else 0
        frequencies.append(avg_freq)
    frequencies = np.array(frequencies)
    weights = frequencies / np.sum(frequencies)
    
    # Check if weights sum to 1
    if not np.isclose(np.sum(weights), 1.0):
        raise ValueError("Sum of weights does not equal 1. Check frequency data and normalization.")
    
    return mach_points, cl_points, weights, surrounding_bins_list


def get_bin_sizes(df):
    """Return the bin size for Mach and CL from the CSV file."""
    mach_bins = np.sort(df["Mach_bin"].unique())
    cl_bins = np.sort(df["CL_bin"].unique())
    mach_bin_size = np.min(np.diff(mach_bins))
    cl_bin_size = np.min(np.diff(cl_bins))
    return mach_bin_size, cl_bin_size



# Load the dataset
file_path = "data/E170_histogram_data.csv"
df = load_histogram_data(file_path)

# Get bin sizes
mach_bin_size, cl_bin_size = get_bin_sizes(df)

# Define the bounding box limits
mach_min, mach_max = 0.70, 0.82
cl_min, cl_max = 0.48, 0.54

# Compute the center of the bounding box
mach_center = (mach_min + mach_max) / 2
cl_center = (cl_min + cl_max) / 2

# Define 5 integration points (center + 4 surrounding points)
mach_points_5 = np.array([mach_min, mach_max, mach_center, mach_min, mach_max])
cl_points_5 = np.array([cl_min, cl_min, cl_center, cl_max, cl_max])

# Plot the 2D histogram with integration points
plot_2d_histogram(df, mach_min, mach_max, cl_min, cl_max, mach_points_5, cl_points_5)

# Compute frequencies and weights using bin sizes
mach_points_5, cl_points_5, weights_5 = compute_frequencies_and_weights(mach_points_5, cl_points_5, df, mach_bin_size, cl_bin_size)

# Display the numerical weights in tabular format
print("\nNumerical Weights for Integration Points:")
print(pd.DataFrame({"Mach": mach_points_5, "CL": cl_points_5, "Weight": weights_5}))
