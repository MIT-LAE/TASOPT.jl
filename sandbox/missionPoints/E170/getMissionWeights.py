import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from numpy.polynomial.hermite import hermgauss

def compute_statistics(df):
    """Compute and print the mean and variance of Mach and CL bins."""
    mean_mach = np.mean(df["Mach_bin"])
    var_mach = np.var(df["Mach_bin"])
    mean_cl = np.mean(df["CL_bin"])
    var_cl = np.var(df["CL_bin"])
    
    print(f"Mean Mach: {mean_mach}, Variance Mach: {var_mach}")
    print(f"Mean CL: {mean_cl}, Variance CL: {var_cl}")
    
    return mean_mach, var_mach, mean_cl, var_cl


def load_histogram_data(file_path):
    """Load the CSV file and return the dataframe."""
    return pd.read_csv(file_path)


def plot_2d_histogram(df, mach_min, mach_max, cl_min, cl_max, mach_points_5, cl_points_5, surrounding_bins_list, weights_5,plot_box, plot_points):
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
    
    if plot_box and plot_points == False:     
    # Add the bounding box
        bbox_x = [mach_min, mach_max, mach_max, mach_min, mach_min]
        bbox_y = [cl_min, cl_min, cl_max, cl_max, cl_min]
        plt.plot(bbox_x, bbox_y, color='gray', linewidth=3, linestyle='dashed')
    
    # Add the 5 integration points
    #plt.scatter(mach_points_5, cl_points_5, color='C0', marker='o', edgecolor='black', s=100, label="Integration Points")
    
    # Add a vertical dashed line at Mach = 0.80 for all Reynolds numbers
    ax1.axvline(x=MMo, color='black', linewidth=2.5, linestyle=':')
    
    
    if plot_points and plot_box == False:
        
        # Add the bounding box
        bbox_x = [mach_min, mach_max, mach_max, mach_min, mach_min]
        bbox_y = [cl_min, cl_min, cl_max, cl_max, cl_min]
        plt.plot(bbox_x, bbox_y, color='gray', linewidth=2, linestyle='dashed')
        
       # Highlight surrounding bins in gray
       # for surrounding_bins in surrounding_bins_list:
        #    plt.scatter(surrounding_bins["Mach_bin"], surrounding_bins["CL_bin"], color='gray', marker='s', s=50, alpha=0.5)
       
        # Add the 5 integration points with size proportional to weight
        integration_point_sizes = weights_5 * 1000  # Scale weights for better visibility
        plt.scatter(mach_points_5, cl_points_5, color='C0', marker='P', edgecolor='black', s=integration_point_sizes)    

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

    if plot_points:
        plt.savefig("Plots/E170_Mach_CL_hist_weights.png")
    if plot_box:
        plt.savefig("Plots/E170_Mach_CL_hist_box.png")
    else:
        plt.savefig("Plots/E170_Mach_CL_hist.png")

def compute_frequencies_and_weights(mach_points, cl_points, df, mach_bin_size, cl_bin_size, filter_radius):
    """Compute average frequencies from surrounding bins and normalize weights. Ensure sum of weights is 1."""
    frequencies = []
    surrounding_bins_list = []
    
    total_frequency = df["Frequency"].sum()  # Total frequency in the histogram
    
    for mach, cl in zip(mach_points, cl_points):
        surrounding_bins = df[(df["Mach_bin"] >= mach - filter_radius * mach_bin_size) & (df["Mach_bin"] <= mach + filter_radius * mach_bin_size) &
                              (df["CL_bin"] >= cl - filter_radius * cl_bin_size) & (df["CL_bin"] <= cl + filter_radius * cl_bin_size)]
        surrounding_bins_list.append(surrounding_bins)
        #avg_freq = surrounding_bins["Frequency"].mean() if not surrounding_bins.empty else 0
        relative_freq = surrounding_bins["Frequency"].sum() / total_frequency if not surrounding_bins.empty else 0
        frequencies.append(relative_freq)
        #frequencies.append(avg_freq)
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


def plot_bivariate_normal_distribution(mu_M, sigma_M, mu_CL, sigma_CL):
    """Plot the filled contour of the Bivariate Normal distribution using mean and variance."""
    plt.figure(figsize=(8, 6))
    
    # Generate a grid for contour plot
    mach_vals = np.linspace(mu_M - 3*sigma_M, mu_M + 3*sigma_M, 100)
    cl_vals = np.linspace(mu_CL - 3*sigma_CL, mu_CL + 3*sigma_CL, 100)
    M_grid, CL_grid = np.meshgrid(mach_vals, cl_vals)
    
    # Compute the bivariate normal distribution
    rv = stats.multivariate_normal(mean=[mu_M, mu_CL], cov=[[sigma_M**2, 0], [0, sigma_CL**2]])
    Z = rv.pdf(np.dstack((M_grid, CL_grid)))
    
    # Overlay filled contour plot
    plt.contourf(M_grid, CL_grid, Z, levels=20, cmap='coolwarm', alpha=0.7)
    plt.colorbar(label="Probability Density")
    
    # Labels and title
    plt.xlabel("Mach Number")
    plt.ylabel("Lift Coefficient (C_L)")
    plt.title("Bivariate Normal Distribution of Mach-CL")
    plt.show()


# Load the dataset
file_path = "data/E170_histogram_data.csv"
df = load_histogram_data(file_path)

# Get bin sizes
mach_bin_size, cl_bin_size = get_bin_sizes(df)

# Define the bounding box limits
mach_min, mach_max = 0.70, 0.82
cl_min, cl_max = 0.48, 0.54
filter_radius = 4

# Compute the center of the bounding box
mach_center = (mach_min + mach_max) / 2
cl_center = (cl_min + cl_max) / 2

# Define 5 integration points (center + 4 surrounding points)
mach_points_5 = np.array([mach_min, mach_max, mach_center, mach_min, mach_max])
cl_points_5 = np.array([cl_min, cl_min, cl_center, cl_max, cl_max])

# Compute frequencies and weights using bin sizes
mach_points_5, cl_points_5, weights_5 , surrounding_bins_list = compute_frequencies_and_weights(mach_points_5, cl_points_5, df, mach_bin_size, cl_bin_size, filter_radius)

plot_box = False
plot_points = True


# Plot the 2D histogram with integration points
plot_2d_histogram(df, mach_min, mach_max, cl_min, cl_max, mach_points_5, cl_points_5, surrounding_bins_list, weights_5, plot_box, plot_points)


# Display the numerical weights in tabular format
print("\nNumerical Weights for Integration Points:")
print(pd.DataFrame({"Mach": mach_points_5, "CL": cl_points_5, "Weight": np.round(weights_5, 2)}))


# Compute and print statistics
#mean_mach, var_mach, mean_cl, var_cl = compute_statistics(df)

# Plot the Bivariate Normal Distribution
#plot_bivariate_normal_distribution(mean_mach, np.sqrt(var_mach), mean_cl, np.sqrt(var_cl))