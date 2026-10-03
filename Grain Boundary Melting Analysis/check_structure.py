import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_dat_file(filename, delimiter=None, x_label='X', y_label='Y', 
                  title='Plot from .dat file', marker='o'):
    """
    Read a .dat file with x and y columns and create a scatter plot.
    The x and y axes have equal unit lengths.
    """
    try:
        data = np.loadtxt(filename, delimiter=delimiter, comments='#')
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if data.ndim == 1 or data.shape[1] < 2:
        print("File must contain at least two columns (x and y).")
        return

    x = data[:, 0]
    y = data[:, 1]

    # Create the figure and axes
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Scatter plot (points only)
    ax.scatter(x, y, marker=marker, color='b', label='Data')
    
    # Set equal scaling: 1 unit in x = 1 unit in y
    ax.set_aspect('equal')
    
    # Optional: If you want the plot box to be exactly square based on data range,
    # you can uncomment the next line:
    # ax.axis('equal')  # This also adjusts the limits to make the box square
    
    # Labels, title, grid, legend
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()
    
    # Improve spacing and display
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = input("Enter the path to the .dat file: ")
    plot_dat_file(filename)