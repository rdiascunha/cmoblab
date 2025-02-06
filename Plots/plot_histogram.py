import argparse
import numpy as np
import matplotlib.pyplot as plt

def plot_histograms(files, legend_labels):
    if len(files) != len(legend_labels):
        print("Error: Number of files and legend labels must be the same.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the figsize parameter as needed

    for i, file in enumerate(files):
        # Load data from the input file
        data = np.loadtxt(file)

        # Calculate bin width
        bin_width = data[1, 0] - data[0, 0]

        # Plot the histogram curve with a different color and legend label for each file
        plt.plot(data[:, 0], data[:, 1], label=legend_labels[i], linewidth=2.5)

    # Add the experimental indication
    plt.axvline(x=30, color='red', linestyle='--', linewidth=1.5)

    # Adjust legend size and font weight
    leg = plt.legend(loc='upper right', prop={'size': 14, 'weight': 'bold'})
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)  # Adjust the legend line width

    # Adjust axis font weight
    plt.xlabel('Radius of Gyration (Å)', fontweight='bold', fontsize=14)
    plt.ylabel('Normalized histogram', fontweight='bold', fontsize=12)
    plt.title('', fontweight='bold')
    
    # Remove 'Y' numbers
    plt.yticks([])

    # Adjust tick label font weight
    plt.xticks(fontweight='bold')
    plt.yticks(fontweight='bold')

    plt.grid(False)

    # Save and show the plot
    plt.savefig('histogram_complete.png', dpi=300, format='png')
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Plot multiple histogram curves.')
    parser.add_argument('-f', '--files', nargs='+', required=True, help='Input histogram data files')
    parser.add_argument('-l', '--legend-labels', nargs='+', required=True, help='Legend labels for each curve')
#    parser.add_argument('--line-width', type=float, default=1.5, help='Width of the dotted horizontal line')

    args = parser.parse_args()

    plot_histograms(args.files, args.legend_labels)

if __name__ == "__main__":
    main()

