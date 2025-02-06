import argparse
import numpy as np
import matplotlib.pyplot as plt

def normalize_histogram(input_file, output_file):
    # Load data from the input file
    data = np.loadtxt(input_file)

    # Calculate bin width
    bin_width = data[1, 0] - data[0, 0]

    # Normalize the histogram by area
    normalized_hist = data[:, 1] / (np.sum(data[:, 1]) * bin_width)

    # Save normalized histogram to the output file
    np.savetxt(output_file, np.column_stack((data[:, 0], normalized_hist)), fmt='%.5f %.5f')

    # Plot the original and normalized histograms
    plt.bar(data[:, 0], data[:, 1], width=bin_width, alpha=0.5, label='Original Histogram')
    plt.bar(data[:, 0], normalized_hist, width=bin_width, alpha=0.5, label='Normalized Histogram')
    plt.xlabel('bin center')
    plt.ylabel('frequency/area')
    plt.title('Original and Normalized Histograms')
    plt.legend()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Normalize histogram data.')
    parser.add_argument('-i', '--input-file', required=True, help='Input histogram data file name')
    parser.add_argument('-o', '--output-file', required=True, help='Output file name for normalized histogram data')

    args = parser.parse_args()

    normalize_histogram(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

