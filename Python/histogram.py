import argparse
import numpy as np

def histogram(data, min_val, max_val, num_bins, output_file):
    # Calculate bin edges
    bin_edges = np.linspace(min_val, max_val, num_bins + 1)

    # Bin the data
    hist, _ = np.histogram(data, bins=bin_edges)

    # Print results to screen
    print("bin(centre) = freq")
    for i in range(len(hist)):
        bin_center = (bin_edges[i] + bin_edges[i + 1]) / 2
        print(f"bin({bin_center:.2f}) = {hist[i]}")

    # Save results to a file for plotting
    np.savetxt(output_file, np.column_stack((bin_edges[:-1], hist)), fmt='%.5f %.5f')

def main():
    parser = argparse.ArgumentParser(description='Generate histogram from data file.')
    parser.add_argument('-f', '--input-file', required=True, help='Input data file name')
    parser.add_argument('-o', '--output-file', required=True, help='Output file name for histogram data')
    parser.add_argument('--min', type=float, default=0, help='Minimum value of x axis')
    parser.add_argument('--max', type=float, default=10, help='Maximum value of x axis')
    parser.add_argument('-b', '--num-bins', type=int, default=5, help='Number of bins')

    args = parser.parse_args()

    # Load data from the input file
    data = np.loadtxt(args.input_file)

    histogram(data, args.min, args.max, args.num_bins, args.output_file)

if __name__ == "__main__":
    main()
