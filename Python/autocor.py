import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def read_data(file):
    """Reads data from a file and returns frames and values."""
    data = np.loadtxt(file)
    frames = data[:, 0]
    values = data[:, 1]
    return frames, values

def calculate_autocorrelation(data):
    """Calculates the normalized autocorrelation of the data."""
    n = len(data)
    mean = np.mean(data)
    variance = np.var(data)
    autocorr = np.correlate(data - mean, data - mean, mode='full') / (n * variance)
    return autocorr[n - 1:]  # Return the second half (positive lags)

def save_autocorrelation(time, autocorrelation, output_file):
    """Saves the autocorrelation values to a .dat file."""
    with open(output_file, 'w') as f:
        f.write("# Time (ns)   Autocorrelation\n")
        for t, ac in zip(time, autocorrelation):
            f.write(f"{t:.6f}   {ac:.6f}\n")

def plot_autocorrelation(time, autocorrelations, labels, output_file, decay_factor=1.0):
    """Plots the autocorrelation along with a control exponential decay."""
    plt.figure(figsize=(10, 6))

    # Plot control exponential decay
    perfect_decay = np.exp(-decay_factor * time)
    plt.plot(time, perfect_decay, 'r-', label='$y = e^{-\lambda t}$', linewidth=2)

    # Colors for each file
    colors = plt.cm.viridis(np.linspace(0, 1, len(autocorrelations)))

    # Plot autocorrelation curves with different colors
    for autocorr, label, color in zip(autocorrelations, labels, colors):
        plt.plot(time, autocorr, label=label, linewidth=2.0, color=color)

    # Customize plot appearance
    plt.xlabel('Time (ns)', fontsize=14, fontweight='bold')
    plt.ylabel('Correlation', fontsize=14, fontweight='bold')
    plt.legend(fontsize=14, loc='upper right', frameon=False)
    plt.xlim(0, time[-1])  # Set x-axis limits from 0 to the last time value

    # Make axes and borders bold
    ax = plt.gca()
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(axis='both', which='major', labelsize=14, width=2)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Plot autocorrelation of data from multiple files.")
    parser.add_argument('-f', '--files', nargs='+', required=True, help="Input data files.")
    parser.add_argument('-o', '--output', default='autocorrelation_plot.png', help="Output file name for the plot.")
    parser.add_argument('--frame_to_ns', type=float, default=5, help="Number of frames corresponding to 1 ns.")
    parser.add_argument('--decay_factor', type=float, default=0.001, help="Decay factor for the ideal exponential function.")
    args = parser.parse_args()

    autocorrelations = []
    labels = []
    time = None

    for idx, file in enumerate(args.files):
        # Read data
        frames, values = read_data(file)
        labels.append(f"Sim-{idx + 1:02}")

        # Convert frames to time in ns
        if time is None:
            time = frames / args.frame_to_ns

        # Calculate autocorrelation
        autocorr = calculate_autocorrelation(values)
        autocorrelations.append(autocorr[:len(time)])  # Truncate to match time length

        # Save autocorrelation to a .dat file
        output_dat_file = os.path.splitext(file)[0] + "_autocorrelation.dat"
        save_autocorrelation(time, autocorr[:len(time)], output_dat_file)

    # Plot results
    plot_autocorrelation(time, autocorrelations, labels, args.output, args.decay_factor)

if __name__ == '__main__':
    main()

