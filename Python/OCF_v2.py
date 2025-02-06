import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

def orientational_correlation(unit_vectors):
    num_vectors = len(unit_vectors)
    correlation = np.zeros(num_vectors - 1)
    
    for s in range(1, num_vectors):
        dot_products = np.array([np.dot(unit_vectors[i], unit_vectors[i + s]) for i in range(num_vectors - s)])
        correlation[s - 1] = np.mean(dot_products)
    
    return correlation

def exponential_decay(s, k):
    return np.exp(-s / k)

def process_pdb(file_path):
    # Load the trajectory
    traj = md.load(file_path)
    
    # Extract the indices of the Cα atoms
    ca_indices = traj.topology.select('name CA')
    
    # Extract the coordinates of the Cα atoms
    ca_coords = traj.atom_slice(ca_indices).xyz[0]
    
    # Calculate unit vectors connecting consecutive Cα atoms
    vectors = ca_coords[1:] - ca_coords[:-1]
    unit_vectors = vectors / np.linalg.norm(vectors, axis=1)[:, None]
    
    # Calculate the orientational correlation function C(s)
    correlation = orientational_correlation(unit_vectors)
    
    # Fit the correlation data to the exponential decay function to estimate k
    s_values = np.arange(1, len(correlation) + 1)
    popt, _ = curve_fit(exponential_decay, s_values, correlation)
    
    # Extract the fitted value of k
    k = popt[0]
    
    # Calculate the persistence length l_p
    lp = k * 3.8  # Angstroms
    
    return k, lp, correlation

def main():
    # Get the list of PDB files in the current directory and sort them
    pdb_files = sorted([f for f in os.listdir('.') if f.endswith('.pdb')])
    
    # Initialize variables for averaging correlation functions
    all_correlations = []
    max_length = 0
    
    # Open the log file
    with open('OCF.log', 'w') as log_file:
        log_file.write("Filename\tk\tl_p\n")
        
        # Process each PDB file
        for pdb_file in pdb_files:
            try:
                k, lp, correlation = process_pdb(pdb_file)
                log_file.write(f"{pdb_file}\t{k:.2f}\t{lp:.2f}\n")
                print(f"Processed {pdb_file}: k = {k:.2f}, l_p = {lp:.2f} Å")
                
                # Store the correlation for averaging
                all_correlations.append(correlation)
                max_length = max(max_length, len(correlation))
            except Exception as e:
                print(f"Error processing {pdb_file}: {e}")
    
    # Calculate the average correlation function
    avg_correlation = np.zeros(max_length)
    count = np.zeros(max_length)
    
    for correlation in all_correlations:
        avg_correlation[:len(correlation)] += correlation
        count[:len(correlation)] += 1
    
    avg_correlation /= count
    
    # Define s_values outside the try block
    s_values = np.arange(1, max_length + 1)
    
    # Save the average correlation data to a file
    np.savetxt('graphic_OCF.dat', np.column_stack((s_values, avg_correlation)), header='s_values avg_correlation', comments='')
    
    # Fit the average correlation function to the exponential decay model
    try:
        popt, _ = curve_fit(exponential_decay, s_values, avg_correlation)
        k_avg = popt[0]
        
        # Set global font properties
        plt.rcParams.update({'font.size': 12, 'font.weight': 'bold'})
        
        # Plot the average orientational correlation function and the fitted curve
        plt.figure(figsize=(10, 5))
        plt.plot(s_values, avg_correlation, label='Average Orientational Correlation Function', linewidth=2)
        plt.plot(s_values, exponential_decay(s_values, k_avg), label=f'Exponential Fit (k = {k_avg:.2f})', linewidth=2)
        
        # Define font properties
        font = {'size': 14, 'weight': 'bold'}
        
        plt.xlabel('Pairwise Residue Separation (s)', fontdict=font)
        plt.ylabel('Orientational Correlation (C(s))', fontdict=font)
#        plt.title(, fontdict={'size': 16, 'weight': 'bold'})
        plt.legend(fontsize=12)
        plt.xticks(fontsize=12, weight='bold')
        plt.yticks(fontsize=12, weight='bold')
        plt.gca().spines['top'].set_linewidth(2)
        plt.gca().spines['right'].set_linewidth(2)
        plt.gca().spines['left'].set_linewidth(2)
        plt.gca().spines['bottom'].set_linewidth(2)
        plt.show()
    except Exception as e:
        print(f"Error fitting and plotting average correlation: {e}")

if __name__ == "__main__":
    main()

