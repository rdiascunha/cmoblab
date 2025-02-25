import re
import os
from math import sqrt

def print_welcome_message():
    print("""
Welcome to the Excited States Python Analysis of Gaussian QM/MMPol calculations
Code developed by: Renato D. Cunha, PhD
Version: 5.1
""")

def print_usage():
    print("""
Usage:
    Run this script in the same directory as your Gaussian log files.
    The results will be automatically saved in a 'results' folder.
    Ensure that all required dependencies are installed.
""")

def extract_excited_states(log_lines):
    excited_states = []
    pattern = re.compile(r"Excited State\s+(\d+):\s+.*?([\d\.]+) eV.*?f=([\d\.]+)")
    for line in log_lines:
        match = pattern.search(line)
        if match:
            state = int(match.group(1))
            energy = float(match.group(2))
            osc_strength = float(match.group(3))
            excited_states.append((state, energy, osc_strength))
    return excited_states

def extract_transition_dipoles(log_lines):
    dipoles = {}
    capture = False
    
    for line in log_lines:
        if "Ground to excited state transition electric dipole moments (Au):" in line:
            capture = True
            continue
        if capture:
            if line.strip() == "" or "velocity" in line:
                break
            parts = line.split()
            if len(parts) >= 6 and parts[0].isdigit():
                state = int(parts[0])
                dipoles[state] = {
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'z': float(parts[3]),
                    'dip': sqrt(float(parts[4]))
                }
    
    return dipoles

def extract_orbital_info(log_lines):
    orbitals = []
    capture_orbitals = False
    pattern = re.compile(r"Alpha (occ|vir) (\d+) OE=([-\d\.]+) is (.*)")
    orbital_count = {'occ': 0, 'vir': 0}

    for line in log_lines:
        if "Atomic contributions to Alpha molecular orbitals:" in line:
            capture_orbitals = True
            continue

        if capture_orbitals:
            if line.strip() == "":
                break
            match = pattern.search(line)
            if match:
                orbital_type = match.group(1)
                orbital_number = int(match.group(2))
                orbital_energy = float(match.group(3))
                atom_contributions = match.group(4)
                cleaned_atoms = re.sub(r"-[^,]+", "", atom_contributions).strip()

                atom_num = int(re.search(r'\d+', cleaned_atoms).group())
                chromophore = "1" if atom_num < 100 else "2"

                if orbital_type == 'occ':
                    orbital_count['occ'] += 1
                    name = "HOMO-1" if orbital_count['occ'] <= 2 else "HOMO"
                else:
                    orbital_count['vir'] += 1
                    name = "LUMO" if orbital_count['vir'] <= 2 else "LUMO+1"

                orbitals.append((name, chromophore, orbital_number, orbital_energy, cleaned_atoms))

    return sorted(orbitals, key=lambda x: x[3])

def analyze_logs(directory_path=".", output_directory="results"):
    os.makedirs(output_directory, exist_ok=True)
    excited_states_data = []
    state_files = {}
    
    with open(os.path.join(output_directory, "orbital.info"), "w") as f:
        f.write("Frame\tOrbital Name\tChromophore\tOrbital Number\tOrbital Energy\tAtoms\n")

    states_list = ['CT12', 'CT21', 'Qx_1', 'Qx_2', 'Qy_1', 'Qy_2']
    for state in states_list:
        file_path = os.path.join(output_directory, f"chromo{state.lower()}.dat")
        state_files[state] = open(file_path, "w")
        state_files[state].write("Frame\tState\tEnergy(eV)\tOsc.Str.\tDipX\tDipY\tDipZ\tDipTot\n")
    
    for filename in sorted(os.listdir(directory_path)):
        if filename.endswith(".log"):
            with open(os.path.join(directory_path, filename), "r") as file:
                log_lines = file.readlines()
            
            match = re.search(r"_DIMER_(\d+)_ex_", filename)
            frame = match.group(1) if match else filename
            
            excited_states = extract_excited_states(log_lines)
            dipoles = extract_transition_dipoles(log_lines)
            orbital_info = extract_orbital_info(log_lines)
            
            for state, energy, osc_strength in excited_states:
                dip = dipoles.get(state, {'x': 0.0, 'y': 0.0, 'z': 0.0, 'dip': 0.0})
                excited_states_data.append((frame, state, energy, osc_strength, 
                                         dip['x'], dip['y'], dip['z'], dip['dip']))
                
                if osc_strength < 0.01:
                    state_type = "CT12" if state % 2 == 0 else "CT21"
                else:
                    if state % 4 == 0:
                        state_type = "Qx_1"
                    elif state % 4 == 1:
                        state_type = "Qx_2"
                    elif state % 4 == 2:
                        state_type = "Qy_1"
                    else:
                        state_type = "Qy_2"
                
                if state_type in state_files:
                    state_files[state_type].write(f"{frame}\t{state}\t{energy:.4f}\t{osc_strength:.4f}\t"
                        f"{dip['x']:.4f}\t{dip['y']:.4f}\t{dip['z']:.4f}\t{dip['dip']:.4f}\n")
   

            with open(os.path.join(output_directory, "orbital.info"), "a") as f:
                for orb in orbital_info:
                    f.write(f"{frame}\t{orb[0]}\t{orb[1]}\t{orb[2]}\t{orb[3]:.3f}\t{orb[4]}\n")


    with open(os.path.join(output_directory, "excited_states.dat"), "w") as f:
        f.write("Frame\tState\tEnergy(eV)\tOsc.Str.\tDipX\tDipY\tDipZ\tDipTot\n")
        for data in excited_states_data:
            f.write(f"{data[0]}\t{data[1]}\t{data[2]:.3f}\t{data[3]:.4f}\t"
                   f"{data[4]:.4f}\t{data[5]:.4f}\t{data[6]:.4f}\t{data[7]:.4f}\n")
    
    for f in state_files.values():
        f.close()

if __name__ == "__main__":
    print_welcome_message()
    print_usage()
    analyze_logs()

