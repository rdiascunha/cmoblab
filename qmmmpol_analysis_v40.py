import re
import os

def print_welcome_message():
    print("""
Welcome to the Excited States Python Analysis of Gaussian QM/MMPol calculations
Code developed by: Renato D. Cunha, PhD
Version: 4.0
""")

def print_usage():
    print("""
Usage:
    Run this script in the same directory as your Gaussian log files.
    The results will be automatically saved in a 'results' folder.
    Ensure that all required dependencies are installed.
""")

def extract_excited_states(log_lines):
    """Extracts excited state energies and oscillator strengths."""
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

def extract_transitions(log_lines):
    """Extracts transition contributions."""
    transitions = {}
    state = None
    
    for line in log_lines:
        match = re.search(r"Excited State\s+(\d+):", line)
        if match:
            state = int(match.group(1))
            transitions[state] = []
        elif "->" in line and state is not None:
            transitions[state].append(line.strip())
    
#    print(transitions)
    return transitions


def extract_orbital_info(log_lines):
    """Extracts atomic contributions to molecular orbitals."""
    orbitals = []
    capture_orbitals = False
    pattern = re.compile(r"Alpha (occ|vir) (\d+) OE=([-\d\.]+) is (.*)")
    
    for line in log_lines:
        if "Atomic contributions to Alpha molecular orbitals:" in line:
            capture_orbitals = True
            continue
        
        if capture_orbitals:
            if line.strip() == "":
                break  # Stop at empty line
            match = pattern.search(line)
            if match:
#                orbital_type = match.group(1)
                orbital_type = "HOMO" if match.group(1) == "occ" else "LUMO"
                orbital_number = int(match.group(2))
                orbital_energy = float(match.group(3))

                atom_contributions = match.group(4)  # Divide os átomos por vírgula
#                print(atom_contributions)
                cleaned_atoms = re.sub(r"-[^,]+", "", atom_contributions).strip()
#                print(cleaned_atoms)

                orbitals.append((orbital_type, orbital_number, orbital_energy, cleaned_atoms))
    
    return orbitals

def analyze_logs(directory_path=".", output_directory="results"):
    """Analyzes log files and generates required output files."""
    os.makedirs(output_directory, exist_ok=True)
    summary_data = []
    excited_states_data = []
    orbital_info_data = []
    
    for filename in sorted(os.listdir(directory_path)):
        if filename.endswith(".log"):
            log_file_path = os.path.join(directory_path, filename)
            try:
                with open(log_file_path, "r") as file:
                    log_lines = file.readlines()
            except FileNotFoundError:
                print(f"Error: Log file '{log_file_path}' not found.")
                continue
           
            # Extract frame number from filename
            match = re.search(r"_DIMER_(\d+)_ex_", filename)
            frame = match.group(1) if match else filename  # Use full filename if no match
#            print(frame)

            excited_states = extract_excited_states(log_lines)
            transitions = extract_transitions(log_lines)
            orbitals = extract_orbital_info(log_lines)
            
            for state, energy, osc_strength in excited_states:
                transition_lines = transitions.get(state, [])
                excited_states_data.append((frame, state, energy, osc_strength))
                for transition in transition_lines:
                    summary_data.append((frame, state, energy, osc_strength, transition))
                if not transition_lines:
                    summary_data.append((frame, energy, osc_strength, ""))
                
            for orbital_type, orbital_number, orbital_energy, cleaned_atoms in orbitals:
                orbital_info_data.append((frame, orbital_type, orbital_number, orbital_energy, cleaned_atoms))
    
    with open(os.path.join(output_directory, "summary.dat"), "w") as f:
        f.write("Frame\tState\tEnergy(eV)\tOsc.Str.\tTransitions\n")
        for frame, state, energy, osc_strength, transition in summary_data:
            f.write(f"{frame}\t{state}\t{energy:.3f}\t{osc_strength:.4f}\t{transition}\n")
    
    with open(os.path.join(output_directory, "excited_states.dat"), "w") as f:
        f.write("Frame\tState\tEnergy (eV)\tOsc. Str.\n")
        for frame, state, energy, osc_strength in excited_states_data:
            f.write(f"{frame}\t{state}\t{energy:.3f}\t{osc_strength:.4f}\n")
    
    with open(os.path.join(output_directory, "orbital.info"), "w") as f:
        f.write("Frame\tOrbital Num.\tType\tOrbital Energy (OE)\tRep. Atom Orbital\n")
        for frame, orbital_type,  orbital_number, orbital_energy, cleaned_atoms in orbital_info_data:
            f.write(f"{frame}\t{orbital_type}\t{orbital_number}\t{orbital_energy:.3f}\t{cleaned_atoms}\n")
    
    print(f"Analysis complete. Results saved in '{output_directory}'.")


def save_to_dat(filename, data, output_directory="results"):
    """Salva os dados em formato de tabela dentro da pasta 'results'."""
    os.makedirs(output_directory, exist_ok=True)
    filepath = os.path.join(output_directory, filename)
    
    try:
        with open(filepath, "w") as f:
            f.write("Frame\tState\tEnergy(eV)\tOsc.Str.\tTransitions\n")
            for entry in data:
                if len(entry) != 6:
                    print(f"Skipping invalid data entry: {entry}")
                    continue
                frame, frame, state, energy, osc_str, transitions = entry
                transition_str = " ".join(f"{t[0]}->{t[1]}" for t in transitions)
                f.write(f"{frame}\t{state}\t{energy:.4f}\t{osc_str:.4f}\t{transition_str}\n")
        print(f"File {filename} saved successfully in {output_directory}.")
    except Exception as e:
        print(f"Error saving {filename}: {e}")

def save_to_dat_2(filename, data, output_directory="results"):
    """Salva os dados em formato de tabela dentro da pasta 'results'."""
    os.makedirs(output_directory, exist_ok=True)
    filepath = os.path.join(output_directory, filename)
    
    try:
        with open(filepath, "w") as f:
            f.write("Frame\tState\tEnergy(eV)\tOsc.Str.\tTransitions\n")
            for entry in data:
                if len(entry) != 5:
                    print(f"Skipping invalid data entry: {entry}")
                    continue
                frame, state, energy, osc_str, transitions = entry
                transition_str = " ".join(f"{t[0]}->{t[1]}" for t in transitions)
                f.write(f"{frame}\t{state}\t{energy:.4f}\t{osc_str:.4f}\t{transition_str}\n")
        print(f"File {filename} saved successfully in {output_directory}.")
    except Exception as e:
        print(f"Error saving {filename}: {e}")

def parse_orbital_info(file_path):
    orbital_atoms = {}
    with open(file_path, "r") as f:
        next(f)  # Skip header
        for line in f:
            parts = line.split()
            frame, orbital_type, orbital_num, _, atom = parts
            frame = int(frame)
            orbital_num = int(orbital_num)
            atom_num = int(re.search(r'\d+', atom).group())
            if frame not in orbital_atoms:
                orbital_atoms[frame] = {}
            orbital_atoms[frame][orbital_num] = atom_num
#    print(orbital_atoms)
    return orbital_atoms

def parse_summary(file_path):
    excited_states = {}
    with open(file_path, "r") as f:
        next(f)  # Skip header
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue  # Skip malformed lines
            frame, state, energy, osc_str, *transitions = parts
            frame = int(frame)
            state = int(state)
            energy = float(energy)
            osc_str = float(osc_str)

            parsed_transitions = []
            for i in range(len(transitions) - 2):  # Ensuring we don't pick coefficients
                if transitions[i+1] == '->':  # Look for correct transition format
                    try:
                        donor = int(transitions[i])
                        acceptor = int(transitions[i+2])
                        parsed_transitions.append((donor, acceptor))
                    except ValueError:
                        continue  # Skip any malformed data

            if frame not in excited_states:
                excited_states[frame] = []
            excited_states[frame].append((frame, state, energy, osc_str, parsed_transitions))
#    print(excited_states)
    return excited_states

def analyze_states(orbital_file, summary_file):
    orbital_data = parse_orbital_info(orbital_file)
    excited_data = parse_summary(summary_file)
    
    chromo1_Qx, chromo1_Qy = [], []
    chromo2_Qx, chromo2_Qy = [], []
    CT_states = []
    
    for frame, states in excited_data.items():
        chromo1_states = []
        chromo2_states = []
        
        for frame, state, energy, osc_str, parsed_transitions in states:
            if energy > 3.500 or not parsed_transitions:
                continue
#            print(orbital_data)
            same_chromophore = all(
                (orbital_data[frame].get(t[0], 999) < 100 and orbital_data[frame].get(t[1], 999) < 100) or
                (orbital_data[frame].get(t[0], 999) > 100 and orbital_data[frame].get(t[1], 999) > 100)
                for t in parsed_transitions
            )
            
            if same_chromophore:
                if orbital_data[frame].get(parsed_transitions[0][0], 0) < 100:
                    chromo1_states.append((frame, state, energy, osc_str, parsed_transitions))
                else:
                    chromo2_states.append((frame, state, energy, osc_str, parsed_transitions))
            
            elif osc_str < 0.01:
                CT_states.append((frame, state, energy, osc_str, parsed_transitions))
       
            for t in parsed_transitions:
                atom1 = orbital_data[frame].get(t[0], 999)
                atom2 = orbital_data[frame].get(t[1], 999)
#                print(f"Frame {frame}, State {state}: Transition {t[0]} -> {t[1]} | Atoms: {atom1} -> {atom2} | Same Chromophore: {same_chromophore}")


        if chromo1_states:
            # Filtrar apenas estados com energia < 2.500 antes de ordenar
            chromo1_states = [s for s in chromo1_states if s[2] < 2.500]
    
            if chromo1_states:  # Certifique-se de que ainda há estados na lista
                chromo1_states.sort(key=lambda x: x[2])  # Ordenar por energia
                chromo1_Qy.append((frame, *chromo1_states[0]))  # O menor em energia é Qy

                # Procurar Qy com energia próxima de Qx, mas distinta
                Qy_energy = chromo1_states[0][2]
                for state in chromo1_states[1:]:
                    if abs(state[2] - Qy_energy) > 0.05:  # Diferença mínima para evitar estados idênticos
                        chromo1_Qx.append((frame, *state))
                        break  # Para ao encontrar o Qy


        if chromo2_states:
            chromo2_states = [s for s in chromo2_states if s[2] < 2.500]

            if chromo2_states:
                chromo2_states.sort(key=lambda x: x[2])
                chromo2_Qy.append((frame, *chromo2_states[0]))

                Qy_energy = chromo2_states[0][2]
                for state in chromo2_states[1:]:
                    if abs(state[2] - Qy_energy) > 0.05:
                        chromo2_Qx.append((frame, *state))
                        break

#    print(energy, osc_str, transitions)
#    print("Chromo1 Qx:", chromo1_Qx)
#    print("Chromo1 Qy:", chromo1_Qy)
#    print("Chromo2 Qx:", chromo2_Qx)
#    print("Chromo2 Qy:", chromo2_Qy)
#    print("CT States:", CT_states)

    save_to_dat("chromo1_Qx.dat", chromo1_Qx)
    save_to_dat("chromo1_Qy.dat", chromo1_Qy)
    save_to_dat("chromo2_Qx.dat", chromo2_Qx)
    save_to_dat("chromo2_Qy.dat", chromo2_Qy)
    save_to_dat_2("CT_states.dat", CT_states)

if __name__ == "__main__":
    print_welcome_message()
    print_usage()
    
    output_directory = "results"
    analyze_logs(output_directory=output_directory)
    analyze_states(
        os.path.join(output_directory, "orbital.info"),
        os.path.join(output_directory, "summary.dat")
    )

