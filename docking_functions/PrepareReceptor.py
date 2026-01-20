'''
# RECEPTOR PDBQT PREPARATION
=================================
This module prepares protein receptors for molecular docking by converting various file 
formats to PDBQT format with customizable OpenBabel parameters.

## Main processes:

1. Format Conversion:
   - Converts multiple file formats (ENT, XYZ, PQR, MCIF, MMCIF, MOL2) to PDB
   - Converts cleaned PDB to PDBQT format with configurable docking parameters

2. Chain Selection and Deduplication:
   - Extracts specific protein chains from multi-protein PDB files
   - Identifies and removes duplicate/redundant chains based on sequence similarity
   - Filters chains by protein name from PDB metadata

3. Structure Cleaning:
   - Removes non-protein residues and heteroatoms
   - Eliminates unknown amino acids (UNK, UNL) from both PDB and PDBQT files
   - Preserves structural integrity through careful residue filtering

## Configurable OpenBabel Parameters:

   - pH: Protonation pH value (default 7.4)
   - Charge Method: Partial charge calculation ('eem', 'gasteiger', 'mmff94')
   - Polar Hydrogens: Remove polar hydrogens for reduced complexity (default True)

## Output:

   - Cleaned PDB file (protein-only structure)
   - PDBQT file (ready for docking simulations)
'''

#  _     _ _                    _           
# | |   (_) |__  _ __ __ _ _ __(_) ___  ___ 
# | |   | | '_ \| '__/ _` | '__| |/ _ \/ __|
# | |___| | |_) | | | (_| | |  | |  __/\__ \
# |_____|_|_.__/|_|  \__,_|_|  |_|\___||___/

# File handling and terminal subprocesses
import os
import tempfile
import subprocess
# Protein structure handling
from Bio import SeqIO, pairwise2
from Bio.SeqUtils import ProtParam
from Bio.PDB import PDBParser, Select, PDBIO
# OpenBabel run in terminal
# !sudo apt install openbabel

#   ____ _                         
#  / ___| | __ _ ___ ___  ___  ___ 
# | |   | |/ _` / __/ __|/ _ \/ __|
# | |___| | (_| \__ \__ \  __/\__ \
#  \____|_|\__,_|___/___/\___||___/

class ProteinSelect(Select):
    """Select only standard amino acid residues (protein residues)"""
    def accept_residue(self, residue):
        # Only accept standard amino acid residues (protein residues)
        return residue.id[0] == ' '
class ProteinNameSelect(Select):
    """Select only standard amino acid residues from specific chains"""
    def __init__(self, chain_ids):
        """
        Parameters:
            chain_ids (list): List of chain identifiers to keep (e.g., ['A', 'B'])
        """
        self.chain_ids = chain_ids
    
    def accept_chain(self, chain):
        return chain.id in self.chain_ids
    
    def accept_residue(self, residue):
        # Only accept standard amino acid residues (protein residues)
        return residue.id[0] == ' '

#     _              _ _ _                  
#    / \  _   ___  _(_) (_) __ _ _ __ _   _ 
#   / _ \| | | \ \/ / | | |/ _` | '__| | | |
#  / ___ \ |_| |>  <| | | | (_| | |  | |_| |
# /_/   \_\__,_/_/\_\_|_|_|\__,_|_|   \__, |
#                                     |___/ 

# 5
def extract_chain_sequence(
        pdb_file: str, 
        chain_id: str
) -> str:
    """
    Extract amino acid sequence from a specific chain in PDB file.
    
    Parameters:
        - pdb_file (str): Path to PDB file
        - chain_id (str): Chain identifier
    
    Return:
        str: Amino acid sequence (single letter code)
    """
    # Open PDB file and retrieve structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    # Amino acids Dictionary    
    amino_acids = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    # Create sequence
    sequence = ''
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[0] == ' ':  # Standard residue
                        res_name = residue.resname.strip()
                        if res_name in amino_acids:
                            sequence += amino_acids[res_name]
                        else:
                            sequence += 'X'  # Unknown residue

    return sequence

# 4
def calculate_sequence_similarity(
        seq1: str,
        seq2: str
) -> tuple:
    """
    Calculate sequence similarity between two sequences using Needleman-Wunsch.
    
    Parameters:
        - seq1 (str): First sequence
        - seq2 (str): Second sequence
    
    Return:
        tuple: (similarity_percentage, alignment1, alignment2)
    """
    if not seq1 or not seq2:
        return 0.0, '', ''
    
    # Perform alignment
    alignments = pairwise2.align.globalxx(seq1, seq2)
    
    if not alignments:
        return 0.0, '', ''
    
    best = alignments[0]
    aligned_seq1, aligned_seq2, score, begin, end = best
    
    # Calculate similarity percentage
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    max_len = max(len(seq1), len(seq2))
    similarity = (matches / max_len * 100) if max_len > 0 else 0.0
    
    return similarity, aligned_seq1, aligned_seq2

# 3
def count_unknown_residues(
        sequence: str
) -> int:
    """Count number of unknown residues (X) in sequence"""
    return sequence.count('X')

# 2
def select_chains_to_keep(
        pdb_file: str,
        chain_ids: list,
        similarity_threshold: float = 80.0
) -> list:
    """
    Analyze chains and select which ones to keep based on similarity.
    Removes duplicate/redundant chains.
    
    Parameters:
        - pdb_file (str): Path to PDB file
        - chain_ids (list): List of chain identifiers to analyze
        - similarity_threshold (float): Similarity percentage threshold (default 80%)
    
    Return:
        list: Chain IDs to keep
    """
    if len(chain_ids) <= 1:
        return chain_ids
    
    # Extract sequences for all chains
    sequences = {}
    for chain_id in chain_ids:
        seq = extract_chain_sequence(pdb_file, chain_id)
        sequences[chain_id] = seq
        unknown_count = count_unknown_residues(seq)
        print(f"Chain {chain_id}: Length={len(seq)}, Unknown residues={unknown_count}")
    
    chains_to_remove = set()
    
    # Compare all pairs of chains
    for i, chain1 in enumerate(chain_ids):
        if chain1 in chains_to_remove:
            continue
        
        for chain2 in chain_ids[i+1:]:
            if chain2 in chains_to_remove:
                continue
            
            seq1 = sequences[chain1]
            seq2 = sequences[chain2]
            
            similarity, aligned1, aligned2 = calculate_sequence_similarity(seq1, seq2)
            
            print(f"\nComparing Chain {chain1} vs Chain {chain2}:")
            print(f"\t- Similarity: {similarity:.1f}%")
            
            if similarity >= similarity_threshold:
                print(f"\t- Chains are {similarity:.1f}% similar (≥{similarity_threshold}%)")
                
                unknown1 = count_unknown_residues(seq1)
                unknown2 = count_unknown_residues(seq2)

                print(f"\t- Chain {chain1}: {unknown1} unknown residues")
                print(f"\t- Chain {chain2}: {unknown2} unknown residues")

                # If one has unknown residues and the other doesn't
                if unknown1 > 0 and unknown2 == 0:
                    chains_to_remove.add(chain1)
                    print(f"\t- Discarding Chain {chain1} (has unknown residues)")
                elif unknown2 > 0 and unknown1 == 0:
                    chains_to_remove.add(chain2)
                    print(f"\t- Discarding Chain {chain2} (has unknown residues)")
                elif unknown1 == unknown2:
                    # Both have same unknown count, keep first and discard second
                    chains_to_remove.add(chain2)
                    print(f"\t- Discarding Chain {chain2} (redundant)")
                else:
                    # Discard the one with more unknown residues
                    if unknown1 > unknown2:
                        chains_to_remove.add(chain1)
                        print(f"\t- Discarding Chain {chain1} (more unknown residues)")
                    else:
                        chains_to_remove.add(chain2)
                        print(f"\t- Discarding Chain {chain2} (more unknown residues)")
    
    chains_to_keep = [c for c in chain_ids if c not in chains_to_remove]
    print(f"\n{'='*60}")
    print(f"Chains to keep: {chains_to_keep}")
    print(f"{'='*60}\n")
    
    return chains_to_keep

# 1
def get_protein_chains(
        pdb_file: str,
        protein_name: str
) -> list:
    """
    Extract chain identifiers for a specific protein from PDB metadata.
    
    Parameters:
        - pdb_file (str): Path to PDB file
        - protein_name (str): Name or partial name of the protein to find
    
    Return:
        list: Chain identifiers belonging to the protein
    """
    chains = []
    current_protein = None
    protein_chains_dict = {}
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('COMPND'):
                    content = line[10:].strip()
                    
                    if 'MOL_ID:' in content:
                        mol_id = content.split('MOL_ID:')[1].strip().split(';')[0]
                        current_protein = mol_id
                        protein_chains_dict[current_protein] = {
                            'name': '',
                            'chains': []
                        }
                    
                    elif 'MOLECULE:' in content and current_protein:
                        mol_name = content.split('MOLECULE:')[1].strip().rstrip(';')
                        protein_chains_dict[current_protein]['name'] = mol_name
                    
                    elif 'CHAIN:' in content and current_protein:
                        chain_str = content.split('CHAIN:')[1].strip().rstrip(';')
                        chain_list = [c.strip() for c in chain_str.split(',')]
                        protein_chains_dict[current_protein]['chains'] = chain_list
    
    except Exception as e:
        print(f"Error reading PDB file: {e}")
        return []
    
    protein_name_lower = protein_name.lower()
    for mol_id, info in protein_chains_dict.items():
        if protein_name_lower in info['name'].lower():
            chains = info['chains']
            print(f"Found protein '{info['name']}' with chains: {chains}")
            return chains
    
    if not chains:
        print(f"Protein '{protein_name}' not found in the PDB file.")
        print("Available proteins:")
        for mol_id, info in protein_chains_dict.items():
            print(f"\t- {info['name']}")
            print(f"\t\tChains: {info['chains']}")
    
    return chains

#  _____                 _   _                 
# |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

# 2
def convert_to_pdb(
        input_file: str, 
        output_pdb: str
):
    '''
    Convert various file formats to PDB using Open Babel.
    
    Parameters:
        - input_file (str): a file with extension: ENT, XYZ, PQR, MCIF, MMCIF
        - output_pdb (str): output PDB filename
    
    Return:
        None. It will store a file in current directory.
    '''

    available_files = ('.ent', '.xyz', '.pqr', '.mcif', '.mmcif', '.pdbqt')

    # Determine the input format based on the file extension
    extension = os.path.splitext(input_file)[1].lower()
        
    if extension not in available_files:
        raise ValueError(f'Unsupported file format: {extension}. Please use an available format from {available_files}')
    else:
        # Remove dot
        if extension.startswith('.'):
            extension = extension[1:]
    
    # Use Open Babel to convert to PDB
    command = f'obabel -i{extension} "{input_file}" -opdb -O "{output_pdb}"'

    try:
            subprocess.run(command, shell=True, check=True)
            print(f'Conversion successful: {input_file} -> {output_pdb}')
    except subprocess.CalledProcessError as e:
        print(f'Error during conversion: {e}')
        raise

# 1
def clean_and_convert_pdb_to_pdbqt(
        input_receptor: str, 
        protein_name: str = None, 
        deduplicate: bool = True, 
        similarity_threshold: float = 80.0,
        ph: float = 7.4,
        remove_polar_hydrogens: bool = True,
        charge_method: str = 'eem'
):
    """
    Clean the PDB file by removing non-protein residues and convert it to PDBQT format.
    Optionally filters to keep only specific protein chains and removes duplicate chains.

    Parameters:
        - input_receptor (str): Initial PDB file
        - protein_name (str): Optional protein name to filter chains
        - deduplicate (bool): Whether to remove redundant chains (default True)
        - similarity_threshold (float): Similarity threshold for deduplication (default 80%)
        - ph (float): pH value for protonation (default 7.4)
        - remove_polar_hydrogens (bool): Whether to remove polar hydrogens (default True)
        - charge_method (str): Partial charge calculation method: 'eem', 'gasteiger', or 'mmff94' (default 'eem')
    
    Return:
        tuple: (cleaned_pdb, output_pdbqt) file paths
    """
    input_extension = os.path.splitext(input_receptor)[1].lower()

    if input_extension == '.pdb':
        input_pdb = input_receptor
    elif input_extension in ('.ent', '.xyz', '.gpr', '.cif', '.mmcif', '.mol2'):
        input_pdb = os.path.splitext(input_receptor)[0] + '.pdb'
        try:
            convert_to_pdb(input_receptor, input_pdb)
        except Exception as e:
            print(e)
            raise
    else:
        raise ValueError("Input receptor must be a PDB file or a valid convertible format.")

    # If protein_name is specified, extract chains for that protein
    if protein_name:
        target_chains = get_protein_chains(input_pdb, protein_name)
        if not target_chains:
            raise ValueError(f"Could not find chains for protein '{protein_name}'")
    else:
        # Get all chains from the structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', input_pdb)
        target_chains = []
        for model in structure:
            for chain in model:
                target_chains.append(chain.id)
        target_chains = list(set(target_chains))

    # Deduplicate chains if enabled
    if deduplicate and len(target_chains) > 1:
        print(f"\nDeduplicating chains with {similarity_threshold}% similarity threshold...")
        target_chains = select_chains_to_keep(input_pdb, target_chains, similarity_threshold)

    print(f"Final chains to keep: {target_chains}")
    selector = ProteinNameSelect(target_chains)

    # Define filenames for cleaned PDB and output PDBQT
    base_name = os.path.splitext(input_pdb)[0]
    if protein_name:
        cleaned_pdb = f"{base_name}_{protein_name.replace(' ', '_')}_cleaned.pdb"
        output_pdbqt = f"{base_name}_{protein_name.replace(' ', '_')}.pdbqt"
    else:
        cleaned_pdb = f"{base_name}_cleaned.pdb"
        output_pdbqt = f"{base_name}.pdbqt"

    # Parse the PDB file and remove non-protein residues
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    # Save the cleaned protein-only PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(cleaned_pdb, select=selector)

    # Remove UNK and UNL residues from the cleaned PDB file
    temp_pdb = cleaned_pdb.replace('.pdb', '_temp.pdb')
    with open(cleaned_pdb, 'r') as infile, open(temp_pdb, 'w') as outfile:
        for line in infile:
            # Skip lines with UNK or UNL residues
            if line.startswith('ATOM') or line.startswith('HETATM'):
                res_name = line[17:20].strip()
                if res_name not in ('UNK', 'UNL'):
                    outfile.write(line)
            else:
                outfile.write(line)

    os.replace(temp_pdb, cleaned_pdb)
    print(f"Removed UNK/UNL residues from {cleaned_pdb}")

    # Build OpenBabel command with parameters
    obabel_options = '-xr' if remove_polar_hydrogens else ''
    command = f'obabel -ipdb "{cleaned_pdb}" -opdbqt -O "{output_pdbqt}" {obabel_options} -p {ph} --partialcharge {charge_method}'
    
    # Execute the command using subprocess
    with open(os.devnull, 'w') as devnull:
        try:
            subprocess.run(command, shell=True, check=True, stderr=devnull)
            print(f'Conversion successful: {cleaned_pdb} -> {output_pdbqt}')
        except subprocess.CalledProcessError as e:
            raise ValueError(f"Error during conversion: {e}")
    
    # Remove UNK/UNL residues from the generated PDBQT file
    temp_pdbqt = output_pdbqt.replace('.pdbqt', '_temp.pdbqt')
    residues_removed = 0
    
    with open(output_pdbqt, 'r') as infile, open(temp_pdbqt, 'w') as outfile:
        current_residue_id = None
        skip_residue = False
        atom_count = 0
        
        for line in infile:
            # Check if this is an ATOM/HETATM line
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Extract residue name (positions 17-20 in PDBQT format)
                if len(line) >= 20:
                    res_name = line[17:20].strip()
                    
                    # Extract residue number (positions 23-26)
                    if len(line) >= 26:
                        residue_id = line[22:26].strip()
                        
                        # Check if we're starting a new residue
                        if current_residue_id != residue_id:
                            # If previous residue was UNK/UNL, we've already been skipping it
                            # Now check if the new residue is UNK/UNL
                            current_residue_id = residue_id
                            skip_residue = res_name in ('UNK', 'UNL')
                            
                            if skip_residue:
                                residues_removed += 1
                    
                    # Skip this atom if it belongs to UNK/UNL residue
                    if skip_residue:
                        continue
                    
                    # Write the atom line if not UNK/UNL
                    outfile.write(line)
                    atom_count += 1
            else:
                # Write non-atom lines (REMARK, TER, etc.)
                outfile.write(line)
    
    # Replace the original PDBQT with the cleaned version
    os.replace(temp_pdbqt, output_pdbqt)
    
    if residues_removed > 0:
        print(f"Removed {residues_removed} UNK/UNL residues from {output_pdbqt}")
    
    return cleaned_pdb, output_pdbqt

#  _____                     _   _             
# | ____|_  _____  ___ _   _| |_(_) ___  _ __  
# |  _| \ \/ / _ \/ __| | | | __| |/ _ \| '_ \ 
# | |___ >  <  __/ (__| |_| | |_| | (_) | | | |
# |_____/_/\_\___|\___|\__,_|\__|_|\___/|_| |_|

if __name__ == "__main__":
    input_filename = '4TPS.pdb'
    input_receptor = os.path.join('docking_functions', input_filename)
    protein_to_extract = 'DnaA'

    try:
        print(f"\n{'='*60}")
        print(f"Processing: {input_filename}")
        print(f"Extracting protein: {protein_to_extract}")
        print(f"{'='*60}\n")
        
        cleaned_pdb, output_pdbqt = clean_and_convert_pdb_to_pdbqt(
            input_receptor, 
            protein_name=protein_to_extract,
            deduplicate=True,
            similarity_threshold=80.0,
            ph=7.4,
            charge_method='gasteiger',
            remove_polar_hydrogens=True
        )
        print(f'Cleaned PDB file generated: {cleaned_pdb}')
        print(f'PDBQT file generated: {output_pdbqt}')
    except Exception as e:
        print(f"Error: {e}")

# END