'''
# LIGAND PDBQT PREPARATION
=================================
This module prepares ligands for molecular docking by converting various chemical 
file formats and SMILES strings to PDBQT format with customizable OpenBabel parameters.

## Main processes:

1. Format Conversion:
   - Converts multiple ligand file formats (CAN, MDL, MOL, MOL2, PDB, SD, SDF, SMI, SMILES, XYZ) to PDBQT
   - Direct file conversion for existing molecular structures

2. SMILES to PDBQT Conversion:
   - Converts SMILES notation to 3D molecular structures
   - Generates 3D coordinates automatically
   - Handles intermediate file management (temporary PDB files)

3. 3D Coordinate Generation:
   - Generates 3D coordinates when not present in input files
   - Uses OpenBabel's geometry optimization for proper spatial arrangement
   - Ensures proper molecular conformation for docking

## Configurable OpenBabel Parameters:

   - pH: Protonation pH value (default 7.4)
   - Charge Method: Partial charge calculation ('gasteiger', 'eem', 'mmff94')
   - Generate 3D: Create 3D coordinates if missing (default True)

## Output:

   - PDBQT file ready for molecular docking simulations
   - Automatic cleanup of intermediate files
   - Maintains ligand nomenclature in output filenames
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
from openbabel import pybel
# OpenBabel run in terminal
# !sudo apt install openbabel

#  _____                 _   _                 
# |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

def prepare_ligand(
        input_file: str,
        output_dir: str = 'Examples',
        ph: float = 7.4,
        charge_method: str = 'gasteiger',
        generate_3d: bool = True
):
    """
    Prepare a ligand in PDBQT format from various input types using Open Babel.

    Parameters:
        - input_file (str): Ligand filename (PDBQT, CAN, MDL, MOL, MOL2, PDB, SD, SDF, SMI, SMILES, XYZ)
        - output_dir (str): Output directory for prepared ligand (default 'Examples')
        - ph (float): pH value for protonation (default 7.4)
        - charge_method (str): Partial charge calculation method: 'gasteiger', 'eem', or 'mmff94' (default 'gasteiger')
        - generate_3d (bool): Generate 3D coordinates if not present (default True)
    
    Returns:
        str: Path to the prepared PDBQT ligand file
    """
    os.makedirs(output_dir, exist_ok=True)
    
    file_types = ('.pdbqt', '.can', '.mdl', '.mol', '.mol2', '.pdb', '.sd', '.sdf', '.smi', '.smiles', '.xyz')
    extension = os.path.splitext(input_file)[1].lower().strip()
    
    if extension not in file_types:
        raise ValueError(f"Invalid file format: {extension}. Try using any format from {file_types}")
    
    # Get output filename
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}.pdbqt")
    
    # Build OpenBabel command with parameters
    gen3d_flag = '--gen3d' if generate_3d else ''
    
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt") as tmp_pdbqt:
            tmp_pdbqt.close()
            # Use Open Babel command to convert input file to PDBQT
            command = f'obabel "{input_file}" -O "{tmp_pdbqt.name}" {gen3d_flag} --partialcharge {charge_method} -p {ph}'
            subprocess.run(command, shell=True, check=True)
            # Move the temporary PDBQT file to the desired output
            os.replace(tmp_pdbqt.name, output_file)
            print(f"Ligand prepared successfully: {output_file}")
        
        return output_file
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Error preparing ligand: {e}")

def prepare_ligand_from_smiles(
        smiles: str,
        output_filename: str,
        output_dir: str = "Examples",
        ph: float = 7.4,
        charge_method: str = 'gasteiger'
):
    """
    Convert SMILES to PDBQT format for docking.
    
    Parameters:
        - smiles (str): SMILES string of the ligand
        - output_filename (str): Desired output filename (with .pdbqt extension)
        - output_dir (str): Directory to save prepared ligand (default 'Examples')
        - ph (float): pH value for protonation (default 7.4)
        - charge_method (str): Partial charge calculation method: 'gasteiger', 'eem', or 'mmff94' (default 'gasteiger')
    
    Returns:
        str: Path to the prepared PDBQT ligand file, or None if preparation failed
    """
    if pybel is None:
        print("Error: OpenBabel Python bindings required for SMILES preparation")
        return None
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Validate output filename has .pdbqt extension
    if not output_filename.lower().endswith('.pdbqt'):
        print(f"Warning: Output filename '{output_filename}' should have .pdbqt extension")
        # If not .pdbqt, append it
        output_filename = os.path.splitext(output_filename)[0] + '.pdbqt'
        print(f"Using '{output_filename}' instead")
    
    # Extract base name (without extension and directory)
    base_name = os.path.splitext(os.path.basename(output_filename))[0]
    
    # Create temporary SMILES file with base name
    temp_smiles = os.path.join(output_dir, f"{base_name}.smi")
    with open(temp_smiles, 'w') as f:
        f.write(smiles)
    
    try:
        # Read SMILES file and convert to molecule
        mol = next(pybel.readfile("smi", temp_smiles))
        
        # Generate 3D coordinates
        mol.make3D()
        
        # Save intermediate PDB with base name
        ligand_pdb = os.path.join(output_dir, f"{base_name}.pdb")
        mol.write("pdb", ligand_pdb, overwrite=True)
        
        # Create final PDBQT path with the exact output filename
        ligand_pdbqt = os.path.join(output_dir, output_filename)

        # Convert PDB to PDBQT using openbabel command with parameters
        command = f'obabel -ipdb "{ligand_pdb}" -opdbqt -O "{ligand_pdbqt}" --partialcharge {charge_method} -p {ph}'

        # Execute the conversion
        subprocess.run(command, shell=True, check=True)
        
        # Clean up intermediate files
        if os.path.exists(temp_smiles):
            os.remove(temp_smiles)
        if os.path.exists(ligand_pdb):
            os.remove(ligand_pdb)
        
        print(f"Ligand prepared and saved as: {ligand_pdbqt}")
        return ligand_pdbqt
    except Exception as e:
        print(f"Error preparing ligand from SMILES: {e}")
        # Clean up temporary files on error
        if os.path.exists(temp_smiles):
            os.remove(temp_smiles)
        if 'ligand_pdb' in locals() and os.path.exists(ligand_pdb):
            os.remove(ligand_pdb)
        return None

if __name__ == "__main__":
    # Input ligand
    input_filename = "caffeine.sdf"
    input_ligand = os.path.join('docking_functions', input_filename)

    try:
        # Clean file and convert to PDBQT
        ligand_pdbqt = prepare_ligand(
             input_file=input_ligand,
             output_dir=os.path.join('tools', 'Examples'),
             ph=7.4,
             charge_method='gasteiger'
        )
    except Exception as e:
        print(e)

    # Input ligand in SMILES
    smiles_string = "CC(=O)OC1=CC=CC=C1C(=O)O"

    try:
        # Prepare ligand from SMILES
        ligand_pdbqt = prepare_ligand_from_smiles(
            smiles=smiles_string,
            output_filename='aspirin.pdbqt',
            output_dir=os.path.join('tools', 'Examples'),
            ph=7.4,
            charge_method='gasteiger'
        )
    except Exception as e:
        print(e)