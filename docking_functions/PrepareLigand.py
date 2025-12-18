import os
import tempfile
import subprocess
import numpy as np
from openbabel import pybel

# =============================
# FUNCTIONS
# =============================

def prepare_ligand(input_file:str, output_dir:str):
    """
    Prepare a ligand in PDBQT format from various input types using Open Babel.

    Parameters:
        - input_file: filename (PDBQT, CAN, MDL, MOL, MOL2, 
          PDB, SD, SDF, SMI, SMILES, XYZ)
        - output_file: The output PDBQT ligand filename.

    Returns:
        - None. A file will be stored in the current directory.
    """

    # Get output filename
    output_file = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(input_file))[0]}.pdbqt")

    file_types = ('.pdbqt', '.can', '.mdl', '.mol', '.mol2', '.pdb', '.sd', '.sdf', '.smi', '.smiles', '.xyz')

    # Use Open Babel to convert the input file directly to PDBQT
    extension = os.path.splitext(input_file)[1].lower().strip()

    if extension in file_types:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdbqt") as tmp_pdbqt:
            tmp_pdbqt.close()
            # Use Open Babel command to convert input file to PDBQT
            subprocess.run([
                'obabel', input_file, '-O', tmp_pdbqt.name,
                '--gen3d', '--partialcharge', 'gasteiger',
                '--protonate3D', '7.4'  # Protonate at pH 7.4
            ], check=True)
            # Move the temporary PDBQT file to the desired output
            os.replace(tmp_pdbqt.name, output_file)
            print("Ligand prepared successfully.")
        
        return output_file
    else:
        raise ValueError(f"Invalid file. Try using any format from {file_types}")

def calculate_ligand_center_of_mass(ligand_pdbqt):
    """
    Calculate the center of mass for a ligand using Pybel and return it as a list of three values.
    """
    mol = next(pybel.readfile("pdbqt", ligand_pdbqt))
    atom_positions = np.array([atom.coords for atom in mol.atoms])
    center_of_mass = np.mean(atom_positions, axis=0)
    return center_of_mass.tolist()

# =============================
# EXAMPLE USAGE
# =============================

def main():
    # Input ligand
    input_filename = "caffeine.sdf"
    output_dir = "docking_functions/"
    input_ligand = os.path.join(output_dir, input_filename)

    try:
        # Clean file and convert to PDBQT
        output_ligand = prepare_ligand(input_ligand, output_dir)
        print(f"Output ligand: {output_ligand}")
        # Calculate center of mass
        com = calculate_ligand_center_of_mass(output_ligand)
        print(f"Center of mass: {com}")
    except Exception as e:
        print(e)

if __name__ == "__main__":
    main()