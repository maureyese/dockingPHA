import os
import subprocess
from Bio.PDB import PDBParser, Select, PDBIO

# =============================
# FUNCTIONS
# =============================

def convert_to_pdb(input_file:str, output_pdb:str):
    '''
    Convert various file formats to PDB using Open Babel.
    
    Parameters:
        - input_file (str): a file with extension: ENT, XYZ, PQR, MCIF, MMCIF
        - output_file (str): output PDB filename
    
    Return:
        None. It will store a file in current directory.
    '''

    available_files = ('.ent', '.xyz', '.pqr', '.mcif', '.mmcif', 'pdbqt')

    # Determine the input format based on the file extension
    extension = os.path.splittext(input_file)[1].lower()
        
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

class ProteinSelect(Select):
    def accept_residue(self, residue):
        # Only accept standard amino acid residues (protein residues)
        return residue.id[0] == ' '

def clean_and_convert_pdb_to_pdbqt(input_receptor:str):
    """
    Clean the PDB file by removing non-protein residues and convert it to PDBQT format.

    Parameters:
    - input_pdb (str): Initial PDB file
    - cleaned_pdb (str): Temporal PDB file prior conversion to PDBQT
    - output_pdbqt (str): Final PDBQT file

    Return:
        None. A file will be stored in current directory
    """
    # Get the extension
    input_extension = os.path.splitext(input_receptor)[1].lower()

    # Verify which type of extension file
    if input_extension == '.pdb':
        # It's a PDB file, use it directly
        input_pdb = input_receptor
        #print('PDB file detected. Starting pre-processing.')
    elif input_extension == ('.ent', '.xyz', '.gpr', '.cif', '.mmcif', '.mol2'):
        input_pdb = os.path.splitext(input_receptor)[0] + '.pdb'
        try:
            convert_to_pdb(input_receptor, input_receptor)
        except Exception as e:
            print(e)
    else:
        raise ValueError("Input receptor must be a PDB ID or a valid file.")
    
    # Define filenames for cleaned PDB and output PDBQT
    cleaned_pdb = f"{os.path.splitext(input_pdb)[0]}_cleaned.pdb"
    output_pdbqt = f"{os.path.splitext(input_pdb)[0]}.pdbqt"
    
    # Parse the PDB file and remove non-protein residues
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)

    # Save the cleaned protein-only PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(cleaned_pdb, select=ProteinSelect()) # Save only protein residues

    # Command to convert the cleaned PDB to PDBQT using Open Babel
    print("Running OpenBabel command. It may take some time...")
    command = f'obabel -ipdb "{cleaned_pdb}" -opdbqt -O "{output_pdbqt}" -xr -p 7.4 --partialcharge eem'
    
    # Execute the command using subprocess
    # Redirect stderr to suppress Open Babel warnings
    with open(os.devnull, 'w') as devnull:
        try:
            subprocess.run(command, shell=True, check=True, stderr=devnull)
            print(f'Conversion successful: {cleaned_pdb} -> {output_pdbqt}')
        except subprocess.CalledProcessError as e:
            raise ValueError(f"Error during conversion: {e}")
    
    return cleaned_pdb, output_pdbqt

# =============================
# EXAMPLE USAGE
# =============================

def main():
    # Input receptor
    input_filename = '2BXP.pdb'
    output_dir = "docking_functions/"
    input_receptor = os.path.join(output_dir, input_filename)

    try:
        # Clean PDB file and convert to PDBQT
        cleaned_pdb, output_pdbqt = clean_and_convert_pdb_to_pdbqt(input_receptor)
        print(f'Cleaned PDB file generated: {cleaned_pdb}')
        print(f'PDBQT file generated: {output_pdbqt}')
    except Exception as e:
        print(e)
    
if __name__ == "__main__":
    main()