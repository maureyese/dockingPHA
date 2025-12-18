import warnings
import numpy as np
import os
import sys
import time
import subprocess
import secrets
import zipfile

from vina import Vina
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import PDBExceptions
from Bio.PDB.Polypeptide import is_aa
from openbabel import pybel
import MDAnalysis as mda

from PrepareLigand import calculate_ligand_center_of_mass, prepare_ligand
from PrepareReceptor import clean_and_convert_pdb_to_pdbqt

# Suppress specific warnings from Biopython
warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)
warnings.simplefilter('ignore', BiopythonWarning)

def create_zip(ruta_carpeta, nombre_zip):
    """
    Crea un archivo ZIP con el contenido de la carpeta especificada.

    :param ruta_carpeta: Ruta de la carpeta que se desea comprimir.
    :param nombre_zip: Nombre (o ruta con nombre) que tendrá el ZIP resultante.
    """
    # 'w' indica que se abrirá en modo escritura (si existe, se sobrescribe).
    with zipfile.ZipFile(nombre_zip, 'w', zipfile.ZIP_DEFLATED) as archivo_zip:
        # Recorremos todo el árbol de directorios.
        for carpeta_raiz, subcarpetas, ficheros in os.walk(ruta_carpeta):
            for fichero in ficheros:
                # Creamos la ruta absoluta al fichero.
                ruta_fichero = os.path.join(carpeta_raiz, fichero)
                # Queremos que dentro del ZIP aparezca la ruta relativa
                # en lugar de la ruta absoluta, para mantener la estructura.
                ruta_relativa = os.path.relpath(ruta_fichero, ruta_carpeta)

                # Agregamos el fichero al ZIP con la ruta relativa.
                archivo_zip.write(ruta_fichero, arcname=ruta_relativa)

def compute_bounding_box(pdb_file):
    '''
    Function to create box size
    '''
    u = mda.Universe(pdb_file)
    positions = u.atoms.positions
    min_coords = np.min(positions, axis=0)
    max_coords = np.max(positions, axis=0)
    size = max_coords - min_coords
    return size

def adjust_box_size(size, padding=5.0):
    '''
    This function increases the size of the bounding box by adding extra space (padding)
    to each dimension of the box, ensuring the ligands have enough room to bind.
    The padding parameter is the extra space added (default is 5.0 Å).
    '''
    adjusted_size = size + padding
    return [float(dim) for dim in adjusted_size]

def calculate_center_of_mass(cleaned_receptor_pdb):
    """
    Function to calculate center of mass (receptor)
    """
    parser = PDBParser()
    receptor_structure = parser.get_structure('receptor', cleaned_receptor_pdb)

    atoms = [atom for atom in receptor_structure.get_atoms()]
    center = [sum(coord) / len(atoms) for coord in zip(*[atom.coord for atom in atoms])]
    return [float(coord) for coord in center]

def align_ligand_receptor(ligand_pdbqt, receptor_center):
    # Load the ligand from the PDBQT file using Pybel
    pybel_mol = next(pybel.readfile("pdbqt", ligand_pdbqt))
    
    # Calculate the ligand's center of mass
    ligand_center = calculate_ligand_center_of_mass(ligand_pdbqt)
    
    # Calculate translation vector to align ligand center with receptor center
    translation_vector = np.array(receptor_center) - np.array(ligand_center)
    
    # Translate the ligand by updating each atom's position
    for atom in pybel_mol.atoms:
        new_coords = np.add(atom.coords, translation_vector)
        atom.OBAtom.SetVector(*new_coords)
    
    # Save the translated ligand back to PDBQT format
    pybel_mol.write("pdbqt", ligand_pdbqt, overwrite=True)
    print(f'Ligand centered to receptor and saved as: {ligand_pdbqt}')

def run_vina_blind_docking(receptor_pdbqt, 
                           ligand_pdbqt,
                           output_folder, 
                           center, 
                           size,
                           ligand_minimized,
                           exhaustiveness = 32,
                           num_modes = 9,
                           output_filename = "output"):
    """
    Runs blind docking using the vina Python library with dynamic grid box resizing.
    """

    # Dictionary to store the results
    vina_results = {}

    # Start vina
    v = Vina(sf_name='vina', cpu=-1)

    # Set receptor and ligand
    # We will suppose user uploaded prepared files for docking
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)

    # Set the score function and search space
    successful_grid = False
    max_attempts = 10
    attempt = 0

    print("\nValidation step to ensure box size covers receptor and ligand. Analyzing...")

    # Ensure box covers both receptor and ligand
    while not successful_grid and attempt < max_attempts:
        try:
            print(f"\nAttempt {attempt+1}: Box size {size}")
            v.compute_vina_maps(center=center, box_size=size)
            
            # Try scoring to confirm ligand is inside the box
            energy = v.score()  
            successful_grid = True  # If scoring is successful, the ligand is inside the box
        except RuntimeError as e:
            print(f"\nError setting grid box or ligand outside grid: {e}")
            # Expand the grid box size and reattempt if the ligand is outside the box
            size = [dim + 15.0 for dim in size]  # Increase the box by 15 Å each time
            print(f"\nNew box size: {size}")
            attempt += 1
    
    if not successful_grid:
        print(f"\nFailed to set the grid box after {max_attempts} attempts.")
        # TODO check if this is the correct way to return an error
        return  f"\nFailed to set the grid box after {max_attempts} attempts." # Exit if unable to find a suitable grid box size

    # Store the final grid box size
    vina_results['grid_box'] = size
    # Store the center of the grid box
    vina_results['grid_center'] = center

    vina_results['minimization'] = {}

    if ligand_minimized == False:
        # Score the current pose before minimization
        energy = v.score()
        print('\nScore before minimization: %.3f (kcal/mol)' % energy[0])
        vina_results['minimization']['before_min'] = energy[0]

        # Minimize locally the current pose of the ligand
        energy_minimized = v.optimize()
        print('\nScore after minimization : %.3f (kcal/mol)' % energy_minimized[0])
        vina_results['minimization']['after_min'] = energy_minimized[0]
        v.write_pose(ligand_pdbqt, overwrite=True)
    else:
        vina_results['minimization']['before_min'] = None
        vina_results['minimization']['after_min'] = None

    # Run docking
    start_time = time.time()

    try:
        v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes)
    except RuntimeError as e:
        print(f"Docking failed: {e}")
        return

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"\nTime taken: {elapsed_time:.2f} seconds")
    
    vina_results["elapsed_time"] = elapsed_time # Store time taken

    # Get and print the scores
    energies = v.energies(n_poses=num_modes)
    print("\nDocking completed.")
    print("Scores for the best poses:")

    # Store energy results
    vina_results['results'] = []

    # Get and print the scores
    energies = v.energies(n_poses=num_modes)
    print("\nDocking completed.")
    print("Scores for the best poses:")
    for i, energy in enumerate(energies):
        print(f"Pose {i+1}: Affinity = {energy[0]:.2f} kcal/mol")
        vina_results['results'].append({"pose": i+1, "affinity": energy[0], "rmsd_lb":None, "rmsd_ub":None})

    # Store poses and save them as zip file
    os.makedirs(output_folder, exist_ok=True)
    
    pose_files = []
    combined_pose_files = []

    for i in range(1, num_modes + 1):
        # Write out each pose by itself
        pose_file = os.path.join(output_folder, f"ligand_pose_{i}.pdbqt")
        v.write_poses(pose_file, n_poses=i, overwrite=True)
        pose_files.append(pose_file)

        # Now combine each pose with the receptor in a single PDBQT
        combined_file = os.path.join(output_folder, f"complex_pose_{i}.pdbqt")
        with open(combined_file, 'w') as out_f:

            # First write the receptor
            with open(receptor_pdbqt, 'r') as rec_f:
                for line in rec_f:
                    out_f.write(line)

            # Then write the pose (skipping any remarks if you prefer)
            with open(pose_file, 'r') as pose_f:
                for line in pose_f:
                    # If you want to keep REMARK lines, remove the check below
                    if not line.startswith('REMARK'):
                        out_f.write(line)

        combined_pose_files.append(combined_file)

    # Zip the pose files
    # TODO añadir al zip los archivos que se llamen complex pose
    # junto con ligand pose
    create_zip("output_folder", f"{output_folder}/{output_filename}.zip") #TODO Cambiar rutas en produccion
    # Save zip name
    vina_results['zip_file'] = f"{output_filename}.zip"

    #TODO Delete the individual pose files after zipping
    # for pose_file in pose_files:
    #     os.remove(pose_file)
    # TODO update Mau: Guardar en la base de datos, tanto archivos individuales como zip
    
    print("\nBest poses saved in 'output' folder and zipped as 'poses.zip'.")

    return vina_results

def main():
    #region Obtenemos variables
    if len(sys.argv) != 5: # TODO repair this
        # Ejemplo python BlindDocking.py ligand.pdbqt receptor.pdbqt blind ligand_minimized.pdbqt
        print("Faltan parametros")
        sys.exit(1)

    # For blind docking
    # Add parameters to the command:
    #   - ligand_filename
    #   - rigid_receptor_filename
    #   - ligand_minimized
    #   - output_folder
        
    LIGAND_FILENAME = sys.argv[1]
    RECEPTOR_FILENAME = sys.argv[2]
    LIGAND_MINIMIZED = sys.argv[3]
    OUTPUT_FILENAME = sys.argv[4]
    print("Ligand: ", LIGAND_FILENAME)
    print("Receptor: ", RECEPTOR_FILENAME)
    print("Ligand minimized: ", LIGAND_MINIMIZED)
    print("Output filename: ", OUTPUT_FILENAME)
    #endregion

    #! Generamos un nombre unico para el archivo ligando
    #TODO Pendiente
    # output_filename = secrets.token_hex(16) + ".zip" 
    # while crud.docking.get_by_ligand_filename(db, ligand_filename):
    #     ligand_filename = secrets.token_hex(16)
    
    # Input ligand
    # ligand_filename = "caffeine.pdbqt"
    # receptor = os.path.join('Examples', ligand_filename)
    ligand = os.path.join('input', LIGAND_FILENAME) # TODO change
    #ligand = "Examples/caffeine.sdf"

    # Input receptor
    # receptor_filename = "2BXP.pdbqt"
    # receptor = os.path.join('Examples', receptor_filename)
    receptor = os.path.join('input', RECEPTOR_FILENAME) # TODO change
    #receptor = "Examples/2BXP.pdb" 

    # Run blind docking
    try:
        vina_results = run_vina_blind_docking(receptor_pdbqt=receptor,
                           ligand_pdbqt=ligand,
                           output_folder = f"output/{OUTPUT_FILENAME}",
                           center=calculate_center_of_mass(receptor),
                           size=adjust_box_size(compute_bounding_box(receptor)),
                           ligand_minimized=LIGAND_MINIMIZED,
                           exhaustiveness=8,
                           num_modes=9,
                           output_filename=OUTPUT_FILENAME)

        print(vina_results)
    except Exception as e:
        print(f"Error running blind docking: {e}") # TODO check if this is the correct way to return an error

if __name__ == "__main__":
    main()