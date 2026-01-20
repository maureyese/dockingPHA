#  _     _ _                    _           
# | |   (_) |__  _ __ __ _ _ __(_) ___  ___ 
# | |   | | '_ \| '__/ _` | '__| |/ _ \/ __|
# | |___| | |_) | | | (_| | |  | |  __/\__ \
# |_____|_|_.__/|_|  \__,_|_|  |_|\___||___/

import warnings
import os
import time
import numpy as np
import csv
import zipfile
import sys
from io import StringIO

# Molecular modeling and bioinformatics libraries
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, PDBIO, Select, PDBExceptions
from openbabel import pybel
import MDAnalysis as mda
from vina import Vina

# Suppress specific warnings
warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)
warnings.simplefilter('ignore', BiopythonWarning)

# Custom functions
from docking_functions.PrepareLigand import prepare_ligand_from_smiles
from docking_functions.PrepareReceptor import clean_and_convert_pdb_to_pdbqt

#   ____ _                         
#  / ___| | __ _ ___ ___  ___  ___ 
# | |   | |/ _` / __/ __|/ _ \/ __|
# | |___| | (_| \__ \__ \  __/\__ \
#  \____|_|\__,_|___/___/\___||___/

class LogCapture:
    """Capture console output to both stdout and a log file"""
    def __init__(self, filepath):
        self.terminal = sys.stdout
        self.log = open(filepath, "w")
   
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()
   
    def flush(self):
        self.terminal.flush()
        self.log.flush()
   
    def close(self):
        self.log.close()

#     _              _ _ _                  
#    / \  _   ___  _(_) (_) __ _ _ __ _   _ 
#   / _ \| | | \ \/ / | | |/ _` | '__| | | |
#  / ___ \ |_| |>  <| | | | (_| | |  | |_| |
# /_/   \_\__,_/_/\_\_|_|_|\__,_|_|   \__, |
#                                     |___/ 

def calculate_receptor_center_of_mass(pdb_file):
    """
    Calculate center of mass of a protein from PDB file
    
    Parameters:
        pdb_file: Path to PDB file
        
    Returns:
        List of [x, y, z] coordinates of center of mass
    """
    parser = PDBParser(QUIET=True)
    receptor_structure = parser.get_structure('receptor', pdb_file)
    atoms = [atom for atom in receptor_structure.get_atoms()]
    center = [sum(coord) / len(atoms) for coord in zip(*[atom.coord for atom in atoms])]
    return [float(coord) for coord in center]

def compute_bounding_box(pdb_file):
    """
    Calculate the bounding box of a protein
    
    Parameters:
        pdb_file: Path to PDB file
        
    Returns:
        numpy array of [size_x, size_y, size_z]
    """
    u = mda.Universe(pdb_file)
    positions = u.atoms.positions
    min_coords = np.min(positions, axis=0)
    max_coords = np.max(positions, axis=0)
    size = max_coords - min_coords
    return size

def adjust_box_size(size, padding=5.0):
    """
    Increase the size of the bounding box by adding padding
    
    Parameters:
        size: Original box size [x, y, z]
        padding: Padding to add to each dimension
        
    Returns:
        Adjusted box size [x, y, z]
    """
    adjusted_size = size + padding
    return [float(dim) for dim in adjusted_size]

def calculate_ligand_center_of_mass(ligand_pdbqt):
    """
    Calculate the center of mass for a ligand using Pybel and return it as a list of three values.
    """
    mol = next(pybel.readfile("pdbqt", ligand_pdbqt))
    atom_positions = np.array([atom.coords for atom in mol.atoms])
    center_of_mass = np.mean(atom_positions, axis=0)
    return center_of_mass.tolist()

def align_ligand_receptor(ligand_pdbqt, receptor_center):
    """
    Align ligand to receptor center of mass
    
    Parameters:
        - ligand_pdbqt (str): Path to ligand PDBQT file
        - receptor_center (list): [x, y, z] coordinates of receptor center
        
    Returns:
        None (ligand file is modified in place)
    """
    try:
        pybel_mol = next(pybel.readfile("pdbqt", ligand_pdbqt))
        atom_positions = np.array([atom.coords for atom in pybel_mol.atoms])
        ligand_center = np.mean(atom_positions, axis=0)
        translation_vector = np.array(receptor_center) - ligand_center
        
        for atom in pybel_mol.atoms:
            new_coords = np.add(atom.coords, translation_vector)
            atom.OBAtom.SetVector(*new_coords)
        
        pybel_mol.write("pdbqt", ligand_pdbqt, overwrite=True)
        print(f'  Ligand centered to receptor at {[round(c, 2) for c in receptor_center]}')
    except Exception as e:
        print(f"\t- Warning: Could not align ligand: {e}")

#  _____                 _   _                 
# |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

def run_vina_blind_docking(
        receptor_pdbqt, 
        ligand_pdbqt,
        output_folder, 
        center,
        size,
        output_filename: str = "docking_result",
        ligand_minimized: bool = False,
        exhaustiveness: int = 32,
        num_modes: int = 8,
        cpu: int = -1,
        align_ligand: bool = True,
        seed: int = 1234,
        save_poses: str = "first"
):
    """
    Runs blind docking using the vina Python library.
    Creates organized output in a subfolder named after output_filename.
    
    Parameters:
        - receptor_pdbqt (str): Path to receptor PDBQT file
        - ligand_pdbqt (str): Path to ligand PDBQT file
        - output_folder (str): Parent directory to save results
        - center (list): Grid center [x, y, z]
        - size (list or numpy array): Grid size [x, y, z]
        - output_filename (str): Base name for output files and subfolder (default 'docking_result')
        - ligand_minimized (bool): Whether ligand is already minimized (default False)
        - exhaustiveness (int): Exhaustiveness parameter for Vina (default 32)
        - num_modes (int): Number of poses to generate (default 8)
        - cpu (int): Number of CPUs for computation (default -1)
        - align_ligand (bool): Align ligand to receptor center before docking (default True)
        - seed (int): Random seed for reproducibility (default 1234)
        - save_poses (str): How to save poses:
            * "none": Don't save any poses
            * "first": Save only the best (first) pose (default)
            * "all": Save all generated poses
            * "1-5" or "1-n": Save poses in range (1-indexed)
        
    Returns:
        tuple: (vina_results, error_message, output_dir) - Dictionary with docking results, error message if any, and path to output directory
    """
    # Create organized output directory structure
    base_name = os.path.splitext(output_filename)[0]
    docking_output_dir = os.path.join(output_folder, base_name)
    os.makedirs(docking_output_dir, exist_ok=True)
    
    # Setup logging
    log_file = os.path.join(docking_output_dir, f"{base_name}_docking.log")
    log_capture = LogCapture(log_file)
    sys.stdout = log_capture
    
    try:
        # Validate center
        if center is None or len(center) == 0:
            print("Please provide a center of receptor.")
            return None, "No center provided", docking_output_dir

        vina_results = {}
        
        # Align ligand to receptor center before docking
        if align_ligand:
            print("\nAligning ligand to receptor center...")
            align_ligand_receptor(ligand_pdbqt, center)
        
        # Set vina 
        v = Vina(sf_name='vina', cpu=cpu, seed=seed)

        # Set receptor and ligand
        v.set_receptor(receptor_pdbqt)
        v.set_ligand_from_file(ligand_pdbqt)

        # Convert size to list if it's a numpy array
        if size is not None:
            size = [float(dim) for dim in size]
            print(f"\t- Grid box size: {[round(s, 1) for s in size]}")

        # Set docking grid to cover all protein
        if size is None or len(size) == 0:
            successful_grid = False
            max_attempts = 10
            attempt = 0
            size = [30.0, 30.0, 30.0]  # Default initial size

            print("\nSetting up docking grid...")

            while not successful_grid and attempt < max_attempts:
                try:
                    print(f"\t- Attempt {attempt+1}: Box size {[round(s, 1) for s in size]}")
                    v.compute_vina_maps(center=center, box_size=size)
                    energy = v.score()  
                    successful_grid = True
                except RuntimeError as e:
                    print(f"\t- Error: {str(e)[:100]}...")
                    size = [dim + 15.0 for dim in size]
                    attempt += 1

            if not successful_grid:
                return None, f"Failed to set the grid box after {max_attempts} attempts.", docking_output_dir
        else:
            # Size was provided, compute maps directly
            print("\nSetting up docking grid...")
            try:
                v.compute_vina_maps(center=center, box_size=size)
                energy = v.score()
                print("\t- Grid setup successful")
            except RuntimeError as e:
                error_msg = str(e)
                print(f"Error computing Vina maps: {error_msg}")
                
                # If ligand is outside grid, try to align and retry
                if "outside the grid box" in error_msg and align_ligand:
                    print("\n  Ligand was outside grid. Attempting to realign and retry...")
                    align_ligand_receptor(ligand_pdbqt, center)
                    v.set_ligand_from_file(ligand_pdbqt)
                    try:
                        v.compute_vina_maps(center=center, box_size=size)
                        energy = v.score()
                        print("\t- Grid setup successful after realignment")
                    except RuntimeError as retry_error:
                        print(f"\t- Realignment failed: {retry_error}")
                        return None, str(retry_error), docking_output_dir
                else:
                    return None, error_msg, docking_output_dir

        # Store values
        vina_results['grid_box'] = size
        vina_results['grid_center'] = center
        vina_results['minimization'] = {}

        # Energy-minimize ligand if not already minimized
        if not ligand_minimized:
            print("\nPerforming ligand minimization...")
            energy = v.score()
            print(f'  Score before minimization: {energy[0]:.3f} kcal/mol')
            vina_results['minimization']['before_min'] = float(energy[0])

            energy_minimized = v.optimize()
            print(f'  Score after minimization: {energy_minimized[0]:.3f} kcal/mol')
            vina_results['minimization']['after_min'] = float(energy_minimized[0])
            v.write_pose(ligand_pdbqt, overwrite=True)
        else:
            vina_results['minimization']['before_min'] = None
            vina_results['minimization']['after_min'] = None

        # Perform docking
        print(f"\nPerforming blind docking with exhaustiveness={exhaustiveness}...")
        start_time = time.time()
        try:
            v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes)
        except RuntimeError as e:
            print(f"Docking failed: {e}")
            return None, str(e), docking_output_dir
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"\t- Docking time: {elapsed_time:.2f} seconds")
        
        vina_results["elapsed_time"] = elapsed_time
        vina_results['results'] = []

        # Get all poses by writing them to file first, then parsing
        try:
            temp_all_poses = os.path.join(docking_output_dir, f"{base_name}_temp_all_poses.pdbqt")
            v.write_poses(temp_all_poses, n_poses=num_modes, overwrite=True)
            
            # Parse poses and extract affinity from REMARK lines
            all_poses_content = []
            current_pose = []
            pose_affinities = []
            
            with open(temp_all_poses, 'r') as f:
                for line in f:
                    if line.startswith('MODEL'):
                        if current_pose:
                            all_poses_content.append(''.join(current_pose))
                        current_pose = [line]
                    elif line.startswith('REMARK VINA RESULT:'):
                        # Parse affinity from REMARK line
                        try:
                            parts = line.split()
                            affinity = float(parts[3])
                            rmsd_lb = float(parts[4])
                            rmsd_ub = float(parts[5])
                            pose_affinities.append([affinity, rmsd_lb, rmsd_ub])
                        except (IndexError, ValueError):
                            pass
                        current_pose.append(line)
                    else:
                        current_pose.append(line)
                if current_pose:
                    all_poses_content.append(''.join(current_pose))
            
            # Build results from parsed poses
            actual_num_poses = len(all_poses_content)
            print(f"\n  Docking completed. Retrieved {actual_num_poses} poses:")
            
            if actual_num_poses > 0:
                for i in range(actual_num_poses):
                    if i < len(pose_affinities):
                        affinity, rmsd_lb, rmsd_ub = pose_affinities[i]
                        print(f"\t- Pose {i+1}: Affinity = {affinity:.2f} kcal/mol, RMSD_LB = {rmsd_lb:.2f}, RMSD_UB = {rmsd_ub:.2f}")
                        vina_results['results'].append({
                            "pose": i+1,
                            "affinity": affinity,
                            "rmsd_lb": rmsd_lb,
                            "rmsd_ub": rmsd_ub
                        })
                
                vina_results['best_affinity'] = pose_affinities[0][0] if pose_affinities else None
            else:
                print("\t- No poses generated")
                vina_results['best_affinity'] = None
                
        except Exception as e:
            print(f"\t- Error getting energies: {e}")
            vina_results['best_affinity'] = None
            all_poses_content = []

        # Determine which poses to save
        poses_to_save = []
        if save_poses.lower() == "none":
            poses_to_save = []
        elif save_poses.lower() == "first":
            poses_to_save = [1] if vina_results['results'] else []
        elif save_poses.lower() == "all":
            poses_to_save = list(range(1, len(vina_results['results']) + 1))
        elif "-" in save_poses:
            # Parse range like "1-5"
            try:
                start, end = map(int, save_poses.split("-"))
                poses_to_save = list(range(start, min(end + 1, len(vina_results['results']) + 1)))
            except:
                print(f"Invalid pose range format: {save_poses}. Using 'first'.")
                poses_to_save = [1] if vina_results['results'] else []
        else:
            poses_to_save = [1] if vina_results['results'] else []

        print(f"\n  Saving poses: {poses_to_save if poses_to_save else 'None'}")

        # Save poses
        vina_results['pose_files'] = []
        vina_results['complex_files'] = []

        if poses_to_save and vina_results.get('best_affinity') is not None and all_poses_content:
            # Save individual poses
            for pose_num in poses_to_save:
                try:
                    if pose_num <= len(all_poses_content):
                        # Save individual pose
                        pose_file = os.path.join(docking_output_dir, f"{base_name}_pose_{pose_num}.pdbqt")
                        with open(pose_file, 'w') as f:
                            f.write(all_poses_content[pose_num - 1])
                        vina_results['pose_files'].append(pose_file)
                        
                        # Save complex (receptor + pose)
                        complex_file = os.path.join(docking_output_dir, f"{base_name}_complex_pose_{pose_num}.pdbqt")
                        with open(complex_file, 'w') as out_f:
                            # Write receptor
                            with open(receptor_pdbqt, 'r') as rec_f:
                                out_f.write(rec_f.read())
                            
                            # Write pose (skip REMARK lines)
                            for line in all_poses_content[pose_num - 1].split('\n'):
                                if not line.startswith('REMARK'):
                                    out_f.write(line + '\n')
                        
                        vina_results['complex_files'].append(complex_file)
                        print(f"\t- Saved pose {pose_num}: {pose_file}")
                        print(f"\t- Saved complex {pose_num}: {complex_file}")
                    else:
                        print(f"\t- Pose {pose_num} not available (only {len(all_poses_content)} poses generated)")
                except Exception as e:
                    print(f"\t- Error saving pose {pose_num}: {e}")
            
            # Clean up temporary file
            if os.path.exists(temp_all_poses):
                os.remove(temp_all_poses)

            # Create zip file if multiple poses
            if len(poses_to_save) > 1 or save_poses.lower() == "all":
                try:
                    zip_file = os.path.join(docking_output_dir, f"{base_name}_all_poses.zip")
                    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
                        # Add all pose files
                        for pose_file in vina_results['pose_files']:
                            arcname = os.path.basename(pose_file)
                            zipf.write(pose_file, arcname)
                        
                        # Add all complex files
                        for complex_file in vina_results['complex_files']:
                            arcname = os.path.basename(complex_file)
                            zipf.write(complex_file, arcname)
                        
                        # Add log file
                        arcname = os.path.basename(log_file)
                        zipf.write(log_file, arcname)
                    
                    vina_results['zip_file'] = zip_file
                    print(f"\n  Created zip file: {zip_file}")
                except Exception as e:
                    print(f"\t- Error creating zip file: {e}")
                    vina_results['zip_file'] = None
        else:
            vina_results['zip_file'] = None

        # Write summary to CSV file
        try:
            summary_file = os.path.join(docking_output_dir, f"{base_name}_summary.csv")
            with open(summary_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write header
                writer.writerow(['Parameter', 'Value'])
                writer.writerow([''])
                
                # Write grid information
                writer.writerow(['Grid Center (X)', f"{vina_results['grid_center'][0]:.2f}"])
                writer.writerow(['Grid Center (Y)', f"{vina_results['grid_center'][1]:.2f}"])
                writer.writerow(['Grid Center (Z)', f"{vina_results['grid_center'][2]:.2f}"])
                writer.writerow(['Grid Box Size (X)', f"{vina_results['grid_box'][0]:.2f}"])
                writer.writerow(['Grid Box Size (Y)', f"{vina_results['grid_box'][1]:.2f}"])
                writer.writerow(['Grid Box Size (Z)', f"{vina_results['grid_box'][2]:.2f}"])
                writer.writerow([''])
                
                # Write minimization information
                writer.writerow(['Ligand Minimization', 'Score (kcal/mol)'])
                if vina_results['minimization']['before_min'] is not None:
                    writer.writerow(['Before Minimization', f"{vina_results['minimization']['before_min']:.3f}"])
                    writer.writerow(['After Minimization', f"{vina_results['minimization']['after_min']:.3f}"])
                    improvement = vina_results['minimization']['before_min'] - vina_results['minimization']['after_min']
                    writer.writerow(['Improvement', f"{improvement:.3f}"])
                else:
                    writer.writerow(['Minimization', 'Not performed'])
                writer.writerow([''])
                
                # Write docking parameters
                writer.writerow(['Docking Parameter', 'Value'])
                writer.writerow(['Exhaustiveness', exhaustiveness])
                writer.writerow(['Number of Modes', num_modes])
                writer.writerow(['CPU', cpu])
                writer.writerow(['Seed', seed])
                writer.writerow(['Elapsed Time (seconds)', f"{vina_results['elapsed_time']:.2f}"])
                writer.writerow([''])
                
                # Write results
                writer.writerow(['Pose', 'Affinity (kcal/mol)', 'RMSD_LB', 'RMSD_UB'])
                for result in vina_results['results']:
                    writer.writerow([
                        result['pose'],
                        f"{result['affinity']:.2f}",
                        f"{result['rmsd_lb']:.2f}" if result['rmsd_lb'] is not None else 'N/A',
                        f"{result['rmsd_ub']:.2f}" if result['rmsd_ub'] is not None else 'N/A'
                    ])
            
            vina_results['summary_file'] = summary_file
            print(f"\t- Summary table saved: {summary_file}")
        except Exception as e:
            print(f"\t- Error creating summary table: {e}")
            vina_results['summary_file'] = None

        print("\n  Docking results saved")
        print(f"\t- Output directory: {docking_output_dir}")
        print(f"\t- Log file: {log_file}")

        return vina_results, None, docking_output_dir

    finally:
        # Close log file
        sys.stdout = log_capture.terminal
        log_capture.close()

# Example usage
if __name__ == "__main__":
    # Example parameters
    estradiol_string = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O"
    receptor_filename = "5TOA.pdb"
    protein_to_extract = "Estrogen receptor beta"
    input_receptor = os.path.join('docking_functions', receptor_filename)
    output_dir = os.path.join('docking_functions')

    # 1. Convert Ligand SMILES to PDBQT
    try:
        # Prepare ligand from SMILES
        ligand_pdbqt = prepare_ligand_from_smiles(
            smiles=estradiol_string,
            output_filename='estradiol.pdbqt',
            output_dir=output_dir,
            ph=7.4,
            charge_method='gasteiger'
        )
    except Exception as e:
        print(e)

    # 2. Convert Receptor PDB file to PDBQT
    try:
        print(f"\n{'='*60}")
        print(f"Processing: {receptor_filename}")
        print(f"Extracting protein: {protein_to_extract}")
        print(f"{'='*60}\n")
        
        cleaned_pdb, receptor_pdbqt = clean_and_convert_pdb_to_pdbqt(
            input_receptor, 
            protein_name=protein_to_extract,
            deduplicate=True,
            similarity_threshold=80.0,
            ph=7.4,
            charge_method='gasteiger',
            remove_polar_hydrogens=True
        )
        print(f'Cleaned PDB file generated: {cleaned_pdb}')
        print(f'PDBQT file generated: {receptor_pdbqt}')
    except Exception as e:
        print(f"Error: {e}")

    try:
        # Get size
        grid_size = compute_bounding_box(receptor_pdbqt)
        # Get center
        grid_center = calculate_receptor_center_of_mass(receptor_pdbqt)

        results, error, log_file = run_vina_blind_docking(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            align_ligand=True,
            output_folder=output_dir,
            center=grid_center,
            size=grid_size,
            output_filename="estradiol_docking",
            ligand_minimized=False,
            exhaustiveness=64,
            num_modes=8,
            cpu=-1,
            seed=1234,
            save_poses="all" # Options: "none", "first", "all", "1-5", etc.
        )

        if error:
            print(f"Error: {error}")
        else:
            print(f"\nDocking Results:")
            print(f"\t- Best affinity: {results['best_affinity']:.2f} kcal/mol")
            print(f"\t- Number of poses saved: {len(results['pose_files'])}")
            print(f"\t- Log file: {log_file}")
            if results.get('zip_file'):
                print(f"\t- Zip file: {results['zip_file']}")

    except Exception as e:
        print(f"Exception: {e}")