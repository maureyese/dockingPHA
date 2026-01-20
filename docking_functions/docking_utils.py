import pandas as pd
import os

def extract_pose_results(vina_results, protein_file, ligand_file, output_name=None, poses='first'):
    """
    Extract and format pose results from Vina docking results.
    
    Parameters:
    -----------
    vina_results : dict
        Dictionary containing Vina docking results
    protein_file : str
        Path to the protein PDBQT file
    ligand_file : str
        Path to the ligand PDBQT file
    output_name : str, optional
        Name of the output folder/experiment
    poses : str, optional
        Which poses to extract: 'first' (default), 'all'
        
    Returns:
    --------
    dict or list of dict
        Dictionary with extracted pose results, or list of dictionaries if poses='all'
        Returns None if no valid results found
    """
    
    if not vina_results or 'results' not in vina_results or len(vina_results['results']) == 0:
        print(f"No valid docking results found for {os.path.basename(protein_file)}")
        return None
    
    protein_filename = os.path.basename(protein_file)
    ligand_filename = os.path.basename(ligand_file)
    
    if output_name is None:
        output_name = f"{os.path.splitext(protein_filename)[0]}_{os.path.splitext(ligand_filename)[0]}"
    
    grid_center = vina_results.get('grid_center', [None, None, None])
    grid_box = vina_results.get('grid_box', [None, None, None])
    
    def create_result_dict(pose_data):
        return {
            'protein_file': protein_filename,
            'ligand_file': ligand_filename,
            'output_name': output_name,
            'pose_number': pose_data['pose'],
            'affinity_kcal_mol': float(pose_data['affinity']),
            'rmsd_lb': pose_data.get('rmsd_lb'),
            'rmsd_ub': pose_data.get('rmsd_ub'),
            'grid_center_x': grid_center[0],
            'grid_center_y': grid_center[1],
            'grid_center_z': grid_center[2],
            'grid_size_x': grid_box[0],
            'grid_size_y': grid_box[1],
            'grid_size_z': grid_box[2],
            'elapsed_time_seconds': vina_results.get('elapsed_time'),
            'minimization_before': vina_results.get('minimization', {}).get('before_min'),
            'minimization_after': vina_results.get('minimization', {}).get('after_min'),
            'zip_file': vina_results.get('zip_file'),
            'output_folder': f"output/docking/{output_name}"
        }
    
    if poses == 'first':
        return create_result_dict(vina_results['results'][0])
    elif poses == 'all':
        return [create_result_dict(pose) for pose in vina_results['results']]
    else:
        print(f"Invalid poses parameter: {poses}. Use 'first' or 'all'")
        return None


def process_docking_results(vina_results_list, protein_files, ligand_files, output_names=None):
    """
    Process multiple docking results and extract poses from each.
    
    Parameters:
    -----------
    vina_results_list : list
        List of Vina docking result dictionaries
    protein_files : list
        List of protein file paths
    ligand_files : list
        List of ligand file paths
    output_names : list, optional
        List of output names for each docking
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing extracted pose results for all dockings
    """
    
    if len(vina_results_list) != len(protein_files) or len(vina_results_list) != len(ligand_files):
        raise ValueError("All input lists must have the same length")
    
    if output_names is not None and len(output_names) != len(vina_results_list):
        raise ValueError("output_names list must have the same length as other lists")
    
    docking_results = []
    
    for i, vina_results in enumerate(vina_results_list):
        protein_file = protein_files[i]
        ligand_file = ligand_files[i]
        output_name = output_names[i] if output_names else None
        
        result_dict = extract_pose_results(
            vina_results, 
            protein_file, 
            ligand_file, 
            output_name,
            poses='all'
        )
        
        if result_dict:
            if isinstance(result_dict, list):
                docking_results.extend(result_dict)
            else:
                docking_results.append(result_dict)
            print(f"- Successfully extracted results for {os.path.basename(protein_file)}")
    
    if docking_results:
        return pd.DataFrame(docking_results)
    else:
        return pd.DataFrame()


def export_docking_results(df_results, output_filename="docking_results.csv", 
                          include_timestamp=True, format="csv"):
    """
    Export docking results to a file.
    
    Parameters:
    -----------
    df_results : pandas.DataFrame
        DataFrame containing docking results
    output_filename : str, optional
        Base name for output file
    include_timestamp : bool, optional
        Whether to include timestamp in filename
    format : str, optional
        Output format: "csv", "excel", or "json"
        
    Returns:
    --------
    str
        Path to the exported file
    """
    
    if df_results.empty:
        print("No results to export")
        return None
    
    if include_timestamp:
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base, ext = os.path.splitext(output_filename)
        output_filename = f"{base}_{timestamp}{ext}"
    
    if format.lower() == "csv":
        df_results.to_csv(output_filename, index=False)
    elif format.lower() == "excel":
        df_results.to_excel(output_filename, index=False)
    elif format.lower() == "json":
        df_results.to_json(output_filename, orient='records', indent=2)
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'csv', 'excel', or 'json'.")
    
    print(f"Results exported to: {output_filename}")
    return output_filename