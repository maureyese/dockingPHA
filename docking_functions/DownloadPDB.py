import os
import time
from typing import Optional
import requests

# =============================
# FUNCTIONS
# =============================

def download_pdb_file(pdb_id: str, output_dir: str) -> Optional[str]:
    '''
    Download PDB file from RCSB PDB
    
    :param pdb_id: PDB ID (e.g., '1ABC')
    :type pdb_id: str
    :param output_dir: Directory to save PDB file
    :type output_dir: str
    :return: Path to downloaded file or None
    :rtype: str or None
    '''
    # Try multiple URL formats
    url_formats = [
        f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb",
        f"https://files.rcsb.org/view/{pdb_id}.pdb",
        f"https://www.rcsb.org/structure/{pdb_id}",  # This might redirect to download
    ]
    
    for url in url_formats:
        time.sleep(0.5)
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                # Check if the response is actually a PDB file
                content = response.text
                if content and ('HEADER' in content or 'ATOM' in content or 'REMARK' in content):
                    pdb_file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
                    
                    with open(pdb_file_path, "w") as f:
                        f.write(content)
                    
                    print(f"  Downloaded PDB file: {pdb_file_path}")
                    return pdb_file_path
                
        except requests.exceptions.RequestException:
            continue
    
    print(f"  Failed to download PDB {pdb_id}")
    return None

# =============================
# EXAMPLE USAGE
# =============================

def main():
    pdb_input = "2BXP" # Add ID from PDB
    output_dir = "docking_functions/" # Write output folder

    # Run function
    try:
        pdb_file_path = download_pdb_file(pdb_id=pdb_input, output_dir=output_dir)

        print(f"PDB file {pdb_input} stored at {pdb_file_path}")
    except Exception as e:
        print(e)

if __name__ == "__main__":
    main()