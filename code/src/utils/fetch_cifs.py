"""
Script to fetch a list of CIF files from the Protein Data Bank

Usage: python fetch_pdbs.py $PDB_IDS $OUT_PATH

$PDB_IDS: file containing PDB IDs. Each ID is expected to be separated by a
newline
$OUT_PATH: output path; directory should NOT exist
"""

from sys import argv
from pathlib import Path
from time import sleep
import subprocess
from tqdm import tqdm


def query_ebi(pdb_id, cif_out_path):

    assert len(pdb_id) == 4

    cif_url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{pdb_id.lower()}_updated.cif"
    cif_file = cif_out_path / f"{pdb_id}.cif"
    cif_call = ["wget", cif_url, "-O", str(cif_file)]

    # Skip examples we have already downloaded
    skipped = 0
    downloaded = 0
    if not cif_file.is_file():
        try:
            subprocess.check_call(cif_call)
            sleep(0.1)
            downloaded = 1
        except subprocess.CalledProcessError as error:
            raise ValueError(f"wget failed: {error}")
        except OSError:
            raise OSError("wget binary not found in PATH")
    else:
        skipped = 1
    
    return skipped, downloaded
    

def fetch_cifs(pdb_ids_f, out_path):
    with open(pdb_ids_f) as f:
        pdb_ids = f.readlines()
        pdb_ids = [x.strip() for x in pdb_ids]

    skipped = 0
    downloaded = 0
    for pdb_id in tqdm(pdb_ids):
        was_skipped, was_downloaded = query_ebi(pdb_id, out_path)
        skipped += was_skipped
        downloaded += was_downloaded
        sleep(0.1)

    print(f'Skipped {skipped} files')
    print(f'Downloaded {downloaded} files')



if __name__ == "__main__":

    pdb_ids_f = argv[1]
    out_path = Path(argv[2])
    out_path.mkdir(exist_ok=True)

    fetch_cifs(pdb_ids_f, out_path)

    print("Et voila!")




