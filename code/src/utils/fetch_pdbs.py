"""
Script to fetch a list of PDB files from the Protein Data Bank

Usage: python fetch_pdbs.py $PDB_IDS $OUT_PATH

$PDB_IDS: file containing PDB IDs. Each ID is expected to be separated by a
newline
$OUT_PATH: output path; directory should NOT exist
"""

from sys import argv
from pathlib import Path
from time import sleep
import prody as prd
from tqdm import tqdm

def fetch_pdbs(pdb_ids_f, out_path):
    with open(pdb_ids_f) as f:
        pdb_ids = f.readlines()
        pdb_ids = [x.strip() for x in pdb_ids]

    # TODO: add try/except
    # ... for what?
    for pdb_id in tqdm(pdb_ids):
        prd.fetchPDB(pdb_id, compressed=False, folder=str(out_path))
        sleep(0.1)


if __name__ == "__main__":

    pdb_ids_f = argv[1]
    out_path = Path(argv[2])
    out_path.mkdir(exist_ok=False)

    fetch_pdbs(pdb_ids_f, out_path)

    print("Et voila!")
