"""
Auxiliary script to create a separate directory with the structures for 
each protein (data/processed/pdb_pairs/chains_by_protein) to facilitate
downstream analysis.
"""

import pandas as pd
import shutil
from pathlib import Path
from collections import defaultdict


if __name__ == "__main__":

    data_path = Path("../../data/processed/pdb_pairs")
    source_path = data_path / "extracted_chains" / "pdb"
    dest_path = data_path / "chains_by_protein"
    dest_path.mkdir(exist_ok=True)

    # Read dataframe
    df = pd.read_csv(data_path / "filtered_df.csv")
    protein_to_chains = defaultdict(set)
    for idx, row in df.iterrows():
        protein_id = row["UNIPROT"]
        chain_a = row["RENAMED_ENTITY_ONE"]
        chain_b = row["RENAMED_ENTITY_TWO"]
        protein_to_chains[protein_id].add(chain_a)
        protein_to_chains[protein_id].add(chain_b)
    
    for protein_id, chains in protein_to_chains.items():
        protein_path = dest_path / protein_id
        protein_path.mkdir(exist_ok=True)

        for chain in chains:
            orig_file = source_path / f"{chain}.pdb"
            shutil.copy(orig_file, protein_path)

    print("Done!")

