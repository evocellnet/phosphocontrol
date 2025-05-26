
#!/usr/bin/env python3

"""
Auxiliary script to create a separate directory with the structures for 
each protein  to facilitate
downstream analysis.

This auxiliary script creates a new directory structure where PDB chain files
are grouped by UniProt identifiers (data/processed/pdb_pairs/chains_by_protein).
This is intended to facilitate downstream analyses that operate on a per-protein
basis.

Usage:
    python organize_chains_by_protein.py

Inputs
------
    filtered_df: CSV file, dataset of paired phosphorylated and non-phosphorylated
                 structures
    ../../data/processed/pdb_pairs/extracted_chains/pdb/: directory with PDB files
                 for individual chains
Output
------
    - ../../data/processed/pdb_pairs/chains_by_protein/: directory where each
                  subdirectory is named after a UniProt ID, containing all relevant
                  PDB chain files

Author: Miguel Correa Marrero
"""

import pandas as pd
import shutil
from pathlib import Path
from collections import defaultdict


if __name__ == "__main__":

    # Define IO paths
    data_path = Path("../../data/processed/pdb_pairs")
    source_path = data_path / "extracted_chains" / "pdb"
    dest_path = data_path / "chains_by_protein"
    dest_path.mkdir(exist_ok=True)

    # Read dataframe of paired structures
    df = pd.read_csv(data_path / "filtered_df.csv")
    protein_to_chains = defaultdict(set)
    for idx, row in df.iterrows():
        protein_id = row["UNIPROT"]
        chain_a = row["RENAMED_ENTITY_ONE"]
        chain_b = row["RENAMED_ENTITY_TWO"]
        protein_to_chains[protein_id].add(chain_a)
        protein_to_chains[protein_id].add(chain_b)
    
    # Copy chain files into directories grouped by protein ID
    for protein_id, chains in protein_to_chains.items():
        protein_path = dest_path / protein_id
        protein_path.mkdir(exist_ok=True)

        for chain in chains:
            orig_file = source_path / f"{chain}.pdb"
            shutil.copy(orig_file, protein_path)

    print("Done!")

