"""
Create dataframe with annotation on the phospho state of each protein.

This is used as input for bio3D's PCA function
"""
from sys import argv
from pathlib import Path
import pandas as pd

from math import isnan

if __name__ == "__main__":

    pdb_input_dir = Path(argv[1])
    cluster_df = Path(argv[2])

    # Sort IDs alphabetically
    pdb_files = sorted(list(pdb_input_dir.glob("*.pdb")))
    # Expand filenames to the full path
    pdb_files = [pdb_file.resolve() for pdb_file in pdb_files]

    # Read input file containing info about the phospho state
    cluster_df = pd.read_csv(argv[2], delimiter='\t')
    cols = ["pdb_id", "chain_id"]
    cluster_df['combined_id'] = cluster_df[cols].apply(
        lambda row: '_'.join(row.values.astype(str)), axis=1)

    header = ["pdb_id", "file_path", "state", "color"]
    annot_df = []

    for pdb_file in pdb_files:

        pdb_id = pdb_file.stem

        match = cluster_df.loc[cluster_df["combined_id"] == pdb_id]
        if len(match) == 0:
            prot_id, chain = pdb_id.split("_")
            prot_id = prot_id.lower()
            pdb_id = f"{prot_id}_{chain}"
            match = cluster_df.loc[cluster_df["combined_id"] == pdb_id]

        phospho_state = match["phospho_idxs"]

        if isinstance(phospho_state, str):
            # Multiple phosphorylations
            annot_df.append([pdb_id, pdb_file, "Phosphorylated", "orange"])
        elif isnan(phospho_state):
            # Non-phosphorylated
            annot_df.append([pdb_id, pdb_file, "Non-phosphorylated", "blue"])
        else:
            # One phosphorylation
            annot_df.append([pdb_id, pdb_file, "Phosphorylated", "orange"])

    annot_df = pd.DataFrame(annot_df, columns=header)
    annot_df.to_csv("annotation_sulfatase.csv", index=False)
    print("All done here!")
