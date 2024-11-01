"""
Script to perform quality control of Pfam domain structures before creating
structural embeddings
"""

import os
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
import pandas as pd
from tqdm import tqdm

from utils.structure import find_abnormal_residues
from warnings import warn

global RESOLUTION_CUTOFF
RESOLUTION_CUTOFF = 3
global COVERAGE_CUTOFF
COVERAGE_CUTOFF = 0.75
global ONLY_CALPHA_CUTOFF
ONLY_CALPHA_CUTOFF = 0.1

def write_structure(pdb_id, structure, out_file):
    try:
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_file)
    except AttributeError:
        # biopython bug handling disordered atoms
        # Simply remove empty file and continue
        print(f"Error on {pdb_id}")
        os.remove(out_file)


if __name__ == "__main__":

    #########
    # Setup #
    #########

    base_path = Path("../../data/processed/pfam_structures")
    pfam_df = pd.read_csv("../../data/processed/pfam_structures/merged_pfam_data.tsv",delimiter='\t')

    input_path = base_path /"extracted_domains"
    output_path = base_path / "extracted_domains_filtered"
    output_path.mkdir(exist_ok=True)

    #############
    # Filtering #
    #############

    # Filter by resolution and coverage
    pfam_df = pfam_df.loc[(pfam_df['Coverage']>COVERAGE_CUTOFF) & (pfam_df['Resolution']<RESOLUTION_CUTOFF)]
    # Do not allow domains without strictly defined boundaries in SIFTS
    pfam_df = pfam_df.loc[(pfam_df["Domain end"]!='None') & (pfam_df["Domain end"]!='None')]

    pfam_df['combined_id'] = pfam_df['PDB ID'] + "_" + pfam_df['Chain ID']
    domain_names = list(pfam_df["Domain ID"].unique())

    final_set = []
    columns = ["Uniprot ID","PDB ID","Chain ID","Domain ID","Is phosphorylated","Phosphosites"]

    parser = PDBParser(QUIET=True)
    for domain_name in tqdm(domain_names):
        pfam_df_subset = pfam_df.loc[pfam_df["Domain ID"]==domain_name]

        domain_input_path = input_path / domain_name
        domain_output_path = output_path / domain_name
        domain_output_path.mkdir(exist_ok=True)

        for idx, row in tqdm(pfam_df_subset.iterrows()):
            uniprot_id = row["Uniprot ID"]
            is_phosphorylated = row["Is phosphorylated"]
            psites = row["Phosphosites"]
            pdb_id = row['combined_id']
            pdb_code, chain_id = pdb_id.split("_")
            in_file = domain_input_path / f"{pdb_id}.pdb"
            
            structure = parser.get_structure("X",in_file)
            if len(structure) > 1:
                # If there are multiple models, select the first one
                structure = structure[0]
                chain = structure[chain_id]
            else:
                chain = structure[0][chain_id]

            # Check for residues containing only one atom
            # (normally this is an alpha carbon)
            # These cause problems downstream for Geometricus
            fraction_abnormal_residues = find_abnormal_residues(chain)

            if fraction_abnormal_residues < ONLY_CALPHA_CUTOFF:
                out_file = domain_output_path / f"{pdb_id}.pdb" 
                write_structure(pdb_id, structure, str(out_file))
                filt_row = [uniprot_id, pdb_code, chain_id, domain_name, is_phosphorylated, psites]
                final_set.append(filt_row)
            else:
                warn(f"Skipped {pdb_id} (domain {domain_name}, {fraction_abnormal_residues} abnormal residues)",RuntimeWarning)

    final_set = pd.DataFrame(final_set,columns=columns)
    print(f"Final set of structures: {final_set}")
    final_set.to_csv(base_path / "filtered_pfam_structures.tsv",sep='\t')
    print("Done!")


            



